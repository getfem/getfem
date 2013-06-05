/*===========================================================================
 
 Copyright (C) 2004-2012 Yves Renard, Konstantinos Poulios.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/


#include "getfem/getfem_contact_and_friction_nodal.h"
#include "getfem/getfem_contact_and_friction_common.h"
#include "getfem/getfem_assembling.h"

#ifndef GETFEM_HAVE_QHULL_QHULL_H
#include <getfem/bgeot_kdtree.h>
#endif

namespace getfem {

  typedef bgeot::convex<base_node>::dref_convex_pt_ct dref_convex_pt_ct;
  typedef bgeot::basic_mesh::ref_mesh_face_pt_ct ref_mesh_face_pt_ct;


  // Computation of an orthonormal basis to a unit vector.
  static void orthonormal_basis_to_unit_vec(size_type d, const base_node &un,
                                            base_node *ut) {
    size_type n = 0;
    for (size_type k = 0; k <= d && n < d; ++k) {
      gmm::resize(ut[n], d+1);
      gmm::clear(ut[n]);
      ut[n][k] = scalar_type(1);

      ut[n] -= gmm::vect_sp(un, ut[n]) * un;
      for (size_type nn = 0; nn < n; ++nn)
        ut[n] -= gmm::vect_sp(ut[nn], ut[n]) * ut[nn];

      if (gmm::vect_norm2(ut[n]) < 1e-3) continue;
      ut[n] /= gmm::vect_norm2(ut[n]);
      ++n;
    }
    GMM_ASSERT1(n == d, "Gram-Schmidt algorithm to find an "
      "orthonormal basis for the tangential displacement failed");
  }

  // "contact_node" is an object which contains data about nodes expected
  // to participate in a contact condition. A contact node refers to a
  // specific mesh_fem.
  struct contact_node {
    const mesh_fem *mf;          // Pointer to the mesh_fem the contact node is
                                 // associated with
    size_type dof;               // first dof id of the node in the considered mesh_fem
    std::vector<size_type> cvs;  // list of ids of neigbouring convexes
    std::vector<short_type> fcs; // list of local ids of neigbouring faces

    contact_node() : mf(0), cvs(0), fcs(0) {}
    contact_node(const mesh_fem &mf_) {mf = &mf_;}
  };

  // contact_node's pair
  struct contact_node_pair {
    contact_node cn_s, cn_m;  // Slave and master contact_node's
    scalar_type dist2;        // Square of distance between slave and master nodes
    bool is_active;
    contact_node_pair(scalar_type threshold=10.) : cn_s(), cn_m()
      {dist2 = threshold * threshold; is_active = false;}
  };

  // contact_node's pair list
  class contact_node_pair_list : public std::vector<contact_node_pair> {

    void contact_node_list_from_region
      (const mesh_fem &mf, size_type contact_region,
       std::vector<contact_node> &cnl) {

      cnl.clear();
      const mesh &m = mf.linked_mesh();
      size_type qdim = mf.get_qdim();
      std::map<size_type, size_type> dof_to_cnid;
      size_type cnid = 0;
      dal::bit_vector dofs = mf.basic_dof_on_region(contact_region);
      for (dal::bv_visitor dof(dofs); !dof.finished(); ++dof)
        if ( dof % qdim == 0) {
          dof_to_cnid[dof] = cnid++;
          contact_node new_cn(mf);
          new_cn.dof = dof;
          cnl.push_back(new_cn);
        }
      for (mr_visitor face(m.region(contact_region));
           !face.finished(); ++face) {
        assert(face.is_face());
        mesh_fem::ind_dof_face_ct
          face_dofs = mf.ind_basic_dof_of_face_of_element(face.cv(),face.f());
        for (size_type it=0; it < face_dofs.size(); it += qdim ) {
          size_type dof = face_dofs[it];
          cnid = dof_to_cnid[dof];
          cnl[cnid].cvs.push_back(face.cv());
          cnl[cnid].fcs.push_back(face.f());
        } // for:it
      } // for:face
    } // append

  public:
    contact_node_pair_list() : std::vector<contact_node_pair>() {}

    void append_min_dist_cn_pairs(const mesh_fem &mf1, const mesh_fem &mf2,
                                  size_type rg1, size_type rg2,
                                  bool slave1=true, bool slave2=false) {

      std::vector<contact_node> cnl1(0), cnl2(0);
      contact_node_list_from_region(mf1, rg1, cnl1);
      contact_node_list_from_region(mf2, rg2, cnl2);

      // Find minimum distance node pairs
      size_type size0 = this->size();
      size_type size1 = slave1 ? cnl1.size() : 0;
      size_type size2 = slave2 ? cnl2.size() : 0;
      this->resize( size0 + size1 + size2 );
# ifndef GETFEM_HAVE_QHULL_QHULL_H
      bgeot::kdtree tree1, tree2;
      for (size_type i1 = 0; i1 < cnl1.size(); ++i1) {
        contact_node *cn1 = &cnl1[i1];
        tree1.add_point_with_id(cn1->mf->point_of_basic_dof(cn1->dof), i1);
      }
      for (size_type i2 = 0; i2 < cnl2.size(); ++i2) {
        contact_node *cn2 = &cnl2[i2];
        tree2.add_point_with_id(cn2->mf->point_of_basic_dof(cn2->dof), i2);
      }
      if (slave1) {
        size_type ii1=size0;
        for (size_type i1 = 0; i1 < cnl1.size(); ++i1, ++ii1) {
          contact_node *cn1 = &cnl1[i1];
          base_node node1 = cn1->mf->point_of_basic_dof(cn1->dof);
          bgeot::index_node_pair ipt;
          scalar_type dist2 = tree2.nearest_neighbor(ipt, node1);
          if (ipt.i >= 0 && dist2 < (*this)[ii1].dist2) {
            (*this)[ii1].cn_s = *cn1;
            (*this)[ii1].cn_m = cnl2[ipt.i];
            (*this)[ii1].dist2 = dist2;
            (*this)[ii1].is_active = true;
          }
        }
      }
      if (slave2) {
        size_type ii2=size0+size1;
        for (size_type i2 = 0; i2 < cnl2.size(); ++i2, ++ii2) {
          contact_node *cn2 = &cnl2[i2];
          base_node node2 = cn2->mf->point_of_basic_dof(cn2->dof);
          bgeot::index_node_pair ipt;
          scalar_type dist2 = tree1.nearest_neighbor(ipt, node2);
          if (ipt.i >= 0 && dist2 < (*this)[ii2].dist2) {
            (*this)[ii2].cn_s = *cn2;
            (*this)[ii2].cn_m = cnl1[ipt.i];
            (*this)[ii2].dist2 = dist2;
            (*this)[ii2].is_active = true;
          }
        }
      }
# else
      std::vector<base_node> pts;
      for (size_type i1 = 0; i1 < cnl1.size(); ++i1) {
        contact_node *cn1 = &cnl1[i1];
        pts.push_back(cn1->mf->point_of_basic_dof(cn1->dof));
      }
      for (size_type i2 = 0; i2 < cnl2.size(); ++i2) {
        contact_node *cn2 = &cnl2[i2];
        pts.push_back(cn2->mf->point_of_basic_dof(cn2->dof));
      }
      gmm::dense_matrix<size_type> simplexes;

      getfem::delaunay(pts, simplexes);

      size_type nb_vertices = gmm::mat_nrows(simplexes);
      std::vector<size_type> facet_vertices(nb_vertices);
      std::vector< std::vector<size_type> > pt1_neighbours(size1);
      for (size_type i = 0; i < gmm::mat_ncols(simplexes); ++i) {
        gmm::copy(gmm::mat_col(simplexes, i), facet_vertices);
        for (size_type iv1 = 0; iv1 < nb_vertices-1; ++iv1) {
          size_type v1 = facet_vertices[iv1];
          bool v1_on_surface1 = (v1 < size1);
          for (size_type iv2 = iv1 + 1; iv2 < nb_vertices; ++iv2) {
            size_type v2 = facet_vertices[iv2];
            bool v2_on_surface1 = (v2 < size1);
            if (v1_on_surface1 ^ v2_on_surface1) {
              bool already_in = false;
              size_type vv1 = (v1_on_surface1 ? v1 : v2);
              size_type vv2 = (v2_on_surface1 ? v1 : v2);
              for (size_type j = 0; j < pt1_neighbours[vv1].size(); ++j)
                if (pt1_neighbours[vv1][j] == vv2) {
                  already_in = true;
                  break;
                }
              if (!already_in) pt1_neighbours[vv1].push_back(vv2);
            }
          }
        }
      }

      for (size_type i1 = 0; i1 < size1; ++i1)
        for (size_type j = 0; j < pt1_neighbours[i1].size(); ++j) {
          size_type i2 = pt1_neighbours[i1][j] - size1;
          size_type ii1 = size0 + i1;
          size_type ii2 = size0 + size1 + i2;
          contact_node *cn1 = &cnl1[i1];
          base_node node1 = cn1->mf->point_of_basic_dof(cn1->dof);
          contact_node *cn2 = &cnl2[i2];
          base_node node2 = cn2->mf->point_of_basic_dof(cn2->dof);
          scalar_type dist2 = gmm::vect_norm2_sqr(node1-node2);
          if (slave1 && dist2 < (*this)[ii1].dist2) {
            (*this)[ii1].cn_s = *cn1;
            (*this)[ii1].cn_m = *cn2;
            (*this)[ii1].dist2 = dist2;
            (*this)[ii1].is_active = true;
          }
          if (slave2 && dist2 < (*this)[ii2].dist2) {
            (*this)[ii2].cn_s = *cn2;
            (*this)[ii2].cn_m = *cn1;
            (*this)[ii2].dist2 = dist2;
            (*this)[ii2].is_active = true;
          }
        }
#endif
    }

    void append_min_dist_cn_pairs(const mesh_fem &mf,
                                  size_type rg1, size_type rg2,
                                  bool slave1=true, bool slave2=false) {
      append_min_dist_cn_pairs(mf, mf, rg1, rg2, slave1, slave2);
    }
  };

  scalar_type projection_on_convex_face
    (const mesh &m, const size_type cv, const short_type fc,
     const base_node &master_node, const base_node &slave_node,
     base_node &un, base_node &proj_node, base_node &proj_node_ref) {

    bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);

    if (pgt->is_linear()) {  //this condition is practically too strict

      un = m.normal_of_face_of_convex(cv,fc);
      un /= gmm::vect_norm2(un);
      //proj_node = slave_node - [(slave_node-master_node)*n] * n
      gmm::add(master_node, gmm::scaled(slave_node, -1.), proj_node);
      gmm::copy(gmm::scaled(un, gmm::vect_sp(proj_node, un)), proj_node);
      gmm::add(slave_node, proj_node);

      bgeot::geotrans_inv_convex gic;
      gic.init(m.points_of_convex(cv), pgt);
      gic.invert(proj_node, proj_node_ref);
      return pgt->convex_ref()->is_in(proj_node_ref);

    } else {

      size_type N = m.dim();
      size_type P = pgt->structure()->dim();

      size_type nb_pts_cv = pgt->nb_points();
      size_type nb_pts_fc = pgt->structure()->nb_points_of_face(fc);

      bgeot::convex_ind_ct ind_pts_fc = pgt->structure()->ind_points_of_face(fc);
      ref_mesh_face_pt_ct pts_fc = m.points_of_face_of_convex(cv, fc);

      // Local base on reference face
      base_matrix base_ref_fc(P-1,N);
      {
        dref_convex_pt_ct dref_pts_fc = pgt->convex_ref()->dir_points_of_face(fc);
        GMM_ASSERT1( dref_pts_fc.size() == P, "Dimensions mismatch");
        base_node vec(dref_pts_fc[0].size());
        for (size_type i = 0; i < P-1; ++i) {
          vec = dref_pts_fc[i+1] - dref_pts_fc[0];
          gmm::copy(vec,gmm::mat_row(base_ref_fc,i));
        }
      }

      GMM_ASSERT1( slave_node.size() == N, "Dimensions mismatch");
      const base_node &xx = slave_node;
      base_node &xxp = proj_node;  xxp.resize(N);
      base_node &xp = proj_node_ref;  xp.resize(P);
      base_node vres(P);
      scalar_type res= 1.;

      // initial guess
      xp = gmm::mean_value(pgt->convex_ref()->points_of_face(fc));

      gmm::clear(xxp);
      base_vector val(nb_pts_fc);
      pgt->poly_vector_val(xp, ind_pts_fc, val);
      for (size_type l = 0; l < nb_pts_fc; ++l)
        gmm::add(gmm::scaled(pts_fc[l], val[l] ), xxp);

      base_matrix G(N, nb_pts_fc);
      vectors_to_base_matrix(G, pts_fc);

      base_matrix K(N,P-1);

      base_matrix grad_fc(nb_pts_fc, P);
      base_matrix grad_fc1(nb_pts_fc, P-1);
      base_matrix B(N,P-1), BB(N,P), CS(P-1,P-1);

      scalar_type EPS = 10E-12;
      unsigned cnt = 50;
      while (res > EPS && --cnt) {
        // computation of the pseudo inverse matrix B at point xp
        pgt->poly_vector_grad(xp, ind_pts_fc, grad_fc);
        gmm::mult(grad_fc, gmm::transposed(base_ref_fc), grad_fc1);
        gmm::mult(G, grad_fc1, K);
        gmm::mult(gmm::transposed(K), K, CS);
        gmm::lu_inverse(CS);
        gmm::mult(K, CS, B);
        gmm::mult(B, base_ref_fc, BB);

        // Projection onto the face of the convex
        gmm::mult_add(gmm::transposed(BB), xx-xxp, xp);
        gmm::clear(xxp);
        pgt->poly_vector_val(xp, ind_pts_fc, val);
        for (size_type l = 0; l < nb_pts_fc; ++l)
          gmm::add(gmm::scaled(pts_fc[l], val[l]), xxp);

        gmm::mult(gmm::transposed(BB), xx - xxp, vres);
        res = gmm::vect_norm2(vres);
      }
      GMM_ASSERT1( res <= EPS,
                  "Iterative pojection on convex face did not converge");
      { // calculate K at the final point
        pgt->poly_vector_grad(xp, ind_pts_fc, grad_fc);
        gmm::mult(grad_fc, gmm::transposed(base_ref_fc), grad_fc1);
        gmm::mult(G, grad_fc1, K);
      }

      // computation of normal vector
      un.resize(N);
      // un = xx - xxp;
      // gmm::scale(un, 1/gmm::vect_norm2(un));
      gmm::clear(un);
      {
        base_matrix KK(N,P);
        { // calculate KK
          base_matrix grad_cv(nb_pts_cv, P);
          pgt->poly_vector_grad(xp, grad_cv);

          base_matrix GG(N, nb_pts_cv);
          vectors_to_base_matrix(GG, m.points_of_convex(cv));

          gmm::mult(GG, grad_cv, KK);
        }

        base_matrix bases_product(P-1, P);
        gmm::mult(gmm::transposed(K), KK, bases_product);

        for (size_type i = 0; i < P; ++i) {
          std::vector<size_type> ind(0);
          for (size_type j = 0; j < P; ++j)
            if (j != i ) ind.push_back(j);
          scalar_type det = gmm::lu_det(gmm::sub_matrix(bases_product,
                                                        gmm::sub_interval(0, P-1),
                                                        gmm::sub_index(ind)      ) );
          gmm::add(gmm::scaled(gmm::mat_col(KK, i), (i % 2) ? -det : +det ), un);
        }
      }
      // normalizing
      gmm::scale(un, 1/gmm::vect_norm2(un));
      // ensure that normal points outwards
      if (gmm::vect_sp(un, gmm::mean_value(pts_fc) -
                           gmm::mean_value(m.points_of_convex(cv))) < 0)
        gmm::scale(un,scalar_type(-1));

      return pgt->convex_ref()->is_in(proj_node_ref);
    }
  }

  void compute_contact_matrices
         (const mesh_fem &mf_disp1, const mesh_fem &mf_disp2,
          contact_node_pair_list &cnpl, model_real_plain_vector &gap,
          CONTACT_B_MATRIX *BN1, CONTACT_B_MATRIX *BN2 = 0,
          CONTACT_B_MATRIX *BT1 = 0, CONTACT_B_MATRIX *BT2 = 0) {

    GMM_ASSERT1(gmm::vect_size(gap) == cnpl.size(),
                "Wrong number of contact node pairs or wrong size of gap");
    gmm::clear(*BN1);
    GMM_ASSERT1( gmm::mat_nrows(*BN1) == cnpl.size(), "Wrong size of BN1");
    if (BN2) {
      gmm::clear(*BN2);
      GMM_ASSERT1( gmm::mat_nrows(*BN2) == cnpl.size(), "Wrong size of BN2");
    }
    dim_type qdim = mf_disp1.get_qdim();
    size_type d = qdim - 1;
    if (BT1) {
      gmm::clear(*BT1);
      GMM_ASSERT1( gmm::mat_nrows(*BT1) == cnpl.size() * d, "Wrong size of BT1");
    }
    if (BT2) {
      gmm::clear(*BT2);
      GMM_ASSERT1( gmm::mat_nrows(*BT2) == cnpl.size() * d, "Wrong size of BT2");
    }
    gmm::fill(gap, scalar_type(10));  //FIXME: Needs a threshold value
    for (size_type row = 0; row < cnpl.size(); ++row) {
      contact_node_pair *cnp = &cnpl[row];
      if (cnp->is_active) {
        contact_node *cn_s = &cnp->cn_s;  //slave contact node
        contact_node *cn_m = &cnp->cn_m;  //master contact node
        const mesh &mesh_m = cn_m->mf->linked_mesh();
        base_node slave_node = cn_s->mf->point_of_basic_dof(cn_s->dof);
        base_node master_node = cn_m->mf->point_of_basic_dof(cn_m->dof);
        base_node un_sel(3), proj_node_sel(3), proj_node_ref_sel(3);
        scalar_type is_in_min = 1e5;  //FIXME
        size_type cv_sel = 0, fc_sel = 0;
        std::vector<size_type>::iterator cv;
        std::vector<short_type>::iterator fc;
        for (cv = cn_m->cvs.begin(), fc = cn_m->fcs.begin();
             cv != cn_m->cvs.end() && fc != cn_m->fcs.end(); cv++, fc++) {
          base_node un(3), proj_node(3), proj_node_ref(3);
          scalar_type is_in = projection_on_convex_face
            (mesh_m, *cv, *fc, master_node, slave_node, un, proj_node, proj_node_ref);
          if (is_in < is_in_min) {
            is_in_min = is_in;
            cv_sel = *cv;
            fc_sel = *fc;
            un_sel = un;
            proj_node_sel = proj_node;
            proj_node_ref_sel = proj_node_ref;
          }
        }
        if (is_in_min < 0.05) {  //FIXME
          gap[row] = gmm::vect_sp(slave_node-proj_node_sel, un_sel);

          base_node ut[3];
          if (BT1) orthonormal_basis_to_unit_vec(d, un_sel, ut);

          CONTACT_B_MATRIX *BN = 0;
          CONTACT_B_MATRIX *BT = 0;
          if (cn_s->mf == &mf_disp1) {
            BN = BN1;
            BT = BT1;
          } else if (cn_s->mf == &mf_disp2) {
            BN = BN2;
            BT = BT2;
          }
          if (BN)
            for (size_type k = 0; k <= d; ++k)
              (*BN)(row, cn_s->dof + k) -= un_sel[k];
          if (BT)
            for (size_type k = 0; k <= d; ++k)
              for (size_type n = 0; n < d; ++n)
                (*BT)(row * d + n, cn_s->dof + k) -= ut[n][k];

          BN = 0;
          const mesh_fem *mf_disp = 0;
          if (cn_m->mf == &mf_disp1) {
            BN = BN1;
            BT = BT1;
            mf_disp = &mf_disp1;
          } else if (cn_m->mf == &mf_disp2) {
            BN = BN2;
            BT = BT2;
            mf_disp = &mf_disp2;
          }
          if (BN) {
            base_matrix G;
            base_matrix M(qdim, mf_disp->nb_basic_dof_of_element(cv_sel));
            bgeot::vectors_to_base_matrix(G, mesh_m.points_of_convex(cv_sel));
            pfem pf = mf_disp->fem_of_element(cv_sel);
            bgeot::pgeometric_trans pgt = mesh_m.trans_of_convex(cv_sel);
            fem_interpolation_context
              ctx(pgt, pf, proj_node_ref_sel, G, cv_sel, fc_sel);
            pf->interpolation (ctx, M, int(qdim));

            mesh_fem::ind_dof_ct
              master_dofs = mf_disp->ind_basic_dof_of_element(cv_sel);

            model_real_plain_vector MT_u(mf_disp->nb_basic_dof_of_element(cv_sel));
            gmm::mult(gmm::transposed(M), un_sel, MT_u);
            for (size_type j = 0; j < master_dofs.size(); ++j)
              (*BN)(row, master_dofs[j]) += MT_u[j];

            if (BT) {
              for (size_type n = 0; n < d; ++n) {
                gmm::mult(gmm::transposed(M), ut[n], MT_u);
                for (size_type j = 0; j < master_dofs.size(); ++j)
                  (*BT)(row * d + n, master_dofs[j]) += MT_u[j];
              }
            }
          } // BN

        }
      } // if:cnp->cn_s
    } // cnp

  } // compute_contact_matrices



  //=========================================================================
  //
  //  Basic Brick (with given BN, BT, gap) and possibly two bodies
  //
  //=========================================================================

  struct Coulomb_friction_brick : public virtual_brick {

    mutable CONTACT_B_MATRIX BN1, BT1, BN2, BT2;
    mutable CONTACT_B_MATRIX DN, DDN, DT, DDT; // For Hughes stabilization
    mutable CONTACT_B_MATRIX BBN1, BBT1, BBN2, BBT2;
    mutable model_real_plain_vector gap, threshold, friction_coeff, alpha;
    mutable model_real_plain_vector RLN, RLT;
    mutable scalar_type r, gamma;
    mutable bool is_init;
    bool Tresca_version, contact_only;
    bool really_stationary, friction_dynamic_term;
    bool two_variables, Hughes_stabilized;
    int augmentation_version; // 0 for non-symmetric Alart-Curnier version
                              // 1 for symmetric Alart-Curnier version
                              // 2 for new version (augmented multipliers)
                              // 3 for new version with De Saxcé projection

    void init_BBN_BBT(void) const {
      gmm::resize(BBN1, gmm::mat_nrows(BN1), gmm::mat_ncols(BN1));
      gmm::copy(BN1, BBN1);
      if (Hughes_stabilized) {
        gmm::resize(DDN, gmm::mat_nrows(DN), gmm::mat_ncols(DN));
        gmm::copy(DN, DDN);
      }
      if (two_variables) {
        gmm::resize(BBN2, gmm::mat_nrows(BN2), gmm::mat_ncols(BN2));
        gmm::copy(BN2, BBN2);
      }
      if (!contact_only) {
        if (Hughes_stabilized) {
          gmm::resize(DDT, gmm::mat_nrows(DT), gmm::mat_ncols(DT));
          gmm::copy(DT, DDT);
        }
        gmm::resize(BBT1, gmm::mat_nrows(BT1), gmm::mat_ncols(BT1));
        gmm::copy(BT1, BBT1);
        if (two_variables) {
          gmm::resize(BBT2, gmm::mat_nrows(BT2), gmm::mat_ncols(BT2));
          gmm::copy(BT2, BBT2);
        }
      }
      size_type nbc = gmm::mat_nrows(BN1);
      size_type d = gmm::mat_nrows(BT1)/nbc;
      for (size_type i = 0; i < nbc; ++i) {
        gmm::scale(gmm::mat_row(BBN1, i), alpha[i]);
        if (Hughes_stabilized) gmm::scale(gmm::mat_row(DDN, i), alpha[i]);
        if (two_variables)
          gmm::scale(gmm::mat_row(BBN2, i), alpha[i]);
        if (!contact_only)
          for (size_type k = 0; k < d; ++k) {
            if (Hughes_stabilized)
              gmm::scale(gmm::mat_row(DDT, d*i+k), alpha[i]);
            gmm::scale(gmm::mat_row(BBT1, d*i+k), alpha[i]);
            if (two_variables)
              gmm::scale(gmm::mat_row(BBT2, d*i+k), alpha[i]);
          }
      }
      is_init = true;
    }

    void precomp(const model_real_plain_vector &u1,
                 const model_real_plain_vector &u2,
                 const model_real_plain_vector &lambda_n,
                 const model_real_plain_vector &lambda_t,
                 const model_real_plain_vector &wt1,
                 const model_real_plain_vector &wt2) const {
      gmm::copy(gmm::scaled(gap, r), RLN);
      for (size_type i = 0; i < gmm::mat_nrows(BN1); ++i) RLN[i] *= alpha[i];
      gmm::add(lambda_n, RLN);
      gmm::mult_add(BBN1, gmm::scaled(u1, -r), RLN);
      if (Hughes_stabilized)
        gmm::mult_add(DDN, gmm::scaled(lambda_n, -r), RLN);
      if (two_variables) gmm::mult_add(BBN2, gmm::scaled(u2, -r), RLN);
      if (!contact_only) {
        gmm::copy(lambda_t, RLT);
        if (friction_dynamic_term) {
          gmm::mult_add(BBT1, gmm::scaled(wt1, -r*gamma), RLT);
          if (two_variables)
            gmm::mult_add(BBT2, gmm::scaled(wt2, -r*gamma), RLT);
        }
        if (!really_stationary) {
          gmm::mult_add(BBT1, gmm::scaled(u1, -r), RLT);
          if (two_variables) gmm::mult_add(BBT2, gmm::scaled(u2, -r), RLT);
        }
        if (Hughes_stabilized)
          gmm::mult_add(DDT, gmm::scaled(lambda_t, -r), RLT);
      }
    }

    // Common part for all contact with friction bricks
    void basic_asm_real_tangent_terms(const model_real_plain_vector &u1,
                                      const model_real_plain_vector &u2,
                                      const model_real_plain_vector &lambda_n,
                                      const model_real_plain_vector &lambda_t,
                                      const model_real_plain_vector &wt1,
                                      const model_real_plain_vector &wt2,
                                      model::real_matlist &matl,
                                      model::real_veclist &vecl,
                                      build_version version) const {
      size_type nbt = 4 + (contact_only ? 0 : 4) + (two_variables ? 3 : 0)
        + (two_variables && !contact_only ? 2 : 0);
      GMM_ASSERT1(matl.size() == nbt,
                  "Wrong number of terms for the contact brick");

      const scalar_type vt1 = scalar_type(1), vt0 = scalar_type(0);
      size_type nbc = gmm::mat_nrows(BN1);
      size_type d = gmm::mat_nrows(BT1)/nbc;

      // Matrices to be filled
      size_type nt = 0;
      model_real_sparse_matrix &T_u1_u1 = matl[nt++], &T_u2_u2 = matl[nt++];
      if (!two_variables) nt--;
      model_real_sparse_matrix &T_u1_n = matl[nt++], &T_n_u1 = matl[nt++];
      if (!two_variables) nt -= 2;
      model_real_sparse_matrix &T_u2_n = matl[nt++], &T_n_u2 = matl[nt++];
      size_type nvec_lambda_n = nt;
      model_real_sparse_matrix &T_n_n = matl[nt++];
      if (contact_only) nt -= 2;
      model_real_sparse_matrix &T_u1_t = matl[nt++], &T_t_u1 = matl[nt++];
      if (contact_only || !two_variables) nt -= 2;
      model_real_sparse_matrix &T_u2_t = matl[nt++], &T_t_u2 = matl[nt++];
      if (contact_only) nt -= 2;
      size_type nvec_lambda_t = nt;
      model_real_sparse_matrix &T_t_t = matl[nt++], &T_t_n = matl[nt++];

      // Rhs to be filled
      model_real_plain_vector &ru1 = vecl[0];
      model_real_plain_vector &ru2 = vecl[1];
      model_real_plain_vector &rlambda_n = vecl[nvec_lambda_n];
      model_real_plain_vector &rlambda_t = vecl[nvec_lambda_t];

      // pre-computations
      if (!is_init) init_BBN_BBT();
      gmm::resize(RLN, nbc);
      if (!contact_only) gmm::resize(RLT, nbc*d);
      if (augmentation_version <= 2)
        precomp(u1, u2, lambda_n, lambda_t, wt1, wt2);

      if (version & model::BUILD_MATRIX) {
        base_matrix pg(d, d);
        base_vector vg(d);

        gmm::clear(T_n_n); gmm::clear(T_n_u1);
        gmm::clear(T_u1_n); gmm::clear(T_u1_u1);
        if (two_variables)
          { gmm::clear(T_u2_u2); gmm::clear(T_n_u2); gmm::clear(T_u2_n); }
        if (!contact_only) {
          gmm::clear(T_u1_t); gmm::clear(T_t_n); gmm::clear(T_t_t);
          if (two_variables) gmm::clear(T_u2_t);
        }

        switch (augmentation_version) {
        case 1: case 2:
          gmm::copy(gmm::scaled(gmm::transposed(BN1), -vt1), T_u1_n);
          if (two_variables)
            gmm::copy(gmm::scaled(gmm::transposed(BN2), -vt1), T_u2_n);
          for (size_type i=0; i < nbc; ++i) {
            if (RLN[i] > vt0) {
              gmm::clear(gmm::mat_col(T_u1_n, i));
              if (two_variables) gmm::clear(gmm::mat_col(T_u2_n, i));
              T_n_n(i, i) = -vt1/(r*alpha[i]);
            }
            if (Hughes_stabilized && RLN[i] <= vt0)
              gmm::copy(gmm::scaled(gmm::mat_row(DN, i), -vt1),
                        gmm::mat_col(T_n_n, i));
          }
          if (Hughes_stabilized) {
            model_real_sparse_matrix aux(nbc, nbc);
            gmm::copy(gmm::transposed(T_n_n), aux);
            gmm::copy(aux, T_n_n);
          }
          gmm::copy(gmm::transposed(T_u1_n), T_n_u1);
          if (two_variables) gmm::copy(gmm::transposed(T_u2_n), T_n_u2);
          if (!contact_only) {
            for (size_type i=0; i < nbc; ++i) {
              gmm::sub_interval SUBI(i*d, d);
              scalar_type th = Tresca_version ? threshold[i]
                : - (std::min(vt0, RLN[i])) * friction_coeff[i];
              ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
              if (!really_stationary)
                for (size_type k1 = 0; k1 < d; ++k1)
                  for (size_type k2 = 0; k2 < d; ++k2) {
                    gmm::add(gmm::scaled(gmm::mat_row(BT1,i*d+k1),-pg(k2,k1)),
                             gmm::mat_col(T_u1_t, i*d+k2));
                    if (two_variables)
                      gmm::add(gmm::scaled(gmm::mat_row(BT2,i*d+k1),
                                           -pg(k2,k1)),
                               gmm::mat_col(T_u2_t, i*d+k2));
                  }

              if (!Tresca_version) {
                if (RLN[i] <= vt0) {
                  ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
                  for (size_type k = 0; k < d; ++k)
                    T_t_n(i*d+k, i) = - friction_coeff[i] * vg[k]/(r*alpha[i]);
                }
              }
              for (size_type k = 0; k < d; ++k) pg(k,k) -= vt1;
	      gmm::copy(gmm::scaled(pg, vt1/(r*alpha[i])),
			gmm::sub_matrix(T_t_t,SUBI));
              if (Hughes_stabilized) {
                for (size_type k = 0; k < d; ++k)
                  for (size_type l = 0; l < d; ++l) {
                    gmm::add(gmm::scaled(gmm::mat_row(DT, d*i+l), -pg(k,l)),
                             gmm::mat_col(T_t_t, d*i+k));
                  }
              }

            }
            if (Hughes_stabilized) {
              model_real_sparse_matrix aux(gmm::mat_nrows(T_t_t),
					   gmm::mat_nrows(T_t_t));
              gmm::copy(gmm::transposed(T_t_t), aux);
              gmm::copy(aux, T_t_t);
            }
            gmm::copy(gmm::transposed(T_u1_t), T_t_u1);
            if (two_variables) gmm::copy(gmm::transposed(T_u2_t), T_t_u2);
          }

          if (augmentation_version == 1) {
            gmm::copy(gmm::scaled(gmm::transposed(BN1), -vt1), T_u1_n);
            if (two_variables)
              gmm::copy(gmm::scaled(gmm::transposed(BN2), -vt1), T_u2_n);
            if (!contact_only) {
              gmm::copy(gmm::scaled(gmm::transposed(BT1), -vt1), T_u1_t);
              if (two_variables)
                gmm::copy(gmm::scaled(gmm::transposed(BT2), -vt1), T_u2_t);
            }
          } else {
            model_real_sparse_matrix tmp1(gmm::mat_ncols(BN1),
                                          gmm::mat_ncols(BN1));
            gmm::mult(gmm::transposed(gmm::scaled(BBN1,-r)), T_n_u1, tmp1);
            gmm::add(tmp1, T_u1_u1);
            if (two_variables) {
              gmm::mult(gmm::transposed(gmm::scaled(BBN2,-r)), T_n_u2, tmp1);
              gmm::add(tmp1, T_u2_u2);
            }

            if (!contact_only) {
              gmm::mult(gmm::transposed(gmm::scaled(BBT1,-r)), T_t_u1, tmp1);
              gmm::add(tmp1, T_u1_u1);
              if (two_variables) {
                gmm::mult(gmm::transposed(gmm::scaled(BBT2,-r)), T_t_u2, tmp1);
                gmm::add(tmp1, T_u2_u2);
              }
            }
	  }

	  if (!contact_only && !Tresca_version) {
	    // should be simplified ... !
	    model_real_sparse_matrix tmp5(gmm::mat_ncols(BT1),
					  gmm::mat_ncols(BT1));
	    model_real_sparse_matrix tmp6(gmm::mat_ncols(BT1),
					  gmm::mat_ncols(BT1));
	    model_real_sparse_matrix tmp7(gmm::mat_ncols(BT2),
					  gmm::mat_ncols(BT2));
	    model_real_sparse_matrix tmp8(gmm::mat_ncols(BT2),
					  gmm::mat_ncols(BT2));
	    model_real_sparse_matrix tmp3(gmm::mat_ncols(T_t_u1),
					  gmm::mat_nrows(T_t_u1));
	    model_real_sparse_matrix tmp4(gmm::mat_ncols(T_t_u2),
					  gmm::mat_nrows(T_t_u2));
	    
	    for (size_type i=0; i < nbc; ++i) {
	      gmm::sub_interval SUBI(i*d, d);
	      scalar_type th = - (std::min(vt0, RLN[i])) * friction_coeff[i];
	      if (RLN[i] <= vt0) {
		ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
		for (size_type k = 0; k < d; ++k) {
		  gmm::add(gmm::scaled(gmm::mat_row(BN1,i),
			 vg[k]*friction_coeff[i]), gmm::mat_col(tmp3, i*d+k));
		  if (two_variables)
		    gmm::add(gmm::scaled(gmm::mat_row(BN2,i),
			 vg[k]*friction_coeff[i]), gmm::mat_col(tmp4, i*d+k));

		  if (augmentation_version == 2) {

		    gmm::add(gmm::scaled(gmm::mat_row(BT1,i*d+k),
			   vg[k]*friction_coeff[i]), gmm::mat_col(T_u1_n, i));
		    if (two_variables)
		      gmm::add(gmm::scaled(gmm::mat_row(BT2,i*d+k),
			   vg[k]*friction_coeff[i]), gmm::mat_col(T_u2_n, i));
		  
		    gmm::copy(gmm::scaled(gmm::mat_row(BBT1,i*d+k),
					  -r*friction_coeff[i]*vg[k]),
			      gmm::mat_col(tmp5, i*d+k));
		    gmm::copy(gmm::mat_row(BN1,i),
			      gmm::mat_col(tmp6, i*d+k));
		    if (two_variables) {
		      gmm::copy(gmm::scaled(gmm::mat_row(BBT2,i*d+k),
					    -r*friction_coeff[i]*vg[k]),
				gmm::mat_col(tmp7, i*d+k));
		      gmm::copy(gmm::mat_row(BN2,i),
				gmm::mat_col(tmp8, i*d+k));
		    }
		  }
		}
	      }
	    }

	    gmm::add(gmm::transposed(tmp3), T_t_u1);
	    if (two_variables)
	      gmm::add(gmm::transposed(tmp4), T_t_u2);
	    
	    if (augmentation_version == 2) {
	      model_real_sparse_matrix tmp1(gmm::mat_ncols(BN1),
					    gmm::mat_ncols(BN1));
	      gmm::mult(tmp5, gmm::transposed(tmp6), gmm::transposed(tmp1));
	      gmm::add(gmm::transposed(tmp1), T_u1_u1);
	      if (two_variables) {
		gmm::mult(tmp7, gmm::transposed(tmp8),gmm::transposed(tmp1));
		gmm::add(gmm::transposed(tmp1), T_u2_u2);
	      }
	    }
          }
          break;

        case 3:
          gmm::copy(gmm::scaled(gmm::transposed(BN1), -vt1), T_u1_n);
          if (two_variables)
            gmm::copy(gmm::scaled(gmm::transposed(BN2), -vt1), T_u2_n);
          for (size_type i=0; i < nbc; ++i) {
            if (lambda_n[i] > vt0) {
              gmm::clear(gmm::mat_col(T_u1_n, i));
              if (two_variables) gmm::clear(gmm::mat_col(T_u2_n, i));
              T_n_n(i, i) = -vt1/r;
            }
          }
          gmm::copy(gmm::scaled(BN1, -vt1), T_n_u1);
          if (two_variables) gmm::copy(gmm::scaled(BN2, -r), T_n_u2);
          if (!contact_only) {
            for (size_type i=0; i < nbc; ++i) {
              gmm::sub_interval SUBI(i*d, d);
              scalar_type th = Tresca_version ? threshold[i]
                : gmm::neg(lambda_n[i]) * friction_coeff[i];
              ball_projection_grad(gmm::sub_vector(lambda_t, SUBI), th, pg);
              if (!really_stationary)
                for (size_type k1 = 0; k1 < d; ++k1)
                  for (size_type k2 = 0; k2 < d; ++k2) {
                    gmm::add(gmm::scaled(gmm::mat_row(BT1,i*d+k1),-pg(k2,k1)),
                             gmm::mat_col(T_u1_t, i*d+k2));
                    if (two_variables)
                      gmm::add(gmm::scaled(gmm::mat_row(BT2,i*d+k1),
                                           -pg(k2,k1)),
                               gmm::mat_col(T_u2_t, i*d+k2));
                  }
              if (!Tresca_version) {
                ball_projection_grad_r(gmm::sub_vector(lambda_t, SUBI),th,vg);
                for (size_type k1 = 0; k1 < d; ++k1) {
                  gmm::add(gmm::scaled(gmm::mat_row(BT1,i*d+k1),
                                       friction_coeff[i]*vg[k1]),
                           gmm::mat_col(T_u1_n, i));
                  if (two_variables)
                    gmm::add(gmm::scaled(gmm::mat_row(BT2,i*d+k1),
                                         friction_coeff[i]*vg[k1]),
                             gmm::mat_col(T_u2_n, i));
                  T_t_n(i*d+k1, i) = friction_coeff[i] * vg[k1] / (r*alpha[i]);
                }
              }
              for (size_type k = 0; k < d; ++k) pg(k,k) -= vt1;

              gmm::copy(gmm::scaled(pg, vt1/(r*alpha[i])),
			gmm::sub_matrix(T_t_t, SUBI));

            }
            gmm::copy(gmm::scaled(BT1, -vt1), T_t_u1);
            if (two_variables) gmm::copy(gmm::scaled(BT2, -r), T_t_u2);
          }
          break;

        case 4: // Desaxce projection
          base_small_vector x(d+1), n(d+1), u(d); n[0] = vt1;
          base_matrix g(d+1, d+1);
          model_real_sparse_matrix T_n_u1_transp(gmm::mat_ncols(T_n_u1), nbc);
          model_real_sparse_matrix T_n_u2_transp(gmm::mat_ncols(T_n_u2), nbc);

          gmm::mult(BT1, u1, RLT);
          if (two_variables) gmm::mult_add(BT2, u2, RLT);

          for (size_type i=0; i < nbc; ++i) {
            x[0] = lambda_n[i];
            for (size_type j=0; j < d; ++j) x[1+j] = lambda_t[i*d+j];
            De_Saxce_projection_grad(x, n, friction_coeff[i], g);

            gmm::add(gmm::scaled(gmm::mat_row(BN1, i), -g(0,0)),
                     gmm::mat_col(T_u1_n, i));
            if (two_variables)
              gmm::add(gmm::scaled(gmm::mat_row(BN2, i), -g(0,0)),
                       gmm::mat_col(T_u2_n, i));
            T_n_n(i, i) = (g(0,0) - vt1)/(r*alpha[i]);

            gmm::copy(gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), u);
            scalar_type nu = gmm::vect_norm2(u);
            if (nu != vt0)
              for (size_type j=0; j < d; ++j) {
                gmm::add(gmm::scaled(gmm::mat_row(BT1, i*d+j),
                                     friction_coeff[i] * u[j] / nu),
                         gmm::mat_col(T_n_u1_transp, i));
                if (two_variables)
                  gmm::add(gmm::scaled(gmm::mat_row(BT2, i*d+j),
                                       friction_coeff[i] * u[j] / nu),
                           gmm::mat_col(T_n_u2_transp, i));
              }

            for (size_type j=0; j < d; ++j) {
              gmm::add(gmm::scaled(gmm::mat_row(BT1, i*d+j), -g(0,j+1)),
                     gmm::mat_col(T_u1_n, i));
              if (two_variables)
                gmm::add(gmm::scaled(gmm::mat_row(BT2, i*d+j), -g(0,j+1)),
                         gmm::mat_col(T_u2_n, i));

              gmm::add(gmm::scaled(gmm::mat_row(BN1, i), -g(j+1,0)),
                       gmm::mat_col(T_u1_t, i*d+j));
              if (two_variables)
                gmm::add(gmm::scaled(gmm::mat_row(BN2, i), -g(j+1,0)),
                         gmm::mat_col(T_u2_t, i*d+j));

              for (size_type k=0; k < d; ++k) {
                gmm::add(gmm::scaled(gmm::mat_row(BT1, i*d+k), -g(1+j,1+k)),
                         gmm::mat_col(T_u1_t, i*d+j));
                if (two_variables)
                  gmm::add(gmm::scaled(gmm::mat_row(BT2, i*d+k), -g(1+j,1+k)),
                           gmm::mat_col(T_u2_t, i*d+j));
                T_t_t(i*d+j, i*d+k) = g(1+j, 1+k)/r;
              }
              T_t_t(i*d+j, i*d+j) -= vt1/(r*alpha[i]);
              T_t_n(i*d+j, i) = g(1+j,0)/(r*alpha[i]);
              // T_n_t(i, i*d+j) = g(0,1+j)/(r*alpha[i]);
            }
          }
          gmm::copy(gmm::scaled(BN1, -vt1), T_n_u1);
          if (two_variables) gmm::copy(gmm::scaled(BN2, -vt1), T_n_u2);
          gmm::add(gmm::transposed(T_n_u1_transp), T_n_u1);
	  if (two_variables) gmm::add(gmm::transposed(T_n_u2_transp), T_n_u2);
          gmm::copy(gmm::scaled(BT1, -vt1), T_t_u1);
          if (two_variables) gmm::copy(gmm::scaled(BT2, -vt1), T_t_u2);
          break;
        }
      }

      if (version & model::BUILD_RHS) {

        switch (augmentation_version) {
        case 1: // unsymmetric Alart-Curnier
          for (size_type i=0; i < nbc; ++i) {
            RLN[i] = std::min(scalar_type(0), RLN[i]);
            if (!contact_only) {
              scalar_type radius = Tresca_version ? threshold[i]
                : -friction_coeff[i]*RLN[i];
              ball_projection
                (gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), radius);
            }
          }
          gmm::mult_add(gmm::transposed(BN1), lambda_n, ru1);
          if (two_variables)
            gmm::mult_add(gmm::transposed(BN2), lambda_n, ru2);
          if (!contact_only) {
            gmm::mult_add(gmm::transposed(BT1), lambda_t, ru1);
            if (two_variables)
              gmm::mult_add(gmm::transposed(BT2), lambda_t, ru2);
          }
	  for (size_type i = 0; i < nbc; ++i) {
	    rlambda_n[i] = (lambda_n[i] - RLN[i]) / (r * alpha[i]);
	    if (!contact_only) 
	      for (size_type k = 0; k < d; ++k)
		rlambda_t[i*d+k]
		  = (lambda_t[i*d+k] - RLT[i*d+k]) / (r * alpha[i]);
	  }
          break;
        case 2: // symmetric Alart-Curnier
          for (size_type i=0; i < nbc; ++i) {
            RLN[i] = std::min(vt0, RLN[i]);
            if (!contact_only) {
              scalar_type radius = Tresca_version ? threshold[i]
                : -friction_coeff[i]*RLN[i];
              ball_projection
                (gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), radius);
            }
          }
          gmm::mult_add(gmm::transposed(BN1), RLN, ru1);
          if (two_variables) gmm::mult_add(gmm::transposed(BN2), RLN, ru2);
          if (!contact_only) {
            gmm::mult_add(gmm::transposed(BT1), RLT, ru1);
            if (two_variables) gmm::mult_add(gmm::transposed(BT2), RLT, ru2);
          }
	  for (size_type i = 0; i < nbc; ++i) {
	    rlambda_n[i] = (lambda_n[i] - RLN[i]) / (r * alpha[i]);
	    if (!contact_only) 
	      for (size_type k = 0; k < d; ++k)
		rlambda_t[i*d+k]
		  = (lambda_t[i*d+k] - RLT[i*d+k]) / (r * alpha[i]);
	  }
          break;
        case 3: // New unsymmetric method
          if (!contact_only) gmm::copy(lambda_t, RLT);
          for (size_type i=0; i < nbc; ++i) {
            RLN[i] = -gmm::neg(lambda_n[i]);
            rlambda_n[i] = gmm::pos(lambda_n[i])/r - alpha[i]*gap[i];

            if (!contact_only) {
              scalar_type radius = Tresca_version ? threshold[i]
                : friction_coeff[i]*gmm::neg(lambda_n[i]);
              ball_projection
                (gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), radius);
            }
          }
          gmm::mult(gmm::transposed(BN1), RLN, ru1);
          if (two_variables) gmm::mult(gmm::transposed(BN2), RLN, ru2);
          gmm::mult_add(BBN1, u1, rlambda_n);
          if (two_variables) gmm::mult_add(BBN2, u2, rlambda_n);
          if (!contact_only) {
            gmm::mult_add(gmm::transposed(BT1), RLT, ru1);
            if (two_variables) gmm::mult_add(gmm::transposed(BT2), RLT, ru2);
            gmm::add(gmm::scaled(lambda_t, vt1/r), gmm::scaled(RLT,-vt1/r),
                      rlambda_t);
            gmm::mult_add(BBT1, u1, rlambda_t);
            if (two_variables) gmm::mult_add(BBT2, u2, rlambda_t);
          }
	  for (size_type i = 0; i < nbc; ++i) {
	    rlambda_n[i] /= alpha[i];
	    if (!contact_only) 
	      for (size_type k = 0; k < d; ++k) rlambda_t[i*d+k] /= alpha[i];
	  }
          break;
        case 4:  // New unsymmetric method with De Saxce projection
          base_small_vector x(d+1), n(d+1);
          n[0] = vt1;
          GMM_ASSERT1(!Tresca_version,
               "Augmentation version incompatible with Tresca friction law");
          gmm::mult(BBT1, u1, rlambda_t);
          if (two_variables)
              gmm::mult_add(BBT2, u2, rlambda_t);
          for (size_type i=0; i < nbc; ++i) {
            x[0] = lambda_n[i];
            gmm::copy(gmm::sub_vector(lambda_t, gmm::sub_interval(i*d,d)),
                      gmm::sub_vector(x, gmm::sub_interval(1, d)));
            De_Saxce_projection(x, n, friction_coeff[i]);
            RLN[i] = x[0];
            gmm::copy(gmm::sub_vector(x, gmm::sub_interval(1, d)),
                      gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)));
            rlambda_n[i] = lambda_n[i]/r - x[0]/r - alpha[i]*gap[i]
              - friction_coeff[i] * gmm::vect_norm2(gmm::sub_vector(rlambda_t,
                                                    gmm::sub_interval(i*d,d)));
          }
          gmm::mult_add(gmm::transposed(BT1), RLT, ru1);
          if (two_variables) gmm::mult_add(gmm::transposed(BT2), RLT, ru2);
          gmm::mult_add(gmm::transposed(BN1), RLN, ru1);
          if (two_variables) gmm::mult_add(gmm::transposed(BN2), RLN, ru2);
          gmm::add(gmm::scaled(lambda_t, vt1/r), rlambda_t);
          gmm::add(gmm::scaled(RLT, -vt1/r), rlambda_t);
          gmm::mult_add(BBN1, u1, rlambda_n);
          if (two_variables) gmm::mult_add(BBN2, u2, rlambda_n);
	  for (size_type i = 0; i < nbc; ++i) {
	    rlambda_n[i] /= alpha[i];
	    if (!contact_only) 
	      for (size_type k = 0; k < d; ++k) rlambda_t[i*d+k] /= alpha[i];
	  }
          break;
        }
      }
    }

    // specific part for the basic bricks : BN, BT, gap, r, alpha are given.
    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type /* region */,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 0, "Contact brick need no mesh_im");
      size_type nbvar = 2 + (contact_only ? 0 : 1) + (two_variables ? 1 : 0);
      GMM_ASSERT1(vl.size() == nbvar,
                  "Wrong number of variables for contact brick");
      size_type nbdl = 3 + (contact_only ? 0 : 1) + (Tresca_version ? 1 : 0)
        + (friction_dynamic_term ? 2 : 0);
      GMM_ASSERT1(dl.size() == nbdl, "Wrong number of data for contact brick, "
                  << dl.size() << " should be " << nbdl);

      size_type nbc = gmm::mat_nrows(BN1);

      // Variables
      // Without friction and one displacement  : u1, lambda_n
      // Without friction and two displacements : u1, u2, lambda_n
      // With friction and one displacement     : u1, lambda_n, lambda_t
      // With friction and two displacements    : u1, u2, lambda_n, lambda_t
      size_type nv = 0;
      const model_real_plain_vector &u1 = md.real_variable(vl[nv++]);
      const model_real_plain_vector &u2 = md.real_variable(vl[nv++]);
      if (!two_variables) nv--;
      const model_real_plain_vector &lambda_n = md.real_variable(vl[nv++]);
      if (contact_only) nv--;
      const model_real_plain_vector &lambda_t = md.real_variable(vl[nv]);

      // Parameters
      // (order : r, gap, alpha, friction_coeff, gamma, wt, threshold)
      size_type np = 0, np_wt1 = 0, np_wt2 = 0, np_alpha = 0;
      const model_real_plain_vector &vr = md.real_variable(dl[np++]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      r = vr[0];
      const model_real_plain_vector &vgap = md.real_variable(dl[np++]);
      GMM_ASSERT1(gmm::vect_size(vgap) == 1 || gmm::vect_size(vgap) == nbc,
                  "Parameter gap has a wrong size");
      gmm::resize(gap, nbc);
      if (gmm::vect_size(vgap) == 1)
        gmm::fill(gap, vgap[0]);
      else
        gmm::copy(vgap, gap);
      np_alpha = np++;
      const model_real_plain_vector &valpha = md.real_variable(dl[np_alpha]);
      GMM_ASSERT1(gmm::vect_size(valpha)== 1 || gmm::vect_size(valpha) == nbc,
                  "Parameter alpha has a wrong size");
      gmm::resize(alpha, nbc);
      if (gmm::vect_size(valpha) == 1)
        gmm::fill(alpha, valpha[0]);
      else
        gmm::copy(valpha, alpha);
      if (!contact_only) {
        const model_real_plain_vector &vfr = md.real_variable(dl[np++]);
        GMM_ASSERT1(gmm::vect_size(vfr)==1 || gmm::vect_size(vfr) == nbc,
                    "Parameter friction_coeff has a wrong size");
        gmm::resize(friction_coeff, nbc);
        if (gmm::vect_size(vfr) == 1)
          gmm::fill(friction_coeff, vfr[0]);
        else
          gmm::copy(vfr, friction_coeff);
        if (friction_dynamic_term) {
          const model_real_plain_vector &vg = md.real_variable(dl[np++]);
          GMM_ASSERT1(gmm::vect_size(vg) == 1,
                      "Parameter gamma should be a scalar");
          gamma = vg[0];
          np_wt1 = np++;
          if (two_variables) np_wt2 = np++;
        }
        if (Tresca_version) {
          const model_real_plain_vector &vth = md.real_variable(dl[np++]);
          GMM_ASSERT1(gmm::vect_size(vth) == 1 || gmm::vect_size(vth) == nbc,
                      "Parameter threshold has a wrong size");
          gmm::resize(threshold, nbc);
          if (gmm::vect_size(vth) == 1)
            gmm::fill(threshold, vth[0]);
          else
            gmm::copy(vth, threshold);
        }
      }

      if (md.is_var_newer_than_brick(dl[np_alpha], ib)) is_init = false;

      basic_asm_real_tangent_terms
        (u1, u2, lambda_n, lambda_t, md.real_variable(dl[np_wt1]),
         md.real_variable(dl[np_wt2]), matl, vecl, version);

    }

    Coulomb_friction_brick(int aug_version, bool contact_only_,
                           bool two_variables_=false,
                           bool Tresca_version_=false,
                           bool Hughes_stabilized_=false,
                           bool friction_dynamic_term_=false) {

#if GETFEM_PARA_LEVEL > 1
    if (!getfem::MPI_IS_MASTER()) GMM_WARNING1("Nodal contact bricks don't support GETFEM_PARA_LEVEL > 1 yet!!!");
#endif

      if (aug_version == 4 && contact_only_) aug_version = 3;
      augmentation_version = aug_version;
      GMM_ASSERT1(aug_version >= 1 && aug_version <= 4,
                  "Wrong augmentation version");
      GMM_ASSERT1(!Hughes_stabilized_ || aug_version <= 2,
                  "The Hughes stabilized version is only for Alart-Curnier "
                  "version");
      contact_only = contact_only_;
      is_init = false;
      Tresca_version = Tresca_version_;
      really_stationary = false;   // for future version ...
      friction_dynamic_term = friction_dynamic_term_;
      two_variables = two_variables_;
      Hughes_stabilized = Hughes_stabilized_;
      set_flags("Coulomb friction brick", false /* is linear*/,
                /* is symmetric */
                (augmentation_version == 2) && (contact_only||Tresca_version),
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

    void set_BN1(CONTACT_B_MATRIX &BN1_) {
      gmm::resize(BN1, gmm::mat_nrows(BN1_), gmm::mat_ncols(BN1_));
      gmm::copy(BN1_, BN1);
      is_init = false;
    }

    void set_DN(CONTACT_B_MATRIX &DN_) {
      gmm::resize(DN, gmm::mat_nrows(DN_), gmm::mat_ncols(DN_));
      gmm::copy(DN_, DN);
      is_init = false;
    }

    void set_DT(CONTACT_B_MATRIX &DT_) {
      gmm::resize(DT, gmm::mat_nrows(DT_), gmm::mat_ncols(DT_));
      gmm::copy(DT_, DT);
      is_init = false;
    }

    void set_BT1(CONTACT_B_MATRIX &BT1_) {
      gmm::resize(BT1, gmm::mat_nrows(BT1_), gmm::mat_ncols(BT1_));
      gmm::copy(BT1_, BT1);
      is_init = false;
    }

    CONTACT_B_MATRIX &get_BN1(void) { return BN1; }
    CONTACT_B_MATRIX &get_DN(void) { return DN; }
    CONTACT_B_MATRIX &get_DT(void) { return DT; }
    CONTACT_B_MATRIX &get_BT1(void) { return BT1; }
    const CONTACT_B_MATRIX &get_BN1(void) const { return BN1; }
    const CONTACT_B_MATRIX &get_DN(void) const { return DN; }
    const CONTACT_B_MATRIX &get_BT1(void) const { return BT1; }

  };


  CONTACT_B_MATRIX &contact_brick_set_BN
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    Coulomb_friction_brick *p = dynamic_cast<Coulomb_friction_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->get_BN1();
  }


  CONTACT_B_MATRIX &contact_brick_set_DN
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    Coulomb_friction_brick *p = dynamic_cast<Coulomb_friction_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->get_DN();
  }

  CONTACT_B_MATRIX &contact_brick_set_DT
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    Coulomb_friction_brick *p = dynamic_cast<Coulomb_friction_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->get_DT();
  }


  CONTACT_B_MATRIX &contact_brick_set_BT
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    Coulomb_friction_brick *p = dynamic_cast<Coulomb_friction_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->get_BT1();
  }

  //=========================================================================
  //  Add a frictionless contact condition with BN, r, alpha given.
  //=========================================================================

  size_type add_basic_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &dataname_r, CONTACT_B_MATRIX &BN,
   std::string dataname_gap, std::string dataname_alpha,
   int aug_version, bool Hughes_stabilized) {

    Coulomb_friction_brick *pbr_
      = new Coulomb_friction_brick(aug_version, true, false, false, Hughes_stabilized);
    pbr_->set_BN1(BN);
    pbrick pbr = pbr_;

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    model::varnamelist dl(1, dataname_r);

    if (dataname_gap.size() == 0) {
      dataname_gap = md.new_name("contact_gap_on_" + varname_u);
      md.add_initialized_fixed_size_data
        (dataname_gap, model_real_plain_vector(1, scalar_type(0)));
    }
    dl.push_back(dataname_gap);

    if (dataname_alpha.size() == 0) {
      dataname_alpha = md.new_name("contact_parameter_alpha_on_"+ multname_n);
      md.add_initialized_fixed_size_data
        (dataname_alpha, model_real_plain_vector(1, scalar_type(1)));
    }
    dl.push_back(dataname_alpha);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }


  //=========================================================================
  //  Add a contact with friction condition with BN, r, alpha given.
  //=========================================================================

  size_type add_basic_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT,
   std::string dataname_friction_coeff,
   std::string dataname_gap, std::string dataname_alpha,
   int aug_version, bool Tresca_version, const std::string dataname_threshold,
   std::string dataname_gamma, std::string dataname_wt, bool Hughes_stabilized) {

    bool dynamic_terms = (dataname_gamma.size() > 0);

    Coulomb_friction_brick *pbr_
      = new Coulomb_friction_brick(aug_version,false, false,
                            Tresca_version, Hughes_stabilized, dynamic_terms);
    pbr_->set_BN1(BN);
    pbr_->set_BT1(BT);
    pbrick pbr = pbr_;

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    tl.push_back(model::term_description(varname_u, multname_t, false));
    tl.push_back(model::term_description(multname_t, varname_u, false));
    tl.push_back(model::term_description(multname_t, multname_t, false));
    tl.push_back(model::term_description(multname_t, multname_n,
                                         (aug_version == 4)));
    model::varnamelist dl(1, dataname_r);
    if (dataname_gap.size() == 0) {
      dataname_gap = md.new_name("contact_gap_on_" + varname_u);
      md.add_initialized_fixed_size_data
        (dataname_gap, model_real_plain_vector(1, scalar_type(0)));
    }
    dl.push_back(dataname_gap);

    if (dataname_alpha.size() == 0) {
      dataname_alpha = md.new_name("contact_parameter_alpha_on_"+ multname_n);
      md.add_initialized_fixed_size_data
        (dataname_alpha, model_real_plain_vector(1, scalar_type(1)));
    }
    dl.push_back(dataname_alpha);
    dl.push_back(dataname_friction_coeff);
    if (dataname_gamma.size()) {
      dl.push_back(dataname_gamma);
      dl.push_back(dataname_wt);
    }
    if (Tresca_version)
      dl.push_back(dataname_threshold);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);
    vl.push_back(multname_t);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }


  //=========================================================================
  //
  //  Brick with a given rigid obstacle (one body, build BN, BT, gap, alpha)
  //
  //=========================================================================
  // TODO : add an option for a weak contact condition

  struct Coulomb_friction_brick_rigid_obstacle
    : public Coulomb_friction_brick {

    std::string obstacle; // obstacle given with a signed distance expression.

  public :

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1, "This contact brick needs one mesh_im");
      size_type nbvar = 2 + (contact_only ? 0 : 1);
      GMM_ASSERT1(vl.size() == nbvar,
                  "Wrong number of variables for contact brick: "
                  << vl.size() << " should be " << nbvar);
      size_type nbdl = 1 + (contact_only ? 0 : 1) + (Tresca_version ? 1 : 0)
        + (friction_dynamic_term ? 1 : 0);
      GMM_ASSERT1(dl.size() == nbdl,
                  "Wrong number of data for contact brick: "
                  << dl.size() << " should be " << nbdl);
      GMM_ASSERT1(!two_variables, "internal error");
      const mesh_im &mim = *mims[0];

      // Variables
      // Without friction and one displacement  : u1, lambda_n
      // With friction and one displacement     : u1, lambda_n, lambda_t
      size_type nv = 0;
      const model_real_plain_vector &u1 = md.real_variable(vl[nv++]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector &lambda_n = md.real_variable(vl[nv++]);
      if (contact_only) nv--;
      const model_real_plain_vector &lambda_t = md.real_variable(vl[nv]);


      // Parameters (order : r, friction_coeff, gamma, wt, threshold)
      size_type np = 0, np_wt1 = 0, nbc;
      const model_real_plain_vector &vr = md.real_variable(dl[np++]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      r = vr[0];

      // Computation of BN, BT, gap and alpha
      if (md.is_var_mf_newer_than_brick(vl[0], ib)) {

        // Verification that mf_u1 is a pure Lagrange fem.
        GMM_ASSERT1(!(mf_u1.is_reduced()),
                    "This contact brick works only for pure Lagrange fems");
        dal::bit_vector dofs = mf_u1.basic_dof_on_region(region);
        for (dal::bv_visitor id(dofs); !id.finished(); ++id) {
          size_type cv = mf_u1.first_convex_of_basic_dof(id);
          GMM_ASSERT1(mf_u1.fem_of_element(cv)->is_lagrange(),
                      "This contact brick works only for pure Lagrange fems");
        }
        size_type d = mf_u1.get_qdim() - 1, i = 0, j = 0;
        nbc = dofs.card() / (d+1);

        // computation of alpha vector.
        base_node Pmin, Pmax;
        mf_u1.linked_mesh().bounding_box(Pmin, Pmax);
        scalar_type l = scalar_type(0);
        for (i = 0; i < Pmin.size(); ++i)
          l = std::max(l, gmm::abs(Pmax[i] - Pmin[i]));

        CONTACT_B_MATRIX MM(mf_u1.nb_dof(), mf_u1.nb_dof());
        asm_mass_matrix(MM, mim, mf_u1, region);
        gmm::resize(alpha, nbc);
        i = 0; j = 0;
        for (dal::bv_visitor id(dofs); !id.finished(); ++id, ++i)
          if ((i % (d+1)) == 0) alpha[j++] = MM(id, id) / l;


#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H


        mu::Parser parser;
        parser.SetExpr(obstacle);

        gmm::resize(gap, nbc);
        gmm::resize(BN1, nbc, mf_u1.nb_dof());
        gmm::clear(BN1);
        if (!contact_only) {
          gmm::resize(BT1, d*nbc, mf_u1.nb_dof());
          gmm::clear(BT1);
        }
        base_node pt(d+1), grad(d+1), ut[3];

        static std::string varn[4] = {"x", "y", "z", "w"};
        for (size_type k = 0; k <= d; ++k)
          parser.DefineVar(varn[k], &pt[k]);

        i = 0; j = 0;
        for (dal::bv_visitor id(dofs); !id.finished(); ++id, ++i) {
          if ((i % (d+1)) == 0) {
            gmm::copy(mf_u1.point_of_basic_dof(id), pt);
            try {

              // Computation of gap
              gap[j] = scalar_type(parser.Eval());

              // computation of BN
              size_type cv = mf_u1.first_convex_of_basic_dof(id);
              scalar_type eps
                = mf_u1.linked_mesh().convex_radius_estimate(cv) * 1E-3;
              for (size_type k = 0; k <= d; ++k) {
                pt[k] += eps;
                grad[k] = (scalar_type(parser.Eval()) - gap[j]) / eps;
                pt[k] -= eps;
              }
              // unit normal vector
              base_node un = - grad / gmm::vect_norm2(grad);

              for (size_type k = 0; k <= d; ++k)
                BN1(j, id + k) = un[k];

              // computation of BT
              if (!contact_only) {

                orthonormal_basis_to_unit_vec(d, un, ut);

                for (size_type k = 0; k <= d; ++k)
                  for (size_type nn = 0; nn < d; ++nn)
                    BT1(j*d+nn, id + k) = ut[nn][k];
              }

            } catch (mu::Parser::exception_type &e) {
              std::cerr << "Message  : " << e.GetMsg()   << std::endl;
              std::cerr << "Formula  : " << e.GetExpr()  << std::endl;
              std::cerr << "Token    : " << e.GetToken() << std::endl;
              std::cerr << "Position : " << e.GetPos()   << std::endl;
              std::cerr << "Errc     : " << e.GetCode()  << std::endl;
              GMM_ASSERT1(false, "Error in signed distance expression");
            }
            ++j;
          }

        }

        GMM_ASSERT1(gmm::vect_size(md.real_variable(vl[1])) == nbc,
                    "Wrong size of multiplier for the contact condition");

        if (!contact_only)
          GMM_ASSERT1(gmm::vect_size(md.real_variable(vl[2])) == nbc*d,
                      "Wrong size of multiplier for the friction condition");

#else

        GMM_ASSERT1(false, "Muparser is not installed, "
                    "You cannot use this contact brick");

#endif

        is_init = false;
      }
      else
        nbc = gmm::mat_nrows(BN1);

      if (!contact_only) {
        const model_real_plain_vector &vfr = md.real_variable(dl[np++]);
        GMM_ASSERT1(gmm::vect_size(vfr)==1 || gmm::vect_size(vfr) == nbc,
                    "Parameter friction_coeff has a wrong size");
        gmm::resize(friction_coeff, nbc);
        if (gmm::vect_size(vfr) == 1)
          gmm::fill(friction_coeff, vfr[0]);
        else
          gmm::copy(vfr, friction_coeff);
        if (friction_dynamic_term) {
          const model_real_plain_vector &vg = md.real_variable(dl[np++]);
          GMM_ASSERT1(gmm::vect_size(vg) == 1,
                      "Parameter gamma should be a scalar");
          gamma = vg[0];
          np_wt1 = np++;
        }
        if (Tresca_version) {
          const model_real_plain_vector &vth = md.real_variable(dl[np++]);
          GMM_ASSERT1(gmm::vect_size(vth) == 1 || gmm::vect_size(vth) == nbc,
                      "Parameter threshold has a wrong size");
          gmm::resize(threshold, nbc);
          if (gmm::vect_size(vth) == 1)
            gmm::fill(threshold, vth[0]);
          else
            gmm::copy(vth, threshold);
        }
      }

      basic_asm_real_tangent_terms
        (u1, u1, lambda_n, lambda_t, md.real_variable(dl[np_wt1]),
         md.real_variable(dl[np_wt1]), matl, vecl, version);

    }

    Coulomb_friction_brick_rigid_obstacle
    (int aug_version, bool contact_only_, const std::string &obs)
      : Coulomb_friction_brick(aug_version, contact_only_), obstacle(obs) {}

  };


  //=========================================================================
  //  Add a frictionless contact condition with a rigid obstacle given
  //  by a signed distance.
  //=========================================================================

  size_type add_nodal_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_r,
   size_type region, const std::string &obstacle, int aug_version) {
    pbrick pbr
      = new Coulomb_friction_brick_rigid_obstacle(aug_version, true, obstacle);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    model::varnamelist dl(1, dataname_r);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //  Add a contact with friction condition with a rigid obstacle given
  //  by a signed distance.
  //=========================================================================

  size_type add_nodal_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, const std::string &obstacle, int aug_version) {
    pbrick pbr
      = new Coulomb_friction_brick_rigid_obstacle(aug_version,false,obstacle);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    tl.push_back(model::term_description(varname_u, multname_t, false));
    tl.push_back(model::term_description(multname_t, varname_u, false));
    tl.push_back(model::term_description(multname_t, multname_t, false));
    tl.push_back(model::term_description(multname_t, multname_n,
                                         (aug_version == 4)));
    model::varnamelist dl(1, dataname_r);
    dl.push_back(dataname_friction_coeff);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);
    vl.push_back(multname_t);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //
  //  Brick with elastic bodies (one or two bodies, build BN, BT, Gap, alpha)
  //
  //=========================================================================
  // To be done:
  // - Large deformations: what happens when cnpl and nbc change during
  //   the iterative solution?

  struct Coulomb_friction_brick_nonmatching_meshes
    : public Coulomb_friction_brick {

    std::vector<size_type> rg1, rg2; // ids of mesh regions expected to come in
                                     // contact. For one displacement they refer
                                     // both to u1. For two displacements they
                                     // respectively refer to u1, u2.
    bool slave1, slave2; // if true, then rg1 or respectively rg2 are treated
                         // as slave surfaces (the contact multipliers are
                         // defined on these surfaces)

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {

      GMM_ASSERT1(mims.size() == 2, "This contact brick needs two mesh_im");
      const mesh_im &mim1 = *mims[0];
      const mesh_im &mim2 = *mims[1];

      // Variables
      // Without friction and one displacement  : u1, lambda_n
      // With friction and one displacement     : u1, lambda_n, lambda_t
      // Without friction and two displacements : u1, u2, lambda_n
      // With friction and two displacements    : u1, u2, lambda_n, lambda_t
      size_type nv = 0;
      std::string varname_u1 = vl[nv];
      const model_real_plain_vector &u1 = md.real_variable(varname_u1);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(varname_u1);
      if (two_variables) nv++;
      std::string varname_u2 = vl[nv++];
      const model_real_plain_vector &u2 = md.real_variable(varname_u2);
      const mesh_fem &mf_u2 = md.mesh_fem_of_variable(varname_u2);
      const model_real_plain_vector &lambda_n = md.real_variable(vl[nv]);
      if (!contact_only) nv++;
      const model_real_plain_vector &lambda_t = md.real_variable(vl[nv]);

      size_type nbc = lambda_n.size();

      // Parameters (order: r, friction_coeff)
      const model_real_plain_vector &vr = md.real_variable(dl[0]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      r = vr[0];
      if (!contact_only) {
        const model_real_plain_vector &vfr = md.real_variable(dl[1]);
        GMM_ASSERT1(gmm::vect_size(vfr)==1 || gmm::vect_size(vfr) == nbc,
                    "Parameter friction_coeff has a wrong size");
        gmm::resize(friction_coeff, nbc);
        if (gmm::vect_size(vfr) == 1)
          gmm::fill(friction_coeff, vfr[0]);
        else
          gmm::copy(vfr, friction_coeff);
      }

      // Computation of BN, BT, gap and alpha
      if (  md.is_var_mf_newer_than_brick(varname_u1, ib)
         || md.is_var_mf_newer_than_brick(varname_u2, ib)) {

        for (size_type i = 0; i <= size_type(two_variables ? 1 : 0); i++) {
          const mesh_fem &mf_u = i ? mf_u2 : mf_u1;
          // Verification that mf_u is a pure Lagrange fem.
          GMM_ASSERT1(!(mf_u.is_reduced()),
                "This contact brick works only for pure Lagrange fems");
          dal::bit_vector dofs = mf_u.basic_dof_on_region(region);
          for (dal::bv_visitor id(dofs); !id.finished(); ++id) {
            size_type cv = mf_u.first_convex_of_basic_dof(id);
            GMM_ASSERT1(mf_u.fem_of_element(cv)->is_lagrange(),
                "This contact brick works only for pure Lagrange fems");
          }
        }

        contact_node_pair_list cnpl;
        for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it)
          cnpl.append_min_dist_cn_pairs
                 (mf_u1, mf_u2, rg1[it], rg2[it], slave1, slave2);

        // Computation of gap, BN and BT
        gmm::resize(gap, nbc);
        gmm::resize(BN1, nbc, mf_u1.nb_dof());
        if (contact_only) {
          if (!two_variables) {
            compute_contact_matrices(mf_u1, mf_u2, cnpl, gap, &BN1);
          } else {
            gmm::resize(BN2, nbc, mf_u2.nb_dof());
            compute_contact_matrices(mf_u1, mf_u2, cnpl, gap, &BN1, &BN2);
          }
        } else {
          size_type d = mf_u1.get_qdim() - 1;
          gmm::resize(BT1, nbc * d, mf_u1.nb_dof());
          if (!two_variables) {
            compute_contact_matrices(mf_u1, mf_u2, cnpl, gap, &BN1,    0, &BT1);
          } else {
            // d == mf_u2.get_qdim() - 1;
            gmm::resize(BN2, nbc, mf_u2.nb_dof());
            gmm::resize(BT2, nbc * d, mf_u2.nb_dof());
            compute_contact_matrices(mf_u1, mf_u2, cnpl, gap, &BN1, &BN2, &BT1, &BT2);
          }
        }

        // computation of alpha vector.
        scalar_type l = scalar_type(0);
        for (size_type i = 0; i <= size_type(two_variables ? 1 : 0); i++) {
          const mesh_fem &mf_u = i ? mf_u2 : mf_u1;
          base_node Pmin, Pmax;
          mf_u.linked_mesh().bounding_box(Pmin, Pmax);
          for (size_type j = 0; j < Pmin.size(); ++j)
            l = std::max(l, gmm::abs(Pmax[j] - Pmin[j]));
        }
        gmm::resize(alpha, nbc);
        size_type mult_id = 0;
        for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it) {
          for (size_type swap = 0; swap <= 1; ++swap) {
            if (swap ? slave2 : slave1) {
              size_type rg = swap ? rg2[it] : rg1[it];
              const mesh_fem &mf_u = swap ? mf_u2 : mf_u1;
              const mesh_im &mim = swap ? mim2 : mim1;
              CONTACT_B_MATRIX MM(mf_u.nb_dof(), mf_u.nb_dof());
              asm_mass_matrix(MM, mim, mf_u, rg);
              size_type qdim = mf_u.get_qdim();
              dal::bit_vector rg_dofs = mf_u.basic_dof_on_region(rg);
              for (dal::bv_visitor id(rg_dofs); !id.finished(); ++id)
                if (id % qdim == 0) alpha[mult_id++] = MM(id, id) / l;
            }
          }
        }
      }

      const model_real_plain_vector dummy_wt;
      basic_asm_real_tangent_terms
        (u1, u2, lambda_n, lambda_t, dummy_wt, dummy_wt, matl, vecl, version);
    }

    Coulomb_friction_brick_nonmatching_meshes
      (int aug_version, bool contact_only_, bool two_variables_,
       const std::vector<size_type> &rg1_, const std::vector<size_type> &rg2_,
       bool slave1_=true, bool slave2_=false)
      : Coulomb_friction_brick(aug_version, contact_only_, two_variables_),
        rg1(rg1_), rg2(rg2_), slave1(slave1_), slave2(slave2_) {}

  };


  //=========================================================================
  //  Add a frictionless contact condition between two faces of one or two
  //  elastic bodies.
  //=========================================================================

  size_type add_nodal_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1, bool slave2, int aug_version) {

    bool two_variables = (varname_u1.compare(varname_u2) != 0);

    pbrick pbr = new Coulomb_friction_brick_nonmatching_meshes
           (aug_version, true, two_variables, rg1, rg2, slave1, slave2);

    // Calculate multipliers size
    const mesh_fem &mf_u1 = md.mesh_fem_of_variable(varname_u1);
    const mesh_fem &mf_u2 = md.mesh_fem_of_variable(varname_u2);
    size_type nbc = 0;
    for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it) {
      for (size_type swap = 0; swap <= 1; ++swap) {
        if (swap ? slave2 : slave1) {
          const mesh_fem &mf = swap ? mf_u2 : mf_u1;
          size_type rg = swap ? rg2[it] : rg1[it];
          dal::bit_vector rg_dofs = mf.basic_dof_on_region(rg);
          nbc += rg_dofs.card() / mf.get_qdim();
        }
      }
    }

    if (multname_n.size() == 0)
      multname_n = md.new_name("contact_multiplier");
    else
      GMM_ASSERT1(multname_n.compare(md.new_name(multname_n)) == 0,
                  "The given name for the multiplier is alraedy reserved in the model");
    md.add_fixed_size_variable(multname_n, nbc);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, false));
    if (two_variables) {
      tl.push_back(model::term_description(varname_u2, varname_u2, false));
    }
    tl.push_back(model::term_description(varname_u1, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u1, false));
    if (two_variables) {
      tl.push_back(model::term_description(varname_u2, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u2, false));
    }
    tl.push_back(model::term_description(multname_n, multname_n, false));

    // Variables (order: varname_u, multname_n)
    model::varnamelist vl;
    vl.push_back(varname_u1);
    if (two_variables) vl.push_back(varname_u2);
    vl.push_back(multname_n);

    // Parameters (order: r, ...)
    model::varnamelist dl;
    dl.push_back(dataname_r);

    model::mimlist ml;
    ml.push_back(&mim1);
    ml.push_back(&mim2);

    return md.add_brick(pbr, vl, dl, tl, ml, size_type(-1));
  }


  //=========================================================================
  //  Add a contact with friction condition between two faces of one or two
  //  elastic bodies.
  //=========================================================================

  size_type add_nodal_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1, bool slave2, int aug_version) {

    bool two_variables = (varname_u1.compare(varname_u2) != 0);

    pbrick pbr = new Coulomb_friction_brick_nonmatching_meshes
           (aug_version, false, two_variables, rg1, rg2, slave1, slave2);

    // Calculate multipliers size
    const mesh_fem &mf_u1 = md.mesh_fem_of_variable(varname_u1);
    const mesh_fem &mf_u2 = md.mesh_fem_of_variable(varname_u2);
    size_type nbc = 0;
    for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it) {
      for (size_type swap = 0; swap <= 1; ++swap) {
        if (swap ? slave2 : slave1) {
          const mesh_fem &mf = swap ? mf_u2 : mf_u1;
          size_type rg = swap ? rg2[it] : rg1[it];
          dal::bit_vector rg_dofs = mf.basic_dof_on_region(rg);
          nbc += rg_dofs.card() / mf.get_qdim();
        }
      }
    }

    if (multname_n.size() == 0)
      multname_n = md.new_name("contact_normal_multiplier");
    else
      GMM_ASSERT1(multname_n.compare(md.new_name(multname_n)) == 0,
                  "The given name for the multiplier is alraedy reserved in the model");
    md.add_fixed_size_variable(multname_n, nbc);
    if (multname_t.size() == 0)
      multname_t = md.new_name("contact_tangent_multiplier");
    else
      GMM_ASSERT1(multname_t.compare(md.new_name(multname_t)) == 0,
                  "The given name for the multiplier is alraedy reserved in the model");
    md.add_fixed_size_variable(multname_t, nbc * (mf_u1.get_qdim() - 1) ); // ??

    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, false));
    if (two_variables) {
      tl.push_back(model::term_description(varname_u2, varname_u2, false));
    }

    tl.push_back(model::term_description(varname_u1, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u1, false));
    if (two_variables) {
      tl.push_back(model::term_description(varname_u2, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u2, false));
    }
    tl.push_back(model::term_description(multname_n, multname_n, false));

    tl.push_back(model::term_description(varname_u1, multname_t, false));
    tl.push_back(model::term_description(multname_t, varname_u1, false));
    if (two_variables) {
      tl.push_back(model::term_description(varname_u2, multname_t, false));
      tl.push_back(model::term_description(multname_t, varname_u2, false));
    }
    tl.push_back(model::term_description(multname_t, multname_t, false));
    tl.push_back(model::term_description(multname_t, multname_n,
                                         (aug_version == 4)));

    // Variables (order: varname_u, multname_n, multname_t)
    model::varnamelist vl;
    vl.push_back(varname_u1);
    if (two_variables) vl.push_back(varname_u2);
    vl.push_back(multname_n);
    vl.push_back(multname_t);

    // Parameters (order: r, friction_coeff)
    model::varnamelist dl;
    dl.push_back(dataname_r);
    dl.push_back(dataname_friction_coeff);

    model::mimlist ml;
    ml.push_back(&mim1);
    ml.push_back(&mim2);

    return md.add_brick(pbr, vl, dl, tl, ml, size_type(-1));
  }

}  /* end of namespace getfem.                                             */
