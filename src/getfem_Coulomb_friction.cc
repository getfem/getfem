// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2011 Yves Renard, Konstantinos Poulios.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================


#include "getfem/getfem_Coulomb_friction.h"

#include <getfem/getfem_arch_config.h>
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif

#ifdef GETFEM_HAVE_QHULL_QHULL_H
#include <getfem/getfem_mesher.h>
#else
#include <getfem/bgeot_kdtree.h>
#endif

#ifdef _MSC_VER
#define xor ^
#endif

namespace getfem {

  typedef bgeot::convex<base_node>::dref_convex_pt_ct dref_convex_pt_ct;
  typedef bgeot::basic_mesh::ref_mesh_face_pt_ct ref_mesh_face_pt_ct;

  //=========================================================================
  //
  //  Projection on a ball and gradient of the projection.
  //
  //=========================================================================

  template<typename VEC> static void ball_projection(const VEC &x,
                                                     scalar_type radius) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
    else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a);
  }

  template<typename VEC, typename VECR>
  static void ball_projection_grad_r(const VEC &x, scalar_type radius,
                                     VECR &g) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius > 0 && a >= radius) {
      gmm::copy(x, g); gmm::scale(g, scalar_type(1)/a);
    }
    else gmm::clear(g);
  }

  template <typename VEC, typename MAT>
  static void ball_projection_grad(const VEC &x, double radius, MAT &g) {
    if (radius <= scalar_type(0)) { gmm::clear(g); return; }
    gmm::copy(gmm::identity_matrix(), g);
    scalar_type a = gmm::vect_norm2(x);
    if (a >= radius) {
      gmm::scale(g, radius/a);
      // gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
      for (size_type i = 0; i < x.size(); ++i)
        for (size_type j = 0; j < x.size(); ++j)
          g(i,j) -= radius*x[i]*x[j] / (a*a*a);
    }
  }

  template<typename VEC>
  static void De_Saxce_projection(const VEC &x, const VEC &n, scalar_type f) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    if (xn >= scalar_type(0) && f * nxt <= xn) {
      gmm::clear(const_cast<VEC&>(x));
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      gmm::add(gmm::scaled(n, -xn), const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), -f / nxt);
      gmm::add(n, const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), (xn - f * nxt) / (f*f+scalar_type(1)));
    }
  }


  template<typename VEC, typename MAT>
  static void De_Saxce_projection_grad(const VEC &x, const VEC &n,
				       scalar_type f, MAT &g) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));

    if (xn >= scalar_type(0) && f * nxt <= xn) {
      gmm::clear(g);
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      static VEC xt;
      gmm::resize(xt, gmm::vect_size(x));
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      gmm::copy(gmm::identity_matrix(), g);
      gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
      gmm::rank_one_update(g, gmm::scaled(xt, -scalar_type(1)), xt);
      gmm::scale(xt, f*(f - xn/nxt));

      
      gmm::scale(xt, -f); gmm::add(n, xt);
      gmm::rank_one_update(g, xt, xt);
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    } else {
      gmm::copy(gmm::identity_matrix(), g);
    }
  }


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
          if (ipt.i > 0 && dist2 < (*this)[ii1].dist2) {
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
          if (ipt.i > 0 && dist2 < (*this)[ii2].dist2) {
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
            if (v1_on_surface1 xor v2_on_surface1) {
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
    bool Tresca_version, symmetrized, contact_only;
    bool really_stationary, friction_dynamic_term;
    bool two_variables, Hughes_stabilized;

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
      gmm::resize(RLN, gmm::mat_nrows(BN1));
      if (!contact_only) gmm::resize(RLT, gmm::mat_nrows(BT1));

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

      const scalar_type vt1 = scalar_type(1);
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
      precomp(u1, u2, lambda_n, lambda_t, wt1, wt2);

      if (version & model::BUILD_MATRIX) {
        // Unilateral contact
        gmm::clear(T_n_n); gmm::clear(T_u1_u1);
        if (two_variables) gmm::clear(T_u2_u2);
        gmm::copy(gmm::scaled(gmm::transposed(BBN1), -vt1), T_u1_n);
        if (two_variables)
          gmm::copy(gmm::scaled(gmm::transposed(BBN2), -vt1), T_u2_n);
        for (size_type i=0; i < nbc; ++i) {
          if (RLN[i] >= scalar_type(0)) {
            gmm::clear(gmm::mat_col(T_u1_n, i));
            if (two_variables) gmm::clear(gmm::mat_col(T_u2_n, i));
            T_n_n(i, i) = -vt1/r;
          }
          else {
            if (Hughes_stabilized)
              gmm::copy(gmm::scaled(gmm::mat_row(DDN, i), -vt1), gmm::mat_col(T_n_n, i));
          }
        }
        if (Hughes_stabilized) {
          model_real_sparse_matrix aux(gmm::mat_nrows(T_n_n), gmm::mat_nrows(T_n_n));
          gmm::copy(gmm::transposed(T_n_n), aux);
          gmm::copy(aux, T_n_n);
        }
        gmm::copy(gmm::transposed(T_u1_n), T_n_u1);
        if (two_variables) gmm::copy(gmm::transposed(T_u2_n), T_n_u2);

        // Friction
        if (!contact_only) {
          base_matrix pg(d, d);
          base_vector vg(d);
          gmm::clear(T_u1_t); gmm::clear(T_t_n); gmm::clear(T_t_t);
          if (two_variables) gmm::clear(T_u2_t);

          for (size_type i=0; i < nbc; ++i) {
            gmm::sub_interval SUBI(i*d, d);
            scalar_type th = Tresca_version ? threshold[i]
              : - lambda_n[i] * friction_coeff[i];
            ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
            if (!really_stationary)
              for (size_type k1 = 0; k1 < d; ++k1)
                for (size_type k2 = 0; k2 < d; ++k2) {
                  gmm::add(gmm::scaled(gmm::mat_row(BBT1, i*d+k1), -pg(k2,k1)),
                           gmm::mat_col(T_u1_t, i*d+k2));
                  if (two_variables)
                    gmm::add(gmm::scaled(gmm::mat_row(BBT2, i*d+k1), -pg(k2,k1)),
                             gmm::mat_col(T_u2_t, i*d+k2));
                }

            if (!Tresca_version) {
              ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
              for (size_type k = 0; k < d; ++k)
                T_t_n(i*d+k, i) = - friction_coeff[i] * vg[k] / r;
            }
            for (size_type k = 0; k < d; ++k) pg(k,k) -= vt1;
            if (Hughes_stabilized) {
              for (size_type k = 0; k < d; ++k)
                for (size_type l = 0; l < d; ++l) {
                  gmm::add(gmm::scaled(gmm::mat_row(DDT, d*i+l), -pg(k,l)),
                           gmm::mat_col(T_t_t, d*i+k));
                  T_t_t(d*i+l, d*i+k) += pg(k,l) / r;
                }
            } else
              gmm::copy(gmm::scaled(pg, vt1/r), gmm::sub_matrix(T_t_t, SUBI));

          }
          if (Hughes_stabilized) {
            model_real_sparse_matrix aux(gmm::mat_nrows(T_t_t), gmm::mat_nrows(T_t_t));
            gmm::copy(gmm::transposed(T_t_t), aux);
            gmm::copy(aux, T_t_t);
          }
          gmm::copy(gmm::transposed(T_u1_t), T_t_u1);
          if (two_variables) gmm::copy(gmm::transposed(T_u2_t), T_t_u2);
        }

        if (symmetrized) {
          // gmm::copy(gmm::transposed(T_n_u1), T_u1_n);  // already done
          // gmm::copy(gmm::transposed(T_n_u2), T_u2_n);  // already done
          model_real_sparse_matrix tmp1(gmm::mat_ncols(BN1),
                                        gmm::mat_ncols(BN1));
          model_real_sparse_matrix tmp2(gmm::mat_ncols(BN2),
                                        gmm::mat_ncols(BN2));
          gmm::mult(gmm::transposed(gmm::scaled(BBN1,-r)), T_n_u1, tmp1);
          gmm::add(tmp1, T_u1_u1);
          if (two_variables) {
            gmm::mult(gmm::transposed(gmm::scaled(BBN2,-r)), T_n_u2, tmp2);
            gmm::add(tmp2, T_u2_u2);
          }

          if (!contact_only) {
            // gmm::copy(gmm::transposed(T_t_u1), T_u1_t);  // already done
            // gmm::copy(gmm::transposed(T_t_u2), T_u2_t);  // already done
            gmm::mult(gmm::transposed(gmm::scaled(BBT1,-r)), T_t_u1, tmp1);
            gmm::add(tmp1, T_u1_u1);
            if (two_variables) {
              gmm::mult(gmm::transposed(gmm::scaled(BBT2,-r)), T_t_u2, tmp2);
              gmm::add(tmp2, T_u2_u2);
            }
          }
        }
        else {
          gmm::copy(gmm::scaled(gmm::transposed(BN1), -vt1), T_u1_n);
          if (two_variables)
            gmm::copy(gmm::scaled(gmm::transposed(BN2), -vt1), T_u2_n);
          if (!contact_only) {
            gmm::copy(gmm::scaled(gmm::transposed(BT1), -vt1), T_u1_t);
            if (two_variables)
              gmm::copy(gmm::scaled(gmm::transposed(BT2), -vt1), T_u2_t);
          }
        }
      }

      if (version & model::BUILD_RHS) {
        for (size_type i=0; i < nbc; ++i) { 
            RLN[i] = std::min(scalar_type(0), RLN[i]);
            if (!contact_only) {
              scalar_type radius = Tresca_version ? threshold[i]
                : -friction_coeff[i]*lambda_n[i];
              ball_projection
                (gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), radius);
            }
        }

        if (symmetrized) {
          gmm::mult_add(gmm::transposed(BN1), RLN, ru1);
          if (two_variables) gmm::mult_add(gmm::transposed(BN2), RLN, ru2);
          if (!contact_only) {
            gmm::mult_add(gmm::transposed(BT1), RLT, ru1);
            if (two_variables) gmm::mult_add(gmm::transposed(BT2), RLT, ru2);
          }
        } else {
          gmm::mult_add(gmm::transposed(BN1), lambda_n, ru1);
          if (two_variables)
            gmm::mult_add(gmm::transposed(BN2), lambda_n, ru2);
          if (!contact_only) {
            gmm::mult_add(gmm::transposed(BT1), lambda_t, ru1);
            if (two_variables)
              gmm::mult_add(gmm::transposed(BT2), lambda_t, ru2);
          }
        }

        
	gmm::add(gmm::scaled(lambda_n, vt1/r), gmm::scaled(RLN,-vt1/r),
		 rlambda_n);
	
	if (!contact_only)
	  gmm::add(gmm::scaled(lambda_t, vt1/r), gmm::scaled(RLT,-vt1/r),
		   rlambda_t);
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
        + (friction_dynamic_term ? 1 : 0);
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

    Coulomb_friction_brick(bool symmetrized_, bool contact_only_,
                           bool two_variables_=false,
                           bool Hughes_stabilized_=false) {
      symmetrized = symmetrized_;
      contact_only = contact_only_;
      is_init = false;
      Tresca_version = false;   // for future version ...
      really_stationary = false;   // for future version ...
      friction_dynamic_term = false;  // for future version ...
      two_variables = two_variables_;
      Hughes_stabilized = Hughes_stabilized_;
      set_flags("Coulomb friction brick", false /* is linear*/,
                /* is symmetric */
                symmetrized && (contact_only || Tresca_version),
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
   bool symmetrized, bool Hughes_stabilized) {
    Coulomb_friction_brick *pbr_=new Coulomb_friction_brick(symmetrized, true, false,
                                                            Hughes_stabilized);
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

  size_type add_basic_contact_with_friction_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &multname_t, const std::string &dataname_r,
   CONTACT_B_MATRIX &BN, CONTACT_B_MATRIX &BT,
   std::string dataname_friction_coeff,
   std::string dataname_gap, std::string dataname_alpha,
   bool symmetrized, bool Hughes_stabilized) {
    Coulomb_friction_brick *pbr_
      = new Coulomb_friction_brick(symmetrized,false, false,Hughes_stabilized);
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
    tl.push_back(model::term_description(multname_t, multname_n, false));
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

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);
    vl.push_back(multname_t);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }


  //=========================================================================
  //
  //  Continuous augmented Lagrangian brick (given obstacle, u, lambda).
  //
  //=========================================================================

  template <typename T> inline static T Heav(T a)
  { return (a < T(0)) ? T(0) : T(1); }

  void friction_nonlinear_term::adjust_tensor_size(void) {
    sizes_.resize(1); sizes_[0] = 1;
    switch (option) {
      // one-dimensional tensors [N]
    case RHS_U_V1:       case RHS_U_V2:       case RHS_U_V3: case RHS_U_V4:
    case RHS_U_V5:       case RHS_U_V6:       case RHS_U_V7:
    case RHS_U_FRICT_V1: case RHS_U_FRICT_V2:
    case RHS_U_FRICT_V3: case RHS_U_FRICT_V4:
    case RHS_L_FRICT_V1: case RHS_L_FRICT_V2: case RHS_L_FRICT_V3:
    case K_UL_V1:        case K_UL_V2:        case K_UL_V3:  case K_UL_V4:
    case UZAWA_PROJ_FRICT:
      sizes_[0] = N; break;
      // two-dimensional tensors [N x N]
    case K_UU_V1: case K_UU_V2:
    case K_UL_FRICT_V1: case K_UL_FRICT_V2: case K_UL_FRICT_V3:
    case K_UL_FRICT_V4: case K_UL_FRICT_V5: case K_UL_FRICT_V6:
    case K_LL_FRICT_V1: case K_LL_FRICT_V2: case K_LL_FRICT_V3:
    case K_UU_FRICT_V1: case K_UU_FRICT_V2:
    case K_UU_FRICT_V3: case K_UU_FRICT_V4: 
      sizes_.resize(2); sizes_[0] = sizes_[1] = N;  break;
    } 
  }

  void friction_nonlinear_term::compute
  (fem_interpolation_context & /* ctx */, bgeot::base_tensor &t) {

    t.adjust_sizes(sizes_);
    scalar_type e, augm_ln; dim_type i, j;

    switch (option) {

    // scalar tensors [1]

    case RHS_L_V1:
      t[0] = (ln+gmm::neg(ln-r*(un - g)))/r; break;
    case RHS_L_V2:
      t[0] = r*(un-g) + gmm::pos(ln); break;

    case K_LL_V1:
      t[0] = (Heav(r*(un-g)-ln) - scalar_type(1))/r; break;
    case K_LL_V2:
      t[0] = -Heav(ln); break;

    case UZAWA_PROJ:
      t[0] = -gmm::neg(ln - r*(un - g)); break;

    // one-dimensional tensors [N]

    case RHS_U_V1:
      for (i=0; i<N; ++i) t[i] = ln * no[i]; break;
    case RHS_U_V2:
      e = -gmm::neg(ln-r*(un - g));
      for (i=0; i<N; ++i) t[i] = e * no[i];
      break;
    case RHS_U_V3:
      e = ln - gmm::pos(un-g) * r;
      for (i=0; i<N; ++i) t[i] = e * no[i];
      break;
    case RHS_U_V4:
      e = -gmm::neg(ln);
      for (i=0; i<N; ++i) t[i] = e * no[i];
      break;
    case RHS_U_V5:
      e = - gmm::pos(un-g) * r;
      for (i=0; i<N; ++i) t[i] = e * no[i];
      break;
    case RHS_U_V6:
      e = - gmm::neg(ln-r*(un - g));
      auxN = lt - zt;  ball_projection(auxN, -f_coeff*e );
      for (i=0; i<N; ++i) t[i] = (e*no[i] + auxN[i]);
      break;
    case RHS_U_V7:
      e = - gmm::neg(-r*(un - g));
      auxN = - zt;  ball_projection(auxN, -f_coeff *e );
      for (i=0; i<N; ++i) t[i] = (e*no[i] + auxN[i]);
      break;
    case RHS_U_FRICT_V1:
      for (i=0; i<N; ++i) t[i] = lnt[i]; break;
    case RHS_U_FRICT_V2:
      e = gmm::neg(ln - r*(un-g));
      auxN = lt - zt;  ball_projection(auxN, -f_coeff * ln);
      for (i=0; i<N; ++i) t[i] = auxN[i] - e*no[i];
      break;
    case RHS_U_FRICT_V3:
      e = r*gmm::pos(un-g);
      for (i=0; i<N; ++i) t[i] = lnt[i] - e*no[i];
      break;
    case RHS_U_FRICT_V4:
      e = -gmm::neg(ln);
      auxN = lt;  ball_projection(auxN, f_coeff * gmm::neg(ln));
      for (i=0; i<N; ++i) t[i] = no[i]*e + auxN[i];
      break;
    case RHS_L_FRICT_V1:
      e = ln+gmm::neg(ln-r*(un-g));
      auxN = zt - lt;  ball_projection(auxN, -f_coeff * ln); auxN += lt;
      for (i=0; i<N; ++i) t[i] = (e*no[i] + auxN[i])/ r;
      break;
    case RHS_L_FRICT_V2:
      e = r*(un-g) + gmm::pos(ln);
      auxN = lt;  ball_projection(auxN, f_coeff * gmm::neg(ln));
      for (i=0; i<N; ++i) t[i] = no[i]*e + zt[i] + lt[i] - auxN[i];
      break;
    case RHS_L_FRICT_V3:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = (lnt[i] - auxN[i])/r;
      break;
    case K_UL_V1:
      for (i=0; i<N; ++i) t[i] = -no[i];
      break;
    case K_UL_V2 :
      e = -Heav(-ln); //Heav(ln)-scalar_type(1);
      for (i=0; i<N; ++i) t[i] = e*no[i];
      break;
    case K_UL_V3:
      e = -Heav(r*(un-g)-ln);
      for (i=0; i<N; ++i) t[i] = e*no[i];
      break;
    case K_UL_V4:
      e = -r;
      for (i=0; i<N; ++i) t[i] = e*no[i];
      break;
    case UZAWA_PROJ_FRICT:
      e = -gmm::neg(ln - r*(un - g));
      auxN = lt - zt;  ball_projection(auxN, -f_coeff * e); 
      for (i=0; i<N; ++i) t[i] = e*no[i] + auxN[i];
      break;


    // two-dimensional tensors [N x N]

    case K_UU_V1:
      e = Heav(un - g) * r;
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = e * no[i] * no[j];
      break;
    case K_UU_V2:
      e = r*Heav(r*(un - g)-ln);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = e * no[i] * no[j];
      break;

    case K_UL_FRICT_V1:
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = ((i == j) ? -scalar_type(1) : scalar_type(0));
      break;
    case K_UL_FRICT_V2:
      e = -Heav(-ln); //Heav(ln)-scalar_type(1);
      ball_projection_grad(lt, f_coeff * gmm::neg(ln), GP);
       e += gmm::vect_sp(GP, no, no);
      ball_projection_grad_r(lt, f_coeff * gmm::neg(ln), V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - GP(i,j) + f_coeff*Heav(-ln)*no[i]*V[j];
      break;
    case K_UL_FRICT_V3:
      e = -Heav(r*(un-g)-ln);
      auxN = lt - zt; ball_projection_grad(auxN, -f_coeff * ln, GP);
      e += gmm::vect_sp(GP, no, no);
      ball_projection_grad_r(auxN, -f_coeff * ln, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - GP(i,j) - f_coeff * V[j] * no[i];
      break;
    case K_UL_FRICT_V4:
      e = -Heav(r*(un-g)-ln);
      auxN = lt - zt; ball_projection_grad(auxN, -f_coeff * ln, GP);
      e += alpha * gmm::vect_sp(GP, no, no);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - alpha*GP(i,j);
      break;
    case K_UL_FRICT_V5:
      e = r*(alpha-scalar_type(1));
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - ((i == j) ? r*alpha : scalar_type(0));
      break;
    case K_UL_FRICT_V6:
      {
	scalar_type nzt = gmm::vect_norm2(zt);
	auxN = lnt - (r*(un-g) - f_coeff * nzt) * no - zt;
	base_matrix A(N, N), B(N, N);
	De_Saxce_projection_grad(auxN, no, f_coeff, A);
	gmm::copy(gmm::identity_matrix(), B); gmm::scale(B, alpha);
	gmm::rank_one_update(B, gmm::scaled(no, scalar_type(1)-alpha), no);
	if (nzt != scalar_type(0)) 
	  gmm::rank_one_update(B, gmm::scaled(no, -f_coeff*alpha/nzt), zt);
	gmm::mult(A, B, GP);
	for (i=0; i<N; ++i) for (j=0; j<N; ++j)
	  t[i*N+j] = -GP(i,j);
      }
      break;
    case K_LL_FRICT_V1:
      e = (Heav(r*(un-g)-ln))/r;
      auxN = lt - zt; ball_projection_grad(auxN, -f_coeff * ln, GP);
      e -= gmm::vect_sp(GP, no, no) / r;
      ball_projection_grad_r(auxN, -f_coeff * ln, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e
          - (((i == j) ? scalar_type(1) : scalar_type(0)) - GP(i,j))/r
          - f_coeff * V[j] * no[i] / r;
      break;
    case K_LL_FRICT_V2:
      e = -Heav(ln) + scalar_type(1);
      ball_projection_grad(lt, f_coeff * gmm::neg(ln), GP);
      e -= gmm::vect_sp(GP, no, no);
      ball_projection_grad_r(lt, f_coeff * gmm::neg(ln), V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - ((i == j) ? scalar_type(1) : scalar_type(0)) + GP(i,j) - f_coeff*Heav(-ln)*no[i]*V[j]; // transposer ?
      break;
    case K_LL_FRICT_V3:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection_grad(auxN, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
	t[i*N+j] = (GP(i,j) - scalar_type((i == j) ? 1 : 0))/r;
      break;
    case K_UU_FRICT_V1:
      e = r*Heav(r*(un-g)-ln);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = no[i]*no[j]*e;
      break;
    case K_UU_FRICT_V2:
      e = Heav(r*(un-g)-ln);
      auxN = lt - zt; ball_projection_grad(auxN, -f_coeff * ln, GP);
      e -= alpha*gmm::vect_sp(GP, no, no);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = r*(no[i]*no[j]*e + alpha*GP(i,j));
      break;
    case K_UU_FRICT_V3:
      e = Heav(r*(un-g)-ln);
      augm_ln = -gmm::neg(ln - r*(un-g));
      auxN = lt - zt; ball_projection_grad(auxN, -f_coeff * augm_ln, GP);
      e -= alpha*gmm::vect_sp(GP, no, no);
      ball_projection_grad_r(auxN, -f_coeff * augm_ln, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = r*(no[i]*no[j]*e + alpha*GP(i,j)
		      - f_coeff*Heav(r*(un-g)-ln)*no[i]*V[j]);
      break;
    case K_UU_FRICT_V4:
      e = Heav(r*(un-g));
      augm_ln = -gmm::neg(- r*(un-g));
      auxN = - zt; ball_projection_grad(auxN, -f_coeff * augm_ln, GP);
      e -= alpha*gmm::vect_sp(GP, no, no);
      ball_projection_grad_r(auxN, -f_coeff * augm_ln, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = r*(no[i]*no[j]*e + alpha*GP(i,j)
		      - f_coeff*Heav(r*(un-g))*no[i]*V[j]);
      break;
    default : GMM_ASSERT1(false, "Invalid option");
    }
  }

  void friction_nonlinear_term::prepare
  (fem_interpolation_context& ctx, size_type nb) {
    size_type cv = ctx.convex_num();

    switch (nb) { // last is computed first
    case 1 : // calculate [un] and [zt] interpolating [U],[WT] on [mf_u]
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (U, gmm::sub_index
                 (mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no);
      if (!contact_only) {
        if (gmm::vect_size(WT) == gmm::vect_size(U)) {
          gmm::copy(gmm::sub_vector
                    (WT, gmm::sub_index
                     (mf_u.ind_basic_dof_of_element(cv))), coeff);
          ctx.pf()->interpolation(ctx, coeff, auxN, N);
          auxN -= gmm::vect_sp(auxN, no) * no;
        } else {
          gmm::clear(auxN);
        }
        zt = ((V - un * no) - auxN) * (r * alpha); // zt = r*alpha*(u_T-w_T)
      }
      break;

    case 2 : // calculate [lnt],[ln],[lt] interpolating [lambda] on [mf_lambda]
      coeff.resize(mf_lambda.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (lambda, gmm::sub_index
                 (mf_lambda.ind_basic_dof_of_element(cv))), coeff);
      if (contact_only) {
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        ln = aux1[0];
      } else {
        ctx.pf()->interpolation(ctx, coeff, lnt, N);
        ln = gmm::vect_sp(lnt, no);
        lt = lnt - ln * no;
      }
      break;

    case 3 : // calculate [g] and [no] interpolating [obs] on [mf_obs]
      coeff.resize(mf_obs.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (obs, gmm::sub_index
                 (mf_obs.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, grad, 1);
      gmm::copy(gmm::mat_row(grad, 0), no);
      no /= -gmm::vect_norm2(no);
      ctx.pf()->interpolation(ctx, coeff, aux1, 1);
      g = aux1[0];
      break;
    case 4 :// calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
      GMM_ASSERT1(!contact_only, "Invalid friction option");
      if (mf_coeff) {
        coeff.resize(mf_coeff->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (friction_coeff, gmm::sub_index
                   (mf_coeff->ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
      }
      break;
    default : GMM_ASSERT1(false, "Invalid option");
    }

  }

  template<typename MAT, typename VECT1>
  void asm_frictionless_continuous_tangent_matrix_Alart_Curnier
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu, const mesh_im &mim,
   const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r, int option,
   const mesh_region &rg) {

    size_type subterm1 = (option == 4) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 4) ? K_UL_V4 : K_UL_V3;
    size_type subterm3 = (option == 4) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    friction_nonlinear_term nterm1(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm1);
    friction_nonlinear_term nterm2(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm2);
    friction_nonlinear_term nterm3(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm3);
    friction_nonlinear_term nterm4(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm4);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4:
      assem.set
       ("M$1(#1,#2)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#2))(i,:,i,:); " // UL
        "M$2(#2,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#2).vBase(#1))(i,:,:,i); " // LU
        "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3).Base(#2).Base(#2))(i,:,:)");    // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#2)+=comp(NonLin$2(#1,#1,#2,#3).vBase(#1).Base(#2))(i,:,i,:); "      // UL
        "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3).Base(#2).Base(#2))(i,:,:);"          // LL
        "M$4(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    case 3:
      assem.set
      ("M$1(#1,#2)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#2))(i,:,i,:); "      // UL
       "M$2(#2,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#2).vBase(#1))(i,:,:,i); "      // LU
       "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3).Base(#2).Base(#2))(i,:,:); "         // LL
       "M$4(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Kul);
    assem.push_mat(Klu);
    assem.push_mat(Kll);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }

  template<typename MAT, typename VECT1>
  void asm_Coulomb_friction_continuous_tangent_matrix_Alart_Curnier
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu, const mesh_im &mim,
   const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r,
   scalar_type alpha, const getfem::mesh_fem *mf_coeff, const VECT1 &f_coeff,
   const VECT1 &WT, int option, const mesh_region &rg) {

    size_type subterm1, subterm2, subterm3;
    switch (option) {
    case 1 : subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V4;
      subterm3 = K_LL_FRICT_V1; break;
    case 2 : subterm1 = K_UL_FRICT_V3; subterm2 = K_UL_FRICT_V4;
      subterm3 = K_LL_FRICT_V1; break;
    case 3 : subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V4;
      subterm3 = K_LL_FRICT_V1; break;
    case 4 : subterm1 = K_UL_FRICT_V2; subterm2 = K_UL_FRICT_V5;
      subterm3 = K_LL_FRICT_V2; break;
    case 5 : subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V6;
      subterm3 = K_LL_FRICT_V3; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    size_type subterm4 = (option == 2) ? K_UU_FRICT_V2 : K_UU_FRICT_V1;

    friction_nonlinear_term nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm1, false, alpha, mf_coeff,
                                   &f_coeff, &WT);
    friction_nonlinear_term nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm2, false, alpha, mf_coeff,
                                   &f_coeff, &WT);
    friction_nonlinear_term nterm3(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm3, false, alpha, mf_coeff,
                                   &f_coeff, &WT);
    friction_nonlinear_term nterm4(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm4, false, alpha, mf_coeff,
                                   &f_coeff, &WT);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4: case 5:
      assem.set
       ("M$1(#1,#2)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1).vBase(#2))(i,j,:,i,:,j); " // UL
        "M$2(#2,#1)+=comp(NonLin$2(#1,#1,#2,#3,#4).vBase(#2).vBase(#1))(i,j,:,j,:,i); " // LU
        "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3,#4).vBase(#2).vBase(#2))(i,j,:,i,:,j)"); // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#2)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1).vBase(#2))(i,j,:,i,:,j); " // UL
        "M$2(#2,#1)+=comp(NonLin$2(#1,#1,#2,#3,#4).vBase(#2).vBase(#1))(i,j,:,j,:,i); " // LU
        "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3,#4).vBase(#2).vBase(#2))(i,j,:,i,:,j);"  // LL
        "M$4(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3,#4).vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    case 3:
      assem.set
      ("M$1(#1,#2)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1).vBase(#2))(i,j,:,i,:,j); " // UL
       "M$2(#2,#1)+=comp(NonLin$2(#1,#1,#2,#3,#4).vBase(#2).vBase(#1))(i,j,:,j,:,i); " // LU
       "M$3(#2,#2)+=comp(NonLin$3(#1,#1,#2,#3,#4).vBase(#2).vBase(#2))(i,j,:,i,:,j); " // LL
       "M$4(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3,#4).vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Kul);
    assem.push_mat(Klu);
    assem.push_mat(Kll);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }

  template<typename VECT1>
  void asm_frictionless_continuous_rhs_Alart_Curnier
  (VECT1 &Ru, VECT1 &Rl, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r, int option,
   const mesh_region &rg) {

    size_type subterm1;
    switch (option) {
    case 1 : subterm1 = RHS_U_V1; break;
    case 2 : subterm1 = RHS_U_V2; break;
    case 3 : subterm1 = RHS_U_V3; break;
    case 4 : subterm1 = RHS_U_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }
    size_type subterm2 = (option == 4) ? RHS_L_V2 : RHS_L_V1;

    friction_nonlinear_term nterm1(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm1);
    friction_nonlinear_term nterm2(mf_u, U, mf_lambda, lambda, mf_obs,
                                   obs, r, subterm2);

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$2(#1,#1,#2,#3).Base(#2))(i,:)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(Ru);
    assem.push_vec(Rl);
    assem.assembly(rg);

  }

  template<typename VECT1>
  void asm_Coulomb_friction_continuous_rhs_Alart_Curnier
  (VECT1 &Ru, VECT1 &Rl, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r,
   scalar_type alpha, const getfem::mesh_fem *mf_coeff, const VECT1 &f_coeff,
   const VECT1 &WT, int option,
   const mesh_region &rg) {

    size_type subterm1, subterm2;
    switch (option) {
    case 1 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V1; break;
    case 2 : subterm1 = RHS_U_FRICT_V2; subterm2 = RHS_L_FRICT_V1; break;
    case 3 : subterm1 = RHS_U_FRICT_V3; subterm2 = RHS_L_FRICT_V1; break;
    case 4 : subterm1 = RHS_U_FRICT_V4; subterm2 = RHS_L_FRICT_V2; break;
    case 5 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V3; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    friction_nonlinear_term nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm1, false, alpha, mf_coeff,
                                   &f_coeff, &WT);
    friction_nonlinear_term nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs,
                                   r, subterm2, false, alpha, mf_coeff,
                                   &f_coeff, &WT);

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$2(#1,#1,#2,#3,#4).vBase(#2))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);    assem.push_vec(Ru);
    assem.push_vec(Rl);
    assem.assembly(rg);
  }

  struct Coulomb_friction_continuous_brick : public virtual_brick {

    bool Tresca_version, contact_only;
    int option;

    // option = 1 : Alart-Curnier
    // option = 2 : symmetric Alart-Curnier (almost with friction)
    // option = 3 : Alart Curnier "over-augmented"
    // option = 4 : New method
    // option = 5 : De-Saxcé
    // option = 6 : De-Saxcé symmetric ?

    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1,
                  "Continuous Coulomb friction brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 2,
                  "Continuous Coulomb friction brick need two variables");
      GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 5,
                  "Wrong number of data for continuous Coulomb friction "
                  << "brick, " << dl.size() << " should be between 2 and 5.");
      GMM_ASSERT1(matl.size() == size_type(3 + (option == 3)
                                           + (option == 2 && !contact_only)),
                  "Wrong number of terms for "
                  "continuous Coulomb friction brick");

      // variables : u, lambda. The variable lambda should be scalar in the
      //             frictionless case and vector valued in the case with
      //             friction.
      // data      : obstacle, r for the version without friction
      //           : obstacle, r, friction_coeff, alpha, w_t for the
      //             version with friction. alpha and w_t are optional and
      //             equal to 1 and 0 by default, respectively.

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));
      const model_real_plain_vector &lambda = md.real_variable(vl[1]);
      const mesh_fem &mf_lambda = *(md.pmesh_fem_of_variable(vl[1]));
      GMM_ASSERT1(mf_lambda.get_qdim() == (contact_only ? 1 : mf_u.get_qdim()),
                  "The contact stress has not the right dimension");
      const model_real_plain_vector &obstacle = md.real_variable(dl[0]);
      const mesh_fem &mf_obstacle = *(md.pmesh_fem_of_variable(dl[0]));
      size_type sl = gmm::vect_size(obstacle) * mf_obstacle.get_qdim()
        / mf_obstacle.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vr = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      const mesh_im &mim = *mims[0];

      const model_real_plain_vector &friction_coeff
        = contact_only ? u : md.real_variable(dl[2]);
      const mesh_fem *mf_coeff = contact_only ? 0 : md.pmesh_fem_of_variable(dl[2]);
      sl = gmm::vect_size(friction_coeff);
      if (mf_coeff) { sl *= mf_coeff->get_qdim(); sl /= mf_coeff->nb_dof(); }
      GMM_ASSERT1(sl == 1 || contact_only,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4) {
        alpha = md.real_variable(dl[3])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3])) == 1,
                    "Parameter alpha should be a scalar");
      }

      model_real_plain_vector voidvec;
      const model_real_plain_vector &WT
        = (!contact_only && dl.size()>=5) ? md.real_variable(dl[4]) : voidvec;

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Continuous Coulomb friction tangent term");
        gmm::clear(matl[0]); gmm::clear(matl[1]); gmm::clear(matl[2]);
        if (matl.size() >= 4) gmm::clear(matl[3]);
        if (contact_only) {
          size_type fourthmat = (matl.size() >= 4) ? 3 : 1;
          asm_frictionless_continuous_tangent_matrix_Alart_Curnier
            (matl[0], matl[1], matl[2], matl[fourthmat], mim, mf_u, u,
             mf_lambda, lambda, mf_obstacle, obstacle, vr[0], option, rg);
        }
        else {
          size_type fourthmat = (matl.size() >= 4) ? 3 : 1;
          asm_Coulomb_friction_continuous_tangent_matrix_Alart_Curnier
            (matl[0], matl[1], matl[2], matl[fourthmat], mim, mf_u, u,
             mf_lambda, lambda, mf_obstacle, obstacle, vr[0], alpha, mf_coeff,
             friction_coeff, WT, option, rg);
        }
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]); gmm::clear(vecl[1]); gmm::clear(vecl[2]);
        if (matl.size() >= 4) gmm::clear(vecl[3]);

        if (contact_only)
          asm_frictionless_continuous_rhs_Alart_Curnier
            (vecl[0], vecl[2], mim, mf_u, u, mf_lambda, lambda,
             mf_obstacle, obstacle, vr[0], option, rg);
        else
          asm_Coulomb_friction_continuous_rhs_Alart_Curnier
            (vecl[0], vecl[2], mim, mf_u, u, mf_lambda, lambda,
             mf_obstacle, obstacle, vr[0], alpha, mf_coeff,
             friction_coeff, WT, option, rg);
      }

    }

    Coulomb_friction_continuous_brick(bool contact_only_, int option_) {
      Tresca_version = false;   // for future version ...
      option = option_;
      contact_only = contact_only_;
      set_flags("Continuous Coulomb friction brick", false /* is linear*/,
                (option==2) && contact_only /* is symmetric */,
                false /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };


  //=========================================================================
  //  Add a frictionless contact condition with a rigid obstacle given
  //  by a level set.
  //=========================================================================

  size_type add_continuous_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_obs,
   const std::string &dataname_r, size_type region, int option) {

    pbrick pbr = new Coulomb_friction_continuous_brick(true, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 4 :
      tl.push_back(model::term_description(varname_u, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u, false));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u, multname_n, true));
      tl.push_back(model::term_description(varname_u, varname_u, true));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      break;
    case 3 :
      tl.push_back(model::term_description(varname_u, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u, false));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      tl.push_back(model::term_description(varname_u, varname_u, true));
      break;
    default :GMM_ASSERT1(false,
                         "Incorrect option for continuous contact brick");

    }
    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //  Add a contact condition with Coulomb friction with a rigid obstacle
  //  given by a level set.
  //=========================================================================

  size_type add_continuous_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_obs,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, int option,
   const std::string &dataname_alpha, const std::string &dataname_wt) {

    pbrick pbr
      = new Coulomb_friction_continuous_brick(false, option);

    model::termlist tl;

    switch (option) {
    case 1: case 4: case 5:
      tl.push_back(model::term_description(varname_u, multname, false));
      tl.push_back(model::term_description(multname, varname_u, false));
      tl.push_back(model::term_description(multname, multname, true));
      break;
    case 2:
      tl.push_back(model::term_description(varname_u, multname, false));
      tl.push_back(model::term_description(multname, varname_u, false));
      tl.push_back(model::term_description(multname, multname, true));
      tl.push_back(model::term_description(varname_u, varname_u, true));
      break;
    case 3:
      tl.push_back(model::term_description(varname_u, multname, false));
      tl.push_back(model::term_description(multname, varname_u, false));
      tl.push_back(model::term_description(multname, multname, true));
      tl.push_back(model::term_description(varname_u, varname_u, true));
      break;
    default :GMM_ASSERT1(false,
                         "Incorrect option for continuous contact brick");
    }
    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size()) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size()) dl.push_back(dataname_wt);
    }

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //
  //  Continuous penalized contact with friction (given obstacle, u, lambda).
  //
  //=========================================================================


  template<typename MAT, typename VECT1>
  void asm_frictionless_penalized_tangent_matrix
  (MAT &Kuu, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, int option, const mesh_region &rg) {

    friction_nonlinear_term nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs,
				  r, (option == 1) ? K_UU_V1 : K_UU_V2);

    getfem::generic_assembly assem;
    assem.set
    ("M(#1,#1)+=comp(NonLin(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }


  template<typename VECT1>
  void asm_frictionless_penalized_rhs
  (VECT1 &Ru, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U,  const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, int option, const mesh_region &rg) {

    friction_nonlinear_term nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs,
				  r, (option == 1) ? RHS_U_V5 : RHS_U_V2);
    
    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru);
    assem.assembly(rg);
  }

  template<typename MAT, typename VECT1>
  void asm_Coulomb_friction_penalized_tangent_matrix
  (MAT &Kuu, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, scalar_type alpha, const getfem::mesh_fem *mf_coeff,
   const VECT1 &f_coeff, const VECT1 &WT, int option, const mesh_region &rg) {

    friction_nonlinear_term nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r,
				  (option==1) ? K_UU_FRICT_V4 : K_UU_FRICT_V3,
				  false, alpha, mf_coeff, &f_coeff, &WT);

    getfem::generic_assembly assem;
    assem.set
    ("M(#1,#1)+=comp(NonLin(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }


  template<typename VECT1>
  void asm_Coulomb_friction_penalized_rhs
  (VECT1 &Ru, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U,  const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, scalar_type alpha, const getfem::mesh_fem *mf_coeff,
   const VECT1 &f_coeff, const VECT1 &WT, int option, const mesh_region &rg) {

    friction_nonlinear_term nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r,
				  (option == 1) ? RHS_U_V7 : RHS_U_V6,
				  false, alpha, mf_coeff, &f_coeff, &WT);
    
    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru);
    assem.assembly(rg);
  }


  struct penalized_Coulomb_friction_brick : public virtual_brick {

    bool Tresca_version, contact_only;
    int option;
    
    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1,
                  "Penalized Coulomb friction brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 1,
                  "Penalized Coulomb friction brick need two variables");
      size_type nb_data_1 = ((option == 1) ? 2 : 3), nb_data_2 = nb_data_1;
      if (!contact_only) nb_data_2 += 3;
      GMM_ASSERT1(dl.size() >= nb_data_1 && dl.size() <= nb_data_2,
                  "Wrong number of data for penalized Coulomb friction "
                  << "brick, " << dl.size() << " should be between "
		  << nb_data_1 << " and " << nb_data_2 << ".");
      GMM_ASSERT1(matl.size() == 1, "Wrong number of terms for "
                  "penalized Coulomb friction brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));
      const model_real_plain_vector &obstacle = md.real_variable(dl[0]);
      const mesh_fem &mf_obstacle = *(md.pmesh_fem_of_variable(dl[0]));
      size_type sl = gmm::vect_size(obstacle) * mf_obstacle.get_qdim()
        / mf_obstacle.nb_dof();
      size_type N = mf_u.linked_mesh().dim();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vr = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      const mesh_im &mim = *mims[0];


      const model_real_plain_vector &lambda
        = (option == 1) ? obstacle : md.real_variable(dl[2]);
      const mesh_fem *mf_lambda = (option == 1) ? &mf_obstacle : md.pmesh_fem_of_variable(dl[2]);
      sl = gmm::vect_size(lambda);
      sl *= mf_lambda->get_qdim(); sl/=mf_lambda->nb_dof();
      GMM_ASSERT1(sl == (contact_only ? 1 : N) || option == 1,
                  "the data corresponding to the contact stress "
                  "has not the right format");

      size_type shift = ((option == 1) ? 0 : 1);
      const model_real_plain_vector &friction_coeff
        = contact_only ? u : md.real_variable(dl[2+shift]);
      const mesh_fem *mf_coeff = contact_only ? 0 : md.pmesh_fem_of_variable(dl[2+shift]);
      sl = gmm::vect_size(friction_coeff);
      if (mf_coeff) { sl *= mf_coeff->get_qdim(); sl /= mf_coeff->nb_dof(); }
      GMM_ASSERT1(sl == 1 || contact_only,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4+shift) {
        alpha = md.real_variable(dl[3+shift])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3+shift])) == 1,
                    "Parameter alpha should be a scalar");
      }

      model_real_plain_vector voidvec;
      const model_real_plain_vector &WT
        = (!contact_only && dl.size()>=5+shift) ?
	md.real_variable(dl[4+shift]) : voidvec;

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
	GMM_TRACE2("Penalized Coulomb friction tangent term");
	gmm::clear(matl[0]);
        if (contact_only)
          asm_frictionless_penalized_tangent_matrix
            (matl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
	     vr[0], option, rg);
        else
          asm_Coulomb_friction_penalized_tangent_matrix
            (matl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
	     vr[0], alpha, mf_coeff, friction_coeff, WT, option, rg);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        if (contact_only)
          asm_frictionless_penalized_rhs
            (vecl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
	     vr[0], option, rg);
        else
	  asm_Coulomb_friction_penalized_rhs
            (vecl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
	     vr[0], alpha, mf_coeff, friction_coeff, WT, option, rg);
      }

    }

    penalized_Coulomb_friction_brick(bool contact_only_, int option_) {
      Tresca_version = false;   // for future version ...
      contact_only = contact_only_;
      option = option_;
      set_flags(contact_only ? "Continuous penalized contact brick"
		: "Continuous penalized Coulomb friction brick",
		false /* is linear*/, contact_only /* is symmetric */,
                true /* is coercive */, true /* is real */,
		false /* is complex */);
    }

  };


  //=========================================================================
  //  Add a frictionless contact condition with a rigid obstacle given
  //  by a level set.
  //=========================================================================

  size_type add_penalized_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   size_type region, int option, const std::string &dataname_n) {

    pbrick pbr = new penalized_Coulomb_friction_brick(true, option);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, true));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    switch (option) {
    case 1: break;
    case 2: dl.push_back(dataname_n); break;
    default: GMM_ASSERT1(false, "Penalized contact brick : invalid option");
    }

    model::varnamelist vl(1, varname_u);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  //=========================================================================
  //  Add a contact condition with friction with a rigid obstacle given
  //  by a level set.
  //=========================================================================

  size_type add_penalized_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   const std::string &dataname_friction_coeff, 
   size_type region, int option, const std::string &dataname_lambda,
   const std::string &dataname_alpha, const std::string &dataname_wt) {

    pbrick pbr = new penalized_Coulomb_friction_brick(false, option);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, true));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    switch (option) {
    case 1: break;
    case 2: dl.push_back(dataname_lambda); break;
    default: GMM_ASSERT1(false, "Penalized contact brick : invalid option");
    }
    dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size() > 0) dl.push_back(dataname_wt);
    }
    
    model::varnamelist vl(1, varname_u);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
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
    (bool symmetrized_, bool contact_only_, const std::string &obs)
      : Coulomb_friction_brick(symmetrized_, contact_only_), obstacle(obs) {}

  };


  //=========================================================================
  //  Add a frictionless contact condition with a rigid obstacle given
  //  by a signed distance.
  //=========================================================================

  size_type add_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_r,
   size_type region, const std::string &obstacle, bool symmetrized) {
    pbrick pbr
      = new Coulomb_friction_brick_rigid_obstacle(symmetrized, true, obstacle);

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

  size_type add_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, const std::string &obstacle, bool symmetrized) {
    pbrick pbr
      = new Coulomb_friction_brick_rigid_obstacle(symmetrized,false,obstacle);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    tl.push_back(model::term_description(varname_u, multname_t, false));
    tl.push_back(model::term_description(multname_t, varname_u, false));
    tl.push_back(model::term_description(multname_t, multname_t, false));
    tl.push_back(model::term_description(multname_t, multname_n, false));
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
      (bool symmetrized_, bool contact_only_, bool two_variables_,
       const std::vector<size_type> &rg1_, const std::vector<size_type> &rg2_,
       bool slave1_=true, bool slave2_=false)
      : Coulomb_friction_brick(symmetrized_, contact_only_, two_variables_),
        rg1(rg1_), rg2(rg2_), slave1(slave1_), slave2(slave2_) {}

  };


  //=========================================================================
  //  Add a frictionless contact condition between two faces of one or two
  //  elastic bodies.
  //=========================================================================

  size_type add_nonmatching_meshes_contact_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1, bool slave2, bool symmetrized) {

    bool two_variables = (varname_u1.compare(varname_u2) != 0);

    pbrick pbr = new Coulomb_friction_brick_nonmatching_meshes
           (symmetrized, true, two_variables, rg1, rg2, slave1, slave2);

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

  size_type add_nonmatching_meshes_contact_with_friction_brick
  (model &md, const mesh_im &mim1, const mesh_im &mim2,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, std::string &multname_t,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::vector<size_type> &rg1, const std::vector<size_type> &rg2,
   bool slave1, bool slave2, bool symmetrized) {

    bool two_variables = (varname_u1.compare(varname_u2) != 0);

    pbrick pbr = new Coulomb_friction_brick_nonmatching_meshes
           (symmetrized, false, two_variables, rg1, rg2, slave1, slave2);

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
    tl.push_back(model::term_description(multname_t, multname_n, false));

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
