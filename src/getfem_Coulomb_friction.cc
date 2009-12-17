// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2010 Yves Renard, Konstantinos Poulios.
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

namespace getfem {

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
  
  template<class VEC, class VECR>
  static void ball_projection_grad_r(const VEC &x, scalar_type radius,
                                     VECR &g) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius > 0 && a >= radius)
      gmm::copy(gmm::scaled(x, scalar_type(1)/a), g);
    else gmm::clear(g);
  }
  
  template <class VEC, class MAT>
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

  // contact_node is an object which contains data about nodes expected to
  // participate in a contact condition either. A contact node refers to a
  // specific mesh_fem.
  struct contact_node {
    const mesh_fem *mf;          // Pointer to the mesh_fem the contact node is
                                 // associated with
    size_type dof;               // first dof id of the node in the considered mesh_fem
    std::vector<size_type> cvs;  // list of ids of neigbouring convexes
    std::vector<short_type> fcs; // list of local ids of neigbouring faces

    contact_node(const mesh_fem &mf_) {mf = &mf_;}
    contact_node() : mf(0) {}
  };

  // The object contact_node_list holds a list of all nodes which participate
  // in the definition of a contact.
  struct contact_node_list : public std::vector<contact_node> {

    contact_node_list() : std::vector<contact_node>() {}

    void append(const mesh_fem &mf, size_type contact_region,
                dal::bit_vector &bv_cn) {

      const mesh &m = mf.linked_mesh();
      size_type qdim = mf.get_qdim();
      std::map<size_type, size_type> dof_to_cnid;
      size_type cnid = this->size();
      dal::bit_vector dofs = mf.basic_dof_on_region(contact_region);
      for (dal::bv_visitor dof(dofs); !dof.finished(); ++dof) {
        if ( dof % qdim == 0) {
          bv_cn.add(cnid);
          dof_to_cnid[dof] = cnid++;
          contact_node new_cn(mf);
          new_cn.dof = dof;
          this->push_back(new_cn);
        }
      }
      for (mr_visitor face(m.region(contact_region));
           !face.finished(); ++face) {
        assert(face.is_face());
        mesh_fem::ind_dof_face_ct
          face_dofs = mf.ind_basic_dof_of_face_of_element(face.cv(),face.f());
        for (size_type it=0; it < face_dofs.size(); it += qdim ) {
          size_type dof = face_dofs[it];
          cnid = dof_to_cnid[dof];
          (*this)[cnid].cvs.push_back(face.cv());
          (*this)[cnid].fcs.push_back(face.f());
        } // for:it
      } // for:face
    } // append

  }; // struct contact_node_list

  // contact_node's pair
  struct contact_node_pair {
    contact_node *cn_s, *cn_m;  // Pointers to the slave and master contact_node's
    scalar_type dist;           // Distance between slave and master nodes
    contact_node_pair(scalar_type threshold=10.)
      {cn_s = 0, cn_m = 0, dist = threshold;}
  };

  // contact_node's pair
  struct contact_node_pair_list : public std::vector<contact_node_pair> {
    contact_node_list cnl;

    void append_min_dist_cn_pairs(const mesh_fem &mf, 
                                  size_type rg1, size_type rg2,
                                  bool slave1=true, bool slave2=false) {

      dal::bit_vector bv_cn1, bv_cn2;
      cnl.append(mf, rg1, bv_cn1);
      cnl.append(mf, rg2, bv_cn2);

      // Find minimum distance node pairs, FIXME: to be optimized
      size_type size0 = this->size();
      size_type size1 = slave1 ? bv_cn1.card() : 0;
      size_type size2 = slave2 ? bv_cn2.card() : 0;
      this->resize( size0 + size1 + size2 );
      size_type i=size0;
      for (dal::bv_visitor i1(bv_cn1); !i1.finished(); ++i1, ++i) {
        contact_node *cn1 = &cnl[i1];
        base_node node1 = cn1->mf->point_of_basic_dof(cn1->dof);
        size_type j=size0+size1;
        for (dal::bv_visitor i2(bv_cn2); !i2.finished(); ++i2, ++j) {
          contact_node *cn2 = &cnl[i2];
          base_node node2 = cn2->mf->point_of_basic_dof(cn2->dof);
          scalar_type dist = gmm::vect_norm2(node1-node2);
          if (slave1 && dist < (*this)[i].dist) {
            (*this)[i].cn_s = cn1;
            (*this)[i].cn_m = cn2;
            (*this)[i].dist = dist;
          }
          if (slave2 && dist < (*this)[j].dist) {
            (*this)[j].cn_s = cn2;
            (*this)[j].cn_m = cn1;
            (*this)[j].dist = dist;
          }
        } // cn2
      } // cn1
    }
  };


  void calculate_contact_matrices(const mesh_fem &mf_disp,
                                  contact_node_pair_list &cnpl,
                                  model_real_plain_vector &gap,
                                  CONTACT_B_MATRIX &BN) {

    dim_type qdim = mf_disp.get_qdim(); //FIXME: Check if qdim and N
 // size_type N = mesh.dim();           //       are considered correctly

    GMM_ASSERT1(gmm::mat_nrows(BN) == cnpl.size(),
                "Wrong number of contact node pairs");
    gmm::clear(BN);
    gmm::fill(gap, scalar_type(10));  //FIXME: Needs a threshold value
    for (size_type row = 0; row < cnpl.size(); ++row) {
      contact_node_pair *cnp = &cnpl[row];
      if (cnp->cn_s) {
        contact_node *cn_s = cnp->cn_s;  //slave contact node
        contact_node *cn_m = cnp->cn_m;  //master contact node
        const mesh &mesh_m = cn_m->mf->linked_mesh();
        base_node slave_node = cn_s->mf->point_of_basic_dof(cn_s->dof);
        base_node master_node = cn_m->mf->point_of_basic_dof(cn_m->dof);
        base_node un_sel(3), proj_node_sel(3), proj_node_ref_sel(3);
        scalar_type is_in_min = 1e5;  //FIXME
        size_type cv_sel, fc_sel;
        std::vector<size_type>::iterator cv;
        std::vector<short_type>::iterator fc;
        for (cv = cn_m->cvs.begin(), fc = cn_m->fcs.begin();
             cv != cn_m->cvs.end() && fc != cn_m->fcs.end(); cv++, fc++) {
          base_node un = mesh_m.normal_of_face_of_convex(*cv,*fc); //FIXME: this normal is just an approximation
          un /= gmm::vect_norm2(un);
          base_node proj_node(3), proj_node_ref(3);
          //FIXME: the following projection calculation is accurate only for linear triangle elements
          //proj_node = slave_node - [(slave_node-master_node)*n] * n
          gmm::add(master_node, gmm::scaled(slave_node, -1.), proj_node);
          gmm::copy(gmm::scaled(un, gmm::vect_sp(proj_node, un)), proj_node);
          gmm::add(slave_node, proj_node);

          bgeot::pgeometric_trans pgt = mesh_m.trans_of_convex(*cv);
          bgeot::geotrans_inv_convex gic;
          gic.init(mesh_m.points_of_convex(*cv), pgt);
          gic.invert(proj_node, proj_node_ref);
          scalar_type is_in = pgt->convex_ref()->is_in(proj_node_ref);
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

          if (cn_s->mf == &mf_disp) {
            BN(row, cn_s->dof)   += -un_sel[0];
            BN(row, cn_s->dof+1) += -un_sel[1];
            BN(row, cn_s->dof+2) += -un_sel[2];
          }

          if (cn_m->mf == &mf_disp) {
            base_matrix G;
            base_matrix M(qdim, mf_disp.nb_basic_dof_of_element(cv_sel));
            bgeot::vectors_to_base_matrix(G, mesh_m.points_of_convex(cv_sel));
            pfem pf = mf_disp.fem_of_element(cv_sel);
            bgeot::pgeometric_trans pgt = mesh_m.trans_of_convex(cv_sel);
            fem_interpolation_context
              ctx(pgt, pf, proj_node_ref_sel, G, cv_sel, fc_sel);
            pf->interpolation (ctx, M, qdim);

            model_real_plain_vector tmpvec(mf_disp.nb_basic_dof_of_element(cv_sel));
            gmm::mult(gmm::transposed(M), un_sel, tmpvec);

            mesh_fem::ind_dof_ct master_dof = mf_disp.ind_basic_dof_of_element(cv_sel);
            size_type j = 0;
            for (mesh_fem::ind_dof_ct::const_iterator it = master_dof.begin();
                 it != master_dof.end(); it++, j++) {
              BN(row, *it) += tmpvec[j];
            }
          }
        }
      } // if:cnp->cn_s
    } // cnp

  } // calculate_contact_matrices



  //=========================================================================
  //
  //  Basic Brick (with given BN, BT, gap) and possibly two bodies
  //
  //=========================================================================

  struct Coulomb_friction_brick : public virtual_brick {

    mutable CONTACT_B_MATRIX BN1, BT1;
    mutable CONTACT_B_MATRIX BN2, BT2;
    mutable CONTACT_B_MATRIX BBN1, BBT1;
    mutable CONTACT_B_MATRIX BBN2, BBT2;
    mutable model_real_plain_vector gap, threshold, friction_coeff, alpha;
    mutable model_real_plain_vector RLN, RLT; 
    mutable scalar_type r, gamma;
    mutable bool is_init;
    bool Tresca_version, symmetrized, contact_only;
    bool really_stationary, friction_dynamic_term;
    bool two_variables;

    void init_BBN_BBT(void) const {
      gmm::resize(BBN1, gmm::mat_nrows(BN1), gmm::mat_ncols(BN1));
      gmm::copy(BN1, BBN1);
      if (two_variables) {
        gmm::resize(BBN2, gmm::mat_nrows(BN2), gmm::mat_ncols(BN2));
        gmm::copy(BN2, BBN2);
      }
      if (!contact_only) {
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
        if (two_variables)
          gmm::scale(gmm::mat_row(BBN2, i), alpha[i]);
        if (!contact_only)
          for (size_type k = 0; k < d; ++k) {
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

      const scalar_type vt1(1);
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
            if (!really_stationary) {
              for (size_type k = 0; k < d; ++k) {
                gmm::copy(gmm::scaled(gmm::mat_row(BBT1, i*d+k), -pg[k]),
                          gmm::mat_col(T_u1_t, i*d+k));
                if (two_variables)
                  gmm::copy(gmm::scaled(gmm::mat_row(BBT2, i*d+k), -pg[k]),
                            gmm::mat_col(T_u2_t, i*d+k));
              }
            }
            
            if (!Tresca_version) {
              ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
              for (size_type k = 0; k < d; ++k)
                T_t_n(i*d+k, i) = - friction_coeff[i] * vg[k] / r;
            }
            for (size_type k = 0; k < d; ++k) pg(k,k) -= vt1;
            gmm::copy(gmm::scaled(pg, vt1/r), gmm::sub_matrix(T_t_t, SUBI));
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
        
        gmm::add(gmm::scaled(lambda_n, vt1/r), gmm::scaled(RLN, -vt1/r),
                 rlambda_n);

        if (!contact_only)
          gmm::add(gmm::scaled(lambda_t, vt1/r), gmm::scaled(RLT, -vt1/r),
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

    Coulomb_friction_brick(bool symmetrized_, bool contact_only_) {
      symmetrized = symmetrized_;
      contact_only = contact_only_;
      is_init = false;
      Tresca_version = false;   // for future version ...
      really_stationary = false;   // for future version ...
      friction_dynamic_term = false;  // for future version ...
      two_variables = false;  // for future version ...
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

    void set_BT1(CONTACT_B_MATRIX &BT1_) {
      gmm::resize(BT1, gmm::mat_nrows(BT1_), gmm::mat_ncols(BT1_));
      gmm::copy(BT1_, BT1);
      is_init = false;
    }

    CONTACT_B_MATRIX &get_BN1(void) { return BN1; }
    CONTACT_B_MATRIX &get_BT1(void) { return BT1; }
    const CONTACT_B_MATRIX &get_BN1(void) const { return BN1; }
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
   bool symmetrized) {
    Coulomb_friction_brick *pbr_=new Coulomb_friction_brick(symmetrized,true);
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
   bool symmetrized) {
    Coulomb_friction_brick *pbr_=new Coulomb_friction_brick(symmetrized,false);
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
                
                // Computation of an orthonormal basis to un.
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








  // Experimental brick including gap, BN, BT calculation
  // To be done:
  // - Calculation of BT
  // - Optimization of minimum distance node pairs search
  // - Check the compatibility of the provided mesh_fem's
  // - Non-Linear projection
  // - Support for two displacement variables
  // - Large deformations: what happens when cnpl and nbc change during
  //   the iterative solution?

  //=========================================================================
  //
  //  Brick with elastic bodies (one or two bodies, build BN, BT, Gap, alpha)
  //
  //=========================================================================


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

      GMM_ASSERT1(mims.size() == 1, "This contact brick needs one mesh_im");
      const mesh_im &mim = *mims[0];

      // Variables
      // Without friction and one displacement  : u1, lambda_n
      // With friction and one displacement     : u1, lambda_n, lambda_t      //FIXME: to be done
      // Without friction and two displacements : u1, u2, lambda_n            //FIXME: to be done
      // With friction and two displacements    : u1, u2, lambda_n, lambda_t  //FIXME: to be done
      const model_real_plain_vector &u1 = md.real_variable(vl[0]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector &lambda_n = md.real_variable(vl[1]);

      size_type nbc = lambda_n.size();

      // Parameters (order: r, ...)
      size_type np = 0;
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

        contact_node_pair_list cnpl;
        for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it)
          cnpl.append_min_dist_cn_pairs(mf_u1, rg1[it], rg2[it], slave1, slave2);

        // Computation of gap and BN
        gmm::resize(gap, nbc);
        gmm::resize(BN1, nbc, mf_u1.nb_dof());
        calculate_contact_matrices(mf_u1, cnpl, gap, BN1);

        // computation of alpha vector.
        base_node Pmin, Pmax;
        mf_u1.linked_mesh().bounding_box(Pmin, Pmax);
        scalar_type l = scalar_type(0);
        for (size_type i = 0; i < Pmin.size(); ++i)
          l = std::max(l, gmm::abs(Pmax[i] - Pmin[i]));

        CONTACT_B_MATRIX MM(mf_u1.nb_dof(), mf_u1.nb_dof());
        gmm::resize(alpha, nbc);
        size_type mult_id = 0;
        for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it) {
          for (size_type swap = 0; swap < 1; ++swap) {
            if (swap ? slave2 : slave1) {
              size_type rg = swap ? rg2[it] : rg1[it];
              asm_mass_matrix(MM, mim, mf_u1, rg);
              size_type qdim = mf_u1.get_qdim();
              dal::bit_vector rg_dofs = mf_u1.basic_dof_on_region(rg);
              for (dal::bv_visitor id(rg_dofs); !id.finished(); ++id)
                if (id % qdim == 0) alpha[mult_id++] = MM(id, id) / l;
            }
          }
        }
      }

      const model_real_plain_vector dummy_lambda_t, dummy_wt;
      basic_asm_real_tangent_terms
        (u1, u1, lambda_n,dummy_lambda_t, dummy_wt, dummy_wt,
         matl, vecl, version);
    }

    Coulomb_friction_brick_nonmatching_meshes
      (bool symmetrized_, bool contact_only_,
       std::vector<size_type> rg1_, std::vector<size_type> rg2_,
       bool slave1_=true, bool slave2_=false)
      : Coulomb_friction_brick(symmetrized_, contact_only_),
        rg1(rg1_), rg2(rg2_), slave1(slave1_), slave2(slave2_) {}

  };


  //=========================================================================
  //  Add a frictionless contact condition between two faces of an elastic
  //  body (one displacement field).
  //=========================================================================
  size_type add_frictionless_contact_brick
  (model &md, const mesh_im &mim,
   const std::string &varname_u1, const std::string &varname_u2,
   std::string &multname_n, const std::string &dataname_r,
   std::vector<size_type> rg1, std::vector<size_type>  rg2,
   bool slave1, bool slave2, bool symmetrized) {

    bool two_variables = (varname_u1.compare(varname_u2) != 0);
    GMM_ASSERT1(!two_variables, "Two variables support not implemented yet");

    pbrick pbr = new Coulomb_friction_brick_nonmatching_meshes
                     (symmetrized, true, rg1, rg2, slave1, slave2);

    // Calculate multipliers size
    const mesh_fem &mf_u1 = md.mesh_fem_of_variable(varname_u1);
    const mesh_fem &mf_u2 = md.mesh_fem_of_variable(varname_u2);
    size_type nbc = 0;
    for (size_type it = 0; it < rg1.size() && it < rg2.size(); ++it) {
      for (size_type swap = 0; swap < 1; ++swap) {
        if (swap ? slave2 : slave1) {
          const mesh_fem &mf = swap ? mf_u2 : mf_u1;
          size_type rg = swap ? rg2[it] : rg1[it];
          dal::bit_vector rg_dofs = mf.basic_dof_on_region(rg);
          nbc += rg_dofs.card() / mf.get_qdim();
        }
      }
    }

    if (multname_n.size() == 0) {
      multname_n = md.new_name("contact_multiplier");
    }
    // FIXME: Assert multname_n is not defined already in md
    md.add_fixed_size_variable(multname_n, nbc);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, false)); // FIXME: second displacement variable
    tl.push_back(model::term_description(varname_u1, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u1, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));

    // Variables (order: varname_u, multname_n)
    model::varnamelist vl;
    vl.push_back(varname_u1);
    if (two_variables) vl.push_back(varname_u2);
    vl.push_back(multname_n);

    // Parameters (order: r, ...)
    model::varnamelist dl;
    dl.push_back(dataname_r);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), size_type(-1));
  }

}  /* end of namespace getfem.                                             */
