/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2013 Yves Renard, Konstantinos Poulios.

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

#include "getfem/bgeot_rtree.h"
#include "getfem/getfem_contact_and_friction_integral.h"
#include "getfem/getfem_contact_and_friction_common.h"
#include "getfem/getfem_assembling.h"
#include "gmm/gmm_condition_number.h"

#include <getfem/getfem_arch_config.h>
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif

namespace getfem {


  //=========================================================================
  // Augmented friction law
  //=========================================================================


#define FRICTION_LAW 1


#if FRICTION_LAW == 1 // Complete law with friction

  template <typename VEC, typename VEC2, typename VECR>
  void aug_friction(const VEC &lambda, scalar_type g, const VEC &Vs,
                    const VEC &n, scalar_type r, const VEC2 &f, VECR &F) {
    scalar_type nn = gmm::vect_norm2(n);
    scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
    scalar_type lambdan_aug = gmm::neg(lambdan + r * g);
    size_type i = gmm::vect_size(f);
    scalar_type tau = ((i >= 3) ? f[2] : scalar_type(0)) + f[0]*lambdan_aug;
    if (i >= 2) tau = std::min(tau, f[1]);

    if (tau > scalar_type(0)) {
      gmm::add(lambda, gmm::scaled(Vs, -r), F);
      scalar_type mu = gmm::vect_sp(F, n)/nn;
      gmm::add(gmm::scaled(n, -mu/nn), F);
      scalar_type norm = gmm::vect_norm2(F);
      if (norm > tau) gmm::scale(F, tau / norm);
    } else { gmm::clear(F); }

    gmm::add(gmm::scaled(n, -lambdan_aug/nn), F);
  }

  template <typename VEC, typename VEC2, typename VECR, typename MAT>
  void aug_friction_grad(const VEC &lambda, scalar_type g, const VEC &Vs,
                         const VEC &n, scalar_type r, const VEC2 &f, VECR &F,
                         MAT &dlambda, VECR &dg, MAT &dn, MAT &dVs) {
    size_type N = gmm::vect_size(lambda);
    scalar_type nn = gmm::vect_norm2(n);
    scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
    scalar_type lambdan_aug = gmm::neg(lambdan + r * g);
    size_type i = gmm::vect_size(f);
    scalar_type tau = ((i >= 3) ? f[2] : scalar_type(0)) + f[0]*lambdan_aug;
    if (i >= 2) tau = std::min(tau, f[1]);
    scalar_type norm(0);

    if (tau > scalar_type(0)) {
      gmm::add(lambda, gmm::scaled(Vs, -r), F);
      scalar_type mu = gmm::vect_sp(F, n)/nn;
      gmm::add(gmm::scaled(n, -mu/nn), F);
      norm = gmm::vect_norm2(F);
      gmm::copy(gmm::identity_matrix(), dn);
      gmm::scale(dn, -mu/nn);
      gmm::rank_one_update(dn, gmm::scaled(n, mu/(nn*nn*nn)), n);
      gmm::rank_one_update(dn, gmm::scaled(n, scalar_type(-1)/(nn*nn)), F);
      gmm::copy(gmm::identity_matrix(), dVs);
      gmm::rank_one_update(dVs, n, gmm::scaled(n, scalar_type(-1)/(nn*nn)));

      if (norm > tau) {
        gmm::rank_one_update(dVs, F,
                             gmm::scaled(F, scalar_type(-1)/(norm*norm)));
        gmm::scale(dVs, tau / norm);
        gmm::copy(gmm::scaled(F, scalar_type(1)/norm), dg);
        gmm::rank_one_update(dn, gmm::scaled(F, mu/(norm*norm*nn)), F);
        gmm::scale(dn, tau / norm);
        gmm::scale(F, tau / norm);
      } else gmm::clear(dg);

    } else { gmm::clear(dg); gmm::clear(dVs); gmm::clear(F); gmm::clear(dn); }
    // At this stage, F = P_{B_T}, dVs = d_v P_{B_T}, dn = d_n P_{B_T}
    // and dg = d_tau P_{B_T}.

    gmm::copy(dVs, dlambda);
    if (norm > tau && ((i <= 1) || tau < f[1]) && ((i <= 2) || tau > f[2])) {
      gmm::rank_one_update(dn, dg, gmm::scaled(lambda, -f[0]/nn));
      gmm::rank_one_update(dn, dg, gmm::scaled(n, f[0]*lambdan/(nn*nn)));
      gmm::rank_one_update(dlambda, dg, gmm::scaled(n, -f[0]/nn));
      gmm::scale(dg, -f[0]*r);
    } else gmm::clear(dg);
    if (lambdan_aug > scalar_type(0)) {
      gmm::add(gmm::scaled(n, r/nn), dg);
      gmm::rank_one_update(dlambda, n, gmm::scaled(n, scalar_type(1)/(nn*nn)));
      gmm::rank_one_update(dn, gmm::scaled(n, scalar_type(1)/(nn*nn)), lambda);
      gmm::rank_one_update(dn,
                           gmm::scaled(n,(lambdan_aug-lambdan)/(nn*nn*nn)), n);
      for (size_type j = 0; j < N; ++j) dn(j,j) -= lambdan_aug/nn;
    }
    gmm::add(gmm::scaled(n, -lambdan_aug/nn), F);

    gmm::scale(dVs, -r);
  }

#elif FRICTION_LAW == 2 // Contact only

  template <typename VEC, typename VEC2, typename VECR>
  void aug_friction(const VEC &lambda, scalar_type g, const VEC &,
                    const VEC &n, scalar_type r, const VEC2 &, VECR &F) {
    scalar_type nn = gmm::vect_norm2(n);
    scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
    scalar_type lambdan_aug = gmm::neg(lambdan + r * g);
    gmm::copy(gmm::scaled(n, -lambdan_aug/nn), F);
  }

  template <typename VEC, typename VEC2, typename VECR, typename MAT>
  void aug_friction_grad(const VEC &lambda, scalar_type g, const VEC &,
                         const VEC &n, scalar_type r, const VEC2 &, VECR &F,
                         MAT &dlambda, VECR &dg, MAT &dn, MAT &dVs) {
    size_type N = gmm::vect_size(lambda);
    scalar_type nn = gmm::vect_norm2(n);
    scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
    scalar_type lambdan_aug = gmm::neg(lambdan + r * g);

    gmm::clear(dg); gmm::clear(dVs); gmm::clear(F);
    gmm::clear(dn); gmm::clear(dlambda);
    // At this stage, F = P_{B_T}, dVs = d_v P_{B_T}, dn = d_n P_{B_T}
    // and dg = d_tau P_{B_T}.

    if (lambdan_aug > scalar_type(0)) {
      gmm::add(gmm::scaled(n, r/nn), dg);
      gmm::rank_one_update(dlambda, n, gmm::scaled(n, scalar_type(1)/(nn*nn)));
      gmm::rank_one_update(dn, gmm::scaled(n, scalar_type(1)/(nn*nn)), lambda);
      gmm::rank_one_update(dn,
                           gmm::scaled(n,(lambdan_aug-lambdan)/(nn*nn*nn)), n);
      for (size_type j = 0; j < N; ++j) dn(j,j) -= lambdan_aug/nn;
    }
    gmm::add(gmm::scaled(n, -lambdan_aug/nn), F);

    gmm::scale(dVs, -r);
  }



#elif FRICTION_LAW == 3 // Dummy law for test

  template <typename VEC, typename VEC2, typename VECR>
  void aug_friction(const VEC &lambda, scalar_type g, const VEC &Vs,
                    const VEC &n, scalar_type r, const VEC2 &f, VECR &F) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]), F); // dummy
    gmm::copy(gmm::scaled(Vs, g*r*f[0]), F);     // dummy

    gmm::copy(n, F);
  }

  template <typename VEC, typename VEC2, typename VECR, typename MAT>
  void aug_friction_grad(const VEC &lambda, scalar_type g, const VEC &Vs,
                         const VEC &n, scalar_type r, const VEC2 &f, VECR &F,
                         MAT &dlambda, VECR &dg, MAT &dn, MAT &dVs) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]), F); // dummy
    gmm::copy(gmm::scaled(Vs, g*r*f[0]), F);     // dummy

    gmm::copy(n, F);
    gmm::clear(dlambda);
    gmm::clear(dg);
    gmm::clear(dVs);
    gmm::copy(gmm::identity_matrix(), dn);
  }

#elif FRICTION_LAW == 4 // Dummy law for test

  template <typename VEC, typename VEC2, typename VECR>
  void aug_friction(const VEC &lambda, scalar_type g, const VEC &Vs,
                    const VEC &n, scalar_type r, const VEC2 &f, VECR &F) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]*n[0]*Vs[0]), F); // dummy
    gmm::copy(lambda, F);
  }

  template <typename VEC, typename VEC2, typename VECR, typename MAT>
  void aug_friction_grad(const VEC &lambda, scalar_type g, const VEC &Vs,
                         const VEC &n, scalar_type r, const VEC2 &f, VECR &F,
                         MAT &dlambda, VECR &dg, MAT &dn, MAT &dVs) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]*n[0]*Vs[0]), F); // dummy
    gmm::clear(dn);
    gmm::clear(dg);
    gmm::clear(dVs);
    gmm::copy(lambda, F);
    gmm::copy(gmm::identity_matrix(), dlambda);
  }

#elif FRICTION_LAW == 5 // Dummy law for test

  template <typename VEC, typename VEC2, typename VECR>
  void aug_friction(const VEC &lambda, scalar_type g, const VEC &Vs,
                    const VEC &n, scalar_type r, const VEC2 &f, VECR &F) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]*n[0]*Vs[0]), F); // dummy
    gmm::clear(F); F[0] = g;
  }

  template <typename VEC, typename VEC2, typename VECR, typename MAT>
  void aug_friction_grad(const VEC &lambda, scalar_type g, const VEC &Vs,
                         const VEC &n, scalar_type r, const VEC2 &f, VECR &F,
                         MAT &dlambda, VECR &dg, MAT &dn, MAT &dVs) {
    gmm::copy(gmm::scaled(lambda, g*r*f[0]*n[0]*Vs[0]), F); // dummy
    gmm::clear(dlambda);
    gmm::clear(dn);
    gmm::clear(dg);
    gmm::clear(dVs);
    gmm::clear(F); F[0] = g;
    dg[0] = 1.;
  }

#endif




  //=========================================================================
  // Some basic assembly functions
  //=========================================================================


  template <typename MAT1, typename MAT2>
  void mat_elem_assembly(const MAT1 &M_, const gmm::sub_interval &I1,
                         const gmm::sub_interval &I2,
                         const MAT2 &Melem,
                         const mesh_fem &mf1, size_type cv1,
                         const mesh_fem &mf2, size_type cv2) {
    MAT1 &M = const_cast<MAT1 &>(M_);
    typedef typename gmm::linalg_traits<MAT1>::value_type T;
    T val;

    mesh_fem::ind_dof_ct cvdof1 = mf1.ind_basic_dof_of_element(cv1);
    mesh_fem::ind_dof_ct cvdof2 = mf2.ind_basic_dof_of_element(cv2);

    GMM_ASSERT1(cvdof1.size() == gmm::mat_nrows(Melem)
                && cvdof2.size() == gmm::mat_ncols(Melem),
                "Dimensions mismatch");

    if (mf1.is_reduced()) {
      if (mf2.is_reduced()) {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (gmm::sub_matrix(M, I1, I2),
                 gmm::mat_row(mf1.extension_matrix(), cvdof1[i]),
                 gmm::mat_row(mf2.extension_matrix(), cvdof2[j]), val);
      } else {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (gmm::sub_matrix(M, I1, I2),
                 gmm::mat_row(mf1.extension_matrix(), cvdof1[i]),
                 cvdof2[j], val);
      }
    } else {
      if (mf2.is_reduced()) {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (gmm::sub_matrix(M, I1, I2), cvdof1[i],
                 gmm::mat_row(mf2.extension_matrix(), cvdof2[j]), val);
      } else {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              M(cvdof1[i]+I1.first(), cvdof2[j]+I2.first()) += val;
      }
    }
  }

  template <typename VEC1, typename VEC2>
  void vec_elem_assembly(const VEC1 &V_, const gmm::sub_interval &I,
                         const VEC2 &Velem, const mesh_fem &mf, size_type cv) {
    VEC1 &V = const_cast<VEC1 &>(V_);
    typedef typename gmm::linalg_traits<VEC1>::value_type T;
    std::vector<size_type> cvdof(mf.ind_basic_dof_of_element(cv).begin(),
                                 mf.ind_basic_dof_of_element(cv).end());

    GMM_ASSERT1(cvdof.size() == gmm::vect_size(Velem), "Dimensions mismatch");

    if (mf.is_reduced()) {
      T val;
      for (size_type i = 0; i < cvdof.size(); ++i)
        if ((val = Velem[i]) != T(0))
          gmm::add(gmm::scaled(gmm::mat_row(mf.extension_matrix(), cvdof[i]),
                               val), gmm::sub_vector(V, I));
    } else {
      for (size_type i = 0; i < cvdof.size(); ++i)
        V[I.first()+cvdof[i]] += Velem[i];
    }
  }

  //=========================================================================
  //
  //  Large sliding brick. Work in progress
  //
  //=========================================================================

  // For the moment, with raytrace detection and integral unsymmetric
  // Alart-Curnier augmented Lagrangian


  struct integral_large_sliding_contact_brick : public virtual_brick {

    multi_contact_frame &mcf;
    bool with_friction;


    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const;

    integral_large_sliding_contact_brick(multi_contact_frame &mcff,
                                         bool with_fric)
      : mcf(mcff), with_friction(with_fric) {
      set_flags("Integral large sliding contact brick",
                false /* is linear*/, false /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

  };

  static void vectorize_base_tensor(const base_tensor &t, base_matrix &vt,
                             size_type ndof, size_type qdim, size_type N) {
    gmm::resize(vt, ndof*N/qdim, N);
    GMM_ASSERT1(qdim == N || qdim == 1, "mixed intrinsic vector and "
                  "tensorised fem is not supported");
    if (qdim == 1) {
      gmm::clear(vt);
      base_tensor::const_iterator it = t.begin();
      for (size_type i = 0; i < ndof; ++i, ++it)
        for (size_type j = 0; j < N; ++j) vt(i*N+j, j) = *it;
    } else if (qdim == N) {
      gmm::copy(t.as_vector(), vt.as_vector());
    }
  }

  static void vectorize_grad_base_tensor(const base_tensor &t, base_tensor &vt,
                                         size_type ndof, size_type qdim,
                                         size_type N) {
    vt.adjust_sizes(bgeot::multi_index(ndof*N/qdim, N, N));
    GMM_ASSERT1(qdim == N || qdim == 1, "mixed intrinsic vector and "
                  "tensorised fem is not supported");
    if (qdim == 1) {
      gmm::clear(vt.as_vector());
      base_tensor::const_iterator it = t.begin();
      for (size_type k = 0; k < N; ++k)
        for (size_type i = 0; i < ndof; ++i, ++it)
          for (size_type j = 0; j < N; ++j) vt(i*N+j, j, k) = *it;
    } else if (qdim == N) {
      gmm::copy(t.as_vector(), vt.as_vector());
    }
  }


  void integral_large_sliding_contact_brick::asm_real_tangent_terms
  (const model &md, size_type /* ib */, const model::varnamelist &vl,
   const model::varnamelist &dl, const model::mimlist &/* mims */,
   model::real_matlist &matl, model::real_veclist &vecl,
   model::real_veclist &, size_type /* region */,
   build_version version) const {

    // Data : r, friction_coeff.
    GMM_ASSERT1((with_friction && dl.size() >= 2 && dl.size() <= 3)
                || (!with_friction && dl.size() >= 1 && dl.size() <= 2),
            "Wrong number of data for integral large sliding contact brick");

    GMM_ASSERT1(vl.size() == mcf.nb_variables() + mcf.nb_multipliers(),
                "For the moment, it is not allowed to add boundaries to "
                "the multi contact frame object after the model brics has "
                "been added.");

    const model_real_plain_vector &vr = md.real_variable(dl[0]);
    GMM_ASSERT1(gmm::vect_size(vr) == 1, "Large sliding contact "
                    "brick: parameter r should be a scalar");
    scalar_type r = vr[0];

    model_real_plain_vector f_coeff;
    if (with_friction) {
      f_coeff = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(f_coeff) <= 3,
                  "Large sliding contact "
                  "brick: the friction law has less than 3 parameters");
    }
    if (gmm::vect_size(f_coeff) == 0) // default: no friction
      { f_coeff.resize(1); f_coeff[0] = scalar_type(0); }

    scalar_type alpha(0);
    size_type ind = with_friction ? 2:1;
    if (dl.size() >= ind+1) {
      GMM_ASSERT1(md.real_variable(dl[ind]).size() == 1,
                  "Large sliding contact "
                  "brick: parameter alpha should be a scalar");
      alpha = md.real_variable(dl[ind])[0];
    }

    GMM_ASSERT1(matl.size() == 1,
                "Large sliding contact brick should have only one term");
    model_real_sparse_matrix &M = matl[0]; gmm::clear(M);
    model_real_plain_vector &V = vecl[0]; gmm::clear(V);

    mcf.set_raytrace(true);
    mcf.set_fem_nodes_mode(0);
    mcf.compute_contact_pairs();
    size_type N = mcf.dim();
    fem_precomp_pool fppool;
    base_vector coeff;
    base_matrix G, Melem, Inx(N, N), Iny(N, N), Inxy(N, N);
    base_matrix dlambdaF(N, N), dnF(N, N), dVsF(N, N);
    base_matrix grad_phix(N, N), grad_phix_inv(N, N);
    base_matrix grad_phiy(N, N), grad_phiy_inv(N, N);
    base_matrix grad_phi_ny(N, N);

    base_tensor base_lx, base_ux, base_uy; // , base_ly;
    base_tensor grad_base_ux, grad_base_uy;
    base_matrix vbase_lx, vbase_ux, vbase_uy;  // , vbase_ly;
    base_tensor vgrad_base_ux, vgrad_base_uy;
    base_small_vector lambda(N), Vs(N), F(N), dgF(N), wx(N), wy(N);
    base_small_vector aux8(N), aux9(N);

    // Stabilization for non-contact zones
    for (size_type i = 0; i < mcf.nb_boundaries(); ++i)
      if (mcf.is_self_contact() || mcf.is_slave_boundary(i)) {
        size_type region = mcf.region_of_boundary(i);
        const std::string &name_lx = mcf.multname_of_boundary(i);
        GMM_ASSERT1(name_lx.size() > 0, "This brick need "
                    "multipliers defined on the multi_contact_frame object");
        const mesh_fem &mflambda = mcf.mflambda_of_boundary(i);
        const mesh_im &mim = mcf.mim_of_boundary(i);
        const gmm::sub_interval &I = md.interval_of_variable(name_lx);

        if (version & model::BUILD_MATRIX) {
          model_real_sparse_matrix M1(mflambda.nb_dof(), mflambda.nb_dof());
          asm_mass_matrix(M1, mim, mflambda, region);
          gmm::add(gmm::scaled(M1,scalar_type(1)/r), gmm::sub_matrix(M, I, I));
        }

        if (version & model::BUILD_RHS) {
          model_real_plain_vector V1(mflambda.nb_dof());
          asm_source_term
            (V1, mim, mflambda, mflambda,
             md.real_variable(mcf.multname_of_boundary(i)), region);
          gmm::add(gmm::scaled(V1, scalar_type(-1)/r), gmm::sub_vector(V, I));
        }
      }


    // Iterations on the contact pairs
    for (size_type icp = 0; icp < mcf.nb_contact_pairs(); ++icp) {
      const multi_contact_frame::contact_pair &cp = mcf.get_contact_pair(icp);

      // Extract information on slave and master points and fems
      // TODO : Could be done by a structure, on demand ... (certain
      // computations are not necessary for rhs only)
      const base_node &x = cp.slave_point;

      bool test_fixed_normal = false;
      base_small_vector nx(N);
      if (test_fixed_normal) {
        nx[0] = 0; nx[1] = -1.;
      } else {
        nx = cp.slave_n;
      }

      gmm::copy(gmm::identity_matrix(), Inx);
      gmm::rank_one_update(Inx, nx, gmm::scaled(nx, scalar_type(-1)));
      size_type irigid = cp.irigid_obstacle;
      bool isrigid = (irigid != size_type(-1));
      size_type ibx = cp.slave_ind_boundary;
      const std::string &name_ux = mcf.varname_of_boundary(ibx);
      const std::string &name_lx = mcf.multname_of_boundary(ibx);
      GMM_ASSERT1(name_lx.size(), "This brick need multipliers defined "
                  "on the multi_contact_frame object");
      const model_real_plain_vector &all_wx = mcf.w_of_boundary(ibx);
      size_type cvx = cp.slave_ind_element;
      size_type fx = cp.slave_ind_face;
      const mesh_fem *mf_ux = &(mcf.mfu_of_boundary(ibx));
      const gmm::sub_interval &I_ux = md.interval_of_variable(name_ux);
      pfem pf_ux = mf_ux->fem_of_element(cvx);
      size_type ndof_ux = pf_ux->nb_dof(cvx);
      size_type qdim_ux = pf_ux->target_dim();
      const mesh_fem *mf_lx = &(mcf.mflambda_of_boundary(ibx));
      const gmm::sub_interval &I_lx = md.interval_of_variable(name_lx);
      pfem pf_lx = mf_lx->fem_of_element(cvx);
      size_type ndof_lx = pf_lx->nb_dof(cvx);
      size_type qdim_lx = pf_lx->target_dim();
      const mesh_im &mim = mcf.mim_of_boundary(ibx);
      pintegration_method pim = mim.int_method_of_element(cvx);
      size_type ii = cp.slave_ind_pt;
      scalar_type weight = pim->approx_method()->coeff(ii);
      const mesh &mx = mim.linked_mesh();
      bgeot::pgeometric_trans pgtx = mx.trans_of_convex(cvx);

      const base_node &y = cp.master_point;
      const base_small_vector &ny = cp.master_n;
      scalar_type nxny = gmm::vect_sp(nx, ny);
      gmm::copy(gmm::identity_matrix(), Iny);
      gmm::rank_one_update(Iny, ny, gmm::scaled(ny, scalar_type(-1)));
      gmm::copy(gmm::identity_matrix(), Inxy);
      gmm::rank_one_update(Inxy, nx, gmm::scaled(ny,scalar_type(-1)/nxny));
      size_type iby = isrigid ? 0 : cp.master_ind_boundary;
      const std::string &name_uy
        = isrigid ? name_ux : mcf.varname_of_boundary(iby);
      const model_real_plain_vector &all_wy
        = isrigid ? all_wx:mcf.w_of_boundary(iby);
      // const std::string &name_ly
      //   = isrigid ? name_ux : mcf.multname_of_boundary(iby);
      const gmm::sub_interval &I_uy
        = isrigid ? I_ux : md.interval_of_variable(name_uy);
      // const gmm::sub_interval &I_ly
      //   = isrigid ? I_lx : md.interval_of_variable(name_ly);
      const mesh_fem *mf_uy = isrigid ? mf_ux : &(mcf.mfu_of_boundary(iby));
      const mesh &my = isrigid ? mx : mf_uy->linked_mesh();
      size_type cvy(0), fy(0), ndof_uy(0), qdim_uy(0);
      pfem pf_uy(0);
      bgeot::pgeometric_trans pgty(0);
      if (!isrigid) {
        cvy = cp.master_ind_element;
        fy = cp.master_ind_face;
        pf_uy = mf_uy->fem_of_element(cvy);
        ndof_uy = pf_uy->nb_dof(cvy);
        qdim_uy = pf_uy->target_dim();
        // const mesh_fem *mf_ly = &(mcf.mflambda_of_boundary(iby));
        // pfem pf_ly = mf_ly->fem_of_element(cvy);
        // size_type ndof_ly = pf_ly->nb_dof(cvy);
        // size_type qdim_ly = pf_ly->target_dim();
        pgty = my.trans_of_convex(cvy);
      }
      const base_node &y_ref = cp.master_point_ref;
      scalar_type g = cp.signed_dist;

      // Fem interpolation context objects
      bgeot::vectors_to_base_matrix(G, mx.points_of_convex(cvx));
      pfem_precomp pfp_ux
        = fppool(pf_ux, &(pim->approx_method()->integration_points()));
      pfem_precomp pfp_lx
        = fppool(pf_lx, &(pim->approx_method()->integration_points()));
      fem_interpolation_context ctx_ux(pgtx, pfp_ux, ii, G, cvx, fx);
      fem_interpolation_context ctx_lx(pgtx, pfp_lx, ii, G, cvx, fx);

      fem_interpolation_context ctx_uy;
      // fem_interpolation_context ctx_ly;
      if (!isrigid) {
        bgeot::vectors_to_base_matrix(G, my.points_of_convex(cvy));
        ctx_uy = fem_interpolation_context(pgty, pf_uy, y_ref, G, cvy, fy);
        // ctx_ly = fem_interpolation_context(pgty, pf_ly, y_ref, G, cvy, fy);
      }

      // Base tensors
      ctx_lx.base_value(base_lx);
      ctx_ux.base_value(base_ux);
      // ctx_ly.base_value(base_ly);
      if (!isrigid) ctx_uy.base_value(base_uy);
      ctx_ux.grad_base_value(grad_base_ux);
      if (!isrigid) ctx_uy.grad_base_value(grad_base_uy);

      // Vectorisation of tensor (to facilitate the treatment)
      vectorize_base_tensor(base_lx, vbase_lx, ndof_lx, qdim_lx, N);
      vectorize_base_tensor(base_ux, vbase_ux, ndof_ux, qdim_ux, N);
      // vectorize_base_tensor(base_ly, vbase_ly, ndof_ly, qdim_ly, N);
      if (!isrigid)
        vectorize_base_tensor(base_uy, vbase_uy, ndof_uy, qdim_uy, N);
      vectorize_grad_base_tensor(grad_base_ux, vgrad_base_ux, ndof_ux,
                                 qdim_ux, N);
      if (!isrigid)
        vectorize_grad_base_tensor(grad_base_uy, vgrad_base_uy, ndof_uy,
                                   qdim_uy, N);
      ndof_lx = ndof_lx*N/qdim_lx;
      ndof_ux = ndof_ux*N/qdim_ux;
      // ndof_ly = ndof_ly*N/qdim_ly;
      if (!isrigid) ndof_uy = ndof_uy*N/qdim_uy;

      // Value of lambda on the slave surface
      slice_vector_on_basic_dof_of_element(*mf_lx, mcf.mult_of_boundary(ibx),
                                           cvx, coeff);
      ctx_lx.pf()->interpolation(ctx_lx, coeff, lambda, dim_type(N));

      // Value of grad phi on the slave surface
      slice_vector_on_basic_dof_of_element(*mf_ux, mcf.disp_of_boundary(ibx),
                                           cvx, coeff);
      ctx_ux.pf()->interpolation_grad(ctx_ux, coeff, grad_phix, dim_type(N));
      gmm::add(gmm::identity_matrix(), grad_phix);
      gmm::copy(grad_phix, grad_phix_inv);
      scalar_type J = gmm::lu_inverse(grad_phix_inv);
      if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);

      // Value of grad phi on the master surface
      if (!isrigid) {
        slice_vector_on_basic_dof_of_element(*mf_uy, mcf.disp_of_boundary(iby),
                                             cvy, coeff);
        ctx_uy.pf()->interpolation_grad(ctx_uy, coeff, grad_phiy, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad_phiy);
        gmm::copy(grad_phiy, grad_phiy_inv);
        J = gmm::lu_inverse(grad_phiy_inv);
        if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
      }

      if (alpha != scalar_type(0)) {
        // Value of w on the slave surface
        if (all_wx.size()) {
          slice_vector_on_basic_dof_of_element(*mf_ux, all_wx, cvx, coeff);
          ctx_ux.pf()->interpolation(ctx_ux, coeff, wx, dim_type(N));
        }  else gmm::clear(wx);
        gmm::add(ctx_ux.xreal(), wx);

        // Value of w on the master surface
        if (!isrigid) {
          if (all_wy.size()) {
            slice_vector_on_basic_dof_of_element(*mf_uy, all_wy, cvy, coeff);
            ctx_uy.pf()->interpolation(ctx_uy, coeff, wy, dim_type(N));
          } else gmm::clear(wy);
          gmm::add(ctx_uy.xreal(), wy);
        } else gmm::copy(y, wy);

        // Value of grad phi(n-1) on the master surface
        if (!isrigid && all_wy.size()) {
          ctx_uy.pf()->interpolation_grad(ctx_uy, coeff, grad_phi_ny,
                                          dim_type(N));
          gmm::add(gmm::identity_matrix(), grad_phi_ny);
        } else gmm::copy(gmm::identity_matrix(), grad_phi_ny);

        // Value of Vs (sliding velocity)
        gmm::add(x, gmm::scaled(y, scalar_type(-1)), Vs);
        gmm::add(gmm::scaled(wx, scalar_type(-1)), Vs);
        gmm::add(wy, Vs);
        gmm::scale(Vs, alpha);
      } else gmm::clear(Vs);

      base_vector aux6(ndof_uy), aux7(ndof_ux), aux12(ndof_lx);


      if (version & model::BUILD_MATRIX) {

// #define DO_TEST_F

#ifdef DO_TEST_F

        scalar_type EPS = 5E-9;
        for (size_type k = 0; k < 100; ++k) {
          base_small_vector lambda_r(N), Vs_r(N), nx_r(N), f_coeff_r(3);
          base_small_vector F2(N), F3(N);
          scalar_type g_r = gmm::random(1.), r_r = gmm::random();
          gmm::fill_random(lambda_r);
          gmm::fill_random(Vs_r);
          gmm::fill_random(nx_r);
          gmm::scale(nx_r, 1./gmm::vect_norm2(nx_r));
          f_coeff_r[0] = gmm::random();
          f_coeff_r[1] = gmm::random();
          f_coeff_r[2] = gmm::random();

          // cout << "lambda_r = " << lambda_r << " Vs_r = " << Vs_r
          //      << " nx_r = " << nx_r << endl << "g_r = " << g_r
          //      << " r_r = " << r_r << " f = " << f_coeff_r << endl;

          aug_friction(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F);
          aug_friction_grad(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F2,
                            dlambdaF, dgF, dnF, dVsF);
          GMM_ASSERT1(gmm::vect_dist2(F2, F) < 1E-7, "bad F");

          base_small_vector dlambda(N);
          gmm::fill_random(dlambda);


          gmm::add(gmm::scaled(dlambda, EPS), nx_r);
          aug_friction(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F2);

          gmm::mult(dnF, gmm::scaled(dlambda, EPS), F, F3);
          if (gmm::vect_dist2(F2, F3)/EPS > 1E-4) {
            cout << "lambda_r = " << lambda_r << " Vs_r = " << Vs_r
                 << " nx_r = " << nx_r << endl << "g_r = " << g_r
                 << " r_r = " << r_r << " f = " << f_coeff_r << endl;
            cout << "diff = " << gmm::vect_dist2(F2, F3)/EPS << endl;
            GMM_ASSERT1(false, "bad n derivative");
          }

          gmm::add(gmm::scaled(dlambda, -EPS), nx_r);


          gmm::add(gmm::scaled(dlambda, EPS), lambda_r);
          aug_friction(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F2);
          gmm::mult(dlambdaF, gmm::scaled(dlambda, EPS), F, F3);
          if (gmm::vect_dist2(F2, F3)/EPS > 1E-6) {
            cout << "diff = " << gmm::vect_dist2(F2, F3)/EPS << endl;
            GMM_ASSERT1(false, "bad lambda derivative");
          }
          gmm::add(gmm::scaled(dlambda, -EPS), lambda_r);


          gmm::add(gmm::scaled(dlambda, EPS), Vs_r);
          aug_friction(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F2);
          gmm::mult(dVsF, gmm::scaled(dlambda, EPS), F, F3);
          if (gmm::vect_dist2(F2, F3)/EPS > 1E-6) {
            cout << "diff = " << gmm::vect_dist2(F2, F3)/EPS << endl;
            GMM_ASSERT1(false, "bad Vs derivative");
          }
          gmm::add(gmm::scaled(dlambda, -EPS), Vs_r);


          g_r += EPS;
          aug_friction(lambda_r, g_r, Vs_r, nx_r, r_r, f_coeff_r, F2);
          gmm::add(gmm::scaled(dgF, EPS), F, F3);
          if (gmm::vect_dist2(F2, F3)/EPS > 1E-6) {
            cout << "diff = " << gmm::vect_dist2(F2, F3)/EPS << endl;
            GMM_ASSERT1(false, "bad g derivative");
          }
          g_r -= EPS;



        }
#endif

        // cout << "pair " << icp << " g = " << g << endl;
        base_matrix aux1(ndof_uy, N), aux4(ndof_uy, ndof_ux);
        base_matrix aux5(ndof_lx, N), aux10(ndof_lx, N);
        base_matrix aux11(ndof_lx, ndof_ux);
        base_matrix aux2(N, N), aux3(N, N);

        aug_friction_grad(lambda, g, Vs, nx, r, f_coeff, F, dlambdaF,
                          dgF, dnF, dVsF);


        base_matrix graddeltaunx(ndof_ux, N);
        for (size_type i = 0; i < ndof_ux; ++i)
          for (size_type j = 0; j < N; ++j)
            for (size_type k = 0; k < N; ++k)
              graddeltaunx(i, j) += nx[k] * vgrad_base_ux(i, k, j);

#define CONSIDER_TERM1
#define CONSIDER_TERM2
#define CONSIDER_TERM3


#ifdef CONSIDER_TERM1
        // Term  -\delta\lambda(X) . \delta v(X) 
        gmm::resize(Melem, ndof_ux, ndof_lx); gmm::clear(Melem);
        gmm::mult(vbase_ux, gmm::transposed(vbase_lx), Melem);
        gmm::scale(Melem, -weight);
        mat_elem_assembly(M, I_ux, I_lx, Melem, *mf_ux, cvx, *mf_lx, cvx);
#endif


#ifdef CONSIDER_TERM2

        if (!isrigid) {
          // Term  \delta\lambda(X) . \delta v(Y)
          gmm::resize(Melem, ndof_uy, ndof_lx); gmm::clear(Melem);
          gmm::mult(vbase_uy, gmm::transposed(vbase_lx), Melem);
          gmm::scale(Melem, weight);
          mat_elem_assembly(M, I_uy, I_lx, Melem, *mf_uy, cvy, *mf_lx, cvx);

          // Term -\lambda(X) . (\nabla \delta v(Y) (\nabla phi)^(-1)\delta y
          gmm::clear(aux1);
          for (size_type i = 0; i < ndof_uy; ++i)
            for (size_type j = 0; j < N; ++j)
              for (size_type k = 0; k < N; ++k)
                aux1(i, j) += lambda[k] * vgrad_base_uy(i, k, j);
          base_matrix lgraddeltavgradphiyinv(ndof_uy, N);
          gmm::mult(aux1, grad_phiy_inv, lgraddeltavgradphiyinv);

          // first sub term
          gmm::resize(Melem, ndof_uy, ndof_uy); gmm::clear(Melem);
          gmm::clear(aux2); 
          gmm::rank_one_update(aux2, nx, gmm::scaled(ny,scalar_type(1)/nxny));
          gmm::mult(lgraddeltavgradphiyinv, aux2, aux1);
          gmm::mult(aux1, gmm::transposed(vbase_uy), Melem);
          gmm::scale(Melem, weight);
          mat_elem_assembly(M, I_uy, I_uy, Melem, *mf_uy, cvy, *mf_uy, cvy);

          // Second sub term
          gmm::resize(Melem, ndof_uy, ndof_ux); gmm::clear(Melem);
          gmm::mult(lgraddeltavgradphiyinv, Inxy, aux1);
          gmm::mult(aux1, gmm::transposed(vbase_ux), Melem);

          // Third sub term
          if (!test_fixed_normal) {
            gmm::mult(Iny, Inx, aux2);
            gmm::mult(aux2, gmm::transposed(grad_phix_inv), aux3);
            gmm::mult(lgraddeltavgradphiyinv, aux3, aux1);
            gmm::mult(aux1, gmm::transposed(graddeltaunx), aux4);
            gmm::scale(aux4, -g/(nxny*nxny));
            gmm::add(aux4, Melem);
          }
          gmm::scale(Melem, weight);
          mat_elem_assembly(M, I_uy, I_ux, Melem, *mf_uy, cvy, *mf_ux, cvx);
        }

#endif

#ifdef CONSIDER_TERM3

        // cout << "dlambdaF = " << dlambdaF << endl;
        // cout << "lambda = " << lambda << endl;

        // Term (1/r)(I-dlambdaF)\delta\lambda\delta\mu
        //   the I of (I-dlambdaF) is skipped because globally added before
        gmm::resize(Melem, ndof_lx, ndof_lx); gmm::clear(Melem);
        gmm::copy(gmm::scaled(dlambdaF, scalar_type(-1)/r), aux2);
        gmm::mult(vbase_lx, aux2, aux5);
        gmm::mult(aux5, gmm::transposed(vbase_lx), Melem);
        gmm::scale(Melem, weight);
        // cout << "Melem1 " << Melem << endl;
        mat_elem_assembly(M, I_lx, I_lx, Melem, *mf_lx, cvx, *mf_lx, cvx);

        // cout << "dnF = " << dnF << endl;
        // cout << "nx = " << nx << endl;
        // cout << "Inx = " << Inx << endl;

        // Term -(1/r)dnF\delta nx\delta\mu
        gmm::resize(Melem, ndof_lx, ndof_ux); gmm::clear(Melem);
        if (!test_fixed_normal) {
          gmm::mult(vbase_lx, dnF, aux5);
          gmm::mult(aux5, Inx, aux10);
          gmm::mult(aux10,  gmm::transposed(grad_phix_inv), aux5);
          gmm::mult(aux5, gmm::transposed(graddeltaunx), Melem);
          gmm::scale(Melem, scalar_type(1)/r);
          // assembly factorized with the next term
        }

        // cout << "dgF = " << dgF << endl;
        // Term -(1/r)dgF\delta g\delta\mu
        base_vector deltamudgF(ndof_lx);
        gmm::mult(vbase_lx, gmm::scaled(dgF, scalar_type(1)/(r*nxny)),
                  deltamudgF);
        // cout << "deltamudgF = " << deltamudgF << endl;

        // first sub term
        gmm::mult(vbase_ux, ny, aux7);

        // second sub term
        if (!test_fixed_normal) {
          gmm::mult(Inx, gmm::scaled(ny, -g), aux8);
          gmm::mult(grad_phix_inv, aux8, aux9);
          gmm::mult_add(graddeltaunx, aux9, aux7);
          gmm::rank_one_update(Melem, deltamudgF,  aux7);
        }
        gmm::scale(Melem, weight);
        mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);

        if (!isrigid) {
          // third sub term
          gmm::resize(Melem, ndof_lx, ndof_uy); gmm::clear(Melem);
          gmm::mult(vbase_uy, ny, aux6);
          gmm::rank_one_update(Melem, deltamudgF,  aux6);
          gmm::scale(Melem, -weight);
          mat_elem_assembly(M, I_lx, I_uy, Melem, *mf_lx, cvx, *mf_uy, cvy);
        }

        if (alpha != scalar_type(0)) {
          // Term -(1/r) d_Vs F \delta Vs\delta\mu
          base_matrix I_gphingphiyinv(N, N);
          if (!isrigid)
            gmm::mult(gmm::scaled(grad_phi_ny, scalar_type(-1)),
                      grad_phiy_inv, I_gphingphiyinv);
          gmm::add(gmm::identity_matrix(), I_gphingphiyinv);

          // first sub term
          gmm::resize(Melem, ndof_lx, ndof_ux); gmm::clear(Melem);
          gmm::mult(I_gphingphiyinv, Inxy, aux2);
          for (size_type j = 0; j < N; ++j) aux2(j,j) -= scalar_type(1);
          gmm::mult(dVsF, aux2, aux3);
          gmm::mult(vbase_lx, gmm::transposed(aux3), aux5);
          gmm::mult(aux5, gmm::transposed(vbase_ux), Melem);

          // second sub term
          if (!test_fixed_normal) {
            gmm::mult(dVsF, I_gphingphiyinv, aux2);
            gmm::mult(aux2, Iny, aux3);
            gmm::mult(aux3, gmm::scaled(Inx, scalar_type(-1)), aux2);
            gmm::mult(aux2, gmm::transposed(grad_phix_inv), aux3);
            gmm::mult(vbase_lx, gmm::transposed(aux3), aux5);
            gmm::mult(aux5, gmm::transposed(graddeltaunx), aux11);
            gmm::scale(aux11, g/(nxny*nxny));
            gmm::add(aux11, Melem);
          }
          gmm::scale(Melem, weight*alpha/r);
          mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);

          if (!isrigid) {
            // third sub term
            gmm::resize(Melem, ndof_lx, ndof_uy); gmm::clear(Melem);
            gmm::mult(dVsF, I_gphingphiyinv, aux2);
            gmm::mult(aux2, nx, aux8);
            gmm::mult(vbase_lx, aux8, aux12);
            gmm::mult(vbase_uy, gmm::scaled(ny, scalar_type(1)/nxny), aux6);
            gmm::rank_one_update(Melem, aux12, aux6);
            gmm::scale(Melem, weight*alpha/r);
            mat_elem_assembly(M, I_lx, I_uy, Melem, *mf_lx, cvx, *mf_uy, cvy);
          }
        }

#endif

      }

      if (version & model::BUILD_RHS) {

        if (!(version & model::BUILD_MATRIX))
          aug_friction(lambda, g, Vs, nx, r, f_coeff, F);

#ifdef CONSIDER_TERM1

        // Term lambda.\delta v(X)
        gmm::mult(vbase_ux, lambda, aux7);
        gmm::scale(aux7, weight);
        vec_elem_assembly(V, I_ux, aux7, *mf_ux, cvx);
#endif

#ifdef CONSIDER_TERM2

        // Term -lambda.\delta v(Y)
        if (!isrigid) {
          gmm::mult(vbase_uy, lambda, aux6);
          gmm::scale(aux6, -weight);
          vec_elem_assembly(V, I_uy, aux6, *mf_uy, cvy);
        }
#endif

#ifdef CONSIDER_TERM3

        // Term -(1/r)(lambda - F).\delta \mu
        // (1/r)(lambda).\delta \mu is skipped because globally added before
        gmm::mult(vbase_lx, gmm::scaled(F, scalar_type(-1)/r), aux12);
        gmm::scale(aux12, -weight);
        vec_elem_assembly(V, I_lx, aux12, *mf_lx, cvx);
#endif
      }

    }
  }



  size_type add_integral_large_sliding_contact_brick
  (model &md, multi_contact_frame &mcf,
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   const std::string &dataname_alpha) {

    bool with_friction = (dataname_friction_coeff.size() > 0);
    integral_large_sliding_contact_brick *pbr
      = new integral_large_sliding_contact_brick(mcf, with_friction);

    model::termlist tl; // A unique global unsymmetric term
    tl.push_back(model::term_description(true, false));

    model::varnamelist dl(1, dataname_r);
    if (with_friction) dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size()) dl.push_back(dataname_alpha);

    model::varnamelist vl;

    bool selfcontact = mcf.is_self_contact();

    dal::bit_vector uvar, mvar;
    for (size_type i = 0; i < mcf.nb_boundaries(); ++i) {
      size_type ind_u = mcf.ind_varname_of_boundary(i);
      if (!(uvar.is_in(ind_u))) {
        vl.push_back(mcf.varname(ind_u));
        uvar.add(ind_u);
      }
      size_type ind_lambda = mcf.ind_multname_of_boundary(i);

      if (selfcontact || mcf.is_slave_boundary(i))
        GMM_ASSERT1(ind_lambda != size_type(-1), "Large sliding contact "
                    "brick: a multiplier should be associated to each slave "
                    "boundary in the multi_contact_frame object.");
      if (ind_lambda != size_type(-1) && !(mvar.is_in(ind_lambda))) {
        vl.push_back(mcf.multname(ind_lambda));
        mvar.add(ind_u);
      }
    }

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }



  //=========================================================================
  //
  //  Large sliding brick with field extension principle. To be adapated with
  //  the new structure for contact pairs.
  //
  //=========================================================================

  //=========================================================================
  // 0)- Some basic assembly functions
  //=========================================================================

  template <typename MAT1, typename MAT2>
  void mat_elem_assembly(const MAT1 &M_, const MAT2 &Melem,
                         const mesh_fem &mf1, size_type cv1,
                         const mesh_fem &mf2, size_type cv2) {
    MAT1 &M = const_cast<MAT1 &>(M_);
    typedef typename gmm::linalg_traits<MAT1>::value_type T;
    T val;
    mesh_fem::ind_dof_ct cvdof1 = mf1.ind_basic_dof_of_element(cv1);
    mesh_fem::ind_dof_ct cvdof2 = mf2.ind_basic_dof_of_element(cv2);

    GMM_ASSERT1(cvdof1.size() == gmm::mat_nrows(Melem)
                && cvdof2.size() == gmm::mat_ncols(Melem),
                "Dimensions mismatch");

    if (mf1.is_reduced()) {
      if (mf2.is_reduced()) {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (M, gmm::mat_row(mf1.extension_matrix(), cvdof1[i]),
                 gmm::mat_row(mf2.extension_matrix(), cvdof2[j]), val);
      } else {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (M, gmm::mat_row(mf1.extension_matrix(), cvdof1[i]),
                 cvdof2[j], val);
      }
    } else {
      if (mf2.is_reduced()) {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              asmrankoneupdate
                (M, cvdof1[i],
                 gmm::mat_row(mf2.extension_matrix(), cvdof2[j]), val);
      } else {
        for (size_type i = 0; i < cvdof1.size(); ++i)
          for (size_type j = 0; j < cvdof2.size(); ++j)
            if ((val = Melem(i,j)) != T(0))
              M(cvdof1[i], cvdof2[j]) += val;
      }
    }
  }


  template <typename VEC1, typename VEC2>
  void vec_elem_assembly(const VEC1 &V_, const VEC2 &Velem,
                         const mesh_fem &mf, size_type cv) {
    VEC1 &V = const_cast<VEC1 &>(V_);
    typedef typename gmm::linalg_traits<VEC1>::value_type T;
    std::vector<size_type> cvdof(mf.ind_basic_dof_of_element(cv).begin(),
                                 mf.ind_basic_dof_of_element(cv).end());

    GMM_ASSERT1(cvdof.size() == gmm::vect_size(Velem), "Dimensions mismatch");

    if (mf.is_reduced()) {
      T val;
      for (size_type i = 0; i < cvdof.size(); ++i)
        if ((val = Velem[i]) != T(0))
          gmm::add(gmm::scaled(gmm::mat_row(mf.extension_matrix(), cvdof[i]),
                               val), V);
    } else {
      for (size_type i = 0; i < cvdof.size(); ++i) V[cvdof[i]] += Velem[i];
    }
  }


  //=========================================================================
  // 1)- Structure which stores the contact boundaries and rigid obstacles
  //=========================================================================

  struct contact_frame {
    bool frictionless;
    size_type N;
    scalar_type friction_coef;
    std::vector<const model_real_plain_vector *> Us;
    std::vector<model_real_plain_vector> ext_Us;
    std::vector<const model_real_plain_vector *> lambdas;
    std::vector<model_real_plain_vector> ext_lambdas;
    struct contact_boundary {
      size_type region;                 // Boundary number
      const getfem::mesh_fem *mfu;      // F.e.m. for the displacement.
      size_type ind_U;                  // Index of displacement.
      const getfem::mesh_fem *mflambda; // F.e.m. for the multiplier.
      size_type ind_lambda;             // Index of multiplier.
    };
    std::vector<contact_boundary> contact_boundaries;

    gmm::dense_matrix< model_real_sparse_matrix * > UU;
    gmm::dense_matrix< model_real_sparse_matrix * > UL;
    gmm::dense_matrix< model_real_sparse_matrix * > LU;
    gmm::dense_matrix< model_real_sparse_matrix * > LL;

    std::vector< model_real_plain_vector *> Urhs;
    std::vector< model_real_plain_vector *> Lrhs;



    std::vector<std::string> coordinates;
    base_node pt_eval;
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
    std::vector<mu::Parser> obstacles_parsers;
#endif
    std::vector<std::string> obstacles;
    std::vector<std::string> obstacles_velocities;

    size_type add_U(const getfem::mesh_fem &mfu,
                    const model_real_plain_vector &U) {
      size_type i = 0;
      for (; i < Us.size(); ++i) if (Us[i] == &U) return i;
      Us.push_back(&U);
      model_real_plain_vector ext_U(mfu.nb_basic_dof()); // means that the structure has to be build each time ... to be changed. ATTENTION : la même variable ne doit pas être étendue dans deux vecteurs différents.
      mfu.extend_vector(U, ext_U);
      ext_Us.push_back(ext_U);
      return i;
    }

    size_type add_lambda(const getfem::mesh_fem &mfl,
                         const model_real_plain_vector &l) {
      size_type i = 0;
      for (; i < lambdas.size(); ++i) if (lambdas[i] == &l) return i;
      lambdas.push_back(&l);
      model_real_plain_vector ext_l(mfl.nb_basic_dof()); // means that the structure has to be build each time ... to be changed. ATTENTION : la même variable ne doit pas être étendue dans deux vecteurs différents.
      mfl.extend_vector(l, ext_l);
      ext_lambdas.push_back(ext_l);
      return i;
    }

    void extend_vectors(void) {
      for (size_type i = 0; i < contact_boundaries.size(); ++i) {
        size_type ind_U = contact_boundaries[i].ind_U;
        contact_boundaries[i].mfu->extend_vector(*(Us[ind_U]), ext_Us[ind_U]);
        size_type ind_lambda = contact_boundaries[i].ind_lambda;
        contact_boundaries[i].mflambda->extend_vector(*(lambdas[ind_lambda]),
                                                      ext_lambdas[ind_lambda]);
      }
    }


    const getfem::mesh_fem &mfu_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mfu); }
    const getfem::mesh_fem &mflambda_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mflambda); }
    const model_real_plain_vector &disp_of_boundary(size_type n) const
    { return ext_Us[contact_boundaries[n].ind_U]; }
    const model_real_plain_vector &lambda_of_boundary(size_type n) const
    { return ext_lambdas[contact_boundaries[n].ind_lambda]; }
    size_type region_of_boundary(size_type n) const
    { return contact_boundaries[n].region; }
    model_real_sparse_matrix &UU_matrix(size_type n, size_type m) const
    { return *(UU(contact_boundaries[n].ind_U, contact_boundaries[m].ind_U)); }
    model_real_sparse_matrix &LU_matrix(size_type n, size_type m) const {
      return *(LU(contact_boundaries[n].ind_lambda,
                  contact_boundaries[m].ind_U));
    }
    model_real_sparse_matrix &UL_matrix(size_type n, size_type m) const {
      return *(UL(contact_boundaries[n].ind_U,
                  contact_boundaries[m].ind_lambda));
    }
    model_real_sparse_matrix &LL_matrix(size_type n, size_type m) const {
      return *(LL(contact_boundaries[n].ind_lambda,
                  contact_boundaries[m].ind_lambda));
    }
    model_real_plain_vector &U_vector(size_type n) const
    { return *(Urhs[contact_boundaries[n].ind_U]); }
    model_real_plain_vector &L_vector(size_type n) const
    { return *(Lrhs[contact_boundaries[n].ind_lambda]); }

    contact_frame(size_type NN) : N(NN), coordinates(N), pt_eval(N) {
      if (N > 0) coordinates[0] = "x";
      if (N > 1) coordinates[1] = "y";
      if (N > 2) coordinates[2] = "z";
      if (N > 3) coordinates[3] = "w";
      GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
                  "dimension greater than 4");
    }

    size_type add_obstacle(const std::string &obs) {
      size_type ind = obstacles.size();
      obstacles.push_back(obs);
      obstacles_velocities.push_back("");
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H

      mu::Parser mu;
      obstacles_parsers.push_back(mu);
      obstacles_parsers[ind].SetExpr(obstacles[ind]);
      for (size_type k = 0; k < N; ++k)
        obstacles_parsers[ind].DefineVar(coordinates[k], &pt_eval[k]);
#else
      GMM_ASSERT1(false, "You have to link muparser with getfem to deal "
                  "with rigid body obstacles");
#endif
      return ind;
    }

    size_type add_boundary(const getfem::mesh_fem &mfu,
                           const model_real_plain_vector &U,
                           const getfem::mesh_fem &mfl,
                           const model_real_plain_vector &l,
                           size_type reg) {
      contact_boundary cb;
      cb.region = reg;
      cb.mfu = &mfu;
      cb.mflambda = &mfl;
      cb.ind_U = add_U(mfu, U);
      cb.ind_lambda = add_lambda(mfl, l);
      size_type ind = contact_boundaries.size();
      contact_boundaries.push_back(cb);
      gmm::resize(UU, ind+1, ind+1);
      gmm::resize(UL, ind+1, ind+1);
      gmm::resize(LU, ind+1, ind+1);
      gmm::resize(LL, ind+1, ind+1);
      gmm::resize(Urhs, ind+1);
      gmm::resize(Lrhs, ind+1);
      return ind;
    }

  };


  //=========================================================================
  // 2)- Structure which computes the contact pairs, rhs and tangent terms
  //=========================================================================

  struct contact_elements {

    contact_frame &cf;   // contact frame description.

    // list des enrichissements pour ses points : y0, d0, element ...
    bgeot::rtree element_boxes;  // influence regions of boundary elements
    // list des enrichissements of boundary elements
    std::vector<size_type> boundary_of_elements;
    std::vector<size_type> ind_of_elements;
    std::vector<size_type> face_of_elements;
    std::vector<base_node> unit_normal_of_elements;

    contact_elements(contact_frame &ccf) : cf(ccf) {}
    void init(void);
    bool add_point_contribution(size_type boundary_num,
                                getfem::fem_interpolation_context &ctxu,
                                getfem::fem_interpolation_context &ctxl,
                                scalar_type weight, scalar_type f_coeff,
                                scalar_type r, model::build_version version);
  };


  void contact_elements::init(void) {
    fem_precomp_pool fppool;
    // compute the influence regions of boundary elements. To be run
    // before the assembly of contact terms.
    element_boxes.clear();
    unit_normal_of_elements.resize(0);
    boundary_of_elements.resize(0);
    ind_of_elements.resize(0);
    face_of_elements.resize(0);

    size_type N = 0;
    base_matrix G;
    model_real_plain_vector coeff;
    cf.extend_vectors();
    for (size_type i = 0; i < cf.contact_boundaries.size(); ++i) {
      size_type bnum = cf.region_of_boundary(i);
      const mesh_fem &mfu = cf.mfu_of_boundary(i);
      const model_real_plain_vector &U = cf.disp_of_boundary(i);
      const mesh &m = mfu.linked_mesh();
      if (i == 0) N = m.dim();
      GMM_ASSERT1(m.dim() == N,
                  "Meshes are of mixed dimensions, cannot deal with that");
      base_node val(N), bmin(N), bmax(N), n0(N), n(N), n_mean(N);
      base_matrix grad(N,N);
      mesh_region region = m.region(bnum);
      GMM_ASSERT1(mfu.get_qdim() == N,
                  "Wrong mesh_fem qdim to compute contact pairs");

      dal::bit_vector points_already_interpolated;
      std::vector<base_node> transformed_points(m.nb_max_points());
      for (getfem::mr_visitor v(region,m); !v.finished(); ++v) {
        size_type cv = v.cv();
        bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
        pfem pf_s = mfu.fem_of_element(cv);
        size_type nbd_t = pgt->nb_points();
        slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
        bgeot::vectors_to_base_matrix
          (G, mfu.linked_mesh().points_of_convex(cv));

        pfem_precomp pfp = fppool(pf_s, &(pgt->geometric_nodes()));
        fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv,
                                      size_type(-1));

        size_type nb_pt_on_face = 0;
        gmm::clear(n_mean);
        for (short_type ip = 0; ip < nbd_t; ++ip) {
          size_type ind = m.ind_points_of_convex(cv)[ip];

          // computation of transformed vertex
          if (!(points_already_interpolated.is_in(ind))) {
            ctx.set_ii(ip);
            pf_s->interpolation(ctx, coeff, val, dim_type(N));
            val += ctx.xreal();
            transformed_points[ind] = val;
            points_already_interpolated.add(ind);
          } else {
            val = transformed_points[ind];
          }
          // computation of unit normal vector if the vertex is on the face
          bool is_on_face = false;
          bgeot::pconvex_structure cvs = pgt->structure();
          for (size_type k = 0; k < cvs->nb_points_of_face(v.f()); ++k)
            if (cvs->ind_points_of_face(v.f())[k] == ip) is_on_face = true;
          if (is_on_face) {
            ctx.set_ii(ip);
            n0 = bgeot::compute_normal(ctx, v.f());
            pf_s->interpolation_grad(ctx, coeff, grad, dim_type(N));
            gmm::add(gmm::identity_matrix(), grad);
            scalar_type J = gmm::lu_inverse(grad);
            if (J <= scalar_type(0)) GMM_WARNING1("Inverted element ! " << J);
            gmm::mult(gmm::transposed(grad), n0, n);
            n /= gmm::vect_norm2(n);
            n_mean += n;
            ++nb_pt_on_face;
          }

          if (ip == 0) // computation of bounding box
            bmin = bmax = val;
          else {
            for (size_type k = 0; k < N; ++k) {
              bmin[k] = std::min(bmin[k], val[k]);
              bmax[k] = std::max(bmax[k], val[k]);
            }
          }
        }

        GMM_ASSERT1(nb_pt_on_face,
                    "This element has not vertex on considered face !");

        // Computation of influence box :
        // offset of the bounding box relatively to its "diameter"
        scalar_type h = bmax[0] - bmin[0];
        for (size_type k = 1; k < N; ++k)
          h = std::max(h, bmax[k] - bmin[k]);
        for (size_type k = 0; k < N; ++k)
          { bmin[k] -= h; bmax[k] += h; }

        // Store the influence box and additional information.
        element_boxes.add_box(bmin, bmax, unit_normal_of_elements.size());
        n_mean /= gmm::vect_norm2(n_mean);
        unit_normal_of_elements.push_back(n_mean);
        boundary_of_elements.push_back(i);
        ind_of_elements.push_back(cv);
        face_of_elements.push_back(v.f());
      }
    }
  }



  bool contact_elements::add_point_contribution
  (size_type boundary_num, getfem::fem_interpolation_context &ctxu,
   getfem::fem_interpolation_context &ctxl, scalar_type weight,
   scalar_type /*f_coeff*/, scalar_type r, model::build_version version) {
    const mesh_fem &mfu = cf.mfu_of_boundary(boundary_num);
    const mesh_fem &mfl = cf.mflambda_of_boundary(boundary_num);
    const model_real_plain_vector &U = cf.disp_of_boundary(boundary_num);
    const model_real_plain_vector &L = cf.lambda_of_boundary(boundary_num);
    size_type N = mfu.get_qdim();
    base_node x0 = ctxu.xreal();
    bool noisy = false;

    // ----------------------------------------------------------
    // Computation of the point coordinates and the unit normal
    // vector in real configuration
    // ----------------------------------------------------------

    base_node n0 = bgeot::compute_normal(ctxu, ctxu.face_num());
    scalar_type face_factor = gmm::vect_norm2(n0);
    size_type cv = ctxu.convex_num();
    base_small_vector n(N), val(N), h(N);
    base_matrix gradinv(N,N), grad(N,N), gradtot(N,N), G;
    size_type cvnbdofu = mfu.nb_basic_dof_of_element(cv);
    size_type cvnbdofl = mfl.nb_basic_dof_of_element(cv);
    base_vector coeff(cvnbdofu);
    slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
    ctxu.pf()->interpolation(ctxu, coeff, val, dim_type(N));
    base_node x = x0 + val;

    ctxu.pf()->interpolation_grad(ctxu, coeff, gradinv, dim_type(N));
    gmm::add(gmm::identity_matrix(), gradinv);
    scalar_type J = gmm::lu_inverse(gradinv); // remplacer par une résolution...
    if (J <= scalar_type(0)) {
      GMM_WARNING1("Inverted element !");

      GMM_ASSERT1(!(version & model::BUILD_MATRIX), "Impossible to build "
                  "tangent matrix for large sliding contact");
      if (version & model::BUILD_RHS) {
        base_vector Velem(cvnbdofl);
        for (size_type i = 0; i < cvnbdofl; ++i) Velem[i] = 1E200;
        vec_elem_assembly(cf.L_vector(boundary_num), Velem, mfl, cv);
        return false;
      }
    }

    gmm::mult(gmm::transposed(gradinv), n0, n);
    n /= gmm::vect_norm2(n);

    // ----------------------------------------------------------
    // Selection of influence boxes
    // ----------------------------------------------------------

    bgeot::rtree::pbox_set bset;
    element_boxes.find_boxes_at_point(x, bset);

    if (noisy) cout << "Number of boxes found : " << bset.size() << endl;

    // ----------------------------------------------------------
    // Eliminates some influence boxes with the mean normal
    // criterion : should at least eliminate the original element.
    // ----------------------------------------------------------

    bgeot::rtree::pbox_set::iterator it = bset.begin(), itnext;
    for (; it != bset.end(); it = itnext) {
      itnext = it; ++itnext;
      if (gmm::vect_sp(unit_normal_of_elements[(*it)->id], n)
          >= -scalar_type(1)/scalar_type(20)) bset.erase(it);
    }

    if (noisy)
      cout << "Number of boxes satisfying the unit normal criterion : "
           << bset.size() << endl;


    // ----------------------------------------------------------
    // For each remaining influence box, compute y0, the corres-
    // ponding unit normal vector and eliminate wrong auto-contact
    // situations with a test on |x0-y0|
    // ----------------------------------------------------------

    it = bset.begin();
    std::vector<base_node> y0s;
    std::vector<base_small_vector> n0_y0s;
    std::vector<scalar_type> d0s;
    std::vector<scalar_type> d1s;
    std::vector<size_type> elt_nums;
    std::vector<fem_interpolation_context> ctx_y0s;
    for (; it != bset.end(); ++it) {
      size_type boundary_num_y0 = boundary_of_elements[(*it)->id];
      size_type cv_y0 = ind_of_elements[(*it)->id];
      short_type face_y0 = short_type(face_of_elements[(*it)->id]);
      const mesh_fem &mfu_y0 = cf.mfu_of_boundary(boundary_num_y0);
      pfem pf_s_y0 = mfu_y0.fem_of_element(cv_y0);
      const model_real_plain_vector &U_y0
        = cf.disp_of_boundary(boundary_num_y0);
      const mesh &m_y0 = mfu_y0.linked_mesh();
      bgeot::pgeometric_trans pgt_y0 = m_y0.trans_of_convex(cv_y0);
      bgeot::pconvex_structure cvs_y0 = pgt_y0->structure();

      // Find an interior point (in order to promote the more interior
      // y0 in case of locally non invertible transformation.
      size_type ind_dep_point = 0;
      for (; ind_dep_point < cvs_y0->nb_points(); ++ind_dep_point) {
        bool is_on_face = false;
        for (size_type k = 0;
             k < cvs_y0->nb_points_of_face(face_y0); ++k)
          if (cvs_y0->ind_points_of_face(face_y0)[k]
              == ind_dep_point) is_on_face = true;
        if (!is_on_face) break;
      }
      GMM_ASSERT1(ind_dep_point < cvs_y0->nb_points(),
                  "No interior point found !");

      base_node y0_ref = pgt_y0->convex_ref()->points()[ind_dep_point];

      slice_vector_on_basic_dof_of_element(mfu_y0, U_y0, cv_y0, coeff);
      // if (pf_s_y0->need_G())
      bgeot::vectors_to_base_matrix(G, m_y0.points_of_convex(cv_y0));

      fem_interpolation_context ctx_y0(pgt_y0, pf_s_y0, y0_ref, G, cv_y0,
                                       size_type(-1));

      size_type newton_iter = 0;
      for(;;) { // Newton algorithm to invert geometric transformation

        pf_s_y0->interpolation(ctx_y0, coeff, val, dim_type(N));
        val += ctx_y0.xreal() - x;
        scalar_type init_res = gmm::vect_norm2(val);

        if (init_res < 1E-12) break;
        if (newton_iter > 100) {
          GMM_WARNING1("Newton has failed to invert transformation"); // il faudrait faire qlq chose d'autre ... !
          GMM_ASSERT1(!(version & model::BUILD_MATRIX), "Impossible to build "
                      "tangent matrix for large sliding contact");
          if (version & model::BUILD_RHS) {
            base_vector Velem(cvnbdofl);
            for (size_type i = 0; i < cvnbdofl; ++i) Velem[i] = 1E200;
            vec_elem_assembly(cf.L_vector(boundary_num), Velem, mfl, cv);
            return false;
          }
        }

        pf_s_y0->interpolation_grad(ctx_y0, coeff, grad, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad);
        gmm::mult(grad, ctx_y0.K(), gradtot);

        std::vector<int> ipvt(N);
        size_type info = gmm::lu_factor(gradtot, ipvt);
        GMM_ASSERT1(!info, "Singular system, pivot = " << info); // il faudrait faire qlq chose d'autre ... perturber par exemple
        gmm::lu_solve(gradtot, ipvt, h, val);

        // line search
        bool ok = false;
        scalar_type alpha;
        for (alpha = 1; alpha >= 1E-5; alpha/=scalar_type(2)) {

          ctx_y0.set_xref(y0_ref - alpha*h);
          pf_s_y0->interpolation(ctx_y0, coeff, val, dim_type(N));
          val += ctx_y0.xreal() - x;

          if (gmm::vect_norm2(val) < init_res) { ok = true; break; }
        }
        if (!ok)
          GMM_WARNING1("Line search has failed to invert transformation");
        y0_ref -= alpha*h;
        ctx_y0.set_xref(y0_ref);
        newton_iter++;
      }

      base_node y0 = ctx_y0.xreal();
      base_node n0_y0 = bgeot::compute_normal(ctx_y0, face_y0);
      scalar_type d0_ref = pgt_y0->convex_ref()->is_in_face(face_y0, y0_ref);
      scalar_type d0 = d0_ref / gmm::vect_norm2(n0_y0);


      scalar_type d1 = d0_ref; // approximatively a distance to the element
      short_type ifd = short_type(-1);

      for (short_type k = 0; k <  pgt_y0->structure()->nb_faces(); ++k) {
        scalar_type dd = pgt_y0->convex_ref()->is_in_face(k, y0_ref);
        if (dd > scalar_type(0) && dd > gmm::abs(d1)) { d1 = dd; ifd = k; }
      }

      if (ifd != short_type(-1)) {
        d1 /= gmm::vect_norm2(bgeot::compute_normal(ctx_y0, ifd));
        if (gmm::abs(d1) < gmm::abs(d0)) d1 = d0;
      } else d1 = d0;

//       size_type iptf = m_y0.ind_points_of_face_of_convex(cv_y0, face_y0)[0];
//       base_node ptf = x0 - m_y0.points()[iptf];
//       scalar_type d2 = gmm::vect_sp(ptf, n0_y0) / gmm::vect_norm2(n0_y0);

      if (noisy) cout << "gmm::vect_norm2(n0_y0) = " << gmm::vect_norm2(n0_y0) << endl;
      // Eliminates wrong auto-contact situations
      if (noisy) cout << "autocontact status : x0 = " << x0 << " y0 = " << y0 << "  " <<  gmm::vect_dist2(y0, x0) << " : " << d0*0.75 << " : " << d1*0.75 << endl;
      if (noisy) cout << "n = " << n << " unit_normal_of_elements[(*it)->id] = " << unit_normal_of_elements[(*it)->id] << endl;

      if (d0 < scalar_type(0)
          && ((&U_y0 == &U
               && (gmm::vect_dist2(y0, x0) < gmm::abs(d1)*scalar_type(3)/scalar_type(4)))
              || gmm::abs(d1) > 0.05)) {
        if (noisy)  cout << "Eliminated x0 = " << x0 << " y0 = " << y0
                         << " d0 = " << d0 << endl;
        continue;
      }

//       if (d0 < scalar_type(0) && &(U_y0) == &U
//           && gmm::vect_dist2(y0, x0) < gmm::abs(d1) * scalar_type(2)
//           && d2 < -ctxu.J() / scalar_type(2)) {
//         if (noisy) cout << "Eliminated x0 = " << x0 << " y0 = " << y0
//                         << " d0 = " << d0 << endl;
//         continue;
//       }

      y0s.push_back(ctx_y0.xreal()); // useful ?
      elt_nums.push_back((*it)->id);
      d0s.push_back(d0);
      d1s.push_back(d1);
      ctx_y0s.push_back(ctx_y0);
      n0_y0 /= gmm::vect_norm2(n0_y0);
      n0_y0s.push_back(n0_y0);

      if (noisy) cout << "dist0 = " << d0 << " dist0 * area = "
                      << pgt_y0->convex_ref()->is_in(y0_ref) << endl;
    }

    // ----------------------------------------------------------
    // Compute the distance to rigid obstacles and selects the
    // nearest boundary/obstacle.
    // ----------------------------------------------------------

    dim_type state = 0;
    scalar_type d0 = 1E100, d1 = 1E100;
    base_small_vector grad_obs(N);

    size_type ibound = size_type(-1);
    for (size_type k = 0; k < y0s.size(); ++k)
      if (d1s[k] < d1) { d0 = d0s[k]; d1 = d1s[k]; ibound = k; state = 1; }


    size_type irigid_obstacle = size_type(-1);
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
    gmm::copy(x, cf.pt_eval);
    for (size_type i = 0; i < cf.obstacles.size(); ++i) {
      scalar_type d0_o = scalar_type(cf.obstacles_parsers[i].Eval());
      if (d0_o < d0) { d0 = d0_o; irigid_obstacle = i; state = 2; }
    }
    if (state == 2) {
      scalar_type EPS = face_factor * 1E-9;
      for (size_type k = 0; k < N; ++k) {
        cf.pt_eval[k] += EPS;
        grad_obs[k] =
          (scalar_type(cf.obstacles_parsers[irigid_obstacle].Eval())-d0)/EPS;
        cf.pt_eval[k] -= EPS;
      }
    }

#else
    if (cf.obstacles.size() > 0)
      GMM_WARNING1("Rigid obstacles are ignored. Recompile with "
                   "muParser to account for rigid obstacles");
#endif


    // ----------------------------------------------------------
    // Print the found contact state ...
    // ----------------------------------------------------------


    if (noisy && state == 1) {
      cout  << "Point : " << x0 << " of boundary " << boundary_num
            << " and element " << cv << " state = " << int(state);
      if (version & model::BUILD_RHS) cout << " RHS";
      if (version & model::BUILD_MATRIX) cout << " MATRIX";
    }
    if (state == 1) {
      size_type boundary_num_y0 = boundary_of_elements[elt_nums[ibound]];
      const mesh_fem &mfu_y0 = cf.mfu_of_boundary(boundary_num_y0);
      const mesh &m_y0 = mfu_y0.linked_mesh();
      size_type cv_y0 = ind_of_elements[elt_nums[ibound]];

      if (noisy) cout << " y0 = " << y0s[ibound] << " of element "
                            << cv_y0  << " of boundary " << boundary_num_y0 << endl;
      for (size_type k = 0; k < m_y0.nb_points_of_convex(cv_y0); ++k)
        if (noisy) cout << "point " << k << " : "
                        << m_y0.points()[m_y0.ind_points_of_convex(cv_y0)[k]] << endl;
      if (boundary_num_y0 == 0 && boundary_num == 0 && d0 < 0.0 && (version & model::BUILD_MATRIX)) GMM_ASSERT1(false, "oups");
    }
    if (noisy) cout << " d0 = " << d0 << endl;

    // ----------------------------------------------------------
    // Add the contributions to the tangent matrices and rhs
    // ----------------------------------------------------------

    GMM_ASSERT1(ctxu.pf()->target_dim() == 1 && ctxl.pf()->target_dim() == 1,
                "Large sliding contact assembly procedure has to be adapted "
                "to intrinsic vectorial elements. To be done.");

    // Éviter les calculs inutiles dans le cas state == 2 ... à voir à la fin
    // regarder aussi si on peut factoriser des mat_elem_assembly ...

    base_matrix Melem;
    base_vector Velem;

    base_tensor tl, tu;
    ctxl.base_value(tl);
    ctxu.base_value(tu);

    base_small_vector lambda(N);
    slice_vector_on_basic_dof_of_element(mfl, L, cv, coeff);
    ctxl.pf()->interpolation(ctxl, coeff, lambda, dim_type(N));
    GMM_ASSERT1(!(isnan(lambda[0])), "internal error");

    // Unstabilized frictionless case for the moment

    // auxiliary variables
    scalar_type aux1, aux2;

    if (state) {

      // zeta = lamda + d0 * r * n
      base_small_vector zeta(N);
      gmm::add(lambda, gmm::scaled(n, r*d0), zeta);

      base_tensor tgradu;
      ctxu.grad_base_value(tgradu);

      // variables for y0
      base_tensor tu_y0;
      size_type boundary_num_y0 = 0, cv_y0 = 0, cvnbdofu_y0 = 0;
      if (state == 1) {
        ctx_y0s[ibound].base_value(tu_y0);
        boundary_num_y0 = boundary_of_elements[elt_nums[ibound]];
        cv_y0 = ind_of_elements[elt_nums[ibound]];
        cvnbdofu_y0 = cf.mfu_of_boundary(boundary_num_y0).nb_basic_dof_of_element(cv_y0);
      }
      const mesh_fem &mfu_y0 = (state == 1) ?
                               cf.mfu_of_boundary(boundary_num_y0) : mfu;

      if (version & model::BUILD_RHS) {
        // Rhs term Lx
        gmm::resize(Velem, cvnbdofl); gmm::clear(Velem);

        // Rhs term Lx: (1/r)\int (\lambda - P(\zeta)).\mu
        base_small_vector vecaux(N);
        gmm::copy(zeta, vecaux);
        De_Saxce_projection(vecaux, n, scalar_type(0));
        gmm::scale(vecaux, -scalar_type(1));
        gmm::add(lambda, vecaux);
        for (size_type i = 0; i < cvnbdofl; ++i)
          Velem[i] = tl[i/N] * vecaux[i%N] * weight/r;
        vec_elem_assembly(cf.L_vector(boundary_num), Velem, mfl, cv);

        // Rhs terms Ux, Uy: \int \lambda.(\psi(x_0) - \psi(y_0))
        gmm::resize(Velem, cvnbdofu); gmm::clear(Velem);
        for (size_type i = 0; i < cvnbdofu; ++i)
          Velem[i] = tu[i/N] * lambda[i%N] * weight;
        vec_elem_assembly(cf.U_vector(boundary_num), Velem, mfu, cv);

        if (state == 1) {
          gmm::resize(Velem, cvnbdofu_y0); gmm::clear(Velem);
          for (size_type i = 0; i < cvnbdofu_y0; ++i)
            Velem[i] = -tu_y0[i/N] * lambda[i%N] * weight;
          vec_elem_assembly(cf.U_vector(boundary_num_y0), Velem, mfu_y0, cv_y0);
        }
      }

      if (version & model::BUILD_MATRIX) {

        base_small_vector gradinv_n(N);
        gmm::mult(gradinv, n, gradinv_n);

        // de Saxce projection gradient and normal gradient at zeta
        base_matrix pgrad(N,N), pgradn(N,N);
        De_Saxce_projection_grad(zeta, n, scalar_type(0), pgrad);
        De_Saxce_projection_gradn(zeta, n, scalar_type(0), pgradn);

        base_small_vector pgrad_n(N), pgradn_n(N);
        gmm::mult(pgrad, n, pgrad_n);
        gmm::mult(pgradn, n, pgradn_n);
        base_matrix gradinv_pgrad(N,N), gradinv_pgradn(N,N);
        gmm::mult(gradinv, gmm::transposed(pgrad), gradinv_pgrad);
        gmm::mult(gradinv, gmm::transposed(pgradn), gradinv_pgradn);

        // Tangent term LxLx
        gmm::resize(Melem, cvnbdofl, cvnbdofl); gmm::clear(Melem);
        // -(1/r) \int \delta\lambda.\mu
        for (size_type i = 0; i < cvnbdofl; i += N) {
          aux1 = -tl[i/N] * weight/r;
          for (size_type j = 0; j < cvnbdofl; j += N) {
            aux2 = aux1 * tl[j/N];
            for (size_type k = 0; k < N; k++) Melem(i+k,j+k) = aux2;
          } // Melem(i+k,j+k) = -tl[i/N] * tl[j/N] * weight/r;
        }
        // (1/r) \int \nabla P(\zeta) (d\zeta/d\lambda)(\delta\lambda) . \mu
        for (size_type i = 0, ii = 0; i < cvnbdofl; ++i, ii = i%N)
          for (size_type j = 0, jj = 0; j < cvnbdofl; ++j, jj = j%N)
            Melem(i,j) += tl[i/N] * tl[j/N] * pgrad(ii,jj) * weight/r;
        mat_elem_assembly(cf.LL_matrix(boundary_num, boundary_num),
                          Melem, mfl, cv, mfl, cv);

        // Tangent term UxLx
        gmm::resize(Melem, cvnbdofu, cvnbdofl); gmm::clear(Melem);
        // \int -\delta\lambda.\psi(x_0)
        for (size_type i = 0; i < cvnbdofu; i += N) {
          aux1 = -tu[i/N] * weight;
          for (size_type j = 0; j < cvnbdofl; j += N) {
            aux2 = aux1 * tl[j/N];
            for (size_type k = 0; k < N; k++) Melem(i+k,j+k) = aux2;
          }
        }
        mat_elem_assembly(cf.UL_matrix(boundary_num, boundary_num),
                          Melem, mfu, cv, mfl, cv);

        // Tangent term LxUx
        if (0) { // DISABLED
        gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);
        // \int d_0(\nabla P(\zeta))(dn/du)(\delta u).\mu
        for (size_type i = 0, ii = 0; i < cvnbdofl; ++i, ii = i%N)
          for (size_type j = 0, jj = 0; j < cvnbdofu; ++j, jj = j%N) {
            aux1 = aux2 = scalar_type(0);
            for (size_type k = 0; k < N; ++k) {
              aux1 += tgradu[j/N+N*k] * gradinv_n[k];
              aux2 += tgradu[j/N+N*k] * gradinv_pgrad(k,ii);
            }
            Melem(i,j) = d0 * tl[i/N] * (pgrad_n[ii] * aux1 - aux2) * n[jj] * weight;
          }

        // (1/r)\int \nabla_n P(zeta) (dn/du)(\delta u) . \mu
        // On peut certainement factoriser d'avantage ce terme avec le
        // précédent. Attendre la version avec frottement.
        for (size_type i = 0, ii = 0; i < cvnbdofl; ++i, ii = i%N)
          for (size_type j = 0, jj = 0; j < cvnbdofu; ++j, jj = j%N) {
            aux1 = aux2 = scalar_type(0);
            for (size_type k = 0; k < N; ++k) {
              aux1 += tgradu[j/N+N*k] * gradinv_n[k];
              aux2 += tgradu[j/N+N*k] * gradinv_pgradn(k,ii);
            }
            Melem(i,j) += tl[i/N] * (pgradn_n[ii] * aux1 - aux2) * n[jj] * weight / r;
          }
        mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num),
                          Melem, mfl, cv, mfu, cv);
        } // DISABLED

        if (state == 1) {

          base_tensor tgradu_y0;
          ctx_y0s[ibound].grad_base_value(tgradu_y0);

          base_matrix gradinv_y0(N,N);
          base_small_vector ntilde_y0(N);
          { // calculate gradinv_y0 and ntilde_y0
            base_matrix grad_y0(N,N);
            base_vector coeff_y0(cvnbdofu_y0);
            const model_real_plain_vector &U_y0
              = cf.disp_of_boundary(boundary_num_y0);
            slice_vector_on_basic_dof_of_element(mfu_y0, U_y0, cv_y0, coeff_y0);
            ctx_y0s[ibound].pf()->interpolation_grad(ctx_y0s[ibound], coeff_y0,
                                                   grad_y0, dim_type(N));
            gmm::add(gmm::identity_matrix(), grad_y0);

            gmm::copy(grad_y0, gradinv_y0);
            gmm::lu_inverse(gradinv_y0); // à proteger contre la non-inversibilité
            gmm::mult(gmm::transposed(gradinv_y0), n0_y0s[ibound], ntilde_y0); // (not unit) normal vector
          }

          // Tangent term UyLx: \int \delta\lambda.\psi(y_0)
          gmm::resize(Melem, cvnbdofu_y0, cvnbdofl); gmm::clear(Melem);
          for (size_type i = 0; i < cvnbdofu_y0; i += N) {
            aux1 = tu_y0[i/N] * weight;
            for (size_type j = 0; j < cvnbdofl; j += N) {
              aux2 = aux1 * tl[j/N];
              for (size_type k = 0; k < N; k++) Melem(i+k,j+k) = aux2;
            }
          }
          mat_elem_assembly(cf.UL_matrix(boundary_num_y0, boundary_num),
                            Melem, mfu_y0, cv_y0, mfl, cv);

          // Tangent terms UyUx, UyUy
          // \int \lambda.((\nabla \psi(y_0))(I+\nabla u(y_0))^{-1}(\delta u(x_0) - \delta u(y_0)))

          // Tangent term UyUx
          gmm::resize(Melem, cvnbdofu_y0, cvnbdofu); gmm::clear(Melem);
          // \int \lambda.((\nabla \psi(y_0))(I+\nabla u(y_0))^{-1}\delta u(x_0))
          for (size_type i = 0, ii = 0; i < cvnbdofu_y0; ++i, ii = i%N)
            for (size_type j = 0, jj = 0; j < cvnbdofu; ++j, jj = j%N) {
              aux1 = scalar_type(0);
              for (size_type k = 0; k < N; ++k)
                aux1 += tgradu_y0[i/N+N*k]* gradinv_y0(k,jj);
              Melem(i,j) = lambda[ii] * aux1 * tu[j/N] * weight;
            }
          mat_elem_assembly(cf.UU_matrix(boundary_num_y0, boundary_num),
                            Melem, mfu_y0, cv_y0, mfu, cv);

          // Tangent term UyUy
          gmm::resize(Melem, cvnbdofu_y0, cvnbdofu_y0); gmm::clear(Melem);
          // -\int \lambda.((\nabla \psi(y_0))(I+\nabla u(y_0))^{-1}\delta u(y_0))
          for (size_type i = 0, ii = 0; i < cvnbdofu_y0; ++i, ii = i%N)
            for (size_type j = 0, jj = 0; j < cvnbdofu_y0; ++j, jj = j%N) {
              aux1 = scalar_type(0);
              for (size_type k = 0; k < N; ++k)
                aux1 += tgradu_y0[i/N+N*k] * gradinv_y0(k,jj);
              Melem(i,j) = - lambda[ii] * aux1 * tu_y0[j/N] * weight;
            }
          mat_elem_assembly(cf.UU_matrix(boundary_num_y0, boundary_num_y0),
                            Melem, mfu_y0, cv_y0, mfu_y0, cv_y0);

          // Tangent term LxUy
          gmm::resize(Melem, cvnbdofl, cvnbdofu_y0); gmm::clear(Melem);
          // -\int (I+\nabla u(y_0))^{-T}\nabla \delta(y_0).\delta u(y_0)(\nabla P(\zeta) n . \mu)
          for (size_type i = 0; i < cvnbdofl; ++i) {
            aux1 = tl[i/N] * pgrad_n[i%N] * weight;
            for (size_type j = 0; j < cvnbdofu_y0; ++j)
              Melem(i,j) = - aux1 * tu_y0[j/N] * ntilde_y0[j%N];
          }
          mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num_y0),
                            Melem, mfl, cv, mfu_y0, cv_y0);

          // Addition to tangent term LxUx
          gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);
          // \int (I+\nabla u(y_0))^{-T}\nabla \delta(y_0).\delta u(x_0)(\nabla P(\zeta) n . \mu)
          for (size_type i = 0; i < cvnbdofl; ++i) {
            aux1 = tl[i/N] * pgrad_n[i%N] * weight;
            for (size_type j = 0; j < cvnbdofu; ++j)
              Melem(i,j) = aux1 * tu[j/N] * ntilde_y0[j%N];
          }
        }
        else {
          // Addition to tangent term LxUx
          gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);
          // \int (I+\nabla u(y_0))^{-T}\nabla \delta(y_0).\delta u(x_0)(\nabla P(\zeta) n . \mu)
          for (size_type i = 0; i < cvnbdofl; ++i) {
            aux1 = tl[i/N] * pgrad_n[i%N] * weight;
            for (size_type j = 0; j < cvnbdofu; ++j)
              Melem(i,j) = aux1 * tu[j/N] * grad_obs[j%N];
          }
        }
        mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num),
                          Melem, mfl, cv, mfu, cv);

      }

    } else { // state == 0

      // Rhs term Lx: (1/r)\int \lambda.\mu
      if (version & model::BUILD_RHS) {
        gmm::resize(Velem, cvnbdofl); gmm::clear(Velem);
        for (size_type i = 0; i < cvnbdofl; ++i)
          Velem[i] = tl[i/N] * lambda[i%N] * weight/r;
        vec_elem_assembly(cf.L_vector(boundary_num), Velem, mfl, cv);
      }

      // Tangent term LxLx: -(1/r)\int \delta\lambda.\mu
      if (version & model::BUILD_MATRIX) {
        gmm::resize(Melem, cvnbdofl, cvnbdofl); gmm::clear(Melem);
        for (size_type i = 0; i < cvnbdofl; i += N) {
          aux1 = -tl[i/N] * weight/r;
          for (size_type j = 0; j < cvnbdofl; j += N) {
            aux2 = aux1 * tl[j/N];
            for (size_type k = 0; k < N; k++) Melem(i+k,j+k) = aux2;
          } // Melem(i+k,j+k) = -tl[i/N] * tl[j/N] * weight/r;
        }
        mat_elem_assembly(cf.LL_matrix(boundary_num, boundary_num),
                          Melem, mfl, cv, mfl, cv);
      }
    }

    return true;
  }

  //=========================================================================
  // 3)- Large sliding contact brick
  //=========================================================================

  struct integral_large_sliding_contact_brick_field_extension : public virtual_brick {


    struct contact_boundary {
      size_type region;
      std::string varname;
      std::string multname;
      const mesh_im *mim;
    };

    std::vector<contact_boundary> boundaries;
    std::vector<std::string> obstacles;

    void add_boundary(const std::string &varn, const std::string &multn,
                      const mesh_im &mim, size_type region) {
      contact_boundary cb;
      cb.region = region; cb.varname = varn; cb.multname = multn; cb.mim=&mim;
      boundaries.push_back(cb);
    }

    void add_obstacle(const std::string &obs)
    { obstacles.push_back(obs); }

    void build_contact_frame(const model &md, contact_frame &cf) const {
      for (size_type i = 0; i < boundaries.size(); ++i) {
        const contact_boundary &cb = boundaries[i];
        cf.add_boundary(md.mesh_fem_of_variable(cb.varname),
                        md.real_variable(cb.varname),
                        md.mesh_fem_of_variable(cb.multname),
                        md.real_variable(cb.multname), cb.region);
      }
      for (size_type i = 0; i < obstacles.size(); ++i)
        cf.add_obstacle(obstacles[i]);
    }


    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const;

    integral_large_sliding_contact_brick_field_extension() {
      set_flags("Integral large sliding contact brick",
                false /* is linear*/, false /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

  };




  void integral_large_sliding_contact_brick_field_extension::asm_real_tangent_terms
  (const model &md, size_type /* ib */, const model::varnamelist &vl,
   const model::varnamelist &dl, const model::mimlist &/* mims */,
   model::real_matlist &matl, model::real_veclist &vecl,
   model::real_veclist &, size_type /* region */,
   build_version version) const {

    fem_precomp_pool fppool;
    base_matrix G;
    size_type N = md.mesh_fem_of_variable(vl[0]).linked_mesh().dim();
    contact_frame cf(N);
    build_contact_frame(md, cf);

    size_type Nvar = vl.size(), Nu = cf.Urhs.size(), Nl = cf.Lrhs.size();
    GMM_ASSERT1(Nvar == Nu+Nl, "Wrong size of variable list for integral "
                "large sliding contact brick");
    GMM_ASSERT1(matl.size() == Nvar*Nvar, "Wrong size of terms for "
                "integral large sliding contact brick");

    if (version & model::BUILD_MATRIX) {
      for (size_type i = 0; i < Nvar; ++i)
        for (size_type j = 0; j < Nvar; ++j) {
          gmm::clear(matl[i*Nvar+j]);
          if (i <  Nu && j <  Nu) cf.UU(i,j)       = &(matl[i*Nvar+j]);
          if (i >= Nu && j <  Nu) cf.LU(i-Nu,j)    = &(matl[i*Nvar+j]);
          if (i <  Nu && j >= Nu) cf.UL(i,j-Nu)    = &(matl[i*Nvar+j]);
          if (i >= Nu && j >= Nu) cf.LL(i-Nu,j-Nu) = &(matl[i*Nvar+j]);
        }
    }
    if (version & model::BUILD_RHS) {
      for (size_type i = 0; i < vl.size(); ++i) {
        if (i < Nu) cf.Urhs[i] = &(vecl[i*Nvar]);
        else cf.Lrhs[i-Nu] = &(vecl[i*Nvar]);
      }
    }

    // Data : r, [friction_coeff,]
    GMM_ASSERT1(dl.size() == 2, "Wrong number of data for integral large "
                "sliding contact brick");

    const model_real_plain_vector &vr = md.real_variable(dl[0]);
    GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");

    const model_real_plain_vector &f_coeff = md.real_variable(dl[1]);
    GMM_ASSERT1(gmm::vect_size(f_coeff) == 1,
                "Friction coefficient should be a scalar");

    contact_elements ce(cf);
    ce.init();

    for (size_type bnum = 0; bnum < boundaries.size(); ++bnum) {
      mesh_region rg(boundaries[bnum].region);
      const mesh_fem &mfu=md.mesh_fem_of_variable(boundaries[bnum].varname);
      const mesh_fem &mfl=md.mesh_fem_of_variable(boundaries[bnum].multname);
      const mesh_im &mim = *(boundaries[bnum].mim);
      const mesh &m = mfu.linked_mesh();
      mfu.linked_mesh().intersect_with_mpi_region(rg);

      for (getfem::mr_visitor v(rg, m); !v.finished(); ++v) {
        // cout << "boundary " << bnum << " element " << v.cv() << endl;
        size_type cv = v.cv();
        bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
        pfem pf_s = mfu.fem_of_element(cv);
        pfem pf_sl = mfl.fem_of_element(cv);
        pintegration_method pim = mim.int_method_of_element(cv);
        bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));

        pfem_precomp pfpu
          = fppool(pf_s,&(pim->approx_method()->integration_points()));
        pfem_precomp pfpl
          = fppool(pf_sl,&(pim->approx_method()->integration_points()));
        fem_interpolation_context ctxu(pgt,pfpu,size_type(-1), G, cv, v.f());
        fem_interpolation_context ctxl(pgt,pfpl,size_type(-1), G, cv, v.f());

        for (size_type k = 0;
             k < pim->approx_method()->nb_points_on_face(v.f()); ++k) {
          size_type ind
            = pim->approx_method()->ind_first_point_on_face(v.f()) + k;
          ctxu.set_ii(ind);
          ctxl.set_ii(ind);
          if (!(ce.add_point_contribution
               (bnum, ctxu, ctxl,pim->approx_method()->coeff(ind),
                f_coeff[0], vr[0], version))) return;
        }
      }
    }
  }


  // r ne peut pas être variable pour le moment.
  // dataname_friction_coeff ne peut pas être variable non plus ...

  size_type add_integral_large_sliding_contact_brick_field_extension
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_r,
   const std::string &dataname_friction_coeff, size_type region) {

    integral_large_sliding_contact_brick_field_extension *pbr
      = new integral_large_sliding_contact_brick_field_extension();

    pbr->add_boundary(varname_u, multname, mim, region);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname,  false));
    tl.push_back(model::term_description(multname,  varname_u, false));
    tl.push_back(model::term_description(multname,  multname,  false));

    model::varnamelist dl(1, dataname_r);
    dl.push_back(dataname_friction_coeff);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  void add_boundary_to_large_sliding_contact_brick
  (model &md, size_type indbrick, const mesh_im &mim,
   const std::string &varname_u, const std::string &multname,
   size_type region) {
    dim_type N = md.mesh_fem_of_variable(varname_u).linked_mesh().dim();
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    integral_large_sliding_contact_brick_field_extension *p
      = dynamic_cast<integral_large_sliding_contact_brick_field_extension *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    p->add_boundary(varname_u, multname, mim, region);
    md.add_mim_to_brick(indbrick, mim);

    contact_frame cf(N);
    p->build_contact_frame(md, cf);

    model::varnamelist vl;
    size_type nvaru = 0;
    for (size_type i = 0; i < cf.contact_boundaries.size(); ++i)
      if (cf.contact_boundaries[i].ind_U >= nvaru)
        { vl.push_back(p->boundaries[i].varname); ++nvaru; }

    size_type nvarl = 0;
    for (size_type i = 0; i < cf.contact_boundaries.size(); ++i)
      if (cf.contact_boundaries[i].ind_lambda >= nvarl)
        { vl.push_back(p->boundaries[i].multname); ++nvarl; }
    md.change_variables_of_brick(indbrick, vl);

    model::termlist tl;
    for (size_type i = 0; i < vl.size(); ++i)
      for (size_type j = 0; j < vl.size(); ++j)
        tl.push_back(model::term_description(vl[i], vl[j], false));

    md.change_terms_of_brick(indbrick, tl);
  }

  void add_rigid_obstacle_to_large_sliding_contact_brick
  (model &md, size_type indbrick, const std::string &obs) { // The velocity field should be added to an (optional) parameter ... (and optionaly represented by a rigid motion only ... the velocity should be modifiable ...
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    integral_large_sliding_contact_brick_field_extension *p
      = dynamic_cast<integral_large_sliding_contact_brick_field_extension *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    p->add_obstacle(obs);
  }

}  /* end of namespace getfem.                                             */
