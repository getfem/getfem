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
#include "getfem/getfem_contact_and_friction_large_sliding.h"
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


  struct gauss_point_precomp {
    size_type N;
    fem_precomp_pool fppool;
    const multi_contact_frame &mcf;
    const model &md;
    const multi_contact_frame::contact_pair *cp;

    const base_node &x(void) const { return cp->slave_point; }
    const base_node &nx(void) const { return cp->slave_n; }
    const base_node &y(void) const { return cp->master_point; }
    const base_node &y_ref(void) const { return cp->master_point_ref; }
    const base_node &ny(void) const { return cp->master_n; }
    scalar_type g(void) const { return cp->signed_dist; }

    base_matrix I_nxnx_;
    bool I_nxnx_computed;
    const base_matrix &I_nxnx(void) {
      if (!I_nxnx_computed) {
        gmm::copy(gmm::identity_matrix(), I_nxnx_);
        gmm::rank_one_update(I_nxnx_, nx(), gmm::scaled(nx(),scalar_type(-1)));
        I_nxnx_computed = true;
      }
      return I_nxnx_;
    }

    base_matrix I_nyny_;
    bool I_nyny_computed;
    const base_matrix &I_nyny(void) {
      if (!I_nyny_computed) {
        gmm::copy(gmm::identity_matrix(), I_nyny_);
        gmm::rank_one_update(I_nyny_, ny(), gmm::scaled(ny(),scalar_type(-1)));
        I_nyny_computed = true;
      }
      return I_nyny_;
    }

    base_matrix I_nxny_;
    bool I_nxny_computed;
    const base_matrix &I_nxny(void) {
      if (!I_nxny_computed) {
        gmm::copy(gmm::identity_matrix(), I_nxny_);
        gmm::rank_one_update(I_nxny_, nx(),
                             gmm::scaled(ny(),scalar_type(-1)/nxny));
        I_nxny_computed = true;
      }
      return I_nxny_;
    }

    scalar_type nxny;
    scalar_type nxdotny(void) const { return nxny; }

    bool isrigid_;
    bool isrigid(void) { return isrigid_; }

    base_tensor base_ux, base_uy, base_lx, base_ly;
    base_matrix vbase_ux_, vbase_uy_, vbase_lx_, vbase_ly_;
    bool vbase_ux_init, vbase_uy_init, vbase_lx_init, vbase_ly_init;
    base_tensor grad_base_ux, grad_base_uy, vgrad_base_ux_, vgrad_base_uy_;
    bool vgrad_base_ux_init, vgrad_base_uy_init;
    bool have_lx, have_ly;

    fem_interpolation_context ctx_ux_, ctx_uy_, ctx_lx_, ctx_ly_;
    bool ctx_ux_init, ctx_uy_init, ctx_lx_init, ctx_ly_init;
    base_matrix Gx, Gy;
    const mesh_fem *mf_ux_, *mf_uy_, *mf_lx_, *mf_ly_;
    gmm::sub_interval I_ux_, I_uy_, I_lx_, I_ly_;
    pfem pf_ux, pf_uy, pf_lx, pf_ly;
    size_type ndof_ux_, qdim_ux, ndof_uy_, qdim_uy, ndof_lx_, qdim_lx;
    size_type ndof_ly_, qdim_ly, cvx_, cvy_, ibx, iby;
    short_type fx, fy;
    bgeot::pgeometric_trans pgtx, pgty;
    const mesh_im *mim;
    pintegration_method pim;
    scalar_type weight_;

    scalar_type weight(void) { return weight_; }

    const mesh &meshx(void) const { return mf_ux_->linked_mesh(); }
    const mesh &meshy(void) const { return mf_uy_->linked_mesh(); }
    const mesh_fem *mf_ux(void) const { return mf_ux_; }
    const mesh_fem *mf_uy(void) const { return mf_uy_; }
    const mesh_fem *mf_lx(void) const { return mf_lx_; }
    const mesh_fem *mf_ly(void) const { return mf_ly_; }
    size_type ndof_ux(void) const { return ndof_ux_; }
    size_type ndof_uy(void) const { return ndof_uy_; }
    size_type ndof_lx(void) const { return ndof_lx_; }
    size_type ndof_ly(void) const { return ndof_ly_; }
    size_type cvx(void) const { return cvx_; }
    size_type cvy(void) const { return cvy_; }
    const gmm::sub_interval I_ux(void) const { return I_ux_; }
    const gmm::sub_interval I_uy(void) const { return I_uy_; }
    const gmm::sub_interval I_lx(void) const { return I_lx_; }
    const gmm::sub_interval I_ly(void) const { return I_ly_; }


    fem_interpolation_context &ctx_ux(void) {
      if (!ctx_ux_init) {
        bgeot::vectors_to_base_matrix(Gx, meshx().points_of_convex(cvx_));
        pfem_precomp pfp_ux
          = fppool(pf_ux, &(pim->approx_method()->integration_points()));
        ctx_ux_ = fem_interpolation_context(pgtx, pfp_ux, cp->slave_ind_pt,
                                            Gx, cvx_, fx);
        ctx_ux_init = true;
      }
      return ctx_ux_;
    }

    fem_interpolation_context &ctx_lx(void) {
      GMM_ASSERT1(have_lx, "No multiplier defined on the slave surface");
      if (!ctx_lx_init) {
        pfem_precomp pfp_lx
          = fppool(pf_lx, &(pim->approx_method()->integration_points()));
        ctx_lx_ = fem_interpolation_context(pgtx, pfp_lx, cp->slave_ind_pt,
                                            ctx_ux().G(), cvx_, fx);
        ctx_lx_init = true;
      }
      return ctx_lx_;
    }

    fem_interpolation_context &ctx_uy(void) {
      GMM_ASSERT1(!isrigid(), "Rigid obstacle master node: no fem defined");
      if (!ctx_uy_init) {
        bgeot::vectors_to_base_matrix(Gy, meshy().points_of_convex(cvy_));
        ctx_uy_ = fem_interpolation_context(pgty, pf_uy, y_ref(), Gy, cvy_, fy);
        ctx_uy_init = true;
      }
      return ctx_uy_;
    }

    fem_interpolation_context &ctx_ly(void) {
      GMM_ASSERT1(have_ly, "No multiplier defined on the master surface");
      if (!ctx_ly_init) {
        ctx_ly_ = fem_interpolation_context(pgty, pf_ly, y_ref(),
                                            ctx_uy().G(), cvy_, fy);
        ctx_ly_init = true;
      }
      return ctx_ly_;
    }

    const base_matrix &vbase_ux(void) {
      if (!vbase_ux_init) {
        ctx_ux().base_value(base_ux);
        vectorize_base_tensor(base_ux, vbase_ux_, ndof_ux_, qdim_ux, N);
        vbase_ux_init = true;
      }
      return vbase_ux_;
    }

    const base_matrix &vbase_uy(void) {
      if (!vbase_uy_init) {
        ctx_uy().base_value(base_uy);
        vectorize_base_tensor(base_uy, vbase_uy_, ndof_uy_, qdim_uy, N);
        vbase_uy_init = true;
      }
      return vbase_uy_;
    }

    const base_matrix &vbase_lx(void) {
      if (!vbase_lx_init) {
        ctx_lx().base_value(base_lx);
        vectorize_base_tensor(base_lx, vbase_lx_, ndof_lx_, qdim_lx, N);
        vbase_lx_init = true;
      }
      return vbase_lx_;
    }

    const base_matrix &vbase_ly(void) {
      if (!vbase_ly_init) {
        ctx_ly().base_value(base_ly);
        vectorize_base_tensor(base_ly, vbase_ly_, ndof_ly_, qdim_ly, N);
        vbase_ly_init = true;
      }
      return vbase_ly_;
    }

    const base_tensor &vgrad_base_ux(void) {
      if (!vgrad_base_ux_init) {
        ctx_ux().grad_base_value(grad_base_ux);
        vectorize_grad_base_tensor(grad_base_ux, vgrad_base_ux_, ndof_ux_,
                                   qdim_ux, N);
        vgrad_base_ux_init = true;
      }
      return vgrad_base_ux_;
    }

    const base_tensor &vgrad_base_uy(void) {
      if (!vgrad_base_uy_init) {
        ctx_uy().grad_base_value(grad_base_uy);
        vectorize_grad_base_tensor(grad_base_uy, vgrad_base_uy_, ndof_uy_,
                                   qdim_uy, N);
        vgrad_base_uy_init = true;
      }
      return vgrad_base_uy_;
    }

    base_small_vector lambda_x_, lambda_y_;
    bool lambda_x_init, lambda_y_init;
    base_vector coeff;

    const base_small_vector &lx(void) {
      if (!lambda_x_init) {
        pfem pf = ctx_lx().pf();
        slice_vector_on_basic_dof_of_element(*mf_lx_,mcf.mult_of_boundary(ibx),
                                             cvx_, coeff);
        pf->interpolation(ctx_lx(), coeff, lambda_x_, dim_type(N));
        lambda_x_init = true;
      }
      return lambda_x_;
    }

    const base_small_vector &ly(void) {
      if (!lambda_y_init) {
        pfem pf = ctx_ly().pf();
        slice_vector_on_basic_dof_of_element(*mf_ly_,mcf.mult_of_boundary(iby),
                                             cvy_, coeff);
        pf->interpolation(ctx_ly(), coeff, lambda_y_, dim_type(N));
        lambda_y_init = true;
      }
      return lambda_y_;
    }

    base_matrix grad_phix_, grad_phix_inv_, grad_phiy_, grad_phiy_inv_;
    bool grad_phix_init, grad_phix_inv_init;
    bool grad_phiy_init, grad_phiy_inv_init;

    const base_matrix &grad_phix(void) {
      if (!grad_phix_init) {
        pfem pf = ctx_ux().pf();
        slice_vector_on_basic_dof_of_element(*mf_ux_,mcf.disp_of_boundary(ibx),
                                             cvx_, coeff);
        pf->interpolation_grad(ctx_ux(), coeff, grad_phix_, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad_phix_);
        grad_phix_init = true;
      }
      return grad_phix_;
    }

    const base_matrix &grad_phix_inv(void) {
      if (!grad_phix_inv_init) {
        gmm::copy(grad_phix(), grad_phix_inv_);
        /* scalar_type J = */ gmm::lu_inverse(grad_phix_inv_);
        // if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
        grad_phix_inv_init = true;
      }
      return grad_phix_inv_;
    }

    const base_matrix &grad_phiy(void) {
      if (!grad_phiy_init) {
        pfem pf = ctx_uy().pf();
        slice_vector_on_basic_dof_of_element(*mf_uy_,mcf.disp_of_boundary(iby),
                                             cvy_, coeff);
        pf->interpolation_grad(ctx_uy(), coeff, grad_phiy_, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad_phiy_);
        grad_phiy_init = true;
      }
      return grad_phiy_;
    }

    const base_matrix &grad_phiy_inv(void) {
      if (!grad_phiy_inv_init) {
        gmm::copy(grad_phiy(), grad_phiy_inv_);
        /* scalar_type J = */ gmm::lu_inverse(grad_phiy_inv_);
        // if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
        grad_phiy_inv_init = true;
      }
      return grad_phiy_inv_;
    }

    scalar_type alpha;
    base_small_vector x0_, y0_, nx0_, Vs_;
    bool x0_init, y0_init, nx0_init, Vs_init;
    base_matrix grad_phiy0_;
    bool grad_phiy0_init;

    const base_small_vector &x0(void) {
      if (!x0_init) {
        const model_real_plain_vector &all_x0 = mcf.disp0_of_boundary(ibx);
        if (all_x0.size()) {
          pfem pf = ctx_ux().pf();
          slice_vector_on_basic_dof_of_element(*mf_ux_, all_x0, cvx_, coeff);
          pf->interpolation(ctx_ux(), coeff, x0_, dim_type(N));
        } else gmm::clear(x0_);
        gmm::add(ctx_ux().xreal(), x0_);
        x0_init = true;
      }
      return x0_;
    }

    const base_small_vector &y0(void) {
      if (!y0_init) {
        if (!isrigid()) {
          const model_real_plain_vector &all_y0 = mcf.disp0_of_boundary(iby);
          if (all_y0.size()) {
            pfem pf = ctx_uy().pf();
            slice_vector_on_basic_dof_of_element(*mf_uy_, all_y0, cvy_, coeff);
            pf->interpolation(ctx_uy(), coeff, y0_, dim_type(N));
          } else gmm::clear(y0_);
          gmm::add(ctx_uy().xreal(), y0_);
        } else gmm::copy(y(), y0_);
        y0_init = true;
      }
      return y0_;
    }

    const base_small_vector &nx0(void) {
      if (!nx0_init) {
        const model_real_plain_vector &all_x0 = mcf.disp0_of_boundary(ibx);
        if (all_x0.size()) {
          pfem pf = ctx_ux().pf();
          slice_vector_on_basic_dof_of_element(*mf_ux_, all_x0, cvx_, coeff);
          base_small_vector nx00_(N);
          base_matrix grad_phix0_(N,N);
          compute_normal(ctx_ux(), fx, false, coeff, nx00_, nx0_, grad_phix0_);
          nx0_ /= gmm::vect_norm2(nx0_);
        } else gmm::clear(nx0_);
        nx0_init = true;
      }
      return nx0_;
    }

    const base_small_vector &Vs(void) { // relative velocity
      if (!Vs_init) {
        if (alpha != scalar_type(0)) {
#ifdef CONSIDER_FRAME_INDIFFERENCE
          if (!isrigid()) {
            gmm::add(y0(), gmm::scaled(x0(), scalar_type(-1)), Vs_);
            gmm::add(gmm::scaled(nx0(), -g()), Vs_);
          } else
#endif
          {
            gmm::add(x(), gmm::scaled(y(), scalar_type(-1)), Vs_);
            gmm::add(gmm::scaled(x0(), scalar_type(-1)), Vs_);
            gmm::add(y0(), Vs_);
          }
          gmm::scale(Vs_, alpha);
        } else gmm::clear(Vs_);
        Vs_init = true;
      }
      return Vs_;
    }

    const base_matrix &grad_phiy0(void) { // grad_phiy of previous time step
      // To be verified ...
      if (!grad_phiy0_init) {
        const model_real_plain_vector &all_y0 = mcf.disp0_of_boundary(iby);
        if (!isrigid() && all_y0.size()) {
          pfem pf = ctx_uy().pf();
          slice_vector_on_basic_dof_of_element(*mf_uy_, all_y0, cvy_, coeff);
          pf->interpolation_grad(ctx_uy(), coeff, grad_phiy0_, dim_type(N));
          gmm::add(gmm::identity_matrix(), grad_phiy0_);
        } else gmm::copy(gmm::identity_matrix(), grad_phiy0_);
        grad_phiy0_init = true;
      }
      return grad_phiy0_;
    }

    base_small_vector un;

    void set_pair(const multi_contact_frame::contact_pair &cp_) {
      cp = &cp_;
      I_nxnx_computed = I_nyny_computed = I_nxny_computed = false;
      ctx_ux_init = ctx_uy_init = ctx_lx_init = ctx_ly_init = false;
      vbase_ux_init = vbase_uy_init = vbase_lx_init = vbase_ly_init = false;
      vgrad_base_ux_init = vgrad_base_uy_init = false;
      lambda_x_init = lambda_y_init = false;
      have_lx = have_ly = false;
      grad_phix_init = grad_phiy_init = false;
      grad_phix_inv_init = grad_phiy_inv_init = false;
      x0_init = y0_init = Vs_init = grad_phiy0_init = false;
      nxny = gmm::vect_sp(nx(), ny());
      isrigid_ = (cp->irigid_obstacle != size_type(-1));

      cvx_ = cp->slave_ind_element;
      ibx = cp->slave_ind_boundary;
      mf_ux_ = &(mcf.mfdisp_of_boundary(ibx));
      pf_ux = mf_ux_->fem_of_element(cvx_);
      qdim_ux = pf_ux->target_dim();
      ndof_ux_ = pf_ux->nb_dof(cvx_) * N / qdim_ux;
      fx = cp->slave_ind_face;
      pgtx = meshx().trans_of_convex(cvx_);
      mim = &(mcf.mim_of_boundary(ibx));
      pim = mim->int_method_of_element(cvx_);
      weight_ = pim->approx_method()->coeff(cp->slave_ind_pt) * ctx_ux().J();
      gmm::mult(ctx_ux().B(), pgtx->normals()[fx], un);
      weight_ *= gmm::vect_norm2(un);
      const std::string &name_ux = mcf.varname_of_boundary(ibx);
      I_ux_ = md.interval_of_variable(name_ux);

      const std::string &name_lx = mcf.multname_of_boundary(ibx);
      have_lx = (name_lx.size() > 0);
      if (have_lx) {
        mf_lx_ = &(mcf.mfmult_of_boundary(ibx));
        I_lx_ = md.interval_of_variable(name_lx);
        pf_lx = mf_lx_->fem_of_element(cvx_);
        qdim_lx = pf_lx->target_dim();
        ndof_lx_ = pf_lx->nb_dof(cvx_) * N / qdim_lx;
      }

      if (!isrigid_) {
        cvy_ = cp->master_ind_element;
        iby = cp->master_ind_boundary;
        fy = cp->master_ind_face;
        mf_uy_ = &(mcf.mfdisp_of_boundary(iby));
        pf_uy = mf_uy_->fem_of_element(cvy_);
        qdim_uy = pf_uy->target_dim();
        ndof_uy_ = pf_uy->nb_dof(cvy_) * N / qdim_uy;
        pgty = meshy().trans_of_convex(cvy_);

        const std::string &name_uy = mcf.varname_of_boundary(iby);
        I_uy_ = md.interval_of_variable(name_uy);
        const std::string &name_ly = mcf.multname_of_boundary(iby);
        have_ly = (name_ly.size() > 0);
        if (have_ly) {
          mf_ly_ = &(mcf.mfmult_of_boundary(iby));
          I_ly_ = md.interval_of_variable(name_ly);
          pf_ly = mf_ly_->fem_of_element(cvy_);
          qdim_ly = pf_ly->target_dim();
          ndof_ly_ = pf_ly->nb_dof(cvy_) * N / qdim_ly;
        }
      }
    }

    gauss_point_precomp(size_type N_, const model &md_,
                        const multi_contact_frame &mcf_, scalar_type alpha_) :
      N(N_), mcf(mcf_), md(md_),
      I_nxnx_(N,N), I_nyny_(N,N), I_nxny_(N,N),
      lambda_x_(N), lambda_y_(N),
      grad_phix_(N, N), grad_phix_inv_(N, N),
      grad_phiy_(N, N), grad_phiy_inv_(N, N), alpha(alpha_),
      x0_(N), y0_(N), nx0_(N), Vs_(N), grad_phiy0_(N, N), un(N) {}

  };

  static void do_test_F(size_type N) {

    base_matrix dlambdaF(N, N), dnF(N, N), dVsF(N, N);
    base_small_vector F(N), dgF(N);

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
                "the multi contact frame object after a model brick has "
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
    mcf.set_nodes_mode(0);
    mcf.compute_contact_pairs();

    size_type N = mcf.dim();
    base_matrix Melem;
    base_matrix dlambdaF(N, N), dnF(N, N), dVsF(N, N);
    base_small_vector F(N), dgF(N);

    scalar_type FMULT = 1.;

    // Stabilization for non-contact zones
    for (size_type i = 0; i < mcf.nb_boundaries(); ++i)
      if (mcf.is_self_contact() || mcf.is_slave_boundary(i)) {
        size_type region = mcf.region_of_boundary(i);
        const std::string &name_lx = mcf.multname_of_boundary(i);
        GMM_ASSERT1(name_lx.size() > 0, "This brick need "
                    "multipliers defined on the multi_contact_frame object");
        const mesh_fem &mflambda = mcf.mfmult_of_boundary(i);
        const mesh_im &mim = mcf.mim_of_boundary(i);
        const gmm::sub_interval &I = md.interval_of_variable(name_lx);

        if (version & model::BUILD_MATRIX) { // LXLX term
          model_real_sparse_matrix M1(mflambda.nb_dof(), mflambda.nb_dof());
          asm_mass_matrix(M1, mim, mflambda, region);
          gmm::add(gmm::scaled(M1, FMULT/r), gmm::sub_matrix(M, I, I));
        }

        if (version & model::BUILD_RHS) { // LX term
          model_real_plain_vector V1(mflambda.nb_dof());
          asm_source_term
            (V1, mim, mflambda, mflambda,
             md.real_variable(mcf.multname_of_boundary(i)), region);
          gmm::add(gmm::scaled(V1, -FMULT/r), gmm::sub_vector(V, I));
        }
      }

    gauss_point_precomp gpp(N, md, mcf, alpha);

    // do_test_F(2); do_test_F(3);

    base_matrix auxNN1(N, N), auxNN2(N, N);
    base_small_vector auxN1(N), auxN2(N);

    // Iterations on the contact pairs
    for (size_type icp = 0; icp < mcf.nb_contact_pairs(); ++icp) {
      const multi_contact_frame::contact_pair &cp = mcf.get_contact_pair(icp);
      gpp.set_pair(cp);
      const base_small_vector &nx = gpp.nx(), &ny = gpp.ny();
      const mesh_fem *mf_ux = gpp.mf_ux(), *mf_lx = gpp.mf_lx(), *mf_uy(0);
      size_type ndof_ux = gpp.ndof_ux(), ndof_uy(0), ndof_lx = gpp.ndof_lx();
      size_type cvx = gpp.cvx(), cvy(0);
      const gmm::sub_interval &I_ux = gpp.I_ux(), &I_lx = gpp.I_lx();
      gmm::sub_interval I_uy;
      bool isrigid = gpp.isrigid();
      if (!isrigid) {
        ndof_uy = gpp.ndof_uy(); I_uy = gpp.I_uy();
        mf_uy = gpp.mf_uy(); cvy =  gpp.cvy();
      }
      scalar_type weight = gpp.weight(), g = gpp.g();
      const base_small_vector &lambda = gpp.lx();

      base_vector auxUX(ndof_ux), auxUY(ndof_uy);

      if (version & model::BUILD_MATRIX) {

        base_matrix auxUYN(ndof_uy, N);
        base_matrix auxLXN1(ndof_lx, N), auxLXN2(ndof_lx, N);

        aug_friction_grad(lambda, g, gpp.Vs(), nx, r, f_coeff, F, dlambdaF,
                          dgF, dnF, dVsF);


        const base_tensor &vgrad_base_ux = gpp.vgrad_base_ux();
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
        gmm::mult(gpp.vbase_ux(), gmm::transposed(gpp.vbase_lx()), Melem);
        gmm::scale(Melem, -weight);
        mat_elem_assembly(M, I_ux, I_lx, Melem, *mf_ux, cvx, *mf_lx, cvx);
#endif

#ifdef CONSIDER_TERM2

        if (!isrigid) {
          // Term  \delta\lambda(X) . \delta v(Y)
          gmm::resize(Melem, ndof_uy, ndof_lx); gmm::clear(Melem);
          gmm::mult(gpp.vbase_uy(), gmm::transposed(gpp.vbase_lx()), Melem);
          gmm::scale(Melem, weight);
          mat_elem_assembly(M, I_uy, I_lx, Melem, *mf_uy, cvy, *mf_lx, cvx);

          // Term \lambda(X) . (\nabla \delta v(Y) (\nabla phi)^(-1)\delta y
          gmm::clear(auxUYN);
          const base_tensor &vgrad_base_uy = gpp.vgrad_base_uy();
          for (size_type i = 0; i < ndof_uy; ++i)
            for (size_type j = 0; j < N; ++j)
              for (size_type k = 0; k < N; ++k)
                auxUYN(i, j) += lambda[k] * vgrad_base_uy(i, k, j);
          base_matrix lgraddeltavgradphiyinv(ndof_uy, N);
          gmm::mult(auxUYN, gpp.grad_phiy_inv(), lgraddeltavgradphiyinv);

          // first sub term
          gmm::resize(Melem, ndof_uy, ndof_uy); gmm::clear(Melem);
          gmm::mult(lgraddeltavgradphiyinv, gpp.I_nxny(), auxUYN);
          gmm::mult(auxUYN, gmm::transposed(gpp.vbase_uy()), Melem);
          // Caution: re-use of auxUYN in second sub term
          gmm::scale(Melem, -weight);
          mat_elem_assembly(M, I_uy, I_uy, Melem, *mf_uy, cvy, *mf_uy, cvy);

          // Second sub term
          gmm::resize(Melem, ndof_uy, ndof_ux); gmm::clear(Melem);
          // Caution: re-use of auxUYN
          // gmm::mult(lgraddeltavgradphiyinv, gpp.I_nxny(), auxUYN);
          gmm::mult(auxUYN, gmm::transposed(gpp.vbase_ux()), Melem);

          // Third sub term
          base_matrix auxUYUX(ndof_uy, ndof_ux);
          gmm::mult(gpp.I_nxny(), gmm::transposed(gpp.grad_phix_inv()), auxNN1);
          gmm::mult(lgraddeltavgradphiyinv, auxNN1, auxUYN);
          gmm::mult(auxUYN, gmm::transposed(graddeltaunx), auxUYUX);
          gmm::scale(auxUYUX, -g);
          gmm::add(auxUYUX, Melem);
          gmm::scale(Melem, weight);
          mat_elem_assembly(M, I_uy, I_ux, Melem, *mf_uy, cvy, *mf_ux, cvx);
        }

#endif


#ifdef CONSIDER_TERM3

        // LXLX term
        // Term (1/r)(I-dlambdaF)\delta\lambda\delta\mu
        //   the I of (I-dlambdaF) is skipped because globally added before
        gmm::resize(Melem, ndof_lx, ndof_lx); gmm::clear(Melem);
        gmm::copy(gmm::scaled(dlambdaF, scalar_type(-1)/r), auxNN1);
        gmm::mult(gpp.vbase_lx(), auxNN1, auxLXN1);
        gmm::mult(auxLXN1, gmm::transposed(gpp.vbase_lx()), Melem);
        gmm::scale(Melem, weight*FMULT);
        mat_elem_assembly(M, I_lx, I_lx, Melem, *mf_lx, cvx, *mf_lx, cvx);

        // LXUX term
        // Term -(1/r)dnF\delta nx\delta\mu
        gmm::resize(Melem, ndof_lx, ndof_ux); gmm::clear(Melem);
        gmm::mult(gpp.vbase_lx(), dnF, auxLXN1);
        gmm::mult(auxLXN1, gpp.I_nxnx(), auxLXN2);
        gmm::mult(auxLXN2,  gmm::transposed(gpp.grad_phix_inv()), auxLXN1);
        gmm::mult(auxLXN1, gmm::transposed(graddeltaunx), Melem);
        gmm::scale(Melem, scalar_type(1)/r);
        // assembly factorized with the next term

        // Term -(1/r)dgF\delta g\delta\mu
        base_vector deltamudgF(ndof_lx);
        gmm::mult(gpp.vbase_lx(),
                  gmm::scaled(dgF, scalar_type(1)/(r*gpp.nxdotny())),
                  deltamudgF);

        // first sub term
        gmm::mult(gpp.vbase_ux(), ny, auxUX);

        // second sub term
        gmm::mult(gpp.I_nxnx(), gmm::scaled(ny, -g), auxN1);
        gmm::mult(gpp.grad_phix_inv(), auxN1, auxN2);
        gmm::mult_add(graddeltaunx, auxN2, auxUX);    // auxUX -> \delta u(X) - g Dn_x
        gmm::rank_one_update(Melem, deltamudgF,  auxUX);
        gmm::scale(Melem, weight*FMULT);
        mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);

        if (!isrigid) {
          // LXUY term
          // third sub term
          gmm::resize(Melem, ndof_lx, ndof_uy); gmm::clear(Melem);
          gmm::mult(gpp.vbase_uy(), ny, auxUY);       // auxUY -> \delta u(Y)
          gmm::rank_one_update(Melem, deltamudgF, auxUY);
          gmm::scale(Melem, -weight*FMULT);
          mat_elem_assembly(M, I_lx, I_uy, Melem, *mf_lx, cvx, *mf_uy, cvy);
        }

        if (alpha != scalar_type(0)) {
          // Term -(1/r) d_Vs F \delta Vs\delta\mu

          if (!isrigid) {
#ifdef CONSIDER_FRAME_INDIFFERENCE
            base_matrix gphiy0gphiyinv(N, N);
            gmm::mult(gpp.grad_phiy0(), gpp.grad_phiy_inv(), gphiy0gphiyinv);
            gmm::mult(gphiy0gphiyinv, gpp.I_nxny(), auxNN1);
            gmm::rank_one_update(auxNN1, gpp.nx0(),
                                 gmm::scaled(gpp.ny(),scalar_type(1)/gpp.nxdotny()));
            gmm::mult(dVsF, auxNN1, auxNN2);
            // Caution: auxNN2 re-used in the second sub term

            // LXUX term
            // first sub term
            gmm::resize(Melem, ndof_lx, ndof_ux); gmm::clear(Melem);
            gmm::mult(gpp.vbase_lx(), gmm::transposed(auxNN2), auxLXN1);
            // Caution: auxLXN1 re-used in the third sub term
            gmm::mult(auxLXN1, gmm::transposed(gpp.vbase_ux()), Melem);

            // second sub term
            base_matrix auxLXUX(ndof_lx, ndof_ux);
            gmm::mult(auxNN2, gpp.I_nxnx(), auxNN1);
            gmm::mult(auxNN1, gmm::transposed(gpp.grad_phix_inv()), auxNN2);
            gmm::mult(gpp.vbase_lx(), gmm::transposed(auxNN2), auxLXN2);
            gmm::mult(auxLXN2, gmm::transposed(graddeltaunx), auxLXUX);
            gmm::scale(auxLXUX, -g);
            gmm::add(auxLXUX, Melem);
            gmm::scale(Melem, -weight*alpha*FMULT/r);
            mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);

            // LXUY term
            // third sub term
            gmm::resize(Melem, ndof_lx, ndof_uy); gmm::clear(Melem);
            // Caution: auxLXN1 re-used
            gmm::mult(auxLXN1, gmm::transposed(gpp.vbase_uy()), Melem);
            gmm::scale(Melem, weight*alpha*FMULT/r);
            mat_elem_assembly(M, I_lx, I_uy, Melem, *mf_lx, cvx, *mf_uy, cvy);
#else
            base_matrix I_gphiy0gphiyinv(N, N);
            gmm::mult(gmm::scaled(gpp.grad_phiy0(), scalar_type(-1)),
                      gpp.grad_phiy_inv(), I_gphiy0gphiyinv);
            gmm::add(gmm::identity_matrix(), I_gphiy0gphiyinv);

            // LXUX term
            // first sub term
            gmm::resize(Melem, ndof_lx, ndof_ux); gmm::clear(Melem);
            gmm::mult(I_gphiy0gphiyinv, gpp.I_nxny(), auxNN1);
            for (size_type j = 0; j < N; ++j) auxNN1(j,j) -= scalar_type(1);
            gmm::mult(dVsF, auxNN1, auxNN2);
            gmm::mult(gpp.vbase_lx(), gmm::transposed(auxNN2), auxLXN2);
            // Caution: auxLXN2 re-used in the third sub term
            gmm::mult(auxLXN2, gmm::transposed(gpp.vbase_ux()), Melem);

            // second sub term
            base_matrix auxLXUX(ndof_lx, ndof_ux);
            gmm::mult(dVsF, I_gphiy0gphiyinv, auxNN1);
            gmm::mult(auxNN1, gpp.I_nxny(), auxNN2);
            gmm::mult(auxNN2, gmm::transposed(gpp.grad_phix_inv()), auxNN1);
            gmm::mult(gpp.vbase_lx(), gmm::transposed(auxNN1), auxLXN1);
            gmm::mult(auxLXN1, gmm::transposed(graddeltaunx), auxLXUX);
            gmm::scale(auxLXUX, -g);
            gmm::add(auxLXUX, Melem);
            gmm::scale(Melem, weight*alpha*FMULT/r);
            mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);

            // LXUY term
            // third sub term
            gmm::resize(Melem, ndof_lx, ndof_uy); gmm::clear(Melem);
            // Caution: auxLXN2 re-used
            gmm::mult(auxLXN2, gmm::transposed(gpp.vbase_uy()), Melem);
            gmm::scale(Melem, -weight*alpha*FMULT/r);
            mat_elem_assembly(M, I_lx, I_uy, Melem, *mf_lx, cvx, *mf_uy, cvy);
#endif
          } else {
            // LXUX term
            gmm::mult(gpp.vbase_lx(), gmm::transposed(dVsF), auxLXN1);
            gmm::mult(auxLXN1, gmm::transposed(gpp.vbase_ux()), Melem);
            gmm::scale(Melem, -weight*alpha*FMULT/r);
            mat_elem_assembly(M, I_lx, I_ux, Melem, *mf_lx, cvx, *mf_ux, cvx);
          }
        }
#endif
      }

      if (version & model::BUILD_RHS) {

        if (!(version & model::BUILD_MATRIX))
          aug_friction(lambda, g, gpp.Vs(), nx, r, f_coeff, F);

#ifdef CONSIDER_TERM1

        // Term lambda.\delta v(X)
        gmm::mult(gpp.vbase_ux(), lambda, auxUX);
        gmm::scale(auxUX, weight);
        vec_elem_assembly(V, I_ux, auxUX, *mf_ux, cvx);
#endif

#ifdef CONSIDER_TERM2

        // Term -lambda.\delta v(Y)
        if (!isrigid) {
          gmm::mult(gpp.vbase_uy(), lambda, auxUY);
          gmm::scale(auxUY, -weight);
          vec_elem_assembly(V, I_uy, auxUY, *mf_uy, cvy);
        }
#endif

#ifdef CONSIDER_TERM3

        // Term -(1/r)(lambda - F).\delta \mu
        // (1/r)(lambda).\delta \mu is skipped because globally added before
        base_vector auxLX(ndof_lx);
        gmm::mult(gpp.vbase_lx(), gmm::scaled(F, weight*FMULT/r), auxLX);
        vec_elem_assembly(V, I_lx, auxLX, *mf_lx, cvx);
#endif
      }

    }
  }


  size_type add_integral_large_sliding_contact_brick_raytrace
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
  //  Large sliding brick with field extension principle.
  //  Deprecated. To be adapated to the high-level generic assembly
  //
  //=========================================================================

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
      model_real_plain_vector ext_U(mfu.nb_basic_dof()); // means that the structure has to be build each time ... to be changed. ATTENTION : la meme variable ne doit pas etre etendue dans deux vecteurs differents.
      mfu.extend_vector(U, ext_U);
      ext_Us.push_back(ext_U);
      return i;
    }

    size_type add_lambda(const getfem::mesh_fem &mfl,
                         const model_real_plain_vector &l) {
      size_type i = 0;
      for (; i < lambdas.size(); ++i) if (lambdas[i] == &l) return i;
      lambdas.push_back(&l);
      model_real_plain_vector ext_l(mfl.nb_basic_dof()); // means that the structure has to be build each time ... to be changed. ATTENTION : la meme variable ne doit pas etre etendue dans deux vecteurs differents.
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
        fem_interpolation_context ctx(pgt, pfp, size_type(-1), G, cv,
                                      short_type(-1));

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
    scalar_type J = gmm::lu_inverse(gradinv); // remplacer par une resolution...
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
                                       short_type(-1));

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

    // Eviter les calculs inutiles dans le cas state == 2 ... a voir a la fin
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
        // precedent. Attendre la version avec frottement.
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
            gmm::lu_inverse(gradinv_y0); // a proteger contre la non-inversibilite
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
        fem_interpolation_context ctxu(pgt, pfpu, size_type(-1), G, cv, v.f());
        fem_interpolation_context ctxl(pgt, pfpl, size_type(-1), G, cv, v.f());

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


  // r ne peut pas etre variable pour le moment.
  // dataname_friction_coeff ne peut pas etre variable non plus ...

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



  // ----------------------------------------------------------------------
  //
  // Brick for large sliding contact with friction using raytracing contact
  // detection and the high-level generic assemly
  //
  // ----------------------------------------------------------------------

  struct intergral_large_sliding_contact_brick_raytracing
    : public virtual_brick {

    struct contact_boundary {
      size_type region;
      std::string varname_u;
      std::string varname_lambda;
      std::string varname_w;
      bool is_master;
      bool is_slave;
      const mesh_im *mim;
      
      std::string expr;
    };

    std::vector<contact_boundary> boundaries;
    std::string transformation_name;
    std::string u_group;
    std::string w_group;
    std::string friction_coeff;
    std::string alpha;
    std::string augmentation_param;
    model::varnamelist vl, dl;
    model::mimlist ml;

    bool sym_version;

    void add_contact_boundary(model &md, const mesh_im &mim, size_type region,
                              bool is_master, bool is_slave,
                              const std::string &u,
                              const std::string &lambda,
                              const std::string &w = "") {
      std::string test_u = "Test_" + sup_previous_to_varname(u);
      std::string test_u_group = "Test_" + sup_previous_to_varname(u_group);
      std::string test_lambda = "Test_" + sup_previous_to_varname(lambda);
      GMM_ASSERT1(is_slave || is_master, "The contact boundary should be "
                  "either master, slave or both");
      const mesh_fem *mf = md.pmesh_fem_of_variable(u);
      GMM_ASSERT1(mf, "The displacement variable should be a f.e.m. one");
      GMM_ASSERT1(&(mf->linked_mesh()) == &(mim.linked_mesh()),
                  "The displacement variable and the integration method "
                  "should share the same mesh");
      if (is_slave) {
        const mesh_fem *mf_l = md.pmesh_fem_of_variable(lambda);
        GMM_ASSERT1(mf, "The multiplier variable should be a f.e.m. one");
        GMM_ASSERT1(&(mf_l->linked_mesh()) == &(mim.linked_mesh()),
                    "The displacement variable and the multiplier one "
                    "should share the same mesh");
      }

      if (w.size()) {
        const mesh_fem *mf2 =  md.pmesh_fem_of_variable(w);
        GMM_ASSERT1(!mf2 || &(mf2->linked_mesh()) == &(mf->linked_mesh()),
                    "The data for the sliding velocity should be defined on "
                    " the same mesh as the displacement variable");
      }

      for (size_type i = 0; i < boundaries.size(); ++i) {
        const contact_boundary &cb = boundaries[i];
        if (&(md.mesh_fem_of_variable(cb.varname_u).linked_mesh())
            == &(mf->linked_mesh()) && cb.region == region)
          GMM_ASSERT1(false, "This contact boundary has already been added");
      }
      if (is_master)
        add_master_contact_boundary_to_raytracing_transformation
          (md, transformation_name, mf->linked_mesh(), u_group, region);
      else
         add_slave_contact_boundary_to_raytracing_transformation
          (md, transformation_name, mf->linked_mesh(), u_group, region);
      
      boundaries.push_back(contact_boundary());
      contact_boundary &cb = boundaries.back();
      cb.region = region;
      cb.varname_u = u;
      if (is_slave) cb.varname_lambda = lambda;
      cb.varname_w = w;
      cb.is_master = is_master;
      cb.is_slave = is_slave;
      cb.mim = &mim;
      if (is_slave) {
        // Coulomb_friction_coupled_projection(lambda,
        //    Transformed_unit_vector(Grad_u, Normal),
        //    (u-Interpolate(ug,trans)-(w-Interpolate(wg,trans)))*alpha,
        //    (Interpolate(X,trans)+Interpolate(ug,trans)-X-u).
        //      Transformed_unit_vector(Grad_u, Normal), f, r)
        std::string coupled_projection_def =
          "Coulomb_friction_coupled_projection("
          + lambda+", Transformed_unit_vector(Grad_"+u+", Normal),"
          + "("+u+"-Interpolate("+u_group+","+transformation_name+")"
          + (w.size() 
             ? ("-"+w+"+Interpolate("+w_group+","+transformation_name+")")
             : "")
          +")*"+alpha+","
          + "(Interpolate(X,"+transformation_name+")+Interpolate("+u_group+","
          + transformation_name+")-X-"+u+").Transformed_unit_vector(Grad_"
          + u+", Normal),"+friction_coeff+","+augmentation_param+")";

        // Coulomb_friction_coupled_projection(lambda,
        //   Transformed_unit_vector(Grad_u, Normal), (u-w)*alpha,
        //   (Interpolate(X,trans)-X-u).Transformed_unit_vector(Grad_u,
        //                                                      Normal), f, r)
        std::string coupled_projection_rig =
          "Coulomb_friction_coupled_projection("
          + lambda+", Transformed_unit_vector(Grad_"+u+", Normal),"
          + "("+u+(w.size() ? ("-"+w):"")+")*"+alpha
          + ",(Interpolate(X,"+transformation_name+")-X-"+u
          + ").Transformed_unit_vector(Grad_"+u+", Normal),"
          + friction_coeff+","+augmentation_param+")";

        cb.expr =
          // -lambda.Test_u for non-symmetric version
          (sym_version ? "" : ("-"+lambda+"." + test_u))
          // -coupled_projection_def.Test_u and -coupled_projection_rig.Test_u
          // for symmetric version
          + (sym_version ? ("+ Interpolate_filter("+transformation_name+",-"
                            +coupled_projection_def+"."+test_u+",1)") : "")
          + (sym_version ? ("+ Interpolate_filter("+transformation_name+",-"
                            +coupled_projection_rig+"."+test_u+",2)") : "")
          // Interpolate_filter(trans,
          //                   lambda.Interpolate(Test_ug, contact_trans), 1)
          // or
          // Interpolate_filter(trans,
          //     coupled_projection_def.Interpolate(Test_ug, contact_trans), 1)
          + "+ Interpolate_filter("+transformation_name+","
          + (sym_version ? coupled_projection_def : lambda)
          + ".Interpolate("+test_u_group+"," + transformation_name+"), 1)"
          // -(1/r)*lambda.Test_lambda
          + "-(1/"+augmentation_param+")*"+lambda+"."+test_lambda
          // Interpolate_filter(trans,
          //   (1/r)*coupled_projection_rig.Test_lambda, 2)
          + "+ Interpolate_filter("+transformation_name+","
          + "(1/"+augmentation_param+")*"+ coupled_projection_rig
          + "."+test_lambda+", 2)"
          // Interpolate_filter(trans,
          //   (1/r)*coupled_projection_def.Test_lambda, 1)
          + "+ Interpolate_filter("+transformation_name+","
          + "(1/"+augmentation_param+")*" + coupled_projection_def + "."
          + test_lambda+", 1)";
      }
    }

    virtual void asm_real_tangent_terms(const model &md, size_type ,
                                        const model::varnamelist &,
                                        const model::varnamelist &,
                                        const model::mimlist &,
                                        model::real_matlist &,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type,
                                        build_version) const {
      // GMM_ASSERT1(mims.size() == 1,
      //            "Generic linear assembly brick needs one and only one "
      //            "mesh_im"); // to be verified ...

      for (size_type i = 0; i < boundaries.size(); ++i) {
        const contact_boundary &cb = boundaries[i];
        if (cb.is_slave)
          md.add_generic_expression(cb.expr, *(cb.mim), cb.region);
      }
    }

    virtual scalar_type asm_real_pseudo_potential(const model &, size_type,
                                                  const model::varnamelist &,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &,
                                                  model::real_veclist &,
                                                  model::real_veclist &,
                                                  size_type) const {
      GMM_WARNING1("Brick " << name << " has no contribution to "
                   << "the pseudo potential !");
      return 0;
    }


    intergral_large_sliding_contact_brick_raytracing
    (const std::string &r,
     const std::string &f_coeff,const std::string &ug,
     const std::string &wg, const std::string &tr,
     const std::string &alpha_ = "1", bool sym_v = false) {
      sym_version = sym_v;
      transformation_name = tr;
      u_group = ug; w_group = wg;
      friction_coeff = f_coeff;
      alpha = alpha_;
      augmentation_param = r;
      
      set_flags("Integral large sliding contact bick raytracing",
                false /* is linear*/,
                false /* is symmetric */, false /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  
  const std::string &transformation_name_of_large_sliding_contact_brick
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    intergral_large_sliding_contact_brick_raytracing *p
      = dynamic_cast<intergral_large_sliding_contact_brick_raytracing *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->transformation_name;
  }

  const std::string &displacement_group_name_of_large_sliding_contact_brick
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    intergral_large_sliding_contact_brick_raytracing *p
      = dynamic_cast<intergral_large_sliding_contact_brick_raytracing *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->u_group;
  }

  const std::string &sliding_data_group_name_of_large_sliding_contact_brick
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    intergral_large_sliding_contact_brick_raytracing *p
      = dynamic_cast<intergral_large_sliding_contact_brick_raytracing *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->w_group;
  }

  void add_rigid_obstacle_to_large_sliding_contact_brick
  (model &md, size_type indbrick, std::string expr, size_type N) {
    pbrick pbr = md.brick_pointer(indbrick);
     intergral_large_sliding_contact_brick_raytracing *p
       = dynamic_cast<intergral_large_sliding_contact_brick_raytracing *>
       (const_cast<virtual_brick *>(pbr.get()));
     GMM_ASSERT1(p, "Wrong type of brick");
     add_rigid_obstacle_to_raytracing_transformation
       (md, p->transformation_name, expr, N);
  }
  
  void add_contact_boundary_to_large_sliding_contact_brick
  (model &md, size_type indbrick, const mesh_im &mim, size_type region,
   bool is_master, bool is_slave, const std::string &u,
   const std::string &lambda, const std::string &w) {
     pbrick pbr = md.brick_pointer(indbrick);
     intergral_large_sliding_contact_brick_raytracing *p
       = dynamic_cast<intergral_large_sliding_contact_brick_raytracing *>
       (const_cast<virtual_brick *>(pbr.get()));
     GMM_ASSERT1(p, "Wrong type of brick");
     
     bool found_u = false, found_lambda = false;
     for (size_type i = 0; i < p->vl.size(); ++i) {
       if (p->vl[i].compare(u) == 0) found_u = true;
       if (p->vl[i].compare(lambda) == 0) found_lambda = true;
     }
     if (!found_u) p->vl.push_back(u);
     GMM_ASSERT1(!is_slave || lambda.size(),
                 "You should define a multiplier on each slave boundary");
     if (is_slave && !found_lambda) p->vl.push_back(lambda);
     if (!found_u || (is_slave && !found_lambda))
       md.change_variables_of_brick(indbrick, p->vl);

     std::vector<std::string> ug = md.variable_group(p->u_group);
     found_u = false;
     for (size_type i = 0; i < ug.size(); ++i)
         if (ug[i].compare(u) == 0) found_u = true;
     if (!found_u) {
       ug.push_back(u);
       md.define_variable_group(p->u_group, ug);
     }

     if (w.size()) {
       bool found_w = false;
       for (size_type i = 0; i < p->dl.size(); ++i)
         if (p->dl[i].compare(w) == 0) found_w = true;
       if (!found_w) { 
         p->dl.push_back(w);
         md.change_data_of_brick(indbrick, p->dl);
       }
       std::vector<std::string> wg = md.variable_group(p->w_group);
       found_w = false;
       for (size_type i = 0; i < wg.size(); ++i)
         if (wg[i].compare(w) == 0) found_w = true;
       if (!found_w) {
         wg.push_back(w);
         md.define_variable_group(p->w_group, wg);
       }
     }
     
     bool found_mim = false;
     for (size_type i = 0; i < p->ml.size(); ++i)
       if (p->ml[i] == &mim) found_mim = true;
     if (!found_mim) {
       p->ml.push_back(&mim);
       md.change_mims_of_brick(indbrick, p->ml);
     }

     p->add_contact_boundary(md, mim, region, is_master, is_slave, u,lambda,w);
  } 

  size_type add_integral_large_sliding_contact_brick_raytracing
  (model &md, const std::string &augm_param,
   scalar_type release_distance, const std::string &f_coeff, 
   const std::string &alpha, bool sym_v) {

    char ugroupname[50], wgroupname[50], transname[50];
    for (int i = 0; i < 10000; ++i) {
      sprintf(ugroupname, "disp__group_raytracing_%d", i);
      if (!(md.variable_group_exists(ugroupname))
          && !(md.variable_exists(ugroupname)))
        break;
    }
    md.define_variable_group(ugroupname, std::vector<std::string>());

    for (int i = 0; i < 10000; ++i) {
      sprintf(wgroupname, "w__group_raytracing_%d", i);
      if (!(md.variable_group_exists(wgroupname))
          && !(md.variable_exists(wgroupname)))
        break;
    }
    md.define_variable_group(wgroupname, std::vector<std::string>());

    for (int i = 0; i < 10000; ++i) {
      sprintf(transname, "trans__raytracing_%d", i);
      if (!(md.interpolate_transformation_exists(transname)))
        break;
    }
    add_raytracing_transformation(md, transname, release_distance);

    model::varnamelist vl, dl;
    if (md.variable_exists(augm_param)) dl.push_back(augm_param);
    if (md.variable_exists(f_coeff)) dl.push_back(f_coeff);
    if (md.variable_exists(alpha)) dl.push_back(alpha);
    
    intergral_large_sliding_contact_brick_raytracing *p
      = new intergral_large_sliding_contact_brick_raytracing
      (augm_param, f_coeff, ugroupname, wgroupname, transname, alpha, sym_v);
    pbrick pbr = p;
    p->dl = dl;

    return md.add_brick(pbr, p->vl, p->dl, model::termlist(),model::mimlist(),
                        size_type(-1));
  }











}  /* end of namespace getfem.                                             */
