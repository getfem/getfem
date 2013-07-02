/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2011-2013 Yves Renard, Konstantinos Poulios.

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
#include "getfem/getfem_projected_fem.h"
#include "gmm/gmm_condition_number.h"

#include <getfem/getfem_arch_config.h>
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif

namespace getfem {



  template <typename T> inline static T Heav(T a)
  { return (a < T(0)) ? T(0) : T(1); }


  //=========================================================================
  //
  //  Basic non linear term used for contact bricks.
  //
  //=========================================================================

  void contact_nonlinear_term::adjust_tensor_size(void) {
    sizes_.resize(1); sizes_[0] = 1;
    switch (option) {
      // one-dimensional tensors [N]
    case RHS_U_V1:       case RHS_U_V2:       case RHS_U_V4:
    case RHS_U_V5:       case RHS_U_FRICT_V6: case RHS_U_FRICT_V7:
    case RHS_U_FRICT_V8: case RHS_U_FRICT_V1:
    case RHS_U_FRICT_V4: case RHS_U_FRICT_V5:
    case RHS_L_FRICT_V1: case RHS_L_FRICT_V2: case RHS_L_FRICT_V4:
    case K_UL_V1:        case K_UL_V2:        case K_UL_V3:
    case UZAWA_PROJ_FRICT: case UZAWA_PROJ_FRICT_SAXCE:
      sizes_[0] = N; break;
      // two-dimensional tensors [N x N]
    case K_UU_V1: case K_UU_V2:
    case K_UL_FRICT_V1: case K_UL_FRICT_V2: case K_UL_FRICT_V3:
    case K_UL_FRICT_V4: case K_UL_FRICT_V5:
    case K_UL_FRICT_V7: case K_UL_FRICT_V8:
    case K_LL_FRICT_V1: case K_LL_FRICT_V2: case K_LL_FRICT_V4:
    case K_UU_FRICT_V1: case K_UU_FRICT_V2:
    case K_UU_FRICT_V3: case K_UU_FRICT_V4: case K_UU_FRICT_V5:
      sizes_.resize(2); sizes_[0] = sizes_[1] = N;  break;
    }

    // adjust temporary variables sizes
    lnt.resize(N); lt.resize(N); zt.resize(N); no.resize(N);
    aux1.resize(1); auxN.resize(N); V.resize(N);
    gmm::resize(GP, N, N);
  }

  void contact_nonlinear_term::friction_law
  (scalar_type p, scalar_type &tau) {
    tau = (p > scalar_type(0)) ? tau_adh + f_coeff * p : scalar_type(0);
    if (tau > tresca_lim) tau = tresca_lim;
  }

  void contact_nonlinear_term::friction_law
  (scalar_type p, scalar_type &tau, scalar_type &tau_grad) {
    if (p <= scalar_type(0)) {
      tau = scalar_type(0);
      tau_grad = scalar_type(0);
    }
    else {
      tau = tau_adh + f_coeff * p;
      if (tau > tresca_lim) {
        tau = tresca_lim;
        tau_grad = scalar_type(0);
      }
      else
        tau_grad = f_coeff;
    }
  }

  void contact_nonlinear_term::compute
  (fem_interpolation_context &/* ctx */, bgeot::base_tensor &t) {

    t.adjust_sizes(sizes_);
    scalar_type e, augm_ln, rho, rho_grad;
    dim_type i, j;
    bool coulomb;

    switch (option) {

    // scalar tensors [1]

    case RHS_L_V1:
      t[0] = (ln+gmm::neg(ln-r*(un - g)))/r; break;
    case RHS_L_V2:
      t[0] = (un-g) + gmm::pos(ln)/r; break;

    case K_LL_V1:
      t[0] = (Heav(r*(un-g)-ln) - scalar_type(1))/r; break;
    case K_LL_V2:
      t[0] = -Heav(ln)/r; break;

    case UZAWA_PROJ:
      t[0] = -gmm::neg(ln - r*(un - g)); break;

    case CONTACT_FLAG:
      // here ln is expected to be a threshold value expressing a penetration
      // (positive value) or separation (negative value) distance
      t[0] = Heav(un-g - ln);  break;
    case CONTACT_PRESSURE:
      t[0] = -ln;  break;

    // one-dimensional tensors [N]

    case RHS_U_V1:
      for (i=0; i<N; ++i) t[i] = ln * no[i]; break;
    case RHS_U_V2:
      e = -gmm::neg(ln-r*(un - g));
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
    case RHS_U_FRICT_V6:
      e = gmm::neg(ln-r*(un - g));
      friction_law(e, rho);
      auxN = lt - zt;  ball_projection(auxN, rho);
      for (i=0; i<N; ++i) t[i] = auxN[i] - e*no[i];
      break;
    case RHS_U_FRICT_V7:
      e = gmm::neg(-r*(un - g));
      friction_law(e, rho);
      auxN = - zt;  ball_projection(auxN, rho);
      for (i=0; i<N; ++i) t[i] = auxN[i] - e*no[i];
      break;
    case RHS_U_FRICT_V8: // ignores friction_law, assumes pure Coulomb friction
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
      break;
    case RHS_U_FRICT_V1:
      for (i=0; i<N; ++i) t[i] = lnt[i]; break;
    case RHS_U_FRICT_V4:
      e = gmm::neg(ln);
      // if (e > 0. && ctx.xreal()[1] > 1.)
      //        cout << "x = " << ctx.xreal() << " e = " << e << endl;
      friction_law(e, rho);
      auxN = lt;  ball_projection(auxN, rho);
      // if (gmm::vect_norm2(auxN) > 0. && ctx.xreal()[1] > 1.)
      //        cout << "x = " << ctx.xreal() << " auxN = " << auxN << endl;
      for (i=0; i<N; ++i) t[i] = auxN[i] - e*no[i];
      break;
    case RHS_U_FRICT_V5: // ignores friction_law, assumes pure Coulomb friction
      auxN = lnt; De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
      break;
    case RHS_L_FRICT_V1:
      e = gmm::neg(ln-r*(un-g));
      friction_law(e, rho);
      auxN = zt - lt;  ball_projection(auxN, rho); auxN += lt;
      for (i=0; i<N; ++i) t[i] = ((e+ln)*no[i] + auxN[i])/ r;
      break;
    case RHS_L_FRICT_V2:
      e = r*(un-g) + gmm::pos(ln);
      friction_law(gmm::neg(ln), rho);
      auxN = lt;  ball_projection(auxN, rho);
      for (i=0; i<N; ++i) t[i] = (no[i]*e + zt[i] + lt[i] - auxN[i])/r;
      break;
    case RHS_L_FRICT_V4: // ignores friction_law, assumes pure Coulomb friction
      auxN = lnt;
      De_Saxce_projection(auxN, no, f_coeff);
      auxN -= lnt + (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no + zt;
      for (i=0; i<N; ++i) t[i] = -auxN[i]/r;
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
    case UZAWA_PROJ_FRICT:
      e = gmm::neg(ln - r*(un - g));
      friction_law(e, rho);
      auxN = lt - zt;  ball_projection(auxN, rho);
      for (i=0; i<N; ++i) t[i] = auxN[i] - e*no[i];
      break;
    case UZAWA_PROJ_FRICT_SAXCE: // ignores friction_law, assumes pure Coulomb friction
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
      break;

    // two-dimensional tensors [N x N]

    case K_UU_V1:
      e = r * Heav(un - g);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = e * no[i] * no[j];
      break;
    case K_UU_V2:
      e = r * Heav(r*(un - g)-ln);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = e * no[i] * no[j];
      break;

    case K_UL_FRICT_V1:
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = ((i == j) ? -scalar_type(1) : scalar_type(0));
      break;
    case K_UL_FRICT_V2:
      friction_law(gmm::neg(ln), rho, rho_grad);
      ball_projection_grad(lt, rho, GP);
      e = gmm::vect_sp(GP, no, no) - Heav(-ln);
      coulomb = (rho_grad > 0) && bool(Heav(-ln));
      if (coulomb) ball_projection_grad_r(lt, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - GP(i,j) +
                   (coulomb ? rho_grad*no[i]*V[j] : scalar_type(0));
      break;
    case K_UL_FRICT_V3:
      augm_ln = ln - r*(un-g);
      friction_law(gmm::neg(augm_ln), rho, rho_grad);
      auxN = lt - zt;
      ball_projection_grad(auxN, rho, GP);
      e = gmm::vect_sp(GP, no, no) - Heav(-augm_ln);
      coulomb = (rho_grad > 0) && bool(Heav(-augm_ln));
      if (coulomb) ball_projection_grad_r(auxN, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - GP(i,j) +
                   (coulomb ? rho_grad*no[i]*V[j] : scalar_type(0));
      break;
    case K_UL_FRICT_V4:
      augm_ln = ln - r*(un-g);
      friction_law(gmm::neg(augm_ln), rho, rho_grad);
      auxN = lt - zt;
      ball_projection_grad(auxN, rho, GP); gmm::scale(GP, alpha);
      e = gmm::vect_sp(GP, no, no) - Heav(-augm_ln);
      coulomb = (rho_grad > 0) && bool(Heav(-augm_ln));
      if (coulomb) ball_projection_grad_r(auxN, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - GP(i,j) +
                   (coulomb ? rho_grad*V[i]*no[j] : scalar_type(0));
      break;
    case K_UL_FRICT_V5:
      e = alpha - scalar_type(1);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - ((i == j) ? alpha : scalar_type(0));
      break;
    case K_UL_FRICT_V7: // ignores friction_law, assumes pure Coulomb friction
      De_Saxce_projection_grad(lnt, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = -GP(j,i);
      break;
    case K_UL_FRICT_V8: // ignores friction_law, assumes pure Coulomb friction
      {
        scalar_type nzt = gmm::vect_norm2(zt);
        gmm::copy(gmm::identity_matrix(), GP); gmm::scale(GP, alpha);
        gmm::rank_one_update(GP, gmm::scaled(no, scalar_type(1)-alpha), no);
        if (nzt != scalar_type(0))
          gmm::rank_one_update(GP, gmm::scaled(no, -f_coeff*alpha/nzt), zt);
        for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = - GP(i,j);
      }
      break;
    case K_LL_FRICT_V1:
      augm_ln = ln - r*(un-g);
      friction_law(gmm::neg(augm_ln), rho, rho_grad);
      auxN = lt - zt;
      ball_projection_grad(auxN, rho, GP);
      e = Heav(-augm_ln) - gmm::vect_sp(GP, no, no);
      coulomb = (rho_grad > 0) && bool(Heav(-augm_ln));
      if (coulomb) ball_projection_grad_r(auxN, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = (no[i]*no[j]*e + GP(i,j)
                    - ((i == j) ? scalar_type(1) : scalar_type(0))
                    - (coulomb ? rho_grad*no[i]*V[j] : scalar_type(0))) / r;
      break;
    case K_LL_FRICT_V2:
      friction_law(gmm::neg(ln), rho, rho_grad);
      ball_projection_grad(lt, rho, GP);
      e = Heav(-ln) - gmm::vect_sp(GP, no, no);
      coulomb = (rho_grad > 0) && bool(Heav(-ln));
      if (coulomb) ball_projection_grad_r(lt, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = (no[i]*no[j]*e + GP(i,j)
                    - ((i == j) ? scalar_type(1) : scalar_type(0))
                    - (coulomb ? rho_grad*no[i]*V[j] : scalar_type(0))) / r;
      break;
    case K_LL_FRICT_V4: // ignores friction_law, assumes pure Coulomb friction
      De_Saxce_projection_grad(lnt, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = (GP(i,j) - ((i == j) ? scalar_type(1) : scalar_type(0)))/r;
      break;
    case K_UU_FRICT_V1:
      e = r * Heav(r*(un-g)-ln);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = no[i]*no[j]*e;
      break;
    case K_UU_FRICT_V2:
      friction_law(-ln, rho, rho_grad);
      auxN = lt - zt;
      ball_projection_grad(auxN, rho, GP); gmm::scale(GP, alpha);
      e = Heav(r*(un-g)-ln) - gmm::vect_sp(GP, no, no);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = r*(no[i]*no[j]*e + GP(i,j));
      break;
    case K_UU_FRICT_V3:
    case K_UU_FRICT_V4:
      augm_ln = (option == K_UU_FRICT_V3) ? ln - r*(un-g) : - r*(un-g);
      auxN = (option == K_UU_FRICT_V3) ? lt - zt : -zt;
      friction_law(gmm::neg(augm_ln), rho, rho_grad);
      ball_projection_grad(auxN, rho, GP); gmm::scale(GP, alpha);
      e = Heav(-augm_ln) - gmm::vect_sp(GP, no, no);
      coulomb = (rho_grad > 0) && bool(Heav(-augm_ln));
      if (coulomb) ball_projection_grad_r(auxN, rho, V);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = r*(no[i]*no[j]*e + GP(i,j)
                      - (coulomb ? rho_grad*no[i]*V[j] : scalar_type(0)));
      break;
    case K_UU_FRICT_V5: // ignores friction_law, assumes pure Coulomb friction
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
        for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = r*GP(j,i);
      }
      break;
    default : GMM_ASSERT1(false, "Invalid option");
    }
  }


  //=========================================================================
  //
  //  Non linear term used for contact with rigid obstacle bricks.
  //
  //=========================================================================

  void contact_rigid_obstacle_nonlinear_term::prepare
  (fem_interpolation_context& ctx, size_type nb) {
    size_type cv = ctx.convex_num();

    switch (nb) { // last is computed first
    case 1 : // calculate [un] and [zt] interpolating [U],[WT],[VT] on [mf_u]
      slice_vector_on_basic_dof_of_element(mf_u, U, cv, coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no);
      if (!contact_only) {
        if (gmm::vect_size(WT) == gmm::vect_size(U)) {
          slice_vector_on_basic_dof_of_element(mf_u, WT, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, auxN, N);
          auxN -= gmm::vect_sp(auxN, no) * no;
          if (gmm::vect_size(VT) == gmm::vect_size(U)) {
            slice_vector_on_basic_dof_of_element(mf_u, VT, cv, coeff);
            ctx.pf()->interpolation(ctx, coeff, vt, N);
            vt -= gmm::vect_sp(vt, no) * no;
            // zt = r*(alpha*(u_T-w_T) + (1-gamma)*v_T)
            zt = (((V - un * no) - auxN) * alpha + vt * (1 - gamma)) * r;
          } else {
            // zt = r*alpha*(u_T-w_T)
            zt = ((V - un * no) - auxN) * (r * alpha);
          }
        } else {
          // zt = r*alpha*u_T
          zt = (V - un * no) * (r * alpha);
        }
      }
      break;

    case 2 : // calculate [g] and [no] interpolating [obs] on [mf_obs]
             // calculate [ln] and [lt] from [lnt] and [no]
      slice_vector_on_basic_dof_of_element(mf_obs, obs, cv, coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, grad, 1);
      gmm::copy(gmm::mat_row(grad, 0), no);
      no /= -gmm::vect_norm2(no);
      ctx.pf()->interpolation(ctx, coeff, aux1, 1);
      g = aux1[0];

      if (!contact_only && pmf_lambda) {
        ln = gmm::vect_sp(lnt, no);
        lt = lnt - ln * no;
      }

      break;

    case 3 : // calculate [ln] or [lnt] interpolating [lambda] on [mf_lambda]
      if (pmf_lambda) {
        slice_vector_on_basic_dof_of_element(*pmf_lambda, lambda, cv, coeff);
        if (contact_only) {
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          ln = aux1[0];
        } else {
          ctx.pf()->interpolation(ctx, coeff, lnt, N);
        }
      }
      break;

    case 4 : // calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
             // calculate [tau_adh] interpolating [tau_adhesion] on [mf_coeff]
             // calculate [tresca_lim] interpolating [tresca_limit] on [mf_coeff]
      GMM_ASSERT1(!contact_only, "Invalid friction option");
      if (pmf_coeff) {
        slice_vector_on_basic_dof_of_element(*pmf_coeff, friction_coeff, cv, coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
        if (gmm::vect_size(tau_adhesion)) {
          slice_vector_on_basic_dof_of_element(*pmf_coeff, tau_adhesion, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          tau_adh = aux1[0];
          if (gmm::vect_size(tresca_limit)) {
            slice_vector_on_basic_dof_of_element(*pmf_coeff, tresca_limit, cv, coeff);
            ctx.pf()->interpolation(ctx, coeff, aux1, 1);
            tresca_lim = aux1[0];
          }
        }
      }
      break;

    default : GMM_ASSERT1(false, "Invalid option");
    }

  }


  //=========================================================================
  //
  //  Non linear term used for contact between non-matching meshes bricks.
  //
  //=========================================================================

  void contact_nonmatching_meshes_nonlinear_term::prepare
  (fem_interpolation_context& ctx, size_type nb) {

    size_type cv = ctx.convex_num();

    // - this method is called for nb=4,3,2,1 corresponding to
    //   NonLin(#i0,#i1,#i2,#i3,#i4) before the compute() method is called
    //   for i0
    // - in each call ctx.pf() corresponds to the fem of the cv convex
    //   on the i4,i3,i2 and i1 mesh_fem respectively
    // - this method expects that i1,i2,i3 and i4 mesh_fems will correspond
    //   to mf_u1, mf_u2, mf_lambda and mf_coeff

    switch (nb) { // last is computed first
    case 1 : // calculate [un] and [zt] interpolating [U1],[WT1] on [mf_u1]
             // and subtracting [un] and [zt] calculated on [mf_u2]
      slice_vector_on_basic_dof_of_element(mf_u1, U1, cv, coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      {
        scalar_type un1 = gmm::vect_sp(V, no);
        if (!contact_only) {
          if (gmm::vect_size(WT1) == gmm::vect_size(U1)) {
            slice_vector_on_basic_dof_of_element(mf_u1, WT1, cv, coeff);
            ctx.pf()->interpolation(ctx, coeff, auxN, N);
            auxN -= gmm::vect_sp(auxN, no) * no;
            zt = ((V - un1 * no) - auxN) * (r * alpha) - zt; // = zt1 - zt2 , with zt = r*alpha*(u_T-w_T)
          } else {
            zt = (V - un1 * no) * (r * alpha) - zt;          // = zt1 - zt2 , with zt = r*alpha*u_T
          }
        }
        un = un1 - un; // = un1 - un2
      }
      break;

    case 2 : // calculate [g] and [no]
             // calculate [ln] and [lt] from [lnt] and [no]
             // calculate [un] and [zt] interpolating [U2],[WT2] on [mf_u2]
      {
        const projected_fem &pfe = dynamic_cast<const projected_fem&>(*ctx.pf());
        pfe.projection_data(ctx, no, g);
        gmm::scale(no, scalar_type(-1)); // pointing outwards from mf_u1
      }

      if (!contact_only && pmf_lambda) {
        ln = gmm::vect_sp(lnt, no);
        lt = lnt - ln * no;
      }

      slice_vector_on_basic_dof_of_element(mf_u2, U2, cv, coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no);
      if (!contact_only) {
        if (gmm::vect_size(WT2) == gmm::vect_size(U2)) {
          slice_vector_on_basic_dof_of_element(mf_u2, WT2, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, auxN, N);
          auxN -= gmm::vect_sp(auxN, no) * no;
          zt = ((V - un * no) - auxN) * (r * alpha); // zt = r*alpha*(u_T-w_T)
        } else {
          zt = (V - un * no) * (r * alpha);          // zt = r*alpha*u_T
        }
      }
      break;

    case 3 : // calculate [ln] or [lnt] interpolating [lambda] on [mf_lambda]
      if (pmf_lambda) {
        slice_vector_on_basic_dof_of_element(*pmf_lambda, lambda, cv, coeff);
        if (contact_only) {
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          ln = aux1[0];
        } else {
          ctx.pf()->interpolation(ctx, coeff, lnt, N);
        }
      }
      break;

    case 4 : // calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
             // calculate [tau_adh] interpolating [tau_adhesion] on [mf_coeff]
             // calculate [tresca_lim] interpolating [tresca_limit] on [mf_coeff]
      GMM_ASSERT1(!contact_only, "Invalid friction option");
      if (pmf_coeff) {
        slice_vector_on_basic_dof_of_element(*pmf_coeff, friction_coeff, cv, coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
        if (gmm::vect_size(tau_adhesion)) {
          slice_vector_on_basic_dof_of_element(*pmf_coeff, tau_adhesion, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          tau_adh = aux1[0];
          if (gmm::vect_size(tresca_limit)) {
            slice_vector_on_basic_dof_of_element(*pmf_coeff, tresca_limit, cv, coeff);
            ctx.pf()->interpolation(ctx, coeff, aux1, 1);
            tresca_lim = aux1[0];
          }
        }
      }
      break;

    default : GMM_ASSERT1(false, "Invalid option");
    }

  }


  //=========================================================================
  //
  //  Integral augmented Lagrangian brick (given obstacle, u, lambda).
  //
  //=========================================================================

  template<typename MAT, typename VECT1>
  void asm_Alart_Curnier_contact_rigid_obstacle_tangent_matrix // frictionless
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    size_type subterm1 = (option == 3) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 3) ? K_UL_V1 : K_UL_V3;
    size_type subterm3 = (option == 3) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm3(subterm3, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm4(subterm4, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 3:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); " // UL
        "M$2(#3,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#1))(i,:,:,i); " // LU
        "M$3(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:)");    // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$2(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); "      // UL
        "M$3(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:);"          // LL
        "M$4(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
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
  void asm_Alart_Curnier_contact_rigid_obstacle_tangent_matrix // with friction
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   scalar_type gamma, const VECT1 *VT,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2, subterm3;
    switch (option) {
    case 1 : subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V4;
      subterm3 = K_LL_FRICT_V1; break;
    case 2 : subterm1 = K_UL_FRICT_V3; subterm2 = K_UL_FRICT_V4;
      subterm3 = K_LL_FRICT_V1; break;
    case 3 : subterm1 = K_UL_FRICT_V2; subterm2 = K_UL_FRICT_V5;
      subterm3 = K_LL_FRICT_V2; break;
    case 4 : subterm1 = K_UL_FRICT_V7; subterm2 = K_UL_FRICT_V8;
      subterm3 = K_LL_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    size_type subterm4 = K_UU_FRICT_V3;

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT),
      nterm3(subterm3, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT),
      nterm4(subterm4, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 3: case 4:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // UL
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU
        "M$3(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j)"); // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // UL
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU
        "M$3(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j);"  // LL
        "M$4(#1,#1)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j)"); // UU
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);
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
  void asm_Alart_Curnier_contact_rigid_obstacle_rhs // frictionless
  (VECT1 &Ru, VECT1 &Rl,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    size_type subterm1;
    switch (option) {
    case 1 : subterm1 = RHS_U_V1; break;
    case 2 : subterm1 = RHS_U_V2; break;
    case 3 : subterm1 = RHS_U_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }
    size_type subterm2 = (option == 3) ? RHS_L_V2 : RHS_L_V1;

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); "
              "V$2(#3)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3))(i,:)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(Ru);
    assem.push_vec(Rl);
    assem.assembly(rg);

  }

  template<typename VECT1>
  void asm_Alart_Curnier_contact_rigid_obstacle_rhs // with friction
  (VECT1 &Ru, VECT1 &Rl,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   scalar_type gamma, const VECT1 *VT,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2;
    switch (option) {
    case 1 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V1; break;
    case 2 : subterm1 = RHS_U_FRICT_V6; subterm2 = RHS_L_FRICT_V1; break;
    case 3 : subterm1 = RHS_U_FRICT_V4; subterm2 = RHS_L_FRICT_V2; break;
    case 4 : subterm1 = RHS_U_FRICT_V5; subterm2 = RHS_L_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT, gamma, VT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); "
              "V$2(#3)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_lambda);
    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(Ru);
    assem.push_vec(Rl);
    assem.assembly(rg);
  }

  struct integral_contact_rigid_obstacle_brick : public virtual_brick {

    bool contact_only;
    int option;

    // option = 1 : Alart-Curnier
    // option = 2 : symmetric Alart-Curnier (with friction, almost symmetric),
    // option = 3 : Unsymmetric method based on augmented multipliers
    // option = 4 : Unsymmetric method based on augmented multipliers
    //              with De-Saxce projection.

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
                  "Integral contact with rigid obstacle bricks need a single mesh_im");
      GMM_ASSERT1(vl.size() == 2,
                  "Integral contact with rigid obstacle bricks need two variables");
      GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 7,
                  "Wrong number of data for integral contact with rigid obstacle "
                  << "brick, " << dl.size() << " should be between 2 and 7.");
      GMM_ASSERT1(matl.size() == size_type(3 + (option == 2 && !contact_only)),
                  "Wrong number of terms for "
                  "integral contact with rigid obstacle brick");

      // variables : u, lambda. The variable lambda should be scalar in the
      //             frictionless case and vector valued in the case with
      //             friction.
      // data      : obstacle, r for the version without friction
      //           : obstacle, r, friction_coeffs, alpha, w_t, gamma, v_t for
      //             the version with friction. alpha, w_t , gamma and v_t
      //             are optional and equal to 1, 0, 1 and 0 by default,
      //             respectively.

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector &lambda = md.real_variable(vl[1]);
      const mesh_fem &mf_lambda = md.mesh_fem_of_variable(vl[1]);
      GMM_ASSERT1(mf_lambda.get_qdim() == (contact_only ? 1 : mf_u.get_qdim()),
                  "The contact stress has not the right dimension");
      const model_real_plain_vector &obstacle = md.real_variable(dl[0]);
      const mesh_fem &mf_obstacle = md.mesh_fem_of_variable(dl[0]);
      size_type sl = gmm::vect_size(obstacle) * mf_obstacle.get_qdim()
        / mf_obstacle.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vr = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      const mesh_im &mim = *mims[0];

      const model_real_plain_vector dummy_vec(0);
      const model_real_plain_vector &friction_coeffs = contact_only
                                                     ? dummy_vec : md.real_variable(dl[2]);
      const mesh_fem *pmf_coeff = contact_only ? 0 : md.pmesh_fem_of_variable(dl[2]);
      sl = gmm::vect_size(friction_coeffs);
      if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
      GMM_ASSERT1(sl == 1 || sl == 2 || sl == 3 || contact_only,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4) {
        alpha = md.real_variable(dl[3])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3])) == 1,
                    "Parameter alpha should be a scalar");
      }

      const model_real_plain_vector *WT = 0;
      if (!contact_only && dl.size() >= 5) {
        if (dl[4].compare(vl[0]) != 0)
          WT = &(md.real_variable(dl[4]));
        else if (md.n_iter_of_variable(vl[0]) > 1)
          WT = &(md.real_variable(vl[0],1));
      }

      scalar_type gamma = 1;
      if (!contact_only && dl.size() >= 6) {
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[5])) == 1,
                    "Parameter gamma should be a scalar");
        gamma = md.real_variable(dl[5])[0];
      }

      const model_real_plain_vector *VT
        = (!contact_only && dl.size()>=7) ? &(md.real_variable(dl[6])) : 0;

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Integral contact with rigid obstacle friction tangent term");
        gmm::clear(matl[0]); gmm::clear(matl[1]); gmm::clear(matl[2]);
        if (matl.size() >= 4) gmm::clear(matl[3]);
        size_type fourthmat = (matl.size() >= 4) ? 3 : 1;
        if (contact_only)
          asm_Alart_Curnier_contact_rigid_obstacle_tangent_matrix
            (matl[0], matl[1], matl[2], matl[fourthmat], mim,
             mf_u, u, mf_obstacle, obstacle, mf_lambda, lambda, vr[0],
             rg, option);
        else
          asm_Alart_Curnier_contact_rigid_obstacle_tangent_matrix
            (matl[0], matl[1], matl[2], matl[fourthmat], mim,
             mf_u, u, mf_obstacle, obstacle, mf_lambda, lambda,
             pmf_coeff, &friction_coeffs, vr[0], alpha, WT, gamma, VT,
             rg, option);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]); gmm::clear(vecl[1]); gmm::clear(vecl[2]);
        if (matl.size() >= 4) gmm::clear(vecl[3]);

        if (contact_only)
          asm_Alart_Curnier_contact_rigid_obstacle_rhs
            (vecl[0], vecl[2], mim,
             mf_u, u, mf_obstacle, obstacle, mf_lambda, lambda, vr[0],
             rg, option);
        else
          asm_Alart_Curnier_contact_rigid_obstacle_rhs
            (vecl[0], vecl[2], mim,
             mf_u, u, mf_obstacle, obstacle, mf_lambda, lambda,
             pmf_coeff, &friction_coeffs, vr[0], alpha, WT, gamma, VT,
             rg, option);
      }

    }

    integral_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
      option = option_;
      contact_only = contact_only_;
      set_flags(contact_only
                ? "Integral contact with rigid obstacle brick"
                : "Integral contact and friction with rigid obstacle brick",
                false /* is linear*/,
                (option==2) && contact_only /* is symmetric */,
                false /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };


  //=========================================================================
  //  Add a frictionless contact condition with a rigid obstacle given
  //  by a level set.
  //=========================================================================

  size_type add_integral_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname_n, const std::string &dataname_obs,
   const std::string &dataname_r, size_type region, int option) {

    pbrick pbr = new integral_contact_rigid_obstacle_brick(true, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 3 :
      tl.push_back(model::term_description(varname_u, multname_n, false)); // UL
      tl.push_back(model::term_description(multname_n, varname_u, false)); // LU
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u, multname_n, true));  // UL
      tl.push_back(model::term_description(varname_u, varname_u, true));   // UU (fourthmat == 1)
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      break;
    default :GMM_ASSERT1(false,
                         "Incorrect option for integral contact brick");

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

  size_type add_integral_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_obs,
   const std::string &dataname_r, const std::string &dataname_friction_coeffs,
   size_type region, int option,
   const std::string &dataname_alpha, const std::string &dataname_wt,
   const std::string &dataname_gamma, const std::string &dataname_vt) {

    pbrick pbr
      = new integral_contact_rigid_obstacle_brick(false, option);

    model::termlist tl;

    switch (option) {
    case 1: case 3: case 4:
      tl.push_back(model::term_description(varname_u, multname, false)); // UL
      tl.push_back(model::term_description(multname, varname_u, false)); // LU
      tl.push_back(model::term_description(multname, multname, true));   // LL
      break;
    case 2:
      tl.push_back(model::term_description(varname_u, multname, false)); // UL
      tl.push_back(model::term_description(multname, varname_u, false)); // LU
      tl.push_back(model::term_description(multname, multname, true));   // LL
      tl.push_back(model::term_description(varname_u, varname_u, true)); // UU
      break;
    default :GMM_ASSERT1(false,
                         "Incorrect option for integral contact brick");
    }
    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    dl.push_back(dataname_friction_coeffs);
    if (dataname_alpha.size()) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size()) {
        dl.push_back(dataname_wt);
        if (dataname_gamma.size()) {
          dl.push_back(dataname_gamma);
          if (dataname_vt.size()) dl.push_back(dataname_vt);
        }
      }
    }

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //
  //  Integral penalized contact with friction (given obstacle, u, lambda).
  //
  //=========================================================================

  template<typename MAT, typename VECT1>
  void asm_penalized_contact_rigid_obstacle_tangent_matrix // frictionless
  (MAT &Kuu,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    contact_rigid_obstacle_nonlinear_term
      nterm((option == 1) ? K_UU_V1 : K_UU_V2, r,
            mf_u, U, mf_obs, obs, pmf_lambda, lambda);

    const std::string aux_fems = pmf_lambda ? "#1,#2,#3": "#1,#2";
    getfem::generic_assembly assem;
    assem.set("M(#1,#1)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }


  template<typename VECT1>
  void asm_penalized_contact_rigid_obstacle_rhs // frictionless
  (VECT1 &Ru,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    contact_rigid_obstacle_nonlinear_term
      nterm((option == 1) ? RHS_U_V5 : RHS_U_V2, r,
            mf_u, U, mf_obs, obs, pmf_lambda, lambda);

    const std::string aux_fems = pmf_lambda ? "#1,#2,#3": "#1,#2";
    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru);
    assem.assembly(rg);
  }

  template<typename MAT, typename VECT1>
  void asm_penalized_contact_rigid_obstacle_tangent_matrix // with friction
  (MAT &Kuu,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  K_UU_FRICT_V4; break;
    case 2 : subterm =  K_UU_FRICT_V3; break;
    case 3 : subterm =  K_UU_FRICT_V5; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(subterm, r, mf_u, U, mf_obs, obs, pmf_lambda, lambda,
            pmf_coeff, f_coeffs, alpha, WT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4"
                                           : (pmf_lambda ? "#1,#2,#3": "#1,#2");
    getfem::generic_assembly assem;
    assem.set("M(#1,#1)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);

    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    else if (pmf_coeff)
      assem.push_mf(*pmf_coeff); // dummy

    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);

    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Kuu);
    assem.assembly(rg);
  }


  template<typename VECT1>
  void asm_penalized_contact_rigid_obstacle_rhs // with friction
  (VECT1 &Ru,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_FRICT_V7; break;
    case 2 : subterm =  RHS_U_FRICT_V6; break;
    case 3 : subterm =  RHS_U_FRICT_V8; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(subterm, r, mf_u, U, mf_obs, obs, pmf_lambda, lambda,
            pmf_coeff, f_coeffs, alpha, WT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4"
                                           : (pmf_lambda ? "#1,#2,#3": "#1,#2");
    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);

    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    else if (pmf_coeff)
      assem.push_mf(*pmf_coeff); // dummy

    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);

    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru);
    assem.assembly(rg);
  }


  struct penalized_contact_rigid_obstacle_brick : public virtual_brick {

    bool contact_only;
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
      // Integration method
      GMM_ASSERT1(mims.size() == 1,
                  "Penalized contact with rigid obstacle bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];

      // Variables : u
      GMM_ASSERT1(vl.size() == 1,
                  "Penalized contact with rigid obstacle bricks need a single variable");
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);

      size_type N = mf_u.linked_mesh().dim();

      // Data : obs, r, [lambda,] [friction_coeffs,] [alpha,] [WT]
      size_type nb_data_1 = ((option == 1) ? 2 : 3) + (contact_only ? 0 : 1);
      size_type nb_data_2 = nb_data_1 + (contact_only ? 0 : 2);
      GMM_ASSERT1(dl.size() >= nb_data_1 && dl.size() <= nb_data_2,
                  "Wrong number of data for penalized contact with rigid obstacle  "
                  << "brick, " << dl.size() << " should be between "
                  << nb_data_1 << " and " << nb_data_2 << ".");

      size_type nd = 0;
      const model_real_plain_vector &obs = md.real_variable(dl[nd]);
      const mesh_fem &mf_obs = md.mesh_fem_of_variable(dl[nd]);
      size_type sl = gmm::vect_size(obs) * mf_obs.get_qdim() / mf_obs.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      nd++;
      const model_real_plain_vector &vr = md.real_variable(dl[nd]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");

      const model_real_plain_vector *lambda = 0;
      const mesh_fem *pmf_lambda = 0;
      if (option != 1) {
        nd++;
        lambda = &(md.real_variable(dl[nd]));
        pmf_lambda = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*lambda) * pmf_lambda->get_qdim() / pmf_lambda->nb_dof();
        GMM_ASSERT1(sl == (contact_only ? 1 : N),
                    "the data corresponding to the contact stress "
                    "has not the right format");
      }

      const model_real_plain_vector *f_coeffs = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      const model_real_plain_vector *WT = 0;
      if (!contact_only) {
        nd++;
        f_coeffs = &(md.real_variable(dl[nd]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*f_coeffs);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1 || sl == 2 || sl == 3,
                    "the data corresponding to the friction coefficient "
                    "has not the right format");
        if (dl.size() > nd+1) {
          nd++;
          alpha = md.real_variable(dl[nd])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[nd])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() > nd+1) {
          nd++;
          if (dl[nd].compare(vl[0]) != 0)
            WT = &(md.real_variable(dl[nd]));
          else if (md.n_iter_of_variable(vl[0]) > 1)
            WT = &(md.real_variable(vl[0],1));
        }
      }

      GMM_ASSERT1(matl.size() == 1, "Wrong number of terms for "
                  "penalized contact with rigid obstacle brick");

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Penalized contact with rigid obstacle tangent term");
        gmm::clear(matl[0]);
        if (contact_only)
          asm_penalized_contact_rigid_obstacle_tangent_matrix
            (matl[0], mim, mf_u, u, mf_obs, obs, pmf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_penalized_contact_rigid_obstacle_tangent_matrix
            (matl[0], mim, mf_u, u, mf_obs, obs, pmf_lambda, lambda,
             pmf_coeff, f_coeffs, vr[0], alpha, WT, rg, option);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        if (contact_only)
          asm_penalized_contact_rigid_obstacle_rhs
            (vecl[0], mim, mf_u, u, mf_obs, obs, pmf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_penalized_contact_rigid_obstacle_rhs
            (vecl[0], mim, mf_u, u, mf_obs, obs, pmf_lambda, lambda,
             pmf_coeff, f_coeffs, vr[0], alpha, WT, rg, option);
      }

    }

    penalized_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
      contact_only = contact_only_;
      option = option_;
      set_flags(contact_only
                ? "Integral penalized contact with rigid obstacle brick"
                : "Integral penalized contact and friction with rigid obstacle brick",
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

    pbrick pbr = new penalized_contact_rigid_obstacle_brick(true, option);

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

  size_type add_penalized_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   const std::string &dataname_friction_coeffs,
   size_type region, int option, const std::string &dataname_lambda,
   const std::string &dataname_alpha, const std::string &dataname_wt) {

    pbrick pbr = new penalized_contact_rigid_obstacle_brick(false, option);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    switch (option) {
    case 1: break;
    case 2: case 3: dl.push_back(dataname_lambda); break;
    default: GMM_ASSERT1(false, "Penalized contact brick : invalid option");
    }
    dl.push_back(dataname_friction_coeffs);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size() > 0) dl.push_back(dataname_wt);
    }

    model::varnamelist vl(1, varname_u);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  //=========================================================================
  //
  //  Integral contact (with friction) between non-matching meshes.
  //
  //=========================================================================

  template<typename MAT, typename VEC>
  void asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix // frictionless
  (MAT &Ku1l, MAT &Klu1, MAT &Ku2l, MAT &Klu2, MAT &Kll,
   MAT &Ku1u1, MAT &Ku2u2, MAT &Ku1u2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VEC &U1,
   const getfem::mesh_fem &mf_u2, const VEC &U2,
   const getfem::mesh_fem &mf_lambda, const VEC &lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    size_type subterm1 = (option == 3) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 3) ? K_UL_V1 : K_UL_V3;
    size_type subterm3 = (option == 3) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm3(subterm3, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm4(subterm4, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 3:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); " // U1L
        "M$2(#3,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#1))(i,:,:,i); " // LU1
        "M$3(#2,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#2).Base(#3))(i,:,i,:); " // U2L
        "M$4(#3,#2)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#2))(i,:,:,i); " // LU2
        "M$5(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:)");    // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$2(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); "      // U1L
        "M$3(#2,#3)+=comp(NonLin$2(#1,#1,#2,#3).vBase(#2).Base(#3))(i,:,i,:); "      // U2L
        "M$5(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:); "         // LL
        "M$6(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j); " // U1U1
        "M$7(#2,#2)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#2).vBase(#2))(i,j,:,i,:,j); " // U2U2
        "M$8(#1,#2)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#2))(i,j,:,i,:,j)"); // U1U2
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Ku1l);
    assem.push_mat(Klu1);
    assem.push_mat(Ku2l);
    assem.push_mat(Klu2);
    assem.push_mat(Kll);
    assem.push_mat(Ku1u1);
    assem.push_mat(Ku2u2);
    assem.push_mat(Ku1u2);
    assem.assembly(rg);

    gmm::scale(Ku2l, scalar_type(-1));
    if (option != 2) // Klu2 was calculated
      gmm::scale(Klu2, scalar_type(-1));
    gmm::scale(Ku1u2, scalar_type(-1));
  }

  template<typename MAT, typename VEC>
  void asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix // with friction
  (MAT &Ku1l, MAT &Klu1, MAT &Ku2l, MAT &Klu2, MAT &Kll,
   MAT &Ku1u1, MAT &Ku2u2, MAT &Ku1u2, MAT &Ku2u1,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VEC &U1,
   const getfem::mesh_fem &mf_u2, const VEC &U2,
   const getfem::mesh_fem &mf_lambda, const VEC &lambda,
   const getfem::mesh_fem *pmf_coeff, const VEC *f_coeffs,
   scalar_type r, scalar_type alpha,
   const VEC *WT1, const VEC *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2, subterm3;
    switch (option) {
    case 1 :
      subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V4; subterm3 = K_LL_FRICT_V1;
      break;
    case 2 :
      subterm1 = K_UL_FRICT_V3; subterm2 = K_UL_FRICT_V4; subterm3 = K_LL_FRICT_V1;
      break;
    case 3 :
      subterm1 = K_UL_FRICT_V2; subterm2 = K_UL_FRICT_V5; subterm3 = K_LL_FRICT_V2;
      break;
    case 4 :
      subterm1 = K_UL_FRICT_V7; subterm2 = K_UL_FRICT_V8; subterm3 = K_LL_FRICT_V4;
      break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    size_type subterm4 = K_UU_FRICT_V3;

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2),
      nterm3(subterm3, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2),
      nterm4(subterm4, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 3: case 4:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // U1L
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU1
        "M$3(#2,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2).vBase(#3))(i,j,:,i,:,j); " // U2L
        "M$4(#3,#2)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#2))(i,j,:,j,:,i); " // LU2
        "M$5(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j)"); // LL
      break;
    case 2:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // U1L
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU1
        "M$3(#2,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2).vBase(#3))(i,j,:,i,:,j); " // U2L
        "M$4(#3,#2)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#2))(i,j,:,j,:,i); " // LU2
        "M$5(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j); " // LL
        "M$6(#1,#1)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j); " // U1U1
        "M$7(#2,#2)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#2).vBase(#2))(i,j,:,i,:,j); " // U2U2
        "M$8(#1,#2)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#1).vBase(#2))(i,j,:,i,:,j); " // U1U2
        "M$9(#2,#1)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#2).vBase(#1))(i,j,:,i,:,j)"); // U2U1
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mf(mf_lambda);
    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Ku1l);
    assem.push_mat(Klu1);
    assem.push_mat(Ku2l);
    assem.push_mat(Klu2);
    assem.push_mat(Kll);
    assem.push_mat(Ku1u1);
    assem.push_mat(Ku2u2);
    assem.push_mat(Ku1u2);
    assem.push_mat(Ku2u1);
    assem.assembly(rg);

    gmm::scale(Ku2l, scalar_type(-1));
    gmm::scale(Klu2, scalar_type(-1));
    gmm::scale(Ku1u2, scalar_type(-1));
  }

  template<typename VECT1>
  void asm_Alart_Curnier_contact_nonmatching_meshes_rhs // frictionless
  (VECT1 &Ru1, VECT1 &Ru2, VECT1 &Rl,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    size_type subterm1;
    switch (option) {
    case 1 : subterm1 = RHS_U_V1; break;
    case 2 : subterm1 = RHS_U_V2; break;
    case 3 : subterm1 = RHS_U_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }
    size_type subterm2 = (option == 3) ? RHS_L_V2 : RHS_L_V1;

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#2))(i,:,i); "
              "V$3(#3)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3))(i,:)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(Ru1);
    assem.push_vec(Ru2);
    assem.push_vec(Rl);
    assem.assembly(rg);

    gmm::scale(Ru2, scalar_type(-1));
  }

  template<typename VECT1>
  void asm_Alart_Curnier_contact_nonmatching_meshes_rhs // with friction
  (VECT1 &Ru1, VECT1 &Ru2, VECT1 &Rl,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs,
   scalar_type r, scalar_type alpha,
   const VECT1 *WT1, const VECT1 *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2;
    switch (option) {
    case 1 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V1; break;
    case 2 : subterm1 = RHS_U_FRICT_V6; subterm2 = RHS_L_FRICT_V1; break;
    case 3 : subterm1 = RHS_U_FRICT_V4; subterm2 = RHS_L_FRICT_V2; break;
    case 4 : subterm1 = RHS_U_FRICT_V5; subterm2 = RHS_L_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeffs, alpha, WT1, WT2);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2))(i,:,i); "
              "V$3(#3)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mf(mf_lambda);
    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(Ru1);
    assem.push_vec(Ru2);
    assem.push_vec(Rl);
    assem.assembly(rg);

    gmm::scale(Ru2, scalar_type(-1));
  }

  struct integral_contact_nonmatching_meshes_brick : public virtual_brick {

    size_type rg1, rg2; // ids of mesh regions on mf_u1 and mf_u2 that are
                        // expected to come in contact.
    mutable getfem::pfem pfem_proj; // cached fem for the projection between nonmatching meshes
    bool contact_only;
    int option;

    // option = 1 : Alart-Curnier
    // option = 2 : symmetric Alart-Curnier (almost symmetric with friction)
    // option = 3 : New method
    // option = 4 : New method based on De-Saxce.

    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      // Integration method
      GMM_ASSERT1(mims.size() == 1,
                  "Integral contact between nonmatching meshes bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];

      // Variables : u1, u2, lambda
      //       the variable lambda should be scalar in the frictionless
      //       case and vector valued in the case with friction.
      GMM_ASSERT1(vl.size() == 3,
                  "Integral contact between nonmatching meshes bricks need three variables");
      const model_real_plain_vector &u1 = md.real_variable(vl[0]);
      const model_real_plain_vector &u2 = md.real_variable(vl[1]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_u2 = md.mesh_fem_of_variable(vl[1]);
      const model_real_plain_vector &lambda = md.real_variable(vl[2]);
      const mesh_fem &mf_lambda = md.mesh_fem_of_variable(vl[2]);
      GMM_ASSERT1(mf_lambda.get_qdim() == (contact_only ? 1 : mf_u1.get_qdim()),
                  "The contact stress variable has not the right dimension");

      // Data : r, [friction_coeffs,] [alpha,] [WT1, WT2]
      //     alpha, WT1, WT2 are optional and equal to 1, 0 and 0 by default respectively.
      if (contact_only) {
        GMM_ASSERT1(dl.size() == 1,
                    "Wrong number of data for integral contact between nonmatching meshes "
                    << "brick, the number of data should be equal to 1 .");
      }
      else {
        GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 5,
                    "Wrong number of data for integral contact between nonmatching meshes "
                    << "brick, it should be between 2 and 5 instead of "  << dl.size() << " .");
      }
      const model_real_plain_vector &vr = md.real_variable(dl[0]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");

      const model_real_plain_vector *f_coeffs = 0, *WT1 = 0, *WT2 = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      if (!contact_only) {
        f_coeffs = &(md.real_variable(dl[1]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[1]);

        size_type sl = gmm::vect_size(*f_coeffs);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1 || sl == 2 || sl ==3,
                    "the data corresponding to the friction coefficient "
                    "has not the right format");

        if (dl.size() >= 3) {
          alpha = md.real_variable(dl[2])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[2])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() >= 4) {
          if (dl[3].compare(vl[0]) != 0)
            WT1 = &(md.real_variable(dl[3]));
          else if (md.n_iter_of_variable(vl[0]) > 1)
            WT1 = &(md.real_variable(vl[0],1));
        }

        if (dl.size() >= 5) {
          if (dl[4].compare(vl[1]) != 0)
            WT2 = &(md.real_variable(dl[4]));
          else if (md.n_iter_of_variable(vl[1]) > 1)
            WT2 = &(md.real_variable(vl[1],1));
        }
      }

      // Matrix terms (T_u1l, T_lu1, T_u2l, T_lu2, T_ll, T_u1u1, T_u2u2, T_u1u2)
      GMM_ASSERT1(matl.size() == size_type(3 +                                // U1L, U2L, LL
                                           2 * !is_symmetric() +              // LU1, LU2
                                           3 * (option == 2) + // U1U1, U2U2, U1U2
                                           1 * (option == 2 && !is_symmetric())), // U2U1
                  "Wrong number of terms for "
                  "integral contact between nonmatching meshes brick");

      mesh_region rg(region);
      mf_u1.linked_mesh().intersect_with_mpi_region(rg);

      size_type N = mf_u1.linked_mesh().dim();

      // projection of the second mesh_fem onto the mesh of the first mesh_fem
      if (!pfem_proj)
        pfem_proj = new_projected_fem(mf_u2, mim, rg2, rg1);

      getfem::mesh_fem mf_u2_proj(mim.linked_mesh(), dim_type(N));
      mf_u2_proj.set_finite_element(mim.linked_mesh().convex_index(), pfem_proj);

      size_type nbdof1 = mf_u1.nb_dof();
      size_type nbdof_lambda = mf_lambda.nb_dof();
      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = mf_u2_proj.nb_basic_dof();

      std::vector<size_type> ind;
      mf_u2_proj.get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      gmm::csc_matrix<scalar_type> Esub(nbsub, nbdof2);
      if (mf_u2.is_reduced())
          gmm::copy(gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    Esub);

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(Esub, u2, u2_proj);
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      model_real_plain_vector WT2_proj(0);
      if (WT2) {
        gmm::resize(WT2_proj, nbsub);
        if (mf_u2.is_reduced())
          gmm::mult(Esub, *WT2, WT2_proj);
        else
          gmm::copy(gmm::sub_vector(*WT2, SUBI), WT2_proj);
      }

      size_type U1L = 0;
      size_type LU1 = is_symmetric() ? size_type(-1) : 1;
      size_type U2L = is_symmetric() ? 1 : 2;
      size_type LU2 = is_symmetric() ? size_type(-1) : 3;
      size_type LL = is_symmetric() ? 2 : 4;
      size_type U1U1 = (option != 2) ? size_type(-1) : (is_symmetric() ? 3 : 5);
      size_type U2U2 = (option != 2) ? size_type(-1) : (is_symmetric() ? 4 : 6);
      size_type U1U2 = (option != 2) ? size_type(-1) : (is_symmetric() ? 5 : 7);
      size_type U2U1 = (option != 2 || is_symmetric()) ? size_type(-1) : 8;

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Integral contact between nonmatching meshes "
                   "tangent term");
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(matl[i]);

        model_real_sparse_matrix dummy_mat(0, 0);
        model_real_sparse_matrix &Klu1 = (LU1 == size_type(-1)) ? dummy_mat : matl[LU1];
        model_real_sparse_matrix &Ku1u1 = (U1U1 == size_type(-1)) ? dummy_mat : matl[U1U1];

        model_real_sparse_matrix Ku2l(nbsub, nbdof_lambda);
        model_real_sparse_matrix Klu2(nbdof_lambda, nbsub);
        model_real_sparse_matrix Ku2u2(nbsub, nbsub);
        model_real_sparse_matrix Ku1u2(nbdof1, nbsub);
        model_real_sparse_matrix Ku2u1(nbsub, nbdof1);

        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[U1L], Klu1, Ku2l, Klu2, matl[LL], Ku1u1, Ku2u2, Ku1u2,
             mim, mf_u1, u1, mf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[U1L], Klu1, Ku2l, Klu2, matl[LL], Ku1u1, Ku2u2, Ku1u2, Ku2u1,
             mim, mf_u1, u1, mf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, f_coeffs, vr[0], alpha, WT1, &WT2_proj, rg, option);

        if (mf_u2.is_reduced()) {
          gmm::mult(gmm::transposed(Esub), Ku2l, matl[U2L]);
          if (LU2 != size_type(-1)) gmm::mult(Klu2, Esub, matl[LU2]);
          if (U2U2 != size_type(-1)) {
            model_real_sparse_matrix tmp(nbsub, nbdof2);
            gmm::mult(Ku2u2, Esub, tmp);
            gmm::mult(gmm::transposed(Esub), tmp, matl[U2U2]);
          }
          if (U1U2 != size_type(-1)) gmm::mult(Ku1u2, Esub, matl[U1U2]);
          if (U2U1 != size_type(-1)) gmm::mult(gmm::transposed(Esub), Ku2u1, matl[U2U1]);
        }
        else {
          gmm::copy(Ku2l, gmm::sub_matrix(matl[U2L], SUBI, gmm::sub_interval(0, nbdof_lambda)));
          if (LU2 != size_type(-1))
            gmm::copy(Klu2, gmm::sub_matrix(matl[LU2], gmm::sub_interval(0, nbdof_lambda), SUBI));
          if (U2U2 != size_type(-1))
            gmm::copy(Ku2u2, gmm::sub_matrix(matl[U2U2], SUBI));
          if (U1U2 != size_type(-1))
            gmm::copy(Ku1u2, gmm::sub_matrix(matl[U1U2], gmm::sub_interval(0, nbdof1), SUBI));
          if (U2U1 != size_type(-1))
            gmm::copy(Ku2u1, gmm::sub_matrix(matl[U2U1], SUBI, gmm::sub_interval(0, nbdof1)));
        }
      }

      if (version & model::BUILD_RHS) {
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(vecl[i]);

        model_real_plain_vector Ru2(nbsub);

        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[U1L], Ru2, vecl[LL], // u1, u2, lambda
             mim, mf_u1, u1, mf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[U1L], Ru2, vecl[LL], // u1, u2, lambda
             mim, mf_u1, u1, mf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, f_coeffs, vr[0], alpha, WT1, &WT2_proj, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(gmm::transposed(Esub), Ru2, vecl[U2L]);
        else
          gmm::copy(Ru2, gmm::sub_vector(vecl[U2L], SUBI));
      }
    }

    integral_contact_nonmatching_meshes_brick(size_type rg1_, size_type rg2_,
                                                bool contact_only_, int option_)
    : rg1(rg1_), rg2(rg2_), pfem_proj(0),
      contact_only(contact_only_), option(option_)
    {
      set_flags(contact_only
                ? "Integral contact between nonmatching meshes brick"
                : "Integral contact and friction between nonmatching "
                  "meshes brick",
                false /* is linear*/,
                (option==2) && contact_only /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

    ~integral_contact_nonmatching_meshes_brick()
    { if (pfem_proj) del_projected_fem(pfem_proj); }

  };


  //=========================================================================
  //  Add a frictionless contact condition between two bodies discretized
  //  with nonmatching meshes.
  //=========================================================================

  size_type add_integral_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &multname_n,
   const std::string &dataname_r,
   size_type region1, size_type region2, int option) {

    pbrick pbr = new integral_contact_nonmatching_meshes_brick
                     (region1, region2, true /* contact_only */, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 3 :
      tl.push_back(model::term_description(varname_u1, multname_n, false)); // U1L
      tl.push_back(model::term_description(multname_n, varname_u1, false)); // LU1
      tl.push_back(model::term_description(varname_u2, multname_n, false)); // U2L
      tl.push_back(model::term_description(multname_n, varname_u2, false)); // LU2
      tl.push_back(model::term_description(multname_n, multname_n, true));  // LL
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u1, multname_n, true)); // U1L
      tl.push_back(model::term_description(varname_u2, multname_n, true)); // U2L
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      tl.push_back(model::term_description(varname_u1, varname_u1, true)); // U1U1
      tl.push_back(model::term_description(varname_u2, varname_u2, true)); // U2U2
      tl.push_back(model::term_description(varname_u1, varname_u2, true)); // U1U2
      break;
    default : GMM_ASSERT1(false,
                          "Incorrect option for integral contact brick");
    }

    model::varnamelist dl(1, dataname_r);

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);
    vl.push_back(multname_n);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region1);
  }


  //=========================================================================
  //  Add a contact condition with Coulomb friction between two bodies
  //  discretized with nonmatching meshes.
  //=========================================================================

  size_type add_integral_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &multname,
   const std::string &dataname_r, const std::string &dataname_friction_coeffs,
   size_type region1, size_type region2, int option,
   const std::string &dataname_alpha,
   const std::string &dataname_wt1, const std::string &dataname_wt2) {

    pbrick pbr = new integral_contact_nonmatching_meshes_brick
                     (region1, region2, false /* contact_only */, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 3 : case 4 :
      tl.push_back(model::term_description(varname_u1, multname, false));  // 0: U1L
      tl.push_back(model::term_description(multname, varname_u1, false));  // 1: LU1
      tl.push_back(model::term_description(varname_u2, multname, false));  // 2: U2L
      tl.push_back(model::term_description(multname, varname_u2, false));  // 3: LU2
      tl.push_back(model::term_description(multname, multname, true));     // 4: LL
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u1, multname, false));  // 0: U1L
      tl.push_back(model::term_description(multname, varname_u1, false));  // 1: LU1
      tl.push_back(model::term_description(varname_u2, multname, false));  // 2: U2L
      tl.push_back(model::term_description(multname, varname_u2, false));  // 3: LU2
      tl.push_back(model::term_description(multname, multname, true));     // 4: LL
      tl.push_back(model::term_description(varname_u1, varname_u1, true)); // 5: U1U1
      tl.push_back(model::term_description(varname_u2, varname_u2, true)); // 6: U2U2
      tl.push_back(model::term_description(varname_u1, varname_u2, true)); // 7: U1U2
      tl.push_back(model::term_description(varname_u2, varname_u1, true)); // 8: U2U1
      break;
    default : GMM_ASSERT1(false,
                          "Incorrect option for integral contact brick");
    }

    model::varnamelist dl(1, dataname_r);   // 0 -> r
    dl.push_back(dataname_friction_coeffs); // 1 -> f_coeff,[tau_adh,tresca_lim]
    if (dataname_alpha.size()) {
      dl.push_back(dataname_alpha);         // 2 -> alpha
      if (dataname_wt1.size()) {
        dl.push_back(dataname_wt1);         // 3 -> WT1
        if (dataname_wt2.size()) {
          dl.push_back(dataname_wt2);       // 4 -> WT2
          // TODO: VT1, VT2, gamma
        }
      }
    }

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);
    vl.push_back(multname);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region1);
  }


  //=========================================================================
  //
  //  Integral penalized contact (with friction) between non-matching meshes.
  //
  //=========================================================================

  template<typename MAT, typename VECT1>
  void asm_penalized_contact_nonmatching_meshes_tangent_matrix // frictionless
  (MAT &Ku1u1, MAT &Ku2u2, MAT &Ku1u2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   const scalar_type r, const mesh_region &rg, int option = 1) {

    contact_nonmatching_meshes_nonlinear_term
      nterm((option == 1) ? K_UU_V1 : K_UU_V2, r,
            mf_u1, U1, mf_u2, U2, pmf_lambda, lambda);

    const std::string aux_fems = pmf_lambda ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly assem;
    assem.set("M$1(#1,#1)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j); "
              "M$2(#2,#2)+=comp(NonLin(#1," + aux_fems + ").vBase(#2).vBase(#2))(i,j,:,i,:,j); "
              "M$3(#1,#2)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#2))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Ku1u1);
    assem.push_mat(Ku2u2);
    assem.push_mat(Ku1u2);
    assem.assembly(rg);

    gmm::scale(Ku1u2, scalar_type(-1));
  }

  template<typename VECT1>
  void asm_penalized_contact_nonmatching_meshes_rhs // frictionless
  (VECT1 &Ru1, VECT1 &Ru2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   scalar_type r, const mesh_region &rg, int option = 1) {

    contact_nonmatching_meshes_nonlinear_term
      nterm((option == 1) ? RHS_U_V5 : RHS_U_V2, r,
            mf_u1, U1, mf_u2, U2, pmf_lambda, lambda);

    const std::string aux_fems = pmf_lambda ? "#1,#2,#3": "#1,#2";
    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru1);
    assem.push_vec(Ru2);
    assem.assembly(rg);

    gmm::scale(Ru2, scalar_type(-1));
  }


  template<typename MAT, typename VECT1>
  void asm_penalized_contact_nonmatching_meshes_tangent_matrix // with friction
  (MAT &Ku1u1, MAT &Ku2u2, MAT &Ku1u2, MAT &Ku2u1,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT1, const VECT1 *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  K_UU_FRICT_V4; break;
    case 2 : subterm =  K_UU_FRICT_V3; break;
    case 3 : subterm =  K_UU_FRICT_V5; break;
    }

    contact_nonmatching_meshes_nonlinear_term
      nterm(subterm, r, mf_u1, U1, mf_u2, U2, pmf_lambda, lambda,
            pmf_coeff, f_coeffs, alpha, WT1, WT2);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4"
                                           : (pmf_lambda ? "#1,#2,#3": "#1,#2");

    getfem::generic_assembly assem;
    assem.set("M$1(#1,#1)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j); "
              "M$2(#2,#2)+=comp(NonLin(#1," + aux_fems + ").vBase(#2).vBase(#2))(i,j,:,i,:,j); "
              "M$3(#1,#2)+=comp(NonLin(#1," + aux_fems + ").vBase(#1).vBase(#2))(i,j,:,i,:,j); "
              "M$4(#2,#1)+=comp(NonLin(#1," + aux_fems + ").vBase(#2).vBase(#1))(i,j,:,i,:,j)");
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);

    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    else if (pmf_coeff)
      assem.push_mf(*pmf_coeff); // dummy

    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);

    assem.push_nonlinear_term(&nterm);
    assem.push_mat(Ku1u1);
    assem.push_mat(Ku2u2);
    assem.push_mat(Ku1u2);
    assem.push_mat(Ku2u1);
    assem.assembly(rg);

    gmm::scale(Ku1u2, scalar_type(-1));
    gmm::scale(Ku2u1, scalar_type(-1));
  }


  template<typename VECT1>
  void asm_penalized_contact_nonmatching_meshes_rhs // with friction
  (VECT1 &Ru1, VECT1 &Ru2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem *pmf_lambda, const VECT1 *lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeffs, scalar_type r,
   scalar_type alpha, const VECT1 *WT1, const VECT1 *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_FRICT_V7; break;
    case 2 : subterm =  RHS_U_FRICT_V6; break;
    case 3 : subterm =  RHS_U_FRICT_V8; break;
    }

    contact_nonmatching_meshes_nonlinear_term
      nterm(subterm, r, mf_u1, U1, mf_u2, U2, pmf_lambda, lambda,
            pmf_coeff, f_coeffs, alpha, WT1, WT2);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4"
                                           : (pmf_lambda ? "#1,#2,#3": "#1,#2");
    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2))(i,:,i)");

    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);

    if (pmf_lambda)
      assem.push_mf(*pmf_lambda);
    else if (pmf_coeff)
      assem.push_mf(*pmf_coeff); // dummy

    if (pmf_coeff)
      assem.push_mf(*pmf_coeff);

    assem.push_nonlinear_term(&nterm);
    assem.push_vec(Ru1);
    assem.push_vec(Ru2);
    assem.assembly(rg);

    gmm::scale(Ru2, scalar_type(-1));
  }


  struct penalized_contact_nonmatching_meshes_brick : public virtual_brick {

    size_type rg1, rg2; // ids of mesh regions on mf_u1 and mf_u2 that are
                        // expected to come in contact.
    mutable getfem::pfem pfem_proj; // cached fem for the projection between nonmatching meshes
    bool contact_only;
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
      // Integration method
      GMM_ASSERT1(mims.size() == 1,
                  "Penalized contact between nonmatching meshes bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];

      // Variables : u1, u2
      GMM_ASSERT1(vl.size() == 2,
                  "Penalized contact between nonmatching meshes bricks need two variables");
      const model_real_plain_vector &u1 = md.real_variable(vl[0]);
      const model_real_plain_vector &u2 = md.real_variable(vl[1]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_u2 = md.mesh_fem_of_variable(vl[1]);

      size_type N = mf_u1.linked_mesh().dim();

      // Data : r, [lambda,] [friction_coeffs,] [alpha,] [WT1, WT2]
      size_type nb_data_1 = ((option == 1) ? 1 : 2) + (contact_only ? 0 : 1);
      size_type nb_data_2 = nb_data_1 + (contact_only ? 0 : 3);
      GMM_ASSERT1(dl.size() >= nb_data_1 && dl.size() <= nb_data_2,
                  "Wrong number of data for penalized contact between nonmatching meshes "
                  << "brick, " << dl.size() << " should be between "
                  << nb_data_1 << " and " << nb_data_2 << ".");

      size_type nd = 0;
      const model_real_plain_vector &vr = md.real_variable(dl[nd]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");

      size_type sl;
      const model_real_plain_vector *lambda = 0;
      const mesh_fem *pmf_lambda = 0;
      if (option != 1) {
        nd++;
        lambda = &(md.real_variable(dl[nd]));
        pmf_lambda = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*lambda) * pmf_lambda->get_qdim() / pmf_lambda->nb_dof();
        GMM_ASSERT1(sl == (contact_only ? 1 : N),
                    "the data corresponding to the contact stress "
                    "has not the right format");
      }

      const model_real_plain_vector *f_coeffs = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      const model_real_plain_vector *WT1 = 0;
      const model_real_plain_vector *WT2 = 0;
      if (!contact_only) {
        nd++;
        f_coeffs = &(md.real_variable(dl[nd]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*f_coeffs);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1 || sl == 2 || sl == 3,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

        if (dl.size() > nd+1) {
          nd++;
          alpha = md.real_variable(dl[nd])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[nd])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() > nd+1) {
          nd++;
          if (dl[nd].compare(vl[0]) != 0)
            WT1 = &(md.real_variable(dl[nd]));
          else if (md.n_iter_of_variable(vl[0]) > 1)
            WT1 = &(md.real_variable(vl[0],1));
        }

        if (dl.size() > nd+1) {
          nd++;
          if (dl[nd].compare(vl[1]) != 0)
            WT2 = &(md.real_variable(dl[nd]));
          else if (md.n_iter_of_variable(vl[1]) > 1)
            WT2 = &(md.real_variable(vl[1],1));
        }
      }

      GMM_ASSERT1(matl.size() == contact_only ? 3 : 4,
                  "Wrong number of terms for penalized contact "
                  "between nonmatching meshes brick");

      mesh_region rg(region);
      mf_u1.linked_mesh().intersect_with_mpi_region(rg); // FIXME: mfu_2?

      // projection of the second mesh_fem onto the mesh of the first mesh_fem
      if (!pfem_proj)
        pfem_proj = new_projected_fem(mf_u2, mim, rg2, rg1);

      getfem::mesh_fem mf_u2_proj(mim.linked_mesh(), dim_type(N));
      mf_u2_proj.set_finite_element(mim.linked_mesh().convex_index(), pfem_proj);

      size_type nbdof1 = mf_u1.nb_dof();
      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = mf_u2_proj.nb_dof();

      std::vector<size_type> ind;
      mf_u2_proj.get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      gmm::csc_matrix<scalar_type> Esub(nbsub, nbdof2);
      if (mf_u2.is_reduced())
          gmm::copy(gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    Esub);

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(Esub, u2, u2_proj);
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      model_real_plain_vector WT2_proj(0);
      if (WT2) {
        gmm::resize(WT2_proj, nbsub);
        if (mf_u2.is_reduced())
          gmm::mult(Esub, *WT2, WT2_proj);
        else
          gmm::copy(gmm::sub_vector(*WT2, SUBI), WT2_proj);
      }

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Penalized contact between nonmatching meshes tangent term");
        gmm::clear(matl[0]);
        gmm::clear(matl[1]);
        gmm::clear(matl[2]);

        model_real_sparse_matrix Ku2u2(nbsub,nbsub);
        model_real_sparse_matrix Ku1u2(nbdof1,nbsub);

        if (contact_only) {
          asm_penalized_contact_nonmatching_meshes_tangent_matrix
            (matl[0], Ku2u2, Ku1u2, mim, mf_u1, u1, mf_u2_proj, u2_proj,
             pmf_lambda, lambda, vr[0], rg, option);
        }
        else {
          gmm::clear(matl[3]);
          model_real_sparse_matrix Ku2u1(nbsub,nbdof1);
          asm_penalized_contact_nonmatching_meshes_tangent_matrix
            (matl[0], Ku2u2, Ku1u2, Ku2u1, mim, mf_u1, u1, mf_u2_proj, u2_proj,
             pmf_lambda, lambda, pmf_coeff, f_coeffs, vr[0], alpha, WT1, &WT2_proj, rg, option);
          if (mf_u2.is_reduced())
            gmm::mult(gmm::transposed(Esub), Ku2u1, matl[3]);
          else
            gmm::copy(Ku2u1, gmm::sub_matrix(matl[3], SUBI, gmm::sub_interval(0, nbdof1)));
        }

        if (mf_u2.is_reduced()) {
          model_real_sparse_matrix tmp(nbsub, nbdof2);
          gmm::mult(Ku2u2, Esub, tmp);
          gmm::mult(gmm::transposed(Esub), tmp, matl[1]);
          gmm::mult(Ku1u2, Esub, matl[2]);
        }
        else {
          gmm::copy(Ku2u2, gmm::sub_matrix(matl[1], SUBI));
          gmm::copy(Ku1u2, gmm::sub_matrix(matl[2], gmm::sub_interval(0, nbdof1), SUBI));
        }
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        gmm::clear(vecl[1]);

        model_real_plain_vector Ru2(nbsub);

        if (contact_only)
          asm_penalized_contact_nonmatching_meshes_rhs
            (vecl[0], Ru2, mim, mf_u1, u1, mf_u2_proj, u2_proj, pmf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_penalized_contact_nonmatching_meshes_rhs
            (vecl[0], Ru2, mim, mf_u1, u1, mf_u2_proj, u2_proj, pmf_lambda, lambda,
             pmf_coeff, f_coeffs, vr[0], alpha, WT1, &WT2_proj, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(gmm::transposed(Esub), Ru2, vecl[1]);
        else
          gmm::copy(Ru2, gmm::sub_vector(vecl[1], SUBI));
      }
    }

    penalized_contact_nonmatching_meshes_brick(size_type rg1_, size_type rg2_,
                                               bool contact_only_, int option_)
    : rg1(rg1_), rg2(rg2_), pfem_proj(0),
      contact_only(contact_only_), option(option_)
    {
      set_flags(contact_only
                ? "Integral penalized contact between nonmatching meshes brick"
                : "Integral penalized contact and friction between nonmatching "
                  "meshes brick",
                false /* is linear*/, contact_only /* is symmetric */,
                true /* is coercive */, true /* is real */,
                false /* is complex */);
    }

    ~penalized_contact_nonmatching_meshes_brick()
    { if (pfem_proj) del_projected_fem(pfem_proj); }

  };


  //=========================================================================
  //  Add a frictionless contact condition between two bodies discretized
  //  with nonmatching meshes.
  //=========================================================================

  size_type add_penalized_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_r,
   size_type region1, size_type region2,
   int option, const std::string &dataname_n) {

    pbrick pbr = new penalized_contact_nonmatching_meshes_brick
                     (region1, region2, true /* contact_only */, option);
    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, true));
    tl.push_back(model::term_description(varname_u2, varname_u2, true));
    tl.push_back(model::term_description(varname_u1, varname_u2, true));

    model::varnamelist dl(1, dataname_r);
    switch (option) {
    case 1: break;
    case 2: dl.push_back(dataname_n); break;
    default: GMM_ASSERT1(false, "Penalized contact brick : invalid option");
    }

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region1);
  }


  //=========================================================================
  //  Add a contact condition with friction between two bodies discretized
  //  with nonmatching meshes.
  //=========================================================================

  size_type add_penalized_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_r,
   const std::string &dataname_friction_coeffs,
   size_type region1, size_type region2, int option,
   const std::string &dataname_lambda, const std::string &dataname_alpha,
   const std::string &dataname_wt1, const std::string &dataname_wt2) {

    pbrick pbr = new penalized_contact_nonmatching_meshes_brick
                     (region1, region2, false /* contact_only */, option);
    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, true)); // 0: U1U1
    tl.push_back(model::term_description(varname_u2, varname_u2, true)); // 1: U2U2
    tl.push_back(model::term_description(varname_u1, varname_u2, true)); // 2: U1U2
    tl.push_back(model::term_description(varname_u2, varname_u1, true)); // 3: U2U1

    model::varnamelist dl(1, dataname_r);
    switch (option) {
    case 1: break;
    case 2: case 3: dl.push_back(dataname_lambda); break;
    default: GMM_ASSERT1(false, "Penalized contact brick : invalid option");
    }
    dl.push_back(dataname_friction_coeffs);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt1.size() > 0) {
        dl.push_back(dataname_wt1);
        if (dataname_wt2.size() > 0)
          dl.push_back(dataname_wt2);
      }
    }

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region1);
  }


  // Computation of contact area and contact forces
  void compute_integral_contact_area_and_force
  (model &md, size_type indbrick, scalar_type &area,
   model_real_plain_vector &F) {

    pbrick pbr = md.brick_pointer(indbrick);
    const model::mimlist &ml = md.mimlist_of_brick(indbrick);
    const model::varnamelist &vl = md.varnamelist_of_brick(indbrick);
    const model::varnamelist &dl = md.datanamelist_of_brick(indbrick);
    size_type reg = md.region_of_brick(indbrick);

    GMM_ASSERT1(ml.size() == 1, "Wrong size");

    if (pbr->brick_name() == "Integral contact with rigid obstacle brick" ||
        pbr->brick_name() == "Integral contact and friction with rigid obstacle brick") {
      integral_contact_rigid_obstacle_brick *p
        = dynamic_cast<integral_contact_rigid_obstacle_brick *>
         (const_cast<virtual_brick *>(pbr.get()));
      GMM_ASSERT1(p, "Wrong type of brick");

      GMM_ASSERT1(vl.size() >= 2, "Wrong size");
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector &lambda = md.real_variable(vl[1]);
      const mesh_fem &mf_lambda = md.mesh_fem_of_variable(vl[1]);
      GMM_ASSERT1(dl.size() >= 1, "Wrong size");
      const model_real_plain_vector &obs = md.real_variable(dl[0]);
      const mesh_fem &mf_obs = md.mesh_fem_of_variable(dl[0]);

      //FIXME: use an adapted integration method
      area = asm_level_set_contact_area(*ml[0], mf_u, u, mf_obs, obs, reg, -1e-3,
                                        &mf_lambda, &lambda, 1e-1);

      gmm::resize(F, mf_u.nb_dof());
      asm_level_set_normal_source_term
        (F, *ml[0], mf_u, mf_obs, obs, mf_lambda, lambda, reg);
    }
    else if (pbr->brick_name() == "Integral penalized contact with rigid obstacle brick" ||
             pbr->brick_name() == "Integral penalized contact and friction with rigid "
                                  "obstacle brick") {
      penalized_contact_rigid_obstacle_brick *p
        = dynamic_cast<penalized_contact_rigid_obstacle_brick *>
         (const_cast<virtual_brick *>(pbr.get()));
      GMM_ASSERT1(p, "Wrong type of brick");
      GMM_ASSERT1(false, "Not implemented yet");
    }
    else if (pbr->brick_name() == "Integral contact between nonmatching meshes brick" ||
             pbr->brick_name() == "Integral contact and friction between nonmatching "
                                  "meshes brick") {
      integral_contact_nonmatching_meshes_brick *p
        = dynamic_cast<integral_contact_nonmatching_meshes_brick *>
         (const_cast<virtual_brick *>(pbr.get()));
      GMM_ASSERT1(p, "Wrong type of brick");

      GMM_ASSERT1(vl.size() == 3, "Wrong size");
      const model_real_plain_vector &u1 = md.real_variable(vl[0]);
      const model_real_plain_vector &u2 = md.real_variable(vl[1]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_u2 = md.mesh_fem_of_variable(vl[1]);
      const model_real_plain_vector &lambda = md.real_variable(vl[2]);
      const mesh_fem &mf_lambda = md.mesh_fem_of_variable(vl[2]);

      getfem::pfem pfem_proj = new_projected_fem(mf_u2, *ml[0], p->rg2, p->rg1);
      getfem::mesh_fem mf_u2_proj(mf_u1.linked_mesh(), mf_u1.linked_mesh().dim());
      mf_u2_proj.set_finite_element(mf_u1.linked_mesh().convex_index(), pfem_proj);

      std::vector<size_type> ind;
      mf_u2_proj.get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = mf_u2_proj.nb_basic_dof();
      model_real_plain_vector u2_proj(nbsub);

      if (mf_u2.is_reduced()) {
        gmm::csc_matrix<scalar_type> Esub(nbsub, nbdof2);
        gmm::copy(gmm::sub_matrix(mf_u2.extension_matrix(),
                                  SUBI, gmm::sub_interval(0, nbdof2)),
                  Esub);
        gmm::mult(Esub, u2, u2_proj);
      }
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      //FIXME: use an adapted integration method
      area = asm_nonmatching_meshes_contact_area
             (*ml[0], mf_u1, u1, mf_u2_proj, u2_proj, reg, -1e-3,
              &mf_lambda, &lambda, 1e-1);

      gmm::resize(F, mf_u1.nb_dof());
      asm_nonmatching_meshes_normal_source_term
        (F, *ml[0], mf_u1, mf_u2_proj, mf_lambda, lambda, reg);

      del_projected_fem(pfem_proj);
    }
    else if (pbr->brick_name() == "Integral penalized contact between nonmatching meshes brick" ||
             pbr->brick_name() == "Integral penalized contact and friction between nonmatching "
                                  "meshes brick") {
      penalized_contact_nonmatching_meshes_brick *p
        = dynamic_cast<penalized_contact_nonmatching_meshes_brick *>
         (const_cast<virtual_brick *>(pbr.get()));
      GMM_ASSERT1(p, "Wrong type of brick");
      GMM_ASSERT1(false, "Not implemented yet");
    }

  }

  //=========================================================================
  //
  //  Contact condition with a rigid obstacle : generic Nitsche's method
  //
  //=========================================================================


  class contact_nitsche_nonlinear_term : public nonlinear_elem_term {
    // Option:
    // 1 : rhs term
    // 2 : tangent term in main unknown (u)
    // 3 : tangent term in auxilliary variable (p)

  protected:
    base_small_vector u;      // tangential relative displacement
    scalar_type un;            // normal relative displacement (positive when
                               //  the first elas. body surface moves outwards)
    base_small_vector no, n;   // surface normal, pointing outwards with
                               // respect to the (first) elastic body
    scalar_type g, f_coeff;    // gap and friction coefficient

    base_small_vector aux1, wt, V, Pr, pgg, zeta;
    base_matrix GPr, grad;
    base_vector coeff;
    const model *md;
    const std::string *varname;
    const std::string *auxvarname;
    const mesh_fem &mf_u;       // mandatory
    const mesh_fem &mf_obs;     // mandatory
    const mesh_fem *pmf_coeff;
    const mesh_fem *mf_p;
    base_vector U, obs, friction_coeff, WT;
    dim_type N;
    size_type option;
    scalar_type gamma, gamma0, theta, alpha;
    base_tensor tG, tp, tpp, tbv, tpaux;
    mutable bgeot::multi_index sizes_;

    void adjust_tensor_size(void) {
      sizes_.resize(1); sizes_[0] = N;
      tG.adjust_sizes(sizes_);
      sizes_.resize(2); sizes_[0] = sizes_[1] = 1;
      switch (option) {
      case 1 : sizes_.resize(1); break;
      case 2 : case 3 :  break;
      }
      gmm::resize(grad, 1, N);
      u.resize(N); no.resize(N); n.resize(N);
      aux1.resize(1); wt.resize(N); V.resize(N); zeta.resize(N);
      gmm::resize(GPr, N, N); gmm::resize(Pr, N); gmm::resize(pgg, N);
    }

  public:
    const bgeot::multi_index &sizes(size_type cv) const {
      if (cv != size_type(-1))
        switch(option) {
        case 1:
          sizes_[0] = short_type(mf_u.nb_basic_dof_of_element(cv));
          break;
        case 2:
          sizes_[0] = sizes_[1]= short_type(mf_u.nb_basic_dof_of_element(cv));
          break;
        case 3:
          sizes_[0] = short_type(mf_u.nb_basic_dof_of_element(cv));
          sizes_[1] = short_type(mf_p->nb_basic_dof_of_element(cv));
          break;
        }
      return sizes_;
    }


    contact_nitsche_nonlinear_term
      (size_type option_, scalar_type gamma0_, scalar_type theta_,
       scalar_type alpha_, const model &md_, const std::string &varname_,
       const mesh_fem &mf_u_, const model_real_plain_vector &U_,
       const mesh_fem &mf_obs_,
       const model_real_plain_vector &obs_,
       const std::string &auxvarname_,
       const mesh_fem *pmf_p_ = 0,
       const mesh_fem *pmf_coeff_ = 0,
       const model_real_plain_vector *f_coeff_ = 0,
       const model_real_plain_vector *WT_ = 0)
      : md(&md_), varname(&varname_), auxvarname(&auxvarname_),
        mf_u(mf_u_), mf_obs(mf_obs_),
        pmf_coeff(pmf_coeff_), mf_p(pmf_p_), U(mf_u.nb_basic_dof()),
        obs(mf_obs.nb_basic_dof()),
        friction_coeff(0), option(option_),
        gamma0(gamma0_), theta(theta_), alpha(alpha_) {
      N = mf_u_.linked_mesh().dim();
      adjust_tensor_size();

      mf_u.extend_vector(U_, U);
      mf_obs.extend_vector(obs_, obs);

      if (!pmf_coeff)
        if (f_coeff_) f_coeff = (*f_coeff_)[0]; else f_coeff = scalar_type(0);
      else {
        friction_coeff.resize(pmf_coeff->nb_basic_dof());
        pmf_coeff->extend_vector(*f_coeff_, friction_coeff);
      }
      if (WT_) {
        WT.resize(mf_u.nb_basic_dof());
        mf_u_.extend_vector(*WT_, WT);
      }
    }


    void compute(fem_interpolation_context &ctx, bgeot::base_tensor &t) {

      md->compute_Neumann_terms(1, *varname, mf_u, U, ctx, n, tG);
      for (size_type i = 0; i < N; ++i)
        zeta[i] = tG[i]
          + ((g-un+alpha*un) * no[i] + alpha*wt[i] - alpha*u[i] ) / gamma;
      if ((option == 1) || (theta != scalar_type(0))) {
        coupled_projection(zeta, no, f_coeff, Pr);
        gmm::add(Pr, gmm::scaled(tG.as_vector(), -scalar_type(1)), pgg);
      }

      switch (option) {
      case 1:
        {
          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          short_type nbdofu = sizes_[0];
          if (theta != scalar_type(0)) {
            sizes_.resize(2);
            sizes_[1] = N;
            tp.adjust_sizes(sizes_);
            sizes_.resize(1);
            md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
          }
          for (size_type i = 0; i < nbdofu; ++i) {
            t[i] = scalar_type(0);
            for (size_type j = 0; j < N; ++j) {
              if (theta != scalar_type(0))
                t[i] -= gamma*pgg[j]*theta*tp(i,j);
              if (qmult == 1) t[i] += Pr[j]*tbv(i,j);
            }
            if (qmult > 1) t[i] += Pr[i%N] * tbv(i/N,0);
          }
        }
        break;

      case 2:
        {
          short_type nbdofu = sizes_[1];
          sizes_[1] = N;
          tp.adjust_sizes(sizes_);
          sizes_[1] = nbdofu;
          md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
          if (theta != scalar_type(0)) {
            sizes_.resize(3); sizes_[2] = N;
            tpp.adjust_sizes(sizes_);
            sizes_.resize(2);
            md->compute_Neumann_terms(3, *varname, mf_u, U, ctx, n, tpp);
          }

          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          coupled_projection_grad(zeta, no, f_coeff, GPr);

          for (size_type i = 0; i < nbdofu; ++i)
            for (size_type j = 0; j < nbdofu; ++j) {
              scalar_type res(0);
              for (size_type k = 0; k < N; ++k) {
                if (theta != scalar_type(0))
                  res -= gamma * theta * tp(i,k) * tp(j,k);
                scalar_type tbvvi(0), tbvvjn(0);
                if (qmult == 1) {
                  tbvvi = tbv(i,k);
                  for (size_type l = 0; l < N; ++l) tbvvjn += no[l]*tbv(j,l);
                } else {
                  tbvvi = ((i%N)==k) ? tbv(i/N,0) : scalar_type(0);
                  tbvvjn = no[j%N]*tbv(j/N,0);
                }
                for (size_type l = 0; l < N; ++l) {
                  scalar_type tbvvj(0);
                  if (qmult == 1)
                    tbvvj = tbv(j,l);
                  else
                    tbvvj=(((j%N)==l) ? tbv(j/N,0):scalar_type(0));
                  res += GPr(k,l)
                    * (gamma*tp(j,l) - alpha*tbvvj
                       - (scalar_type(1)-alpha)*no[l]*tbvvjn)
                    * (theta * tp(i,k) - tbvvi/gamma);
                }

                if (theta != scalar_type(0))
                  res += theta*gamma*pgg[k] * tpp(i,j,k);
              }
              t(i,j) = res;
            }
        }
        break;

      case 3:
        {
          short_type nbdofu = sizes_[0];
          short_type nbdofp = sizes_[1];
          sizes_[0] = nbdofp; sizes_[1] = N;
          tpaux.adjust_sizes(sizes_);
          sizes_[0] = nbdofu; sizes_[1] = nbdofp;
          md->compute_auxilliary_Neumann_terms(2, *varname, mf_u, U,
                                               *auxvarname, ctx, n, tpaux);

          if (theta != scalar_type(0)) {
            sizes_[1] = N;
            tp.adjust_sizes(sizes_);
            sizes_[1] = nbdofp;
            md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
            sizes_.resize(3); sizes_[2] = N;
            tpp.adjust_sizes(sizes_);
            sizes_.resize(2);
            md->compute_auxilliary_Neumann_terms(3, *varname, mf_u, U,
                                                 *auxvarname, ctx, n, tpp);
          }

          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          coupled_projection_grad(zeta, no, f_coeff, GPr);

          for (size_type i = 0; i < nbdofu; ++i)
            for (size_type j = 0; j < nbdofp; ++j) {
              scalar_type res(0);
              for (size_type k = 0; k < N; ++k) {
                if (theta != scalar_type(0))
                  res -= gamma * theta * tp(i,k) * tpaux(j,k);
                scalar_type gttpik(0), tbvvi(0);
                if (theta != scalar_type(0)) gttpik = gamma*theta*tp(i,k);
                if (qmult == 1) tbvvi = tbv(i,k);
                else tbvvi=(((i%N)==k) ? tbv(i/N,0):scalar_type(0));
                for (size_type l = 0; l < N; ++l)
                  res += GPr(k,l) * tpaux(j,l) * (gttpik - tbvvi);
                if (theta != scalar_type(0))
                  res += theta*gamma*pgg[k] * tpp(i,j,k);
              }
              t(i,j) = res;
            }
        }
        break;

      default : GMM_ASSERT1(false, "Invalid option");
      }
    }


    void prepare(fem_interpolation_context& ctx, size_type nb) {

      size_type cv = ctx.convex_num();

      switch (nb) { // last is computed first
      case 1 : // calculate [u] and [un] interpolating [U] and [WT] on [mf_u]
        slice_vector_on_basic_dof_of_element(mf_u, U, cv, coeff);
        ctx.pf()->interpolation(ctx, coeff, u, N);
        un = gmm::vect_sp(u, no);
        if (gmm::vect_size(WT) == gmm::vect_size(U)) {
          slice_vector_on_basic_dof_of_element(mf_u, WT, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, wt, N);
          wt -= gmm::vect_sp(wt, no) * no;
        }
        // computation of h for gamma = gamma0*h
        scalar_type emax, emin; gmm::condition_number(ctx.K(),emax,emin);
        gamma = gamma0 * emax * sqrt(scalar_type(N));
        break;

      case 2 : // calculate [g], [n] and [no] interpolating [obs] on [mf_obs]
        slice_vector_on_basic_dof_of_element(mf_obs, obs, cv, coeff);
        ctx.pf()->interpolation_grad(ctx, coeff, grad, 1);
        gmm::copy(gmm::mat_row(grad, 0), no);
        no /= -gmm::vect_norm2(no);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        g = aux1[0];
        n = bgeot::compute_normal(ctx, ctx.face_num());
        n /= gmm::vect_norm2(n);
        break;

      case 3 : // calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
        if (pmf_coeff) {
          slice_vector_on_basic_dof_of_element(*pmf_coeff, friction_coeff, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          f_coeff = aux1[0];
        }
        break;

      default : GMM_ASSERT1(false, "Invalid option");
      }
    }
  };



  void asm_Nitsche_contact_rigid_obstacle_rhs
    (model_real_plain_vector &R, const mesh_im &mim, const model &md,
     const std::string &varname,
     const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
     const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
     const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
     const model_real_plain_vector *WT,
     scalar_type gamma0, scalar_type theta, scalar_type alpha,
     const mesh_region &rg) {

    contact_nitsche_nonlinear_term
      nterm(1, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, "", 0, pmf_coeff, f_coeff, WT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly assem("V(#1)+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }


  template<typename MAT>
  void asm_Nitsche_contact_rigid_obstacle_tangent_matrix
  (MAT &K, const mesh_im &mim, const model &md, const std::string &varname,
   const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
   const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
   const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
   const model_real_plain_vector *WT,
   scalar_type gamma0, scalar_type theta, scalar_type alpha,
   const mesh_region &rg) {

    contact_nitsche_nonlinear_term
      nterm(2, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, "", 0, pmf_coeff, f_coeff, WT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly
      assem("M(#1,#1)+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.assembly(rg);
  }

  template<typename MAT>
  void asm_Nitsche_contact_rigid_obstacle_tangent_matrix_auxilliary
  (MAT &K, const mesh_im &mim, const model &md, const std::string &varname,
   const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
   const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
   const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
   const model_real_plain_vector *WT,
   scalar_type gamma0, scalar_type theta, scalar_type alpha,
   const std::string &auxvarname, const getfem::mesh_fem &mf_p,
   const mesh_region &rg) {

    contact_nitsche_nonlinear_term
      nterm(3, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, auxvarname, &mf_p, pmf_coeff, f_coeff, WT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";
    const std::string p_fem = pmf_coeff ? "#4" : "#3";

    getfem::generic_assembly
      assem("M(#1,"+p_fem+")+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.assembly(rg);
  }


  struct Nitsche_contact_rigid_obstacle_brick : public virtual_brick {

    scalar_type theta;
    bool contact_only;

    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      // Integration method
      GMM_ASSERT1(mims.size() == 1, "Nitsche contact with rigid obstacle "
                  "bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];


      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);

      // Data : obs, r, theta, [alpha,] [WT]
      GMM_ASSERT1(dl.size() >= (contact_only ? 2:3),
                  "Wrong number of data for Nitsche "
                  "contact with rigid obstacle brick");

      const model_real_plain_vector &obs = md.real_variable(dl[0]);
      const mesh_fem &mf_obs = md.mesh_fem_of_variable(dl[0]);
      size_type sl = gmm::vect_size(obs) * mf_obs.get_qdim() / mf_obs.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vgamma0 = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vgamma0) == 1,
                  "Parameter gamma0 should be a scalar");
      scalar_type gamma0 = vgamma0[0];

      const model_real_plain_vector *f_coeff = 0;
      const mesh_fem *pmf_coeff = 0;

      if (!contact_only) {
        f_coeff = &(md.real_variable(dl[2]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[2]);
        sl = gmm::vect_size(*f_coeff);
        if (pmf_coeff)
          { sl*= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1, "the data corresponding to the friction "
                    "coefficient has not the right format");
      }

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4) {
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3])) == 1,
                    "Parameter alpha should be a scalar");
        alpha = md.real_variable(dl[3])[0];
      }

      const model_real_plain_vector *WT
        = (!contact_only && dl.size()>=5) ? &(md.real_variable(dl[4])) : 0;


      GMM_ASSERT1(matl.size() == vl.size(), "Wrong number of terms for "
                  "Nitsche contact with rigid obstacle brick");


      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Nitsche contact with rigid obstacle tangent term");
        gmm::clear(matl[0]);
        asm_Nitsche_contact_rigid_obstacle_tangent_matrix
          (matl[0], mim, md, vl[0], mf_u, u, mf_obs, obs,  pmf_coeff,
           f_coeff, WT, gamma0, theta, alpha, rg);

        for (size_type i = 1; i < vl.size(); ++i) { // Auxilliary variables
          gmm::clear(matl[i]);
          asm_Nitsche_contact_rigid_obstacle_tangent_matrix_auxilliary
            (matl[i], mim, md, vl[0], mf_u, u, mf_obs, obs, pmf_coeff,
             f_coeff, WT, gamma0, theta, alpha, vl[i],
             md.mesh_fem_of_variable(vl[i]), rg);
        }
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        asm_Nitsche_contact_rigid_obstacle_rhs
          (vecl[0], mim, md, vl[0], mf_u, u, mf_obs, obs,  pmf_coeff,
           f_coeff, WT, gamma0, theta, alpha, rg);
      }

    }

    Nitsche_contact_rigid_obstacle_brick(scalar_type theta_, bool nofriction) {
      theta = theta_;
      contact_only = nofriction;
      bool co = (theta_ == scalar_type(1)) && nofriction;
      set_flags("Integral Nitsche contact and friction with rigid "
                "obstacle brick",
                false /* is linear*/, co /* is symmetric */,
                co /* is coercive */, true /* is real */,
                false /* is complex */, false /* compute each time */,
                false /* has a Neumann term */);
    }

  };


  size_type add_Nitsche_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_gamma0,
   scalar_type theta,
   const std::string &dataname_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt,
   size_type region) {

    bool nofriction = (dataname_friction_coeff.size() == 0);
    pbrick pbr = new Nitsche_contact_rigid_obstacle_brick(theta, nofriction);

    bool co = (theta == scalar_type(1)) && nofriction;
    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, co));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_gamma0);
    if (!nofriction) dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size() > 0) dl.push_back(dataname_wt);
    }

    model::varnamelist vl(1, varname_u);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname_u, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname_u, aux_vars[i], false));
    }

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }




#ifdef EXPERIMENTAL_PURPOSE_ONLY


  //=========================================================================
  //
  //  Contact condition with a rigid obstacle : generic Nitsche's method
  //  Experimental for midpoint scheme
  //
  //=========================================================================


  class contact_nitsche_nonlinear_term_midpoint : public nonlinear_elem_term {
    // Option:
    // 1 : rhs term
    // 2 : tangent term in main unknown (u)
    // 3 : tangent term in auxilliary variable (p)

  protected:
    base_small_vector u;      // tangential relative displacement
    scalar_type un, wn;       // normal relative displacement (positive when
                               //  the first elas. body surface moves outwards)
    base_small_vector no, n;   // surface normal, pointing outwards with
                               // respect to the (first) elastic body
    scalar_type g, f_coeff;    // gap and friction coefficient

    base_small_vector aux1, wt, V, Pr, pgg, zeta;
    base_matrix GPr, grad;
    base_vector coeff;
    const model *md;
    const std::string *varname;
    const std::string *auxvarname;
    const mesh_fem &mf_u;       // mandatory
    const mesh_fem &mf_obs;     // mandatory
    const mesh_fem *pmf_coeff;
    const mesh_fem *mf_p;
    base_vector U, obs, friction_coeff, WT, UPLUSWT;
    dim_type N;
    size_type option;
    scalar_type gamma, gamma0, theta, alpha;
    base_tensor tG, tp, tpp, tbv, tpaux;
    mutable bgeot::multi_index sizes_;
    size_type option_midpoint;

    void adjust_tensor_size(void) {
      sizes_.resize(1); sizes_[0] = N;
      tG.adjust_sizes(sizes_);
      sizes_.resize(2); sizes_[0] = sizes_[1] = 1;
      switch (option) {
      case 1 : sizes_.resize(1); break;
      case 2 : case 3 :  break;
      }
      gmm::resize(grad, 1, N);
      u.resize(N); no.resize(N); n.resize(N);
      aux1.resize(1); wt.resize(N); V.resize(N); zeta.resize(N);
      gmm::resize(GPr, N, N); gmm::resize(Pr, N); gmm::resize(pgg, N);
    }

  public:
    const bgeot::multi_index &sizes(size_type cv) const {
      if (cv != size_type(-1))
        switch(option) {
        case 1:
          sizes_[0] = short_type(mf_u.nb_basic_dof_of_element(cv));
          break;
        case 2:
          sizes_[0] = sizes_[1]= short_type(mf_u.nb_basic_dof_of_element(cv));
          break;
        case 3:
          sizes_[0] = short_type(mf_u.nb_basic_dof_of_element(cv));
          sizes_[1] = short_type(mf_p->nb_basic_dof_of_element(cv));
          break;
        }
      return sizes_;
    }


    contact_nitsche_nonlinear_term_midpoint
      (size_type option_, scalar_type gamma0_, scalar_type theta_,
       scalar_type alpha_, const model &md_, const std::string &varname_,
       const mesh_fem &mf_u_, const model_real_plain_vector &U_,
       const mesh_fem &mf_obs_,
       const model_real_plain_vector &obs_,
       const std::string &auxvarname_,
       const mesh_fem *pmf_p_ = 0,
       const mesh_fem *pmf_coeff_ = 0,
       const model_real_plain_vector *f_coeff_ = 0,
       const model_real_plain_vector *WT_ = 0, size_type option_midpoint_ = 1)
      : md(&md_), varname(&varname_), auxvarname(&auxvarname_),
        mf_u(mf_u_), mf_obs(mf_obs_),
        pmf_coeff(pmf_coeff_), mf_p(pmf_p_), U(mf_u.nb_basic_dof()),
        obs(mf_obs.nb_basic_dof()),
        friction_coeff(0), option(option_),
        gamma0(gamma0_), theta(theta_), alpha(alpha_),
        option_midpoint(option_midpoint_) {
      N = mf_u_.linked_mesh().dim();
      adjust_tensor_size();

      mf_u.extend_vector(U_, U);
      mf_obs.extend_vector(obs_, obs);

      if (!pmf_coeff)
        if (f_coeff_) f_coeff = (*f_coeff_)[0]; else f_coeff = scalar_type(0);
      else {
        friction_coeff.resize(pmf_coeff->nb_basic_dof());
        pmf_coeff->extend_vector(*f_coeff_, friction_coeff);
      }
      if (WT_) {
        WT.resize(mf_u.nb_basic_dof());
        mf_u_.extend_vector(*WT_, WT);
        UPLUSWT.resize(mf_u.nb_basic_dof());
        gmm::add(U, gmm::scaled(WT, -scalar_type(1)/scalar_type(2)), UPLUSWT);
      }
    }


    void compute(fem_interpolation_context &ctx, bgeot::base_tensor &t) {

      md->compute_Neumann_terms(1, *varname, mf_u, WT, ctx, n, tG);

      scalar_type Pw = wn - gamma * gmm::vect_sp(tG.as_vector(), no);
      cout << "Pw = " << Pw << endl;


      if (option_midpoint == 2)
        md->compute_Neumann_terms(1, *varname, mf_u, UPLUSWT, ctx, n, tG);
      else
        md->compute_Neumann_terms(1, *varname, mf_u, U, ctx, n, tG);
      for (size_type i = 0; i < N; ++i)
        if (option_midpoint == 2)
          zeta[i] = tG[i]
            + ((g-un+wn/scalar_type(2)+alpha*(un-wn/scalar_type(2))) * no[i]
               + alpha*wt[i] - alpha*u[i] ) / gamma;
        else
          zeta[i] = tG[i]
            + ((g-un+alpha*un) * no[i] + alpha*wt[i] - alpha*u[i] ) / gamma;

      if (option_midpoint == 2)
        md->compute_Neumann_terms(1, *varname, mf_u, U, ctx, n, tG);

      if ((option == 1) || (theta != scalar_type(0))) {
        coupled_projection(zeta, no, f_coeff, Pr);
        gmm::add(Pr, gmm::scaled(tG.as_vector(), -scalar_type(1)), pgg);
      }

      switch (option) {
      case 1:
        {
          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          short_type nbdofu = sizes_[0];
          if (theta != scalar_type(0)) {
            sizes_.resize(2);
            sizes_[1] = N;
            tp.adjust_sizes(sizes_);
            sizes_.resize(1);
            md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
          }
          for (size_type i = 0; i < nbdofu; ++i) {
            t[i] = scalar_type(0);
            for (size_type j = 0; j < N; ++j) {
              if (theta != scalar_type(0))
                t[i] -= gamma*pgg[j]*theta*tp(i,j);
              if (qmult == 1) t[i] += Pr[j]*tbv(i,j);
            }
            if (qmult > 1) t[i] += Pr[i%N] * tbv(i/N,0);
          }
        }
        break;

      case 2:
        {
          short_type nbdofu = sizes_[1];
          sizes_[1] = N;
          tp.adjust_sizes(sizes_);
          sizes_[1] = nbdofu;
          if (option_midpoint == 2)
            md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
          else
            md->compute_Neumann_terms(2, *varname, mf_u, UPLUSWT, ctx, n, tp);
          if (theta != scalar_type(0)) {
            sizes_.resize(3); sizes_[2] = N;
            tpp.adjust_sizes(sizes_);
            sizes_.resize(2);
            if (option_midpoint == 1)
              md->compute_Neumann_terms(3, *varname, mf_u, UPLUSWT, ctx,n,tpp);
            else
              md->compute_Neumann_terms(3, *varname, mf_u, U, ctx, n, tpp);
          }

          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          coupled_projection_grad(zeta, no, f_coeff, GPr);

          for (size_type i = 0; i < nbdofu; ++i)
            for (size_type j = 0; j < nbdofu; ++j) {
              scalar_type res(0);
              for (size_type k = 0; k < N; ++k) {
                if (theta != scalar_type(0))
                  res -= gamma * theta * tp(i,k) * tp(j,k);
                scalar_type tbvvi(0), tbvvjn(0);
                if (qmult == 1) {
                  tbvvi = tbv(i,k);
                  for (size_type l = 0; l < N; ++l) tbvvjn += no[l]*tbv(j,l);
                } else {
                  tbvvi = ((i%N)==k) ? tbv(i/N,0) : scalar_type(0);
                  tbvvjn = no[j%N]*tbv(j/N,0);
                }
                for (size_type l = 0; l < N; ++l) {
                  scalar_type tbvvj(0);
                  if (qmult == 1)
                    tbvvj = tbv(j,l);
                  else
                    tbvvj=(((j%N)==l) ? tbv(j/N,0):scalar_type(0));
                  res += GPr(k,l)
                    * (gamma*tp(j,l) - alpha*tbvvj
                       - (scalar_type(1)-alpha)*no[l]*tbvvjn)
                    * (theta * tp(i,k) - tbvvi/gamma);
                }

                if (theta != scalar_type(0))
                  res += theta*gamma*pgg[k] * tpp(i,j,k);
              }
              t(i,j) = res;
            }
        }
        break;

      case 3:
        {
          short_type nbdofu = sizes_[0];
          short_type nbdofp = sizes_[1];
          sizes_[0] = nbdofp; sizes_[1] = N;
          tpaux.adjust_sizes(sizes_);
          sizes_[0] = nbdofu; sizes_[1] = nbdofp;
          if (option_midpoint == 2)
            md->compute_auxilliary_Neumann_terms(2, *varname, mf_u, UPLUSWT,
                                                 *auxvarname, ctx, n, tpaux);
          else
            md->compute_auxilliary_Neumann_terms(2, *varname, mf_u, U,
                                                 *auxvarname, ctx, n, tpaux);

          if (theta != scalar_type(0)) {
            sizes_[1] = N;
            tp.adjust_sizes(sizes_);
            sizes_[1] = nbdofp;
            if (option_midpoint == 2)
              md->compute_Neumann_terms(2, *varname, mf_u, UPLUSWT, ctx, n,tp);
            else
              md->compute_Neumann_terms(2, *varname, mf_u, U, ctx, n, tp);
            sizes_.resize(3); sizes_[2] = N;
            tpp.adjust_sizes(sizes_);
            sizes_.resize(2);
            if (option_midpoint == 2)
              md->compute_auxilliary_Neumann_terms(3, *varname, mf_u, UPLUSWT,
                                                   *auxvarname, ctx, n, tpp);
            else
              md->compute_auxilliary_Neumann_terms(3, *varname, mf_u, U,
                                                   *auxvarname, ctx, n, tpp);
          }

          ctx.pf()->real_base_value(ctx, tbv);
          size_type qmult = N / ctx.pf()->target_dim();
          coupled_projection_grad(zeta, no, f_coeff, GPr);

          for (size_type i = 0; i < nbdofu; ++i)
            for (size_type j = 0; j < nbdofp; ++j) {
              scalar_type res(0);
              for (size_type k = 0; k < N; ++k) {
                if (theta != scalar_type(0))
                  res -= gamma * theta * tp(i,k) * tpaux(j,k);
                scalar_type gttpik(0), tbvvi(0);
                if (theta != scalar_type(0)) gttpik = gamma*theta*tp(i,k);
                if (qmult == 1) tbvvi = tbv(i,k);
                else tbvvi=(((i%N)==k) ? tbv(i/N,0):scalar_type(0));
                for (size_type l = 0; l < N; ++l)
                  res += GPr(k,l) * tpaux(j,l) * (gttpik - tbvvi);
                if (theta != scalar_type(0))
                  res += theta*gamma*pgg[k] * tpp(i,j,k);
              }
              t(i,j) = res;
            }
        }
        break;

      default : GMM_ASSERT1(false, "Invalid option");
      }

      switch (option_midpoint) {
      case 1:
        gmm::scale(t.as_vector(), gmm::Heaviside(Pw));
        break;
      case 2:
        gmm::scale(t.as_vector(),
                   (scalar_type(1) - gmm::Heaviside(Pw)));
        break;
      default:
        GMM_ASSERT1(false, "Wrong option");
      }

    }


    void prepare(fem_interpolation_context& ctx, size_type nb) {

      size_type cv = ctx.convex_num();

      switch (nb) { // last is computed first
      case 1 : // calculate [u] and [un] interpolating [U] and [WT] on [mf_u]
        slice_vector_on_basic_dof_of_element(mf_u, U, cv, coeff);
        ctx.pf()->interpolation(ctx, coeff, u, N);
        un = gmm::vect_sp(u, no);
        if (gmm::vect_size(WT) == gmm::vect_size(U)) {
          slice_vector_on_basic_dof_of_element(mf_u, WT, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, wt, N);
          wn = gmm::vect_sp(wt, no);
          wt -= gmm::vect_sp(wt, no) * no;
        }
        // computation of h for gamma = gamma0*h
        scalar_type emax, emin; gmm::condition_number(ctx.K(),emax,emin);
        gamma = gamma0 * emax * sqrt(scalar_type(N));
        break;

      case 2 : // calculate [g], [n] and [no] interpolating [obs] on [mf_obs]
        slice_vector_on_basic_dof_of_element(mf_obs, obs, cv, coeff);
        ctx.pf()->interpolation_grad(ctx, coeff, grad, 1);
        gmm::copy(gmm::mat_row(grad, 0), no);
        no /= -gmm::vect_norm2(no);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        g = aux1[0];
        n = bgeot::compute_normal(ctx, ctx.face_num());
        n /= gmm::vect_norm2(n);
        break;

      case 3 : // calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
        if (pmf_coeff) {
          slice_vector_on_basic_dof_of_element(*pmf_coeff, friction_coeff, cv, coeff);
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          f_coeff = aux1[0];
        }
        break;

      default : GMM_ASSERT1(false, "Invalid option");
      }
    }
  };



  void asm_Nitsche_contact_rigid_obstacle_rhs_midpoint
    (model_real_plain_vector &R, const mesh_im &mim, const model &md,
     const std::string &varname,
     const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
     const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
     const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
     const model_real_plain_vector *WT,
     scalar_type gamma0, scalar_type theta, scalar_type alpha,
     const mesh_region &rg, size_type option) {

    contact_nitsche_nonlinear_term_midpoint
      nterm(1, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, "", 0, pmf_coeff, f_coeff, WT, option);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly assem("V(#1)+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }


  template<typename MAT>
  void asm_Nitsche_contact_rigid_obstacle_tangent_matrix_midpoint
  (MAT &K, const mesh_im &mim, const model &md, const std::string &varname,
   const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
   const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
   const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
   const model_real_plain_vector *WT,
   scalar_type gamma0, scalar_type theta, scalar_type alpha,
   const mesh_region &rg, size_type option) {

    contact_nitsche_nonlinear_term_midpoint
      nterm(2, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, "", 0, pmf_coeff, f_coeff, WT, option);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly
      assem("M(#1,#1)+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.assembly(rg);
  }

  template<typename MAT>
  void asm_Nitsche_contact_rigid_obstacle_tangent_matrix_auxilliary_midpoint
  (MAT &K, const mesh_im &mim, const model &md, const std::string &varname,
   const getfem::mesh_fem &mf_u, const model_real_plain_vector &U,
   const getfem::mesh_fem &mf_obs, const model_real_plain_vector &obs,
   const getfem::mesh_fem *pmf_coeff, const model_real_plain_vector *f_coeff,
   const model_real_plain_vector *WT,
   scalar_type gamma0, scalar_type theta, scalar_type alpha,
   const std::string &auxvarname, const getfem::mesh_fem &mf_p,
   const mesh_region &rg, size_type option) {

    contact_nitsche_nonlinear_term_midpoint
      nterm(3, gamma0, theta, alpha, md, varname, mf_u, U, mf_obs,
            obs, auxvarname, &mf_p, pmf_coeff, f_coeff, WT, option);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";
    const std::string p_fem = pmf_coeff ? "#4" : "#3";

    getfem::generic_assembly
      assem("M(#1,"+p_fem+")+=comp(NonLin$1(#1,"+aux_fems+"));");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.assembly(rg);
  }


  struct Nitsche_midpoint_contact_rigid_obstacle_brick : public virtual_brick {

    scalar_type theta;
    bool contact_only;
    size_type option;

    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {

      // Integration method
      GMM_ASSERT1(mims.size() == 1, "Nitsche contact with rigid obstacle "
                  "bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];


      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);

      // Data : obs, r, theta, [alpha,] [WT]
      GMM_ASSERT1(dl.size() >= (contact_only ? 2:3),
                  "Wrong number of data for Nitsche "
                  "contact with rigid obstacle brick");

      const model_real_plain_vector &obs = md.real_variable(dl[0]);
      const mesh_fem &mf_obs = md.mesh_fem_of_variable(dl[0]);
      size_type sl = gmm::vect_size(obs) * mf_obs.get_qdim() / mf_obs.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vgamma0 = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vgamma0) == 1,
                  "Parameter gamma0 should be a scalar");
      scalar_type gamma0 = vgamma0[0];

      const model_real_plain_vector *f_coeff = 0;
      const mesh_fem *pmf_coeff = 0;

      if (!contact_only) {
        f_coeff = &(md.real_variable(dl[2]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[2]);
        sl = gmm::vect_size(*f_coeff);
        if (pmf_coeff)
          { sl*= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1, "the data corresponding to the friction "
                    "coefficient has not the right format");
      }

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4) {
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3])) == 1,
                    "Parameter alpha should be a scalar");
        alpha = md.real_variable(dl[3])[0];
      }

      const model_real_plain_vector *WT
        = (dl.size()>=5) ? &(md.real_variable(dl[4])) : 0;

      GMM_ASSERT1(matl.size() == vl.size(), "Wrong number of terms for "
                  "Nitsche contact with rigid obstacle brick");


      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Nitsche contact with rigid obstacle tangent term");
        gmm::clear(matl[0]);
        asm_Nitsche_contact_rigid_obstacle_tangent_matrix_midpoint
          (matl[0], mim, md, vl[0], mf_u, u, mf_obs, obs,  pmf_coeff,
           f_coeff, WT, gamma0, theta, alpha, rg, option);

        for (size_type i = 1; i < vl.size(); ++i) { // Auxilliary variables
          gmm::clear(matl[i]);
          asm_Nitsche_contact_rigid_obstacle_tangent_matrix_auxilliary_midpoint
            (matl[i], mim, md, vl[0], mf_u, u, mf_obs, obs, pmf_coeff,
             f_coeff, WT, gamma0, theta, alpha, vl[i],
             md.mesh_fem_of_variable(vl[i]), rg, option);
        }
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        asm_Nitsche_contact_rigid_obstacle_rhs_midpoint
          (vecl[0], mim, md, vl[0], mf_u, u, mf_obs, obs,  pmf_coeff,
           f_coeff, WT, gamma0, theta, alpha, rg, option);
      }
    }

    Nitsche_midpoint_contact_rigid_obstacle_brick(scalar_type theta_, bool nofriction, size_type option_) {
      theta = theta_;
      contact_only = nofriction;
      option = option_;
      bool co = (theta_ == scalar_type(1)) && nofriction;
      set_flags("Integral Nitsche contact and friction with rigid "
                "obstacle brick",
                false /* is linear*/, co /* is symmetric */,
                co /* is coercive */, true /* is real */,
                false /* is complex */, false /* compute each time */,
                false /* has a Neumann term */);
    }

  };


  size_type add_Nitsche_midpoint_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_gamma0,
   scalar_type theta,
   const std::string &dataname_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt,
   size_type region, size_type option) {

    bool nofriction = (dataname_friction_coeff.size() == 0);
    pbrick pbr = new Nitsche_midpoint_contact_rigid_obstacle_brick(theta, nofriction, option);

    bool co = (theta == scalar_type(1)) && nofriction;
    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, co));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_gamma0);
    if (!nofriction) dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt.size() > 0) dl.push_back(dataname_wt);
    }

    model::varnamelist vl(1, varname_u);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname_u, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname_u, aux_vars[i], false));
    }

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }




#endif












  //=========================================================================
  //
  //  Fictitious domain contact condition (HPP) : generic Nitsche's method
  //
  //=========================================================================


  struct Nitsche_fictitious_domain_contact_brick : public virtual_brick {

    scalar_type theta;
    bool contact_only;

    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type /* region */,
                                        build_version version) const {
					  
      cout << "begining assembly" << endl;

      // Integration method
      GMM_ASSERT1(mims.size() == 1, "Nitsche fictitious domain contact "
                  "bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];
      const mesh &m = mim.linked_mesh();
      size_type N = m.dim();

      GMM_ASSERT1(vl.size() <= 2, "Auxilliary variable not taken into "
                  "account for the moment");

      GMM_ASSERT1(vl.size() >= 2, "Nitsche fictitious domain contact "
                  "bricks need two variables");

      const model_real_plain_vector &U1 = md.real_variable(vl[0]);
      const mesh_fem &mf_u1 = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector &U2 = md.real_variable(vl[1]);
      const mesh_fem &mf_u2 = md.mesh_fem_of_variable(vl[1]);


      GMM_ASSERT1(dl.size() > 2, "Nitsche fictitious domain contact "
                  "bricks need at least 2 data");

      const model_real_plain_vector &D1 = md.real_variable(dl[0]);
      const mesh_fem &mf_d1 = md.mesh_fem_of_variable(dl[0]);
      const model_real_plain_vector &D2 = md.real_variable(dl[1]);
      const mesh_fem &mf_d2 = md.mesh_fem_of_variable(dl[1]);

      const model_real_plain_vector &GAMMA0 = md.real_variable(dl[2]);
      GMM_ASSERT1(GAMMA0.size() == 1, "Gamma0 should be a scalar parameter");
      scalar_type gamma0 = GAMMA0[0];

      scalar_type f_coeff(0), alpha(0);
      const model_real_plain_vector *WT1 = 0, *WT2 = 0;
      if (dl.size() > 3) {
        const model_real_plain_vector &FRICT = md.real_variable(dl[3]);
        GMM_ASSERT1(FRICT.size() == 1, "The friction coefficient should "
                    "be a scalar parameter");
        f_coeff = FRICT[0];

        if (dl.size() > 4) {
          const model_real_plain_vector &ALPHA = md.real_variable(dl[4]);
          GMM_ASSERT1(ALPHA.size() == 1, "Alpha should be a scalar parameter");
          alpha = ALPHA[0];

          if (dl.size() > 6) {
            WT1 = &(md.real_variable(dl[5]));
            GMM_ASSERT1(&mf_u1 == &(md.mesh_fem_of_variable(dl[5])),
                        "wt1 should be described on the same fem than u1");
            WT2 = &(md.real_variable(dl[6]));
            GMM_ASSERT1(&mf_u2 == &(md.mesh_fem_of_variable(dl[6])),
                        "wt2 should be described on the same fem than u2");
          }
        }
      }




      GMM_ASSERT1(&(mf_u1.linked_mesh()) == &m && &(mf_u2.linked_mesh()) == &m
                  && &(mf_d1.linked_mesh()) == &m
                  && &(mf_d2.linked_mesh()) == &m,
                  "All data and variables should be defined on the same mesh");

      cout << "Computing projection ..." << endl;

      bgeot::rtree tree;


	    
	    
      
      for (dal::bv_visitor cv(mf_d2.convex_index()); !cv.finished(); ++cv) {
	base_node min,max;
//      base_node min = m.points_of_convex(cv)[0], max = min;
//      for (size_type i = 1; i <  m.nb_points_of_convex(cv); // pourquoi ça?
//          ++i) {
//       cout << " cv = " << cv << ", min = " << min << ", max = " << max << endl;
//       for (size_type k = 0; k < N; ++k) {
//            const base_node &x = m.points_of_convex(cv)[k];
//           min[k] = std::min(min[k], x[k]);
//           max[k] = std::max(max[k], x[k]);
//         }
//       cout << " cv = " << cv << ", min = " << min << ", max = " << max << endl;
//               }
//    
//       for (size_type k = 0; k < N; ++k) {
//         min[k] -= (max[k] - min[k]) / 5.;
//         max[k] += (max[k] - min[k]) / 5.;
//          }
       scalar_type EPS = 1E-13;
       bounding_box(min, max, mf_d2.linked_mesh().points_of_convex(cv),
                    mf_d2.linked_mesh().trans_of_convex(cv));
       for (unsigned k=0; k < min.size(); ++k) { min[k]-=EPS; max[k]+=EPS; }

 
 
 
 
 
        tree.add_box(min, max, cv);
//	       cout << " cv = " << cv << ", min = " << min << ", max = " << max << endl;


      }

       cout << "Projection computed." << endl;

      if (version & model::BUILD_MATRIX) {
        gmm::clear(matl[0]);
        gmm::clear(matl[1]);
        gmm::clear(matl[2]);
        gmm::clear(matl[3]);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        gmm::clear(vecl[1]);
        gmm::clear(vecl[2]);
        gmm::clear(vecl[3]);
      }



      base_matrix G1, G2, GPr(N,N);
      base_vector coeff, Velem, wt1(N), wt2(N);
      base_matrix Melem, grad_d2(1, N), grad_d1(1, N);
      base_small_vector d2(1), n1(N), n2(N), Pr(N), zeta(N), u1(N), u2(N);
      base_tensor tG1, tGdu1, tGddu1, tbv1, tbv2;
      scalar_type gap, u1n, u2n, tbv1n(0), tbv2n(0);
      size_type cv2(-1);

      bgeot::multi_index sizes_tGdu1(1), sizes_tGddu1(3);
      sizes_tGdu1[0] = N;
      tG1.adjust_sizes(sizes_tGdu1);
      sizes_tGdu1.push_back(N);
      sizes_tGddu1[2] = N;
      
      cout << "begining gauss points loop" << endl;
      
      for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv) {

        cout << "element " << cv << endl;
     
        pintegration_method pim = mim.int_method_of_element(cv);
        if (pim->type() != IM_APPROX) continue; 


        cout << "pim = " << int(pim->type()) << endl;
        cout << "pim = " << pim->approx_method() << endl;

        bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv));
        
        bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
        pfem pf_u1 = mf_u1.fem_of_element(cv);
        pfem pf_d1 = mf_d1.fem_of_element(cv);
        size_type nbdof1 = mf_u1.nb_basic_dof_of_element(cv);
        sizes_tGddu1[0] = sizes_tGddu1[1] = sizes_tGdu1[0] = nbdof1;
        tGdu1.adjust_sizes(sizes_tGdu1);
        tGddu1.adjust_sizes(sizes_tGddu1);
        
          
        size_type nbpt = pim->approx_method()->nb_points();
        for (size_type ipt = 0; ipt < nbpt; ++ipt) {
          
          const base_node xref = pim->approx_method()->integration_points()[ipt];
          cout << "xref = " << xref << endl;
          
          fem_interpolation_context ctx_u1(pgt, pf_u1, xref, G1, cv);
          base_node x0 = ctx_u1.xreal();


	  
          scalar_type weight = pim->approx_method()->coeff(ipt) * ctx_u1.J();
       
	  cout << "pim->approx_method()=" << pim->approx_method()->coeff(ipt) * ctx_u1.J()<< endl;
	  
          // computation of h for gamma = gamma0*h
          scalar_type gamma(0);
          if (ipt == 0) {
            scalar_type emax, emin;
            gmm::condition_number(ctx_u1.K(),emax,emin);
            gamma = gamma0 * emax * sqrt(scalar_type(N));
          }
          
	  
          // computation of u1, w1, f_friction
          slice_vector_on_basic_dof_of_element(mf_u1, U1, cv, coeff);
          ctx_u1.pf()->interpolation(ctx_u1, coeff, u1, bgeot::dim_type(N));
          if (WT1) {
            slice_vector_on_basic_dof_of_element(mf_u1, *WT1, cv, coeff);
            ctx_u1.pf()->interpolation(ctx_u1, coeff, wt1, bgeot::dim_type(N));
          }
		  

          // Computation of n1
          fem_interpolation_context ctx_d1(pgt, pf_d1, xref, G1, cv);
          slice_vector_on_basic_dof_of_element(mf_d1, D1, cv, coeff);
          ctx_d1.pf()->interpolation_grad(ctx_d1, coeff, grad_d1, 1);
          gmm::copy(grad_d1.as_vector(), n1);
          gmm::scale(n1, 1./gmm::vect_norm2(n1));
	  

          cout << " weight ="  << weight << endl;
          cout << " Element " << cv << " point " << ipt << " elt ref : " <<
          pim->approx_method()->integration_points()[ipt] << " elt reel : " << x0 << endl;

 	  // REMARQUE : on peut faire  ctx_u1.xref() plutôt que  pim->approx_method()->integration_points()[ipt]      

	    
           //Definition de la projection

          pfem pf_d2 = mf_d2.fem_of_element(cv);

          fem_interpolation_context ctx_d2(pgt, pf_d2, xref, G1, cv);

          slice_vector_on_basic_dof_of_element(mf_d2, D2, cv, coeff);

          ctx_d2.pf()->interpolation(ctx_d2, coeff, d2, 1);

          ctx_d2.pf()->interpolation_grad(ctx_d2, coeff, grad_d2, 1);

          cout << "grad_d2 = " << grad_d2 << endl;
          cout << "d2 = " << d2[0] << endl;
 
          base_node y0 = x0, yref(N);
          gmm::add(gmm::scaled(gmm::mat_row(grad_d2, 0),
                   -d2[0] / gmm::vect_norm2_sqr(gmm::mat_row(grad_d2, 0))),
                   y0);
	  
	  
	  
          bgeot::rtree::pbox_set pbs;
	  
          tree.find_boxes_at_point(y0, pbs);
	   
          bgeot::rtree::pbox_set::const_iterator it = pbs.begin();
	  
	  
          bool found = false;
          size_type nbdof2(0);
          
	  for (; it != pbs.end(); ++it) {
	    cv2 = (*it)->id;
            bgeot::pgeometric_trans pgty =  m.trans_of_convex(cv2);
            nbdof2 = mf_u2.nb_basic_dof_of_element(cv2);
            bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cv2));

            bgeot::geotrans_inv_convex gic;
            gic.init(m.points_of_convex(cv2), pgty);

            gic.invert(y0, yref);
            if (pgty->convex_ref()->is_in(yref) < 1E-10)
              { found = true; break; }
          }

          GMM_ASSERT1(found && (cv2 != size_type(-1)),
                      "Projection not found ...");

          cout << "y0 = " << y0 << endl;
          cout << "Found element : " << cv2 << " yref = " << yref << endl;
     
	  
	  
	  
         
	  
	  
          pfem pf_u2 = mf_u2.fem_of_element(cv2);
          fem_interpolation_context ctx_u2(pgt, pf_u2, yref, G2, cv2);

          // computation of u2
          slice_vector_on_basic_dof_of_element(mf_u2, U2, cv2, coeff);
          ctx_u2.pf()->interpolation(ctx_u2, coeff, u2, bgeot::dim_type(N));
          if (WT2) {
            slice_vector_on_basic_dof_of_element(mf_u2, *WT2, cv2, coeff);
            ctx_u2.pf()->interpolation(ctx_u2, coeff, wt2, bgeot::dim_type(N));
          }

	  cout << "1" << endl;

          u1n = gmm::vect_sp(u1, n2); u2n = gmm::vect_sp(u2, n2);

	  cout << "2" << endl;

          
          md.compute_Neumann_terms(1, vl[0], mf_u1, U1, ctx_u1, n1, tG1);
	  cout << " tG1 = " << tG1 << endl; 
          md.compute_Neumann_terms(2, vl[0], mf_u1, U1, ctx_u1, n1, tGdu1);
          md.compute_Neumann_terms(3, vl[0], mf_u1, U1, ctx_u1, n1, tGddu1);
	  cout << "2.0" << endl;
	  cout << " tGdu1 = " << tGdu1 << endl;
	  cout << " tGddu1 = " << tGddu1 << endl;
          ctx_u1.pf()->real_base_value(ctx_u1, tbv1);
  	  cout << "2.1" << endl;	

          ctx_u2.pf()->real_base_value(ctx_u2, tbv2);
  	  cout << "2.2" << endl;	

          for(size_type i=0; i<N; ++i)
            zeta[i] = tG1[i] +
              ( (gap + (alpha-1.)*u1n-(alpha-1.)*u2n)*n2[i]
                + alpha*(wt1[i]-wt2[i]) - alpha*u1[i] + alpha*u2[i]) / gamma;

	  cout << "3" << endl;

          coupled_projection(zeta, n2, f_coeff, Pr);
	  coupled_projection_grad(zeta, n2, f_coeff, GPr);

	  for (size_type k = 0; k < nbdof1; ++k)
	    for (size_type l = 0; l < N; ++l){
	       tbv1n += n2[l]*tbv1(k,l);
	       tbv2n += n2[l]*tbv2(k,l);
	    }


          gap = scalar_type(0);
          for( size_type i=0; i<N; ++i)
            gap += (y0[i]-x0[i])*n2[i];

          cout << "gap = " << gap << endl;

	  cout << "4" << endl;


            // Plan tangent 
          if (version & model::BUILD_MATRIX){
            
	    // Matrice en u1,u1
            gmm::resize(Melem, nbdof1, nbdof1); gmm::clear(Melem);
            for (size_type j = 0; j < nbdof1; ++j)
              for (size_type k = 0; k < nbdof1; ++k){
		scalar_type res(0);
                for (size_type i = 0; i < N; ++i) {
                  if (theta != scalar_type(0)) { 
                    res -= theta*gamma*tGdu1(j, i) * tGdu1(k, i);
                    res += theta*gamma*(Pr[i]-tG1[i])*(tGddu1(j,k,i));
                    for (size_type l =0; l<N; ++l){
                      res += theta*GPr(i,l)*(gamma*tGdu1(k,l)
                                                     -alpha*tbv1(k,l)-(scalar_type(1)-alpha)*n2[l]*tbv1n)*tGdu1(j,i);
		    }
		    }
		  for (size_type l =0; l<N;++l){
                    res -= GPr(i,l)*(tGdu1(k,l)+(-alpha*tbv1(k,l)-(scalar_type(1)-alpha)*n2[l]*tbv1n)/(gamma))*tbv1(j,i);
		  }				   
	      }
	      Melem(j, k)=res;
	      }
	      gmm::scale(Melem,weight);
            mat_elem_assembly(matl[0], Melem, mf_u1, cv, mf_u1, cv);

            // Matrice en fontion de u1,u2
            gmm::resize(Melem, nbdof1, nbdof2); gmm::clear(Melem);
            for (size_type j = 0; j < nbdof1; ++j)
              for (size_type k = 0; k < nbdof2; ++k){
		scalar_type res(0);
                for (size_type i = 0; i < N; ++i) 
		  for (size_type l =0; l<N;++l)
                    res += GPr(i,l)*(tGdu1(k,l)+(-alpha*tbv1(k,l)
                                                         -(scalar_type(1)-alpha)*n2[l]*tbv1n)/(gamma))*tbv2(j,i);
		Melem(j, k)=res;
	      }
	     gmm::scale(Melem,weight);
            mat_elem_assembly(matl[1], Melem, mf_u1, cv, mf_u2, cv2);

            // Matrice en u2,u1
            gmm::resize(Melem, nbdof2, nbdof1); gmm::clear(Melem);
            for (size_type j = 0; j < nbdof2; ++j)
              for (size_type k = 0; k < nbdof1; ++k){
		scalar_type res(0);
                for (size_type i = 0; i < N; ++i) {
                  if (theta != scalar_type(0)) {
                    for (size_type l =0; l<N;++l)
                      res += theta*GPr(i,l)*(alpha*tbv2(k,l)+(scalar_type(1)-alpha)*n2[l]*tbv2n)*tGdu1(j,i);
						}
		  for (size_type l =0; l<N;++l)
                    res -= GPr(i,l)*(alpha*tbv2(k,l)+(scalar_type(1)-alpha)*n2[l]*tbv2n)*tbv1(j,i)/(gamma);				
						   } 
		Melem(j, k)=res;				   
	      }
		gmm::scale(Melem,weight);				   
            mat_elem_assembly(matl[2], Melem, mf_u2, cv2, mf_u1, cv);	    
	    
	    // Matrice en u2,u2
            gmm::resize(Melem, nbdof2, nbdof2); gmm::clear(Melem);
            for (size_type j = 0; j < nbdof2; ++j)
              for (size_type k = 0; k < nbdof2; ++k){
		scalar_type res(0);
                for (size_type i = 0; i < N; ++i) {
                  for (size_type l =0; l<N;++l)
                    res += GPr(i,l)*(alpha*tbv2(k,l)+(scalar_type(1)-alpha)*n2[l]*tbv2n)*tbv2(j,i)/(gamma);			
						   } 
		Melem(j, k)=res;				   
	      }
		gmm::scale(Melem,weight);				   
            mat_elem_assembly(matl[3], Melem, mf_u2, cv2, mf_u2, cv2);	
				             }
				             
		// Matrice du second Membre		             
          if (version & model::BUILD_RHS){
            
            gmm::resize(Velem, nbdof1); gmm::clear(Velem);
            for (size_type j = 0; j < nbdof1; ++j){
	      scalar_type res(0);
              for (size_type i = 0; i < N; ++i){
                if (theta != scalar_type(0)){
                  res += theta*gamma*tG1[i] * tGdu1(j, i); // A supprimer si theta==0
                  res -=theta*gamma*Pr[i]*tGdu1(i,j); // A supprimer si theta==0
                }
                res -=Pr[i]*tbv1(j,i);
              }
              Velem[j]=res;
	    }
              gmm::scale(Velem,weight);
            vec_elem_assembly(vecl[0], Velem, mf_u1, cv);

            gmm::resize(Velem, nbdof2);gmm::clear(Velem);
            for (size_type j = 0; j < nbdof2; ++j){
	      scalar_type res(0);
              for (size_type i = 0; i < N; ++i)
                res -= Pr[i]*tbv2(j,i);
	      Velem[j]=res;
	    }
			        
	    gmm::scale(Velem,weight);
            vec_elem_assembly(vecl[1], Velem, mf_u1, cv);
          }


//           size_type nit = 0;
//           while (gmm::abs(d0) > 1E-10 && ++nit < 1000) {
//             for (size_type k = 0; k < N; ++k) {
//               pt_eval[k] += EPS;
//               d1 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
//               n[k] = (d1 - d0) / EPS;
//               pt_eval[k] -= EPS;
//             }

//             gmm::add(gmm::scaled(n, -d0 / gmm::vect_norm2_sqr(n)), pt_eval);
//             // A simple line search could be added
//             d0 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
//           }
//           GMM_ASSERT1(nit < 1000, "Projection on rigid obstacle did not converge");

//           ct.master_point.resize(N);
//           gmm::copy(pt_eval, ct.master_point);






        }

      }


      cout << "end assembly" << endl;


    }
    
    Nitsche_fictitious_domain_contact_brick(scalar_type theta_,
                                            bool nofriction) {
      theta = theta_;
      contact_only = nofriction;
      set_flags("Integral Nitsche contact and friction with rigid "
                "obstacle brick",
                false /* is linear*/, false /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */, false /* compute each time */,
                false /* has a Neumann term */);
    }

  };



  size_type add_Nitsche_fictitious_domain_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &dataname_d1,
   const std::string &dataname_d2, const std::string &dataname_gamma0,
   scalar_type theta,
   const std::string &dataname_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt1, const std::string &dataname_wt2) {

    bool nofriction = (dataname_friction_coeff.size() == 0);
    pbrick pbr = new Nitsche_fictitious_domain_contact_brick(theta,nofriction);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u1, varname_u1, false));
    tl.push_back(model::term_description(varname_u1, varname_u2, false));
    tl.push_back(model::term_description(varname_u2, varname_u1, false));
    tl.push_back(model::term_description(varname_u2, varname_u2, false));


    model::varnamelist dl(1, dataname_d1);
    dl.push_back(dataname_d2);
    dl.push_back(dataname_gamma0);
    if (!nofriction) dl.push_back(dataname_friction_coeff);
    if (dataname_alpha.size() > 0) {
      dl.push_back(dataname_alpha);
      if (dataname_wt1.size() > 0)
        { dl.push_back(dataname_wt1); dl.push_back(dataname_wt2); }
    }

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname_u1, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname_u1, aux_vars[i], false));
    }
//     aux_vars.resize(0);
//     md.auxilliary_variables_of_Neumann_terms(varname_u2, aux_vars);
//     for (size_type i = 0; i < aux_vars.size(); ++i) {
//       vl.push_back(aux_vars[i]);
//       tl.push_back(model::term_description(varname_u2, aux_vars[i], false));
//     }

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim),
                        size_type(-1));
  }









}  /* end of namespace getfem.                                             */
