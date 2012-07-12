/*===========================================================================
 
 Copyright (C) 2011-2012 Yves Renard, Konstantinos Poulios.
 
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
    case RHS_U_V1:       case RHS_U_V2:       case RHS_U_V3: case RHS_U_V4:
    case RHS_U_V5:       case RHS_U_V6:       case RHS_U_V7: case RHS_U_V8:
    case RHS_U_FRICT_V1: case RHS_U_FRICT_V2:
    case RHS_U_FRICT_V3: case RHS_U_FRICT_V4: case RHS_U_FRICT_V5:
    case RHS_L_FRICT_V1: case RHS_L_FRICT_V2:
    case RHS_L_FRICT_V3: case RHS_L_FRICT_V4:
    case K_UL_V1:        case K_UL_V2:        case K_UL_V3:  case K_UL_V4:
    case UZAWA_PROJ_FRICT: case UZAWA_PROJ_FRICT_SAXCE:
      sizes_[0] = N; break;
      // two-dimensional tensors [N x N]
    case K_UU_V1: case K_UU_V2:
    case K_UL_FRICT_V1: case K_UL_FRICT_V2: case K_UL_FRICT_V3:
    case K_UL_FRICT_V4: case K_UL_FRICT_V5: case K_UL_FRICT_V6:
    case K_UL_FRICT_V7: case K_UL_FRICT_V8:
    case K_LL_FRICT_V1: case K_LL_FRICT_V2:
    case K_LL_FRICT_V3: case K_LL_FRICT_V4:
    case K_UU_FRICT_V1: case K_UU_FRICT_V2:
    case K_UU_FRICT_V3: case K_UU_FRICT_V4: case K_UU_FRICT_V5:
      sizes_.resize(2); sizes_[0] = sizes_[1] = N;  break;
    }

    // adjust temporary variables sizes
    lnt.resize(N); lt.resize(N); zt.resize(N); no.resize(N);
    aux1.resize(1); auxN.resize(N); V.resize(N);
    gmm::resize(GP, N, N);
  }

  void contact_nonlinear_term::compute
  (fem_interpolation_context &/* ctx */, bgeot::base_tensor &t) {

    t.adjust_sizes(sizes_);
    scalar_type e, augm_ln;
    dim_type i, j;

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
      // here ln is expected to be an estimation of the mesh size
      // and r should be a threshold coefficient expressing a penetration
      // or separation distance as percentage of the mesh size 
      t[0] = Heav(un-g - r*ln);  break;

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
    case RHS_U_V8:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
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
      // if (e < 0. && ctx.xreal()[1] > 1.)
      //        cout << "x = " << ctx.xreal() << " e = " << e << endl;
      auxN = lt;  ball_projection(auxN, f_coeff * gmm::neg(ln));
      // if (gmm::vect_norm2(auxN) > 0. && ctx.xreal()[1] > 1.)
      //        cout << "x = " << ctx.xreal() << " auxN = " << auxN << endl;
      for (i=0; i<N; ++i) t[i] = no[i]*e + auxN[i];
      break;
    case RHS_U_FRICT_V5:
      auxN = lnt; De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
      break;
    case RHS_L_FRICT_V1:
      e = ln+gmm::neg(ln-r*(un-g));
      auxN = zt - lt;  ball_projection(auxN, -f_coeff * ln); auxN += lt;
      for (i=0; i<N; ++i) t[i] = (e*no[i] + auxN[i])/ r;
      break;
    case RHS_L_FRICT_V2:
      e = r*(un-g) + gmm::pos(ln);
      auxN = lt;  ball_projection(auxN, f_coeff * gmm::neg(ln));
      for (i=0; i<N; ++i) t[i] = (no[i]*e + zt[i] + lt[i] - auxN[i])/r;
      break;
    case RHS_L_FRICT_V3:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = (lnt[i] - auxN[i])/r;
      break;
    case RHS_L_FRICT_V4:
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
    case K_UL_V4:
      for (i=0; i<N; ++i) t[i] = -no[i];
      break;
    case UZAWA_PROJ_FRICT:
      e = -gmm::neg(ln - r*(un - g));
      auxN = lt - zt;  ball_projection(auxN, -f_coeff * e);
      for (i=0; i<N; ++i) t[i] = e*no[i] + auxN[i];
      break;
    case UZAWA_PROJ_FRICT_SAXCE:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection(auxN, no, f_coeff);
      for (i=0; i<N; ++i) t[i] = auxN[i];
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
      e = (alpha-scalar_type(1));
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = no[i]*no[j]*e - ((i == j) ? alpha : scalar_type(0));
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
        for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = -GP(i,j);
      }
      break;
    case K_UL_FRICT_V7:
      De_Saxce_projection_grad(lnt, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = -GP(j,i);
      break;
    case K_UL_FRICT_V8:
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
        t[i*N+j] = (no[i]*no[j]*e - ((i == j) ? scalar_type(1) : scalar_type(0)) + GP(i,j) - f_coeff*Heav(-ln)*no[i]*V[j])/r;
      break;
    case K_LL_FRICT_V3:
      auxN = lnt - (r*(un-g) - f_coeff * gmm::vect_norm2(zt)) * no - zt;
      De_Saxce_projection_grad(auxN, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = (GP(i,j) - scalar_type((i == j) ? 1 : 0))/r;
      break;
    case K_LL_FRICT_V4:
      De_Saxce_projection_grad(lnt, no, f_coeff, GP);
      for (i=0; i<N; ++i) for (j=0; j<N; ++j)
        t[i*N+j] = (GP(i,j) - ((i == j) ? scalar_type(1) : scalar_type(0)))/r;
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
    case K_UU_FRICT_V5:
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
          if (gmm::vect_size(VT) == gmm::vect_size(U)) {
            gmm::copy(gmm::sub_vector
                      (VT, gmm::sub_index
                       (mf_u.ind_basic_dof_of_element(cv))), coeff);
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
      coeff.resize(mf_obs.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (obs, gmm::sub_index
                 (mf_obs.ind_basic_dof_of_element(cv))), coeff);
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
        coeff.resize(pmf_lambda->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (lambda, gmm::sub_index
                   (pmf_lambda->ind_basic_dof_of_element(cv))), coeff);
        if (contact_only) {
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          ln = aux1[0];
        } else {
          ctx.pf()->interpolation(ctx, coeff, lnt, N);
        }
      }
      break;

    case 4 :// calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
      GMM_ASSERT1(!contact_only, "Invalid friction option");
      if (pmf_coeff) {
        coeff.resize(pmf_coeff->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (friction_coeff, gmm::sub_index
                   (pmf_coeff->ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
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
      coeff.resize(mf_u1.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (U1, gmm::sub_index
                 (mf_u1.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no) - un;
      if (!contact_only) {
        if (gmm::vect_size(WT1) == gmm::vect_size(U1)) {
          gmm::copy(gmm::sub_vector
                    (WT1, gmm::sub_index
                     (mf_u1.ind_basic_dof_of_element(cv))), coeff);
          ctx.pf()->interpolation(ctx, coeff, auxN, N);
          auxN -= gmm::vect_sp(auxN, no) * no;
          zt = ((V - un * no) - auxN) * (r * alpha) - zt; // zt = r*alpha*(u_T-w_T)
        } else {
          zt = (V - un * no) * (r * alpha) - zt;          // zt = r*alpha*u_T
        }
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

      coeff.resize(mf_u2.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (U2, gmm::sub_index
                 (mf_u2.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no);
      if (!contact_only) {
        if (gmm::vect_size(WT2) == gmm::vect_size(U2)) {
          gmm::copy(gmm::sub_vector
                    (WT2, gmm::sub_index
                     (mf_u2.ind_basic_dof_of_element(cv))), coeff);
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
        coeff.resize(pmf_lambda->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (lambda, gmm::sub_index
                   (pmf_lambda->ind_basic_dof_of_element(cv))), coeff);
        if (contact_only) {
          ctx.pf()->interpolation(ctx, coeff, aux1, 1);
          ln = aux1[0];
        } else {
          ctx.pf()->interpolation(ctx, coeff, lnt, N);
        }
      }
      break;

    case 4 :// calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
      GMM_ASSERT1(!contact_only, "Invalid friction option");
      if (pmf_coeff) {
        coeff.resize(pmf_coeff->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (friction_coeff, gmm::sub_index
                   (pmf_coeff->ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
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

    size_type subterm1 = (option == 4) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 4) ? K_UL_V4 : K_UL_V3;
    size_type subterm3 = (option == 4) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm3(subterm3, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda),
      nterm4(subterm4, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4:
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
    case 3:
      assem.set
      ("M$1(#1,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); "      // UL
       "M$2(#3,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#1))(i,:,:,i); "      // LU
       "M$3(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:); "         // LL
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   scalar_type gamma, const VECT1 *VT,
   const mesh_region &rg, int option = 1) {

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
    case 6 : subterm1 = K_UL_FRICT_V7; subterm2 = K_UL_FRICT_V8;
      subterm3 = K_LL_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    size_type subterm4 = (option == 2) ? K_UU_FRICT_V2 : K_UU_FRICT_V1;

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT),
      nterm3(subterm3, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT),
      nterm4(subterm4, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4: case 5: case 6:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // UL
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU
        "M$3(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j)"); // LL
      break;
    case 2: case 3:
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
    case 3 : subterm1 = RHS_U_V3; break;
    case 4 : subterm1 = RHS_U_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }
    size_type subterm2 = (option == 4) ? RHS_L_V2 : RHS_L_V1;

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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   scalar_type gamma, const VECT1 *VT,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2;
    switch (option) {
    case 1 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V1; break;
    case 2 : subterm1 = RHS_U_FRICT_V2; subterm2 = RHS_L_FRICT_V1; break;
    case 3 : subterm1 = RHS_U_FRICT_V3; subterm2 = RHS_L_FRICT_V1; break;
    case 4 : subterm1 = RHS_U_FRICT_V4; subterm2 = RHS_L_FRICT_V2; break;
    case 5 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V3; break;
    case 6 : subterm1 = RHS_U_FRICT_V5; subterm2 = RHS_L_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    contact_rigid_obstacle_nonlinear_term
      nterm1(subterm1, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT, gamma, VT);

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

    bool Tresca_version, contact_only;
    int option;

    // option = 1 : Alart-Curnier
    // option = 2 : symmetric Alart-Curnier (almost with friction)
    // option = 3 : Alart Curnier "over-augmented"
    // option = 4 : New method
    // option = 5 : De-Saxce
    // option = 6 : New method based on De-Saxce.

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
      GMM_ASSERT1(matl.size() == size_type(3 + (option == 3)
                                           + (option == 2 && !contact_only)),
                  "Wrong number of terms for "
                  "integral contact with rigid obstacle brick");

      // variables : u, lambda. The variable lambda should be scalar in the
      //             frictionless case and vector valued in the case with
      //             friction.
      // data      : obstacle, r for the version without friction
      //           : obstacle, r, friction_coeff, alpha, w_t, gamma, v_t for
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

      const model_real_plain_vector &friction_coeff
        = contact_only ? u : md.real_variable(dl[2]);
      const mesh_fem *pmf_coeff = contact_only ? 0 : md.pmesh_fem_of_variable(dl[2]);
      sl = gmm::vect_size(friction_coeff);
      if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
      GMM_ASSERT1(sl == 1 || contact_only,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

      scalar_type alpha = 1;
      if (!contact_only && dl.size() >= 4) {
        alpha = md.real_variable(dl[3])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[3])) == 1,
                    "Parameter alpha should be a scalar");
      }

      const model_real_plain_vector *WT
        = (!contact_only && dl.size()>=5) ? &(md.real_variable(dl[4])) : 0;

      scalar_type gamma = 1;
      if (!contact_only && dl.size() >= 6) {
        gamma = md.real_variable(dl[5])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[5])) == 1,
                    "Parameter gamma should be a scalar");
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
             pmf_coeff, &friction_coeff, vr[0], alpha, WT, gamma, VT,
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
             pmf_coeff, &friction_coeff, vr[0], alpha, WT, gamma, VT,
             rg, option);
      }

    }

    integral_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
      Tresca_version = false;   // for future version ...
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
    case 1 : case 4 :
      tl.push_back(model::term_description(varname_u, multname_n, false)); // UL
      tl.push_back(model::term_description(multname_n, varname_u, false)); // LU
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u, multname_n, true));  // UL
      tl.push_back(model::term_description(varname_u, varname_u, true));   // UU (fourthmat == 1)
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      break;
    case 3 :
      tl.push_back(model::term_description(varname_u, multname_n, false)); // UL
      tl.push_back(model::term_description(multname_n, varname_u, false)); // LU
      tl.push_back(model::term_description(multname_n, multname_n, true)); // LL
      tl.push_back(model::term_description(varname_u, varname_u, true));   // UU
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
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region, int option,
   const std::string &dataname_alpha, const std::string &dataname_wt,
   const std::string &dataname_gamma, const std::string &dataname_vt) {

    pbrick pbr
      = new integral_contact_rigid_obstacle_brick(false, option);

    model::termlist tl;

    switch (option) {
    case 1: case 4: case 5: case 6:
      tl.push_back(model::term_description(varname_u, multname, false)); // UL
      tl.push_back(model::term_description(multname, varname_u, false)); // LU
      tl.push_back(model::term_description(multname, multname, true));   // LL
      break;
    case 2: case 3:
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
    dl.push_back(dataname_friction_coeff);
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
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
            pmf_coeff, f_coeff, alpha, WT);

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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
   scalar_type alpha, const VECT1 *WT,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_V7; break;
    case 2 : subterm =  RHS_U_V6; break;
    case 3 : subterm =  RHS_U_V8; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(subterm, r, mf_u, U, mf_obs, obs, pmf_lambda, lambda,
            pmf_coeff, f_coeff, alpha, WT);

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

      // Data : obs, r, [lambda,] [friction_coeff,] [alpha,] [WT]
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

      const model_real_plain_vector *f_coeff = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      const model_real_plain_vector *WT = 0;
      if (!contact_only) {
        nd++;
        f_coeff = &(md.real_variable(dl[nd]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*f_coeff);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1,
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
          WT = &(md.real_variable(dl[nd]));
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
             pmf_coeff, f_coeff, vr[0], alpha, WT, rg, option);
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
             pmf_coeff, f_coeff, vr[0], alpha, WT, rg, option);
      }

    }

    penalized_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
      Tresca_version = false;   // for future version ...
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
   const std::string &dataname_friction_coeff,
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

    size_type subterm1 = (option == 4) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 4) ? K_UL_V4 : K_UL_V3;
    size_type subterm3 = (option == 4) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm3(subterm3, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda),
      nterm4(subterm4, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4:
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
    case 3:
      assem.set
      ("M$1(#1,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); "      // U1L
       "M$2(#3,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#1))(i,:,:,i); "      // LU1
       "M$3(#2,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#2).Base(#3))(i,:,i,:); "      // U2L
       "M$4(#3,#2)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#2))(i,:,:,i); "      // LU2
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
   MAT &Ku1u1, MAT &Ku2u2, MAT &Ku1u2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VEC &U1,
   const getfem::mesh_fem &mf_u2, const VEC &U2,
   const getfem::mesh_fem &mf_lambda, const VEC &lambda,
   const getfem::mesh_fem *pmf_coeff, const VEC *f_coeff,
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
      subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V4; subterm3 = K_LL_FRICT_V1;
      break;
    case 4 :
      subterm1 = K_UL_FRICT_V2; subterm2 = K_UL_FRICT_V5; subterm3 = K_LL_FRICT_V2;
      break;
    case 5 :
      subterm1 = K_UL_FRICT_V1; subterm2 = K_UL_FRICT_V6; subterm3 = K_LL_FRICT_V3;
      break;
    case 6 :
      subterm1 = K_UL_FRICT_V7; subterm2 = K_UL_FRICT_V8; subterm3 = K_LL_FRICT_V4;
      break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    size_type subterm4 = (option == 2) ? K_UU_FRICT_V2 : K_UU_FRICT_V1;

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2),
      nterm3(subterm3, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2),
      nterm4(subterm4, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3,#4" : "#1,#2,#3";

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4: case 5: case 6:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // U1L
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU1
        "M$3(#2,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2).vBase(#3))(i,j,:,i,:,j); " // U2L
        "M$4(#3,#2)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#2))(i,j,:,j,:,i); " // LU2
        "M$5(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j)"); // LL
      break;
    case 2: case 3:
      assem.set
       ("M$1(#1,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#1).vBase(#3))(i,j,:,i,:,j); " // U1L
        "M$2(#3,#1)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#1))(i,j,:,j,:,i); " // LU1
        "M$3(#2,#3)+=comp(NonLin$1(#1," + aux_fems + ").vBase(#2).vBase(#3))(i,j,:,i,:,j); " // U2L
        "M$4(#3,#2)+=comp(NonLin$2(#1," + aux_fems + ").vBase(#3).vBase(#2))(i,j,:,j,:,i); " // LU2
        "M$5(#3,#3)+=comp(NonLin$3(#1," + aux_fems + ").vBase(#3).vBase(#3))(i,j,:,i,:,j); " // LL
        "M$6(#1,#1)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#1).vBase(#1))(i,j,:,i,:,j); " // U1U1
        "M$7(#2,#2)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#2).vBase(#2))(i,j,:,i,:,j); " // U2U2
        "M$8(#1,#2)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#1).vBase(#2))(i,j,:,i,:,j)"); // U1U2
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
    case 3 : subterm1 = RHS_U_V3; break;
    case 4 : subterm1 = RHS_U_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }
    size_type subterm2 = (option == 4) ? RHS_L_V2 : RHS_L_V1;

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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff,
   scalar_type r, scalar_type alpha,
   const VECT1 *WT1, const VECT1 *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm1, subterm2;
    switch (option) {
    case 1 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V1; break;
    case 2 : subterm1 = RHS_U_FRICT_V2; subterm2 = RHS_L_FRICT_V1; break;
    case 3 : subterm1 = RHS_U_FRICT_V3; subterm2 = RHS_L_FRICT_V1; break;
    case 4 : subterm1 = RHS_U_FRICT_V4; subterm2 = RHS_L_FRICT_V2; break;
    case 5 : subterm1 = RHS_U_FRICT_V1; subterm2 = RHS_L_FRICT_V3; break;
    case 6 : subterm1 = RHS_U_FRICT_V5; subterm2 = RHS_L_FRICT_V4; break;
    default : GMM_ASSERT1(false, "Incorrect option");
    }

    contact_nonmatching_meshes_nonlinear_term
      nterm1(subterm1, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, f_coeff, alpha, WT1, WT2);

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
    mutable getfem::pfem pfem_proj;        // cached fem and mesh_fem for the
    mutable getfem::mesh_fem *pmf_u2_proj; // projection between nonmatching meshes
    bool Tresca_version, contact_only;
    int option;

    // option = 1 : Alart-Curnier
    // option = 2 : symmetric Alart-Curnier (almost with friction)
    // option = 3 : Alart Curnier "over-augmented"
    // option = 4 : New method
    // option = 5 : De-Saxce
    // option = 6 : New method based on De-Saxce.

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

      // Data : r, [friction_coeff,] [alpha,] [WT1, WT2]
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

      const model_real_plain_vector *f_coeff = 0, *WT1 = 0, *WT2 = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      if (!contact_only) {
        f_coeff = &(md.real_variable(dl[1]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[1]);

        size_type sl = gmm::vect_size(*f_coeff);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1,
                    "the data corresponding to the friction coefficient "
                    "has not the right format");

        if (dl.size() >= 3) {
          alpha = md.real_variable(dl[2])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[2])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() >= 4)
          WT1 = &(md.real_variable(dl[3]));

        if (dl.size() >= 5)
          WT2 = &(md.real_variable(dl[4]));
      }

      // Matrix terms (T_u1l, T_lu1, T_u2l, T_lu2, T_ll, T_u1u1, T_u2u2, T_u1u2)
      GMM_ASSERT1(matl.size() == size_type(3 +                                // U1L, U2L, LL
                                           2 * !is_symmetric() +              // LU1, LU2
                                           3 * (option == 3 || option == 2)), // U1U1, U2U2, U1U2
                  "Wrong number of terms for "
                  "integral contact between nonmatching meshes brick");

      mesh_region rg(region);
      mf_u1.linked_mesh().intersect_with_mpi_region(rg);

      size_type N = mf_u1.linked_mesh().dim();

      // projection of the second mesh_fem onto the mesh of the first mesh_fem
      if (!pmf_u2_proj) {
        pmf_u2_proj = new getfem::mesh_fem(mim.linked_mesh(), dim_type(N));
        pfem_proj = new_projected_fem(mf_u2, mim, rg2, rg1);
        pmf_u2_proj->set_finite_element(mim.linked_mesh().convex_index(), pfem_proj);
      }

      size_type nbdof1 = mf_u1.nb_dof();
      size_type nbdof_lambda = mf_lambda.nb_dof();
      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = pmf_u2_proj->nb_basic_dof();

      std::vector<size_type> ind;
      pmf_u2_proj->get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      gmm::csc_matrix<scalar_type> Rsub(nbdof2, nbsub), Esub(nbsub, nbdof2);
      if (mf_u2.is_reduced()) {
          gmm::copy(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    Rsub);
          gmm::copy(gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    Esub);
      }

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(Esub, u2, u2_proj);
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      size_type U1L = 0;
      size_type LU1 = U1L + (is_symmetric() ? 0 : 1);
      size_type U2L = LU1 + 1;
      size_type LU2 = U2L + (is_symmetric() ? 0 : 1);
      size_type LL  = LU2 + 1;
      size_type U1U1 = (option == 1 || option == 4) ? U1L : LL + 1;
      size_type U2U2 = (option == 1 || option == 4) ? U2L : LL + 2;
      size_type U1U2 = (option == 1 || option == 4) ? U1L : LL + 3;

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Integral contact between nonmatching meshes "
                   "tangent term");
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(matl[i]);

        model_real_sparse_matrix Ku2l(nbsub, nbdof_lambda);
        model_real_sparse_matrix Klu2(nbdof_lambda, nbsub);
        model_real_sparse_matrix Ku2u2(nbsub, nbsub);
        model_real_sparse_matrix Ku1u2(nbdof1, nbsub);

        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[U1L], matl[LU1], Ku2l, Klu2, matl[LL], matl[U1U1], Ku2u2, Ku1u2,
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[U1L], matl[LU1], Ku2l, Klu2, matl[LL], matl[U1U1], Ku2u2, Ku1u2,
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, f_coeff, vr[0], alpha, WT1, WT2, rg, option);

        if (mf_u2.is_reduced()) {
          gmm::mult(Rsub, Ku2l, matl[U2L]);
          if (LU2 != U2L) gmm::mult(Klu2, Esub, matl[LU2]);
          if (U2U2 != U2L) {
            model_real_sparse_matrix tmp(nbsub, nbdof2);
            gmm::mult(Ku2u2, Esub, tmp);
            gmm::mult(Rsub, tmp, matl[U2U2]);
            gmm::mult(Ku1u2, Esub, matl[U1U2]);
          }
        }
        else {
          gmm::copy(Ku2l, gmm::sub_matrix(matl[U2L], SUBI, gmm::sub_interval(0, nbdof_lambda)));
          if (LU2 != U2L)
            gmm::copy(Klu2, gmm::sub_matrix(matl[LU2], gmm::sub_interval(0, nbdof_lambda), SUBI));
          if (U2U2 != U2L) {
            gmm::copy(Ku2u2, gmm::sub_matrix(matl[U2U2], SUBI));
            gmm::copy(Ku1u2, gmm::sub_matrix(matl[U1U2], gmm::sub_interval(0, nbdof1), SUBI));
          }
        }
      }

      if (version & model::BUILD_RHS) {
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(vecl[i]);

        model_real_plain_vector Ru2(nbsub);

        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[U1L], Ru2, vecl[LL], // u1, u2, lambda
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[U1L], Ru2, vecl[LL], // u1, u2, lambda
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, f_coeff, vr[0], alpha, WT1, WT2, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(Rsub, Ru2, vecl[U2L]);
        else
          gmm::copy(Ru2, gmm::sub_vector(vecl[U2L], SUBI));
      }

    }

    integral_contact_nonmatching_meshes_brick(size_type rg1_, size_type rg2_,
                                                bool contact_only_, int option_)
    : rg1(rg1_), rg2(rg2_), pfem_proj(0), pmf_u2_proj(0),
      contact_only(contact_only_), option(option_)
    {
      Tresca_version = false;   // for future version ...
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
    { if (pmf_u2_proj) delete pmf_u2_proj; }

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
    case 1 : case 4 :
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
    case 3 :
      tl.push_back(model::term_description(varname_u1, multname_n, false)); // U1L
      tl.push_back(model::term_description(multname_n, varname_u1, false)); // LU1
      tl.push_back(model::term_description(varname_u2, multname_n, false)); // U2L
      tl.push_back(model::term_description(multname_n, varname_u2, false)); // LU2
      tl.push_back(model::term_description(multname_n, multname_n, true));  // LL
      tl.push_back(model::term_description(varname_u1, varname_u1, true));  // U1U1
      tl.push_back(model::term_description(varname_u2, varname_u2, true));  // U2U2
      tl.push_back(model::term_description(varname_u1, varname_u2, true));  // U1U2
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
   const std::string &dataname_r, const std::string &dataname_friction_coeff,
   size_type region1, size_type region2, int option,
   const std::string &dataname_alpha,
   const std::string &dataname_wt1, const std::string &dataname_wt2) {

    pbrick pbr = new integral_contact_nonmatching_meshes_brick
                     (region1, region2, false /* contact_only */, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 4 : case 5 : case 6 :
      tl.push_back(model::term_description(varname_u1, multname, false));  // 0: U1L
      tl.push_back(model::term_description(multname, varname_u1, false));  // 1: LU1
      tl.push_back(model::term_description(varname_u2, multname, false));  // 2: U2L
      tl.push_back(model::term_description(multname, varname_u2, false));  // 3: LU2
      tl.push_back(model::term_description(multname, multname, true));     // 4: LL
      break;
    case 2 : case 3 :
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

    model::varnamelist dl(1, dataname_r);  // 0 -> r
    dl.push_back(dataname_friction_coeff); // 1 -> f_coeff
    if (dataname_alpha.size()) {
      dl.push_back(dataname_alpha);        // 2 -> alpha
      if (dataname_wt1.size()) {
        dl.push_back(dataname_wt1);        // 3 -> WT1
        if (dataname_wt2.size()) {
          dl.push_back(dataname_wt2);      // 4 -> WT2
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
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
            pmf_coeff, f_coeff, alpha, WT1, WT2);

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
   const getfem::mesh_fem *pmf_coeff, const VECT1 *f_coeff, scalar_type r,
   scalar_type alpha, const VECT1 *WT1, const VECT1 *WT2,
   const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_V7; break;
    case 2 : subterm =  RHS_U_V6; break;
    case 3 : subterm =  RHS_U_V8; break;
    }

    contact_nonmatching_meshes_nonlinear_term
      nterm(subterm, r, mf_u1, U1, mf_u2, U2, pmf_lambda, lambda,
            pmf_coeff, f_coeff, alpha, WT1, WT2);

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
    mutable getfem::pfem pfem_proj;        // cached fem and mesh_fem for the
    mutable getfem::mesh_fem *pmf_u2_proj; // projection between nonmatching meshes
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

      // Data : r, [lambda,] [friction_coeff,] [alpha,] [WT1, WT2]
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

      const model_real_plain_vector *f_coeff = 0;
      const mesh_fem *pmf_coeff = 0;
      scalar_type alpha = 1;
      const model_real_plain_vector *WT1 = 0;
      const model_real_plain_vector *WT2 = 0;
      if (!contact_only) {
        nd++;
        f_coeff = &(md.real_variable(dl[nd]));
        pmf_coeff = md.pmesh_fem_of_variable(dl[nd]);
        sl = gmm::vect_size(*f_coeff);
        if (pmf_coeff) { sl *= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
        GMM_ASSERT1(sl == 1,
                  "the data corresponding to the friction coefficient "
                  "has not the right format");

        if (dl.size() > nd) {
          nd++;
          alpha = md.real_variable(dl[nd])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[nd])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() > nd) {
          nd++;
          WT1 = &(md.real_variable(dl[nd]));
        }

        if (dl.size() > nd) {
          nd++;
          WT2 = &(md.real_variable(dl[nd]));
        }
      }

      GMM_ASSERT1(matl.size() == contact_only ? 3 : 4,
                  "Wrong number of terms for penalized contact "
                  "between nonmatching meshes brick");

      mesh_region rg(region);
      mf_u1.linked_mesh().intersect_with_mpi_region(rg); // FIXME: mfu_2?

      // projection of the second mesh_fem onto the mesh of the first mesh_fem
      if (!pmf_u2_proj) {
        pmf_u2_proj = new getfem::mesh_fem(mim.linked_mesh(), dim_type(N));
        pfem_proj = new_projected_fem(mf_u2, mim, rg2, rg1);
        pmf_u2_proj->set_finite_element(mim.linked_mesh().convex_index(), pfem_proj);
      }

      size_type nbdof1 = mf_u1.nb_dof();
      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = pmf_u2_proj->nb_dof();

      std::vector<size_type> ind;
      pmf_u2_proj->get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      gmm::csc_matrix<scalar_type> Rsub(nbdof2, nbsub), Esub(nbsub, nbdof2);
      if (mf_u2.is_reduced()) {
          gmm::copy(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    Rsub);
          gmm::copy(gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    Esub);
      }

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(Esub, u2, u2_proj);
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Penalized contact between nonmatching meshes tangent term");
        gmm::clear(matl[0]);
        gmm::clear(matl[1]);
        gmm::clear(matl[2]);

        model_real_sparse_matrix Ku2u2(nbsub,nbsub);
        model_real_sparse_matrix Ku1u2(nbdof1,nbsub);

        if (contact_only) {
          asm_penalized_contact_nonmatching_meshes_tangent_matrix
            (matl[0], Ku2u2, Ku1u2, mim, mf_u1, u1, *pmf_u2_proj, u2_proj,
             pmf_lambda, lambda, vr[0], rg, option);
        }
        else {
          gmm::clear(matl[3]);
          model_real_sparse_matrix Ku2u1(nbsub,nbdof1);
          asm_penalized_contact_nonmatching_meshes_tangent_matrix
            (matl[0], Ku2u2, Ku1u2, Ku2u1, mim, mf_u1, u1, *pmf_u2_proj, u2_proj,
             pmf_lambda, lambda, pmf_coeff, f_coeff, vr[0], alpha, WT1, WT1, rg, option);
          gmm::copy(Ku2u1, gmm::sub_matrix(matl[3], SUBI, gmm::sub_interval(0, nbdof1)));
        }

        if (mf_u2.is_reduced()) {
          model_real_sparse_matrix tmp(nbsub, nbdof2);
          gmm::mult(Ku2u2, Esub, tmp);
          gmm::mult(Rsub, tmp, matl[1]);
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
            (vecl[0], Ru2, mim, mf_u1, u1, *pmf_u2_proj, u2_proj, pmf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_penalized_contact_nonmatching_meshes_rhs
            (vecl[0], Ru2, mim, mf_u1, u1, *pmf_u2_proj, u2_proj, pmf_lambda, lambda,
             pmf_coeff, f_coeff, vr[0], alpha, WT1, WT2, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(Rsub, Ru2, vecl[1]);
        else
          gmm::copy(Ru2, gmm::sub_vector(vecl[1], SUBI));
      }
    }

    penalized_contact_nonmatching_meshes_brick(size_type rg1_, size_type rg2_,
                                               bool contact_only_, int option_)
    : rg1(rg1_), rg2(rg2_), pfem_proj(0), pmf_u2_proj(0),
      contact_only(contact_only_), option(option_) {
      Tresca_version = false;   // for future version ...
      set_flags(contact_only
                ? "Integral penalized contact between nonmatching meshes brick"
                : "Integral penalized contact and friction between nonmatching "
                  "meshes brick",
                false /* is linear*/, contact_only /* is symmetric */,
                true /* is coercive */, true /* is real */,
                false /* is complex */);
    }

    ~penalized_contact_nonmatching_meshes_brick()
    { if (pmf_u2_proj) delete pmf_u2_proj; }

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
   const std::string &dataname_friction_coeff,
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
    dl.push_back(dataname_friction_coeff);
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

      area = asm_level_set_contact_area(*ml[0], mf_u, u, mf_obs, obs, reg, -1e-3);

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

      std::vector<size_type> ind;
      p->pmf_u2_proj->get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = p->pmf_u2_proj->nb_basic_dof();
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

      area = asm_nonmatching_meshes_contact_area
             (*ml[0], mf_u1, u1, *(p->pmf_u2_proj), u2_proj, reg, -1e-3);

      gmm::resize(F, mf_u1.nb_dof());
      asm_nonmatching_meshes_normal_source_term
        (F, *ml[0], mf_u1, *(p->pmf_u2_proj), mf_lambda, lambda, reg);

    }
    else if (pbr->brick_name() == "Integral penalized contact between nonmatching meshes brick" ||
             pbr->brick_name() == "Integral penalized contact and friction between nonmatching "
                                  "meshes brick") {
      penalized_contact_nonmatching_meshes_brick *p
        = dynamic_cast<penalized_contact_nonmatching_meshes_brick *>
         (const_cast<virtual_brick *>(pbr.get()));
      GMM_ASSERT1(p, "Wrong type of brick");
    }

  }



#ifdef EXPERIMENTAL_PURPOSE_ONLY


  // Experimental implementation of contact condition with Nitsche method.
  // To be deleted when a more general implementation will be designed.


  class contact_nitsche_nonlinear_term : public nonlinear_elem_term {

  protected:
    base_small_vector lnt, lt; // multiplier lambda and its tangential component lambda_t
    scalar_type ln;            // normal component lambda_n of the multiplier
    base_small_vector ut;      // tangential relative displacement
    scalar_type un;            // normal relative displacement (positive when the first
                               // elastic body surface moves outwards)
    base_small_vector no, n;   // surface normal, pointing outwards with respect
                               // to the (first) elastic body
    scalar_type g, f_coeff;    // gap and coefficient of friction values
    scalar_type lambda, mu;    // Lame coefficients

    base_small_vector aux1, auxN, V;
    base_matrix GP, grad;
    base_vector coeff;
    const mesh_fem &mf_u;       // mandatory
    const mesh_fem &mf_obs;     // mandatory
    const mesh_fem *pmf_coeff;
    base_vector U, obs, friction_coeff;

    void adjust_tensor_size(void) {
      sizes_.resize(4); sizes_[0] = sizes_[1] = sizes_[2] = sizes_[3] = N;
      switch (option) {
      case 1 : sizes_.resize(1); break;
      case 2 : case 3 :  sizes_.resize(2); break;
      case 4 : case 5 :  sizes_.resize(3); break;
      }
      gmm::resize(grad, 1, N);
      lnt.resize(N); lt.resize(N); ut.resize(N); no.resize(N); n.resize(N);
      aux1.resize(1); auxN.resize(N); V.resize(N);
      gmm::resize(GP, N, N);
    }

  public:
    dim_type N;
    size_type option;
    scalar_type r;

    bgeot::multi_index sizes_;

    template <typename VECT1>
    contact_nitsche_nonlinear_term(size_type option_, scalar_type r_,
				   scalar_type lambda_, scalar_type mu_,
				   const mesh_fem &mf_u_, const VECT1 &U_,
				   const mesh_fem &mf_obs_, const VECT1 &obs_,
				   const mesh_fem *pmf_coeff_ = 0,
				   const VECT1 *f_coeff_ = 0)		  
      : lambda(lambda_), mu(mu_), mf_u(mf_u_),
      mf_obs(mf_obs_), pmf_coeff(pmf_coeff_), U(mf_u.nb_basic_dof()),
	obs(mf_obs.nb_basic_dof()), friction_coeff(0), option(option_), r(r_) {
      N = mf_u_.linked_mesh().dim();
      adjust_tensor_size();

      mf_u.extend_vector(U_, U);
      mf_obs.extend_vector(obs_, obs);

      if (!pmf_coeff)
	f_coeff = (*f_coeff_)[0];
      else {
	friction_coeff.resize(pmf_coeff->nb_basic_dof());
	pmf_coeff->extend_vector(*f_coeff_, friction_coeff);
      }
    }

    const bgeot::multi_index &sizes() const { return sizes_; }

    virtual void compute(fem_interpolation_context&, bgeot::base_tensor &t);
    virtual void prepare(fem_interpolation_context& /*ctx*/, size_type /*nb*/);

  };

  void contact_nitsche_nonlinear_term::compute
  (fem_interpolation_context &/* ctx */, bgeot::base_tensor &t) {

    t.adjust_sizes(sizes_);
    scalar_type e;
    dim_type i, j, k, l;

    if (option >= 3) { // computation of matrix A
      e = f_coeff*gmm::neg(ln-r*(un-g));
      auxN = lt - r*ut;
      ball_projection_grad(auxN, e, GP);
      ball_projection_grad_r(auxN, e, V);
      e = Heav(r*(un-g) - ln);
      gmm::rank_one_update(GP, no, gmm::scaled(V, -e*f_coeff));
      gmm::rank_one_update(GP, gmm::scaled(no, e-gmm::vect_sp(GP,no,no)), no);
      gmm::scale(GP, 1./r);
    } else { // computation of vector W
      e = gmm::neg(ln-r*(un-g));
      V = lt - r*ut;
      ball_projection(V, f_coeff*e);
      V -= e*no;
    }

    switch (option) {
      // one-dimensional tensors [N]
    case 1:
      for (i=0; i < N; ++i) t[i] = V[i];
      break;

      // two-dimensional tensors [N x N]
    case 2:
      V -= lnt;
      gmm::scale(V, -1./r);
      e = gmm::vect_sp(V, n);
      for (i=0; i < N; ++i)
	for (j=0; j < N; ++j) {
	  t(i,j) = mu*(V[i]*n[j]+V[j]*n[i]);
	  if (i == j) t(i,j) += lambda*e;
	}
      break;

    case 3:
      for (i=0; i < N; ++i)
	for (j=0; j < N; ++j)
	  t(i,j) = r*r*GP(j,i);
      break;

    // three-dimensional tensors [N x N x N]
    case 4:
      gmm::mult(gmm::transposed(GP), n, V);
      for (i=0; i < N; ++i)
	for (j=0; j < N; ++j)
	  for (k=0; k < N; ++k) {
	    t(i,j,k) = -r*mu*(GP(j,i)*n[k] + GP(k,i)*n[j]);
	    if (j == k) t(i,j,k) -= r*lambda*V[i];
	  } 
      break;
       
    case 5:
      gmm::mult(GP, n, V);
      for (i=0; i < N; ++i)
	for (j=0; j < N; ++j)
	  for (k=0; k < N; ++k) {
	    t(i,j,k) = -r*mu*(GP(k,i)*n[j] + GP(k,j)*n[i]);
	    if (i == j) t(i,j,k) -= r*lambda*V[k];
	  } 
      break;
      
    // four-dimensional tensors [N x N x N x N]

    case 6:

      for (i=0; i < N; ++i) GP(i,i) -= 1./r;  // matrix B

      e = gmm::vect_sp(GP, n, n);
      gmm::mult(gmm::transposed(GP), n, auxN);
      gmm::mult(GP, n, V);

      for (i=0; i < N; ++i)
	for (j=0; j < N; ++j)
	  for (k=0; k < N; ++k)
	    for (l=0; l < N; ++l) {
	      t(i,j,k,l) = mu*mu*(n[i]*GP(k,j)*n[l] + n[j]*GP(k,i)*n[l]
				  + n[j]*GP(l,i)*n[k] + n[i]*GP(l,j)*n[k]);
	      if (i == j && k == l) t(i,j,k,l) += lambda*lambda*e;
	      if (i == j) t(i,j,k,l) += lambda*mu*(V[k]*n[l] + V[l]*n[k]);
	      if (k == l) t(i,j,k,l) += lambda*mu*(auxN[j]*n[i]+auxN[i]*n[j]);
	    }

      break;
    default : GMM_ASSERT1(false, "Invalid option");
    }
  }


  void contact_nitsche_nonlinear_term::prepare
  (fem_interpolation_context& ctx, size_type nb) {
    size_type cv = ctx.convex_num();

    switch (nb) { // last is computed first
    case 1 : // calculate [un] and [ut] interpolating [U] on [mf_u]
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index
				(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, V, N);
      un = gmm::vect_sp(V, no);
      ut = V - un * no;
      ctx.pf()->interpolation_grad(ctx, coeff, GP, N);
      lnt = lambda*(gmm::mat_trace(GP))*n;
      gmm::mult_add(GP, gmm::scaled(n, mu), lnt);
      gmm::mult_add(gmm::transposed(GP), gmm::scaled(n, mu), lnt);      
      ln = gmm::vect_sp(lnt, no);
      lt = lnt - ln * no;
      break;

    case 2 : // calculate [g] and [no] interpolating [obs] on [mf_obs]
             // calculate [ln] and [lt] from [lnt] and [no]
      coeff.resize(mf_obs.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (obs, gmm::sub_index
                 (mf_obs.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, grad, 1);
      gmm::copy(gmm::mat_row(grad, 0), no);
      no /= -gmm::vect_norm2(no);
      ctx.pf()->interpolation(ctx, coeff, aux1, 1);
      g = aux1[0];
      n = bgeot::compute_normal(ctx, ctx.face_num());
      n /= gmm::vect_norm2(n);
      break;

    case 3 :// calculate [f_coeff] interpolating [friction_coeff] on [mf_coeff]
      if (pmf_coeff) {
        coeff.resize(pmf_coeff->nb_basic_dof_of_element(cv));
        gmm::copy(gmm::sub_vector
                  (friction_coeff, gmm::sub_index
                   (pmf_coeff->ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        f_coeff = aux1[0];
      }
      break;

    default : GMM_ASSERT1(false, "Invalid option");
    }
  }




  template<typename MAT, typename VECT1>
  void asm_Nitsche_contact_rigid_obstacle_tangent_matrix
  (MAT &K, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   scalar_type gamma, scalar_type lambda, scalar_type mu,
   const mesh_region &rg, int option = 1) {

    contact_nitsche_nonlinear_term
      nterm1(6, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff),
      nterm2(3, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff),
      nterm3(4, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff),
      nterm4(5, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff);

    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";

    getfem::generic_assembly assem;
    std::string as_str 
      = ((option == 0) ? "w1=comp(NonLin$1(#1,"+aux_fems+")(i,j,k,l).vGrad(#1)(:,i,j).vGrad(#1)(:,k,l));" : "")
      + "w2=comp(NonLin$2(#1,"+aux_fems+").vBase(#1).vBase(#1))(i,j,:,i,:,j);"
      + "w3=comp(NonLin$3(#1,"+aux_fems+").vBase(#1).vGrad(#1))(i,j,k,:,i,:,j,k);"
      + ((option == 0) ? "w4=comp(NonLin$4(#1,"+aux_fems+").vGrad(#1).vBase(#1))(i,j,k,:,i,j,:,k);" : "")
      + ((option == 0) ? "M(#1,#1)+=w1+w2+w3+w4;" : "M(#1,#1)+=w2+w3;");

    assem.set(as_str);
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    assem.push_nonlinear_term(&nterm4);
    assem.push_mat(K);
    assem.assembly(rg);
  }



  template<typename VECT1>
  void asm_Nitsche_contact_rigid_obstacle_rhs
  (VECT1 &R, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   scalar_type gamma, scalar_type lambda, scalar_type mu,
   const mesh_region &rg, int option = 1) {

    contact_nitsche_nonlinear_term
      nterm1(1, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff),
      nterm2(2, gamma, lambda, mu, mf_u, U, mf_obs, obs, pmf_coeff, &f_coeff);
    
    const std::string aux_fems = pmf_coeff ? "#1,#2,#3" : "#1,#2";
    
    getfem::generic_assembly assem;
    std::string as_str = 
      "V(#1)+=comp(NonLin$1(#1,"+aux_fems+").vBase(#1))(i,:,i); "
      + ((option == 0) ? "V(#1)+=comp(NonLin$2(#1,"+aux_fems+").vGrad(#1))(i,j,:,i,j)" : "");

    assem.set(as_str);
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_obs);
    if (pmf_coeff) assem.push_mf(*pmf_coeff);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_vec(R);
    assem.assembly(rg);
  }


  struct Nitsche_contact_rigid_obstacle_brick : public virtual_brick {

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

      // Variables : u
      GMM_ASSERT1(vl.size() == 1,
                  "Nitsche contact with rigid obstacle bricks need a "
		  "single variable");
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);

      // Data : obs, r, [lambda,] [friction_coeff,] [alpha,] [WT]
      GMM_ASSERT1(dl.size() == 5, "Wrong number of data for Nitsche "
		  "contact with rigid obstacle brick");

      const model_real_plain_vector &obs = md.real_variable(dl[0]);
      const mesh_fem &mf_obs = md.mesh_fem_of_variable(dl[0]);
      size_type sl = gmm::vect_size(obs) * mf_obs.get_qdim() / mf_obs.nb_dof();
      GMM_ASSERT1(sl == 1, "the data corresponding to the obstacle has not "
                  "the right format");

      const model_real_plain_vector &vr = md.real_variable(dl[1]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");

      const model_real_plain_vector *f_coeff = 0;
      const mesh_fem *pmf_coeff = 0;
      
      f_coeff = &(md.real_variable(dl[2]));
      pmf_coeff = md.pmesh_fem_of_variable(dl[2]);
      sl = gmm::vect_size(*f_coeff);
      if (pmf_coeff) { sl*= pmf_coeff->get_qdim(); sl /= pmf_coeff->nb_dof(); }
      GMM_ASSERT1(sl == 1, "the data corresponding to the friction "
		  "coefficient has not the right format");
      
      const model_real_plain_vector &vlambda = md.real_variable(dl[3]);
      GMM_ASSERT1(gmm::vect_size(vlambda) == 1,
		  "Parameter lambda should be a scalar");
      const model_real_plain_vector &vmu = md.real_variable(dl[4]);
      GMM_ASSERT1(gmm::vect_size(vmu) == 1, "Parameter mu should be a scalar");


      GMM_ASSERT1(matl.size() == 1, "Wrong number of terms for "
                  "Nitsche contact with rigid obstacle brick");

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Nitsche contact with rigid obstacle tangent term");
        gmm::clear(matl[0]);
	asm_Nitsche_contact_rigid_obstacle_tangent_matrix
	  (matl[0], mim, mf_u, u, mf_obs, obs,  pmf_coeff, *f_coeff,
	   vr[0], vlambda[0], vmu[0], rg);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
	asm_Nitsche_contact_rigid_obstacle_rhs
	  (vecl[0], mim, mf_u, u, mf_obs, obs, pmf_coeff, *f_coeff,
	   vr[0], vlambda[0], vmu[0], rg);
      }

    }

    Nitsche_contact_rigid_obstacle_brick(void) {
      set_flags("Integral Nitsche contact and friction with rigid "
                "obstacle brick",
                false /* is linear*/, false /* is symmetric */,
                true /* is coercive */, true /* is real */,
                false /* is complex */);
    }

  };


  size_type add_Nitsche_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   const std::string &dataname_friction_coeff,
   const std::string &dataname_lambda, const std::string &dataname_mu,
   size_type region) {

    pbrick pbr = new Nitsche_contact_rigid_obstacle_brick;

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));

    model::varnamelist dl(1, dataname_obs);
    dl.push_back(dataname_r);
    dl.push_back(dataname_friction_coeff);
    dl.push_back(dataname_lambda);
    dl.push_back(dataname_mu);

    model::varnamelist vl(1, varname_u);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

#endif



  //=========================================================================
  //
  //  Large sliding brick.
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
    std::vector<size_type> cvdof1(mf1.ind_basic_dof_of_element(cv1).begin(),
				  mf1.ind_basic_dof_of_element(cv1).end());
    std::vector<size_type> cvdof2(mf2.ind_basic_dof_of_element(cv2).begin(),
				  mf2.ind_basic_dof_of_element(cv2).end());

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
      if (N > 2) coordinates[3] = "z";
      if (N > 3) coordinates[4] = "w";
      GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
		  "dimension greater than 4");
    }

    size_type add_obstacle(const std::string &obs) {
      size_type ind = obstacles.size();
      obstacles.push_back(obs);
      obstacles_velocities.push_back("");
      mu::Parser mu;
      obstacles_parsers.push_back(mu);
      obstacles_parsers[ind].SetExpr(obstacles[ind]);
      for (size_type k = 0; k < N; ++k)
	obstacles_parsers[ind].DefineVar(coordinates[k], &pt_eval[k]);
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
	size_type cvnbdof = mfu.nb_basic_dof_of_element(cv);
	coeff.resize(cvnbdof);
	mesh_fem::ind_dof_ct::const_iterator
	  itdof = mfu.ind_basic_dof_of_element(cv).begin();
	for (size_type k = 0; k < cvnbdof; ++k, ++itdof) coeff[k]=U[*itdof];
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
   scalar_type f_coeff, scalar_type r, model::build_version version) {
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
    gmm::copy(gmm::sub_vector
	      (U, gmm::sub_index
	       (mfu.ind_basic_dof_of_element(cv))), coeff);
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
    std::vector<base_node> y0s, y0_refs;
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
      pfem pf_s = mfu_y0.fem_of_element(cv_y0);
      const model_real_plain_vector &U_y0
	= cf.disp_of_boundary(boundary_num_y0);
      const mesh &m = mfu_y0.linked_mesh();
      bgeot::pgeometric_trans pgt_y0 = m.trans_of_convex(cv_y0);
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
      
      size_type cvnbdof_y0 = mfu_y0.nb_basic_dof_of_element(cv_y0);
      coeff.resize(cvnbdof_y0);
      mesh_fem::ind_dof_ct::const_iterator
	itdof = mfu_y0.ind_basic_dof_of_element(cv_y0).begin();
      for (size_type k = 0; k < cvnbdof_y0; ++k, ++itdof)
	coeff[k] = U_y0[*itdof];
      // if (pf_s->need_G()) 
      bgeot::vectors_to_base_matrix
	(G, mfu_y0.linked_mesh().points_of_convex(cv_y0));
      
      fem_interpolation_context ctx_y0(pgt_y0, pf_s, y0_ref, G, cv_y0,
				       size_type(-1));
      
      size_type newton_iter = 0;
      for(;;) { // Newton algorithm to invert geometric transformation
	
	pf_s->interpolation(ctx_y0, coeff, val, dim_type(N));
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
	
	pf_s->interpolation_grad(ctx_y0, coeff, grad, dim_type(N));
	
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
	  pf_s->interpolation(ctx_y0, coeff, val, dim_type(N));
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

      
//       size_type iptf = m.ind_points_of_face_of_convex(cv_y0, face_y0)[0];
//       base_node ptf = x0 - m.points()[iptf];
//       scalar_type d2 = gmm::vect_sp(ptf, n0_y0) / gmm::vect_norm2(n0_y0);


      
      if (noisy) cout << "gmm::vect_norm2(n0_y0) = " << gmm::vect_norm2(n0_y0) << endl;
      // Eliminates wrong auto-contact situations
      if (noisy) cout << "autocontact status : x0 = " << x0 << " y0 = " << y0 << "  " <<  gmm::vect_dist2(y0, x0) << " : " << d0*0.75 << " : " << d1*0.75 << endl;
      if (noisy) cout << "n = " << n << " unit_normal_of_elements[(*it)->id] = " << unit_normal_of_elements[(*it)->id] << endl;



      if (d0 < scalar_type(0)
	  && ((&(U_y0) == &U
	       && (gmm::vect_dist2(y0, x0) < gmm::abs(d1)*scalar_type(3)/scalar_type(4)))
	      || gmm::abs(d1) > 0.05)) {
	if (noisy)  cout << "Eliminated x0 = " << x0 << " y0 = " << y0
			<< " d0 = " << d0 << endl;
	continue;
      }


//       if (d0 < scalar_type(0) && &(U_y0) == &U
// 	  && gmm::vect_dist2(y0, x0) < gmm::abs(d1) * scalar_type(2)
// 	  && d2 < -ctxu.J() / scalar_type(2)) {
// 	/*if (noisy) */ cout << "Eliminated x0 = " << x0 << " y0 = " << y0
// 			<< " d0 = " << d0 << endl;
// 	continue;
//       }
      
      y0s.push_back(ctx_y0.xreal()); // Usefull ?
      y0_refs.push_back(y0_ref);
      elt_nums.push_back((*it)->id);
      d0s.push_back(d0);
      d1s.push_back(d1);
      ctx_y0s.push_back(ctx_y0);
      n0_y0 /= gmm::vect_norm2(n0_y0);
      n0_y0s.push_back(n0_y0);
      
      if (noisy) cout << "dist0 = " << d0 << " dist1 = "
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
    for (size_type k = 0; k < y0_refs.size(); ++k)
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
      size_type nbo = boundary_of_elements[elt_nums[ibound]];
      const mesh_fem &mfu_y0 = cf.mfu_of_boundary(nbo);
      const mesh &m = mfu_y0.linked_mesh();
      size_type icv = ind_of_elements[elt_nums[ibound]];
      
      if (noisy) cout << " y0 = " << y0s[ibound] << " of element "
			    << icv  << " of boundary " << nbo << endl;
      for (size_type k = 0; k < m.nb_points_of_convex(icv); ++k)
	if (noisy) cout << "point " << k << " : "
			<< m.points()[m.ind_points_of_convex(icv)[k]] << endl;
      if (nbo == 0 && boundary_num == 0 && d0 < 0.0 && (version & model::BUILD_MATRIX)) GMM_ASSERT1(false, "oups");
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
    base_small_vector lambda(N), zeta(N), vv(N);
    ctxl.base_value(tl);
    ctxu.base_value(tu);
    
    coeff.resize(cvnbdofl);
    gmm::copy(gmm::sub_vector
	      (L, gmm::sub_index
	       (mfl.ind_basic_dof_of_element(cv))), coeff);
    ctxl.pf()->interpolation(ctxl, coeff, lambda, dim_type(N));
    GMM_ASSERT1(!(std::isnan(lambda[0])), "internal error");

    // Tangent term -(1/r)\int \delta\lambda.\mu
    if (version & model::BUILD_MATRIX) {
      gmm::resize(Melem, cvnbdofl, cvnbdofl); gmm::clear(Melem);
      for (size_type i = 0; i < cvnbdofl; ++i)
	for (size_type j = 0; j < cvnbdofl; ++j)
	  if (i%N == j%N) Melem(i,j) = -tl[i/N]*tl[j/N]*weight/r;
      mat_elem_assembly(cf.LL_matrix(boundary_num, boundary_num),
			Melem, mfl, cv, mfl, cv);
    }
    
    // Rhs term (1/r)\int (\lambda - P(\zeta)).\mu
    // Unstabilized frictionless case for the moment
    if (state) gmm::add(lambda, gmm::scaled(n, r*d0), zeta);
    if (version & model::BUILD_RHS) {
      gmm::clear(vv);
      if (state) {
	gmm::copy(zeta, vv);
	De_Saxce_projection(vv, n, scalar_type(0));
	gmm::scale(vv, -scalar_type(1));
	gmm::add(lambda, vv);
      } else gmm::copy(lambda, vv);
      gmm::resize(Velem,  cvnbdofl); gmm::clear(Velem);
      for (size_type i = 0; i < cvnbdofl; ++i)
	Velem[i] = (tl[i/N] * vv[i%N])*weight/r;
      vec_elem_assembly(cf.L_vector(boundary_num), Velem, mfl, cv);
    }

    if (state) {
      base_matrix grad_y0(N, N), gradinv_y0(N, N), gradaux(N,N);
      base_vector coeff_y0;
      base_small_vector vvv(N), ntilde_y0(N);
      base_tensor tgradu, tu_y0, tgradu_y0;
      size_type cv_y0 = 0, cvnbdofu_y0 = 0;
      size_type boundary_num_y0
	= (state == 1) ? boundary_of_elements[elt_nums[ibound]] : 0;
      const mesh_fem &mfu_y0
	= (state == 1) ? cf.mfu_of_boundary(boundary_num_y0) : mfu;
      ctxu.grad_base_value(tgradu);
      
      if (state == 1) {
	cv_y0 = ind_of_elements[elt_nums[ibound]];
	cvnbdofu_y0 = mfu_y0.nb_basic_dof_of_element(cv_y0);
	const model_real_plain_vector &U_y0
	  = cf.disp_of_boundary(boundary_num_y0);
	mesh_fem::ind_dof_ct::const_iterator
	  itdof = mfu_y0.ind_basic_dof_of_element(cv_y0).begin();
	coeff_y0.resize(cvnbdofu_y0);
	gmm::copy(gmm::sub_vector
		  (U_y0, gmm::sub_index
		   (mfu_y0.ind_basic_dof_of_element(cv_y0))), coeff_y0);
	ctx_y0s[ibound].pf()->interpolation_grad(ctx_y0s[ibound], coeff_y0,
						 grad_y0, dim_type(N));
	gmm::add(gmm::identity_matrix(), grad_y0);
	gmm::copy(grad_y0, gradinv_y0);
	gmm::lu_inverse(gradinv_y0);// à proteger contre la non-inversibilité
	ctx_y0s[ibound].base_value(tu_y0);
	ctx_y0s[ibound].grad_base_value(tgradu_y0);
	gmm::mult(gmm::transposed(gradinv_y0), n0_y0s[ibound], ntilde_y0); // (not unit) normal vector
      }
      
      // Rhs term \int \lambda.(\psi(x_0) - \psi(y_0))
      if (version & model::BUILD_RHS) {
	gmm::resize(Velem,  cvnbdofu);gmm::clear(Velem);
	for (size_type i = 0; i < cvnbdofu; ++i)
	  Velem[i] = tu[i/N] * lambda[i%N]*weight;
	vec_elem_assembly(cf.U_vector(boundary_num), Velem, mfu, cv);
	
	if (state == 1) {
	  gmm::resize(Velem,  cvnbdofu_y0); gmm::clear(Velem);
	  for (size_type i = 0; i < cvnbdofu_y0; ++i)
	    Velem[i] = -tu_y0[i/N] * lambda[i%N]*weight;
	  vec_elem_assembly(cf.U_vector(boundary_num_y0), Velem, mfu_y0,cv_y0);
	}
      }	
      
      if (version & model::BUILD_MATRIX) {
	// Tangent term \int (\delta \lambda).(\psi(y_0) - \psi(x_0))
	gmm::resize(Melem, cvnbdofu, cvnbdofl); gmm::clear(Melem);
	for (size_type i = 0; i < cvnbdofu; ++i)
	  for (size_type j = 0; j < cvnbdofl; ++j)
	    if (i%N == j%N) Melem(i,j) = -tu[i/N]*tl[j/N]*weight;
	mat_elem_assembly(cf.UL_matrix(boundary_num, boundary_num),
			  Melem, mfu, cv, mfl, cv);
	
	if (state == 1) {
	  gmm::resize(Melem, cvnbdofu_y0, cvnbdofl); gmm::clear(Melem);
	  for (size_type i = 0; i < cvnbdofu_y0; ++i)
	    for (size_type j = 0; j < cvnbdofl; ++j)
	      if (i%N == j%N) Melem(i,j) = tu_y0[i/N]*tl[j/N]*weight;
	  mat_elem_assembly(cf.UL_matrix(boundary_num_y0, boundary_num),
			    Melem, mfu_y0, cv_y0, mfl, cv);
	}
	
	// Tangent term \int \lambda.((\nabla \psi(y_0))(I+\nabla u(y_0))^{-1}(\delta u(x_0) - \delta u(y_0)))
	if (state == 1) {
	  gmm::resize(Melem, cvnbdofu_y0, cvnbdofu); gmm::clear(Melem);
	  for (size_type i = 0; i < cvnbdofu_y0; ++i)
	    for (size_type j = 0; j < cvnbdofu; ++j)
	      for (size_type k = 0; k < N; ++k)
		Melem(i, j) += lambda[i%N] * tgradu_y0[i-(i%N)+k]
		  * gradinv_y0(k, j%N) * tu[j/N]*weight;
	  mat_elem_assembly(cf.UU_matrix(boundary_num_y0, boundary_num),
			    Melem, mfu_y0, cv_y0, mfu, cv);
	  
	  gmm::resize(Melem, cvnbdofu_y0, cvnbdofu_y0); gmm::clear(Melem);
	  for (size_type i = 0; i < cvnbdofu_y0; ++i)
	    for (size_type j = 0; j < cvnbdofu_y0; ++j)
	      for (size_type k = 0; k < N; ++k)
		Melem(i, j) -= lambda[i%N] * tgradu_y0[i-(i%N)+k]
		  * gradinv_y0(k, j%N) * tu_y0[j/N]*weight;
	  mat_elem_assembly(cf.UU_matrix(boundary_num_y0, boundary_num_y0),
			    Melem, mfu_y0, cv_y0, mfu_y0, cv_y0);
	}
	
	// Tangent term (1/r)\int \nabla P(zeta) (dzeta/dlambda)(\delta lambda) . \mu
	De_Saxce_projection_grad(zeta, n, scalar_type(0), grad);
	gmm::resize(Melem, cvnbdofl, cvnbdofl); gmm::clear(Melem);
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofl; ++j)
	    Melem(i,j) = tl[i/N]*tl[j/N]*grad(i%N,j%N)*weight/r;
	mat_elem_assembly(cf.LL_matrix(boundary_num, boundary_num),
			  Melem, mfl, cv, mfl, cv);

      
	// Tangent term \int (I+\nabla u(y_0))^{-T}\nabla delta(y_0).\delta u(x_0)(\nabla P(zeta) n . \mu)
	gmm::mult(grad, n, vv);
	gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofu; ++j)
	    Melem(i, j) = tl[i/N]*vv[i%N]*tu[j/N]
	      *((state == 1) ? ntilde_y0[j%N] : grad_obs[j%N])*weight;
	mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num),
			  Melem, mfl, cv, mfu, cv);
	
	// Tangent term -\int (I+\nabla u(y_0))^{-T}\nabla delta(y_0).\delta u(y_0)(\nabla P(zeta) n . \mu)
	if (state == 1) {
	  gmm::resize(Melem, cvnbdofl, cvnbdofu_y0); gmm::clear(Melem);
	  for (size_type i = 0; i < cvnbdofl; ++i)
	    for (size_type j = 0; j < cvnbdofu_y0; ++j)
	      Melem(i, j) = -tl[i/N]*vv[i%N]*tu_y0[j/N]*ntilde_y0[j%N]*weight;
	  mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num_y0),
			    Melem, mfl, cv, mfu_y0, cv_y0);
	}
	
	// Tangent term \int d_0(\nabla P)(dn/du)(\delta u).\mu
	gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);
	gmm::mult(grad, n, vv);
	gmm::mult(gradinv, n, vvv);
	gmm::mult(gradinv, gmm::transposed(grad), gradaux);
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofu; ++j)
	    for (size_type k = 0; k < N; ++k)
	      Melem(i,j) += d0*tl[i/N]*vv[i%N]
		*tgradu[j-(j%N)+k]*n[j%N]*vvv[k]*weight;
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofu; ++j)
	    for (size_type k = 0; k < N; ++k)
	      Melem(i,j) -= d0*tl[i/N]*gradaux(k,i%N)*tgradu[j-(j%N)+k]
		*n[j%N]*weight;

	
	
	// Tangent term (1/r)\int \nabla_n P(zeta) (dn/du)(\delta u) . \mu
	// On peut certainement factoriser d'avantage ce terme avec le
	// précédent. Attendre la version avec frottement.
	De_Saxce_projection_gradn(zeta, n, scalar_type(0), grad);
	gmm::mult(gradinv, gmm::transposed(grad), gradaux);
	gmm::mult(grad, n, vv);
	gmm::mult(gradinv, n, vvv);
	// gmm::resize(Melem, cvnbdofl, cvnbdofu); gmm::clear(Melem);factorised
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofu; ++j)
	    for (size_type k = 0; k < N; ++k)
	      Melem(i,j) += tl[i/N]*vv[i%N]
		*tgradu[j-(j%N)+k]*n[j%N]*vvv[k]*weight/r;
	for (size_type i = 0; i < cvnbdofl; ++i)
	  for (size_type j = 0; j < cvnbdofu; ++j)
	    for (size_type k = 0; k < N; ++k)
	      Melem(i,j) -= tl[i/N]*gradaux(k,i%N)*tgradu[j-(j%N)+k]
		*n[j%N]*weight/r;
	mat_elem_assembly(cf.LU_matrix(boundary_num, boundary_num),
			  Melem, mfl, cv, mfu, cv);
      }
    }
    return true;
  }

  //=========================================================================
  // 3)- Large sliding contact brick
  //=========================================================================

  struct integral_large_sliding_contact_brick : public virtual_brick {

    
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

    integral_large_sliding_contact_brick() {
      set_flags("Integral large sliding contact brick",
                false /* is linear*/, false /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

  };




  void integral_large_sliding_contact_brick::asm_real_tangent_terms
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

  size_type add_integral_large_sliding_contact_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &multname, const std::string &dataname_r,
   const std::string &dataname_friction_coeff, size_type region) {

    integral_large_sliding_contact_brick *pbr
      = new integral_large_sliding_contact_brick();
    
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
    integral_large_sliding_contact_brick *p
      = dynamic_cast<integral_large_sliding_contact_brick *>
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
     integral_large_sliding_contact_brick *p
       = dynamic_cast<integral_large_sliding_contact_brick *>
       (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    p->add_obstacle(obs);
  }

}  /* end of namespace getfem.                                             */
