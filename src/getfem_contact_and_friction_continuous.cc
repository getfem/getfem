// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2011 Yves Renard.
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


#include "getfem/getfem_contact_and_friction_continuous.h"
#include "getfem/getfem_projected_fem.h"


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
    // Verified and correct in 2D and 3D.
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    size_type N = gmm::vect_size(x);

    if (xn > scalar_type(0) && f * nxt <= xn) {
      gmm::clear(g);
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      static VEC xt;
      gmm::resize(xt, N);
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      if (N > 2) {
        gmm::copy(gmm::identity_matrix(), g);
        gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
        gmm::rank_one_update(g, gmm::scaled(xt, -scalar_type(1)), xt);
        gmm::scale(g, f*(f - xn/nxt));
      } else {
        gmm::clear(g);
      }

      gmm::scale(xt, -f); gmm::add(n, xt);
      gmm::rank_one_update(g, xt, xt);
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    } else {
      gmm::copy(gmm::identity_matrix(), g);
    }
  }

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
      t[0] = Heav(-ln);  break;

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
  //  Continuous augmented Lagrangian brick (given obstacle, u, lambda).
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
    if (option == 2 || option == 3)
      assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Kul);
    assem.push_mat(Klu);
    assem.push_mat(Kll);
    if (option == 2 || option == 3)
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   const VECT1 *WT, const VECT1 *VT,
   scalar_type r, scalar_type alpha, scalar_type gamma,
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
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma),
      nterm3(subterm3, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma),
      nterm4(subterm4, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma);

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
    if (option == 2 || option == 3)
      assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Kul);
    assem.push_mat(Klu);
    assem.push_mat(Kll);
    if (option == 2 || option == 3)
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   const VECT1 *WT, const VECT1 *VT,
   scalar_type r, scalar_type alpha, scalar_type gamma,
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
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma),
      nterm2(subterm2, r, mf_u, U, mf_obs, obs, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT, alpha, VT, gamma);

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

  struct continuous_contact_rigid_obstacle_brick : public virtual_brick {

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
                  "Continuous contact with rigid obstacle bricks need a single mesh_im");
      GMM_ASSERT1(vl.size() == 2,
                  "Continuous contact with rigid obstacle bricks need two variables");
      GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 7,
                  "Wrong number of data for continuous contact with rigid obstacle "
                  << "brick, " << dl.size() << " should be between 2 and 7.");
      GMM_ASSERT1(matl.size() == size_type(3 + (option == 3)
                                           + (option == 2 && !contact_only)),
                  "Wrong number of terms for "
                  "continuous contact with rigid obstacle brick");

      // variables : u, lambda. The variable lambda should be scalar in the
      //             frictionless case and vector valued in the case with
      //             friction.
      // data      : obstacle, r for the version without friction
      //           : obstacle, r, friction_coeff, alpha, w_t, gamma, v_t for
      //             the version with friction. alpha, w_t , gamma and v_t
      //             are optional and equal to 1, 0, 1 and 0 by default,
      //             respectively.

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
        GMM_TRACE2("Continuous Coulomb friction tangent term");
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
             pmf_coeff, friction_coeff, WT, VT, vr[0], alpha, gamma,
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
             pmf_coeff, friction_coeff, WT, VT, vr[0], alpha, gamma,
             rg, option);
      }

    }

    continuous_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
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

    pbrick pbr = new continuous_contact_rigid_obstacle_brick(true, option);

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
   const std::string &dataname_alpha, const std::string &dataname_wt,
   const std::string &dataname_gamma, const std::string &dataname_vt) {

    pbrick pbr
      = new continuous_contact_rigid_obstacle_brick(false, option);

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
                         "Incorrect option for continuous contact brick");
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
  //  Continuous penalized contact with friction (given obstacle, u, lambda).
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

    getfem::generic_assembly assem;
    if (pmf_lambda)
      assem.set("M(#1,#1)+=comp(NonLin(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)");
    else
      assem.set("M(#1,#1)+=comp(NonLin(#1,#1,#2).vBase(#1).vBase(#1))(i,j,:,i,:,j)");
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

    getfem::generic_assembly assem;
    if (pmf_lambda)
      assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); ");
    else
      assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2).vBase(#1))(i,:,i); ");
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff, const VECT1 *WT,
   scalar_type r, scalar_type alpha, const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  K_UU_FRICT_V4; break;
    case 2 : subterm =  K_UU_FRICT_V3; break;
    case 3 : subterm =  K_UU_FRICT_V5; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(subterm, r, mf_u, U, mf_obs, obs, pmf_lambda, lambda,
            pmf_coeff, &f_coeff, WT, alpha);

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
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff, const VECT1 *WT,
   scalar_type r, scalar_type alpha, const mesh_region &rg, int option = 1) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_V7; break;
    case 2 : subterm =  RHS_U_V6; break;
    case 3 : subterm =  RHS_U_V8; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(subterm, r, mf_u, U, mf_obs, obs, pmf_lambda, lambda,
            pmf_coeff, &f_coeff, WT, alpha);

    getfem::generic_assembly assem;
    if (pmf_coeff)
      assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1))(i,:,i); ");
    else if (pmf_lambda)
      assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); ");
    else
      assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2).vBase(#1))(i,:,i); ");
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
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

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
      const mesh_fem &mf_obs = *(md.pmesh_fem_of_variable(dl[nd]));
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
             pmf_coeff, *f_coeff, WT, vr[0], alpha, rg, option);
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
             pmf_coeff, *f_coeff, WT, vr[0], alpha, rg, option);
      }

    }

    penalized_contact_rigid_obstacle_brick(bool contact_only_, int option_) {
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

  size_type add_penalized_contact_with_friction_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &dataname_obs, const std::string &dataname_r,
   const std::string &dataname_friction_coeff,
   size_type region, int option, const std::string &dataname_lambda,
   const std::string &dataname_alpha, const std::string &dataname_wt) {

    pbrick pbr = new penalized_contact_rigid_obstacle_brick(false, option);

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, true));

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
  //  Continuous contact with friction between non-matching meshes.
  //
  //=========================================================================

  template<typename MAT, typename VECT1>
  void asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix // frictionless
  (MAT &Ku1l, MAT &Klu1, MAT &Ku2l, MAT &Klu2, MAT &Kll, MAT &Ku1u1, MAT &Ku2u2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
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
        "M$7(#2,#2)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#2).vBase(#2))(i,j,:,i,:,j)"); // U2U2
      break;
    case 3:
      assem.set
      ("M$1(#1,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1).Base(#3))(i,:,i,:); "      // U1L
       "M$2(#3,#1)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#1))(i,:,:,i); "      // LU1
       "M$3(#2,#3)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#2).Base(#3))(i,:,i,:); "      // U2L
       "M$4(#3,#2)+=comp(NonLin$2(#1,#1,#2,#3).Base(#3).vBase(#2))(i,:,:,i); "      // LU2
       "M$5(#3,#3)+=comp(NonLin$3(#1,#1,#2,#3).Base(#3).Base(#3))(i,:,:); "         // LL
       "M$6(#1,#1)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#1).vBase(#1))(i,j,:,i,:,j)"   // U1U1
       "M$7(#2,#2)+=comp(NonLin$4(#1,#1,#2,#3).vBase(#2).vBase(#2))(i,j,:,i,:,j)"); // U2U2
      break;
    }
    assem.push_mi(mim);
    assem.push_mf(mf_u1);
    assem.push_mf(mf_u2);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm1);
    assem.push_nonlinear_term(&nterm2);
    assem.push_nonlinear_term(&nterm3);
    if (option == 2 || option == 3)
      assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Ku1l);
    assem.push_mat(Klu1);
    assem.push_mat(Ku2l);
    assem.push_mat(Klu2);
    assem.push_mat(Kll);
    if (option == 2 || option == 3) {
      assem.push_mat(Ku1u1);
      assem.push_mat(Ku2u2);
    }
    assem.assembly(rg);

    gmm::scale(Ku2l, scalar_type(-1));
    if (option != 2) // Klu2 was calculated
      gmm::scale(Klu2, scalar_type(-1));
  }

  template<typename MAT, typename VECT1>
  void asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix // with friction
  (MAT &Ku1l, MAT &Klu1, MAT &Ku2l, MAT &Klu2, MAT &Kll, MAT &Ku1u1, MAT &Ku2u2,
   const mesh_im &mim,
   const getfem::mesh_fem &mf_u1, const VECT1 &U1,
   const getfem::mesh_fem &mf_u2, const VECT1 &U2,
   const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   const VECT1 *WT1, const VECT1 *WT2,
   scalar_type r, scalar_type alpha,
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
             pmf_coeff, &f_coeff, WT1, WT2, alpha),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT1, WT2, alpha),
      nterm3(subterm3, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT1, WT2, alpha),
      nterm4(subterm4, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT1, WT2, alpha);

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
        "M$7(#2,#2)+=comp(NonLin$4(#1," + aux_fems + ").vBase(#2).vBase(#2))(i,j,:,i,:,j)"); // U2U2
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
    if (option == 2 || option == 3)
      assem.push_nonlinear_term(&nterm4);
    assem.push_mat(Ku1l);
    assem.push_mat(Klu1);
    assem.push_mat(Ku2l);
    assem.push_mat(Klu2);
    assem.push_mat(Kll);
    if (option == 2 || option == 3) {
      assem.push_mat(Ku1u1);
      assem.push_mat(Ku2u2);
    }
    assem.assembly(rg);

    gmm::scale(Ku2l, scalar_type(-1));
    gmm::scale(Klu2, scalar_type(-1));
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
   const getfem::mesh_fem *pmf_coeff, const VECT1 &f_coeff,
   const VECT1 *WT1, const VECT1 *WT2,
   scalar_type r, scalar_type alpha,
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
             pmf_coeff, &f_coeff, WT1, WT2, alpha),
      nterm2(subterm2, r, mf_u1, U1, mf_u2, U2, &mf_lambda, &lambda,
             pmf_coeff, &f_coeff, WT1, WT2, alpha);

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

  struct continuous_contact_nonmatching_meshes_brick : public virtual_brick {

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
                  "Continuous contact between nonmatching meshes bricks need a single mesh_im");
      const mesh_im &mim = *mims[0];

      // Variables : u1, u2, lambda
      //       the variable lambda should be scalar in the frictionless
      //       case and vector valued in the case with friction.
      GMM_ASSERT1(vl.size() == 3,
                  "Continuous contact between nonmatching meshes bricks need three variables");
      const model_real_plain_vector &u1 = md.real_variable(vl[0]);
      const model_real_plain_vector &u2 = md.real_variable(vl[1]);
      const mesh_fem &mf_u1 = *(md.pmesh_fem_of_variable(vl[0]));
      const mesh_fem &mf_u2 = *(md.pmesh_fem_of_variable(vl[1]));
      const model_real_plain_vector &lambda = md.real_variable(vl[2]);
      const mesh_fem &mf_lambda = *(md.pmesh_fem_of_variable(vl[2]));
      GMM_ASSERT1(mf_lambda.get_qdim() == (contact_only ? 1 : mf_u1.get_qdim()),
                  "The contact stress variable has not the right dimension");

      // Data : r, [friction_coeff,] [alpha,] [WT1, WT2]
      //     alpha, WT1, WT2 are optional and equal to 1, 0 and 0 by default respectively.
      if (contact_only) {
        GMM_ASSERT1(dl.size() == 1,
                    "Wrong number of data for continuous contact between nonmatching meshes "
                    << "brick, the number of data should be equal to 1 .");
      }
      else {
        GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 5,
                    "Wrong number of data for continuous contact between nonmatching meshes "
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

      // Matrix terms (T_u1l, T_lu1, T_u2l, T_lu2, T_ll, T_u1u1, T_u2u2) // FIXME: T_u1u2, T_u2u1 ???
      GMM_ASSERT1(matl.size() == size_type(5 + 2*(option == 3 ||
                                                  (option == 2 && !contact_only))),
                  "Wrong number of terms for "
                  "continuous contact between nonmatching meshes brick");

      mesh_region rg(region);
      mf_u1.linked_mesh().intersect_with_mpi_region(rg);

      size_type N = mf_u1.linked_mesh().dim();

      // projection of the second mesh_fem onto the mesh of the first mesh_fem
      if (!pmf_u2_proj) {
        pmf_u2_proj = new getfem::mesh_fem(mim.linked_mesh(), dim_type(N));
        pfem_proj = new_projected_fem(mf_u2, mim, rg2, rg1);
        pmf_u2_proj->set_finite_element(mim.linked_mesh().convex_index(), pfem_proj);
      }

      size_type nbdof_lambda = mf_lambda.nb_dof();
      size_type nbdof2 = mf_u2.nb_dof();
      size_type nbsub = pmf_u2_proj->nb_basic_dof();

      std::vector<size_type> ind;
      pmf_u2_proj->get_global_dof_index(ind);
      gmm::unsorted_sub_index SUBI(ind);

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(gmm::sub_matrix(mf_u2.extension_matrix(),
                                  SUBI, gmm::sub_interval(0, nbdof2)),
                  u2,
                  u2_proj);
      else
        gmm::copy(gmm::sub_vector(u2, SUBI), u2_proj);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Continuous contact between nonmatching meshes "
                   "tangent term");
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(matl[i]);

        model_real_sparse_matrix Ku2l(nbsub, nbdof_lambda);
        model_real_sparse_matrix Klu2(nbdof_lambda, nbsub);
        model_real_sparse_matrix Ku2u2(nbsub, nbsub);

        size_type sixthmat = (matl.size() >= 6) ? 5 : 1;
        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[0] /* u1l */, matl[1] /* lu1 */, Ku2l, Klu2, matl[4] /* ll */,
             matl[sixthmat] /* u1u1 */, Ku2u2,
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_tangent_matrix
            (matl[0] /* u1l */, matl[1] /* lu1 */, Ku2l, Klu2, matl[4] /* ll */,
             matl[sixthmat] /* u1u1 */, Ku2u2,
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, *f_coeff, WT1, WT2,
             vr[0], alpha, rg, option);

        if (mf_u2.is_reduced()) {
          gmm::mult(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    Ku2l,
                    matl[2]);
          gmm::mult(Klu2,
                    gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    matl[3]);
          if (matl.size() > 6) {
            model_real_sparse_matrix tmp(nbsub, nbdof2);
            gmm::mult(Ku2u2,
                      gmm::sub_matrix(mf_u2.extension_matrix(),
                                      SUBI, gmm::sub_interval(0, nbdof2)),
                      tmp);
            gmm::mult(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                      gmm::sub_interval(0, nbdof2), SUBI),
                      tmp,
                      matl[6]);
          }
        }
        else {
          gmm::copy(Ku2l, gmm::sub_matrix(matl[2], SUBI, gmm::sub_interval(0, nbdof_lambda)));
          gmm::copy(Klu2, gmm::sub_matrix(matl[3], gmm::sub_interval(0, nbdof_lambda), SUBI));
          if (matl.size() > 6)
            gmm::copy(Ku2u2, gmm::sub_matrix(matl[6], SUBI));
        }
      }

      if (version & model::BUILD_RHS) {
        for (size_type i = 0; i < matl.size(); i++) gmm::clear(vecl[i]);

        model_real_plain_vector Ru2(nbsub);

        if (contact_only)
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[0], Ru2, vecl[4], // u1, u2, lambda
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             vr[0], rg, option);
        else
          asm_Alart_Curnier_contact_nonmatching_meshes_rhs
            (vecl[0], Ru2, vecl[4], // u1, u2, lambda
             mim, mf_u1, u1, *pmf_u2_proj, u2_proj, mf_lambda, lambda,
             pmf_coeff, *f_coeff, WT1, WT2,
             vr[0], alpha, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    Ru2,
                    vecl[2]);
        else
          gmm::copy(Ru2, gmm::sub_vector(vecl[2], SUBI));
      }

    }

    continuous_contact_nonmatching_meshes_brick(size_type rg1_, size_type rg2_,
                                                bool contact_only_, int option_)
    : rg1(rg1_), rg2(rg2_), pfem_proj(0), pmf_u2_proj(0),
      contact_only(contact_only_), option(option_)
    {
      Tresca_version = false;   // for future version ...
      set_flags(contact_only
                ? "Continuous contact between nonmatching meshes brick"
                : "Continuous contact with friction between nonmatching "
                  "meshes brick",
                false /* is linear*/,
                (option==2) && contact_only /* is symmetric */,
                false /* is coercive */, true /* is real */,
                false /* is complex */);
    }

    ~continuous_contact_nonmatching_meshes_brick()
    { if (pmf_u2_proj) delete pmf_u2_proj; }

  };


  //=========================================================================
  //  Add a frictionless contact condition between two bodies discretized
  //  with nonmatching meshes.
  //=========================================================================

  size_type add_continuous_contact_between_nonmatching_meshes_brick
  (model &md, const mesh_im &mim, const std::string &varname_u1,
   const std::string &varname_u2, const std::string &multname_n,
   const std::string &dataname_r,
   size_type region1, size_type region2, int option) {

    pbrick pbr = new continuous_contact_nonmatching_meshes_brick
                     (region1, region2, true /* contact_only */, option);

    model::termlist tl;

    switch (option) {
    case 1 : case 4 :
      tl.push_back(model::term_description(varname_u1, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u1, false));
      tl.push_back(model::term_description(varname_u2, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u2, false));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      break;
    case 2 :
      tl.push_back(model::term_description(varname_u1, multname_n, true));
      tl.push_back(model::term_description(varname_u2, multname_n, true));
      tl.push_back(model::term_description(varname_u1, varname_u1, true));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      break;
    case 3 :
      tl.push_back(model::term_description(varname_u1, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u1, false));
      tl.push_back(model::term_description(varname_u2, multname_n, false));
      tl.push_back(model::term_description(multname_n, varname_u2, false));
      tl.push_back(model::term_description(multname_n, multname_n, true));
      tl.push_back(model::term_description(varname_u1, varname_u1, true));
      tl.push_back(model::term_description(varname_u2, varname_u2, true));
      break;
    default : GMM_ASSERT1(false,
                          "Incorrect option for continuous contact brick");

    }
    model::varnamelist dl(1, dataname_r);

    model::varnamelist vl(1, varname_u1);
    vl.push_back(varname_u2);
    vl.push_back(multname_n);

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region1);
  }


  //=========================================================================
  //
  //  Continuous penalized contact with friction between non-matching meshes.
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

    getfem::generic_assembly assem;
    if (pmf_lambda)
      assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#1))(i,:,i); "
                "V$2(#2)+=comp(NonLin$1(#1,#1,#2,#3).vBase(#2))(i,:,i)");
    else
      assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2).vBase(#1))(i,:,i); "
                "V$2(#2)+=comp(NonLin$1(#1,#1,#2).vBase(#2))(i,:,i)");
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
      const mesh_fem &mf_u1 = *(md.pmesh_fem_of_variable(vl[0]));
      const mesh_fem &mf_u2 = *(md.pmesh_fem_of_variable(vl[1]));

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

        if (dl.size() > nd) {
          nd++;
          alpha = md.real_variable(dl[nd])[0];
          GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[nd])) == 1,
                      "Parameter alpha should be a scalar");
        }

        if (dl.size() > nd) {
          nd++;
          WT = &(md.real_variable(dl[nd]));
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

      model_real_plain_vector u2_proj(nbsub);
      if (mf_u2.is_reduced())
        gmm::mult(gmm::sub_matrix(mf_u2.extension_matrix(),
                                  SUBI, gmm::sub_interval(0, nbdof2)),
                  u2,
                  u2_proj);
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
//        else {
//          gmm::clear(matl[3]);
//          model_real_sparse_matrix Ku2u1(nbsub,nbdof1);
//          asm_penalized_contact_nonmatching_meshes_tangent_matrix
//            (matl[0], Ku2u2, Ku1u2, Ku2u1, mim, mf_u1, u1, *pmf_u2_proj, u2_proj,
//             pmf_lambda, lambda, vr[0], alpha, mf_coeff, friction_coeff, WT, rg, option);
//          gmm::copy(Ku2u1, gmm::sub_matrix(matl[3], SUBI, gmm::sub_interval(0, nbdof1)));
//        }

        if (mf_u2.is_reduced()) {
          model_real_sparse_matrix tmp(nbsub, nbdof2);
          gmm::mult(Ku2u2,
                    gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    tmp);
          gmm::mult(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    tmp,
                    matl[1]);
          gmm::mult(Ku1u2,
                    gmm::sub_matrix(mf_u2.extension_matrix(),
                                    SUBI, gmm::sub_interval(0, nbdof2)),
                    matl[2]);
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
//        else
//          asm_penalized_contact_nonmatching_meshes_rhs
//            (vecl[0], mim, mf_u, u, mf_lambda, lambda, mf_obstacle, obstacle,
//             vr[0], alpha, mf_coeff, friction_coeff, WT, rg, option);

        if (mf_u2.is_reduced())
          gmm::mult(gmm::sub_matrix(mf_u2.reduction_matrix(),
                                    gmm::sub_interval(0, nbdof2), SUBI),
                    Ru2,
                    vecl[1]);
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
                ? "Continuous penalized contact between nonmatching meshes brick"
                : "Continuous penalized contact with friction between nonmatching "
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
                     (region1, region2, /* contact_only */ true, option);
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

}  /* end of namespace getfem.                                             */
