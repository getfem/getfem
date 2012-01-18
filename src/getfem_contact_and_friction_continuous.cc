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

    if (xn >= scalar_type(0) && f * nxt <= xn) {
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
      t[0] = r*(un-g) + gmm::pos(ln); break;

    case K_LL_V1:
      t[0] = (Heav(r*(un-g)-ln) - scalar_type(1))/r; break;
    case K_LL_V2:
      t[0] = -Heav(ln); break;

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
      for (i=0; i<N; ++i) t[i] = no[i]*e + zt[i] + lt[i] - auxN[i];
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
      for (i=0; i<N; ++i) t[i] = -auxN[i];
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
        for (i=0; i<N; ++i) for (j=0; j<N; ++j) t[i*N+j] = -r * GP(i,j);
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
        t[i*N+j] = no[i]*no[j]*e - ((i == j) ? scalar_type(1) : scalar_type(0)) + GP(i,j) - f_coeff*Heav(-ln)*no[i]*V[j];
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
        t[i*N+j] = GP(i,j) - ((i == j) ? scalar_type(1) : scalar_type(0));
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
        } else {
          gmm::clear(auxN);
        }
        if (gmm::vect_size(VT) == gmm::vect_size(U)) {
          gmm::copy(gmm::sub_vector
                    (VT, gmm::sub_index
                     (mf_u.ind_basic_dof_of_element(cv))), coeff);
          ctx.pf()->interpolation(ctx, coeff, vt, N);
          vt -= gmm::vect_sp(vt, no) * no;
        } else {
          gmm::clear(vt);
        }
	// zt = r*(alpha*(u_T-w_T) + (1-gamma)*v_T)
        zt = (((V - un * no) - auxN) * alpha + vt * (1 - gamma)) * r;
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
    // - this method expexts that i1,i2,i3 and i4 mesh_fems will correspond
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
             // calculate [ln] and [lnt] 
             // calculate [un] and [zt] interpolating [U2],[WT2] on [mf_u2]
      {
        const projected_fem &pfe = dynamic_cast<const projected_fem&>(*ctx.pf());
        pfe.projection_data(ctx, no, g);
        gmm::scale(no, scalar_type(-1)); // pointing outwards from mf_u1
      }

      if (!contact_only) {
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
      coeff.resize(mf_lambda.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (lambda, gmm::sub_index
                 (mf_lambda.ind_basic_dof_of_element(cv))), coeff);
      if (contact_only) {
        ctx.pf()->interpolation(ctx, coeff, aux1, 1);
        ln = aux1[0];
      } else {
        ctx.pf()->interpolation(ctx, coeff, lnt, N);
      }
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


  //=========================================================================
  //
  //  Continuous augmented Lagrangian brick (given obstacle, u, lambda).
  //
  //=========================================================================

  template<typename MAT, typename VECT1>
  void asm_continuous_contact_tangent_matrix_Alart_Curnier
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu, const mesh_im &mim,
   const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r, int option,
   const mesh_region &rg) {

    size_type subterm1 = (option == 4) ? K_UL_V2 : K_UL_V1;
    size_type subterm2 = (option == 4) ? K_UL_V4 : K_UL_V3;
    size_type subterm3 = (option == 4) ? K_LL_V2 : K_LL_V1;
    size_type subterm4 = (option == 2) ? K_UU_V2 : K_UU_V1;

    contact_rigid_obstacle_nonlinear_term
      nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm1);
    contact_rigid_obstacle_nonlinear_term
      nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm2);
    contact_rigid_obstacle_nonlinear_term
      nterm3(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm3);
    contact_rigid_obstacle_nonlinear_term
      nterm4(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm4);

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
  void asm_continuous_contact_with_friction_tangent_matrix_Alart_Curnier
  (MAT &Kul, MAT &Klu, MAT &Kll, MAT &Kuu, const mesh_im &mim,
   const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r,
   scalar_type alpha, scalar_type gamma, const getfem::mesh_fem *mf_coeff,
   const VECT1 &f_coeff, const VECT1 &WT, const VECT1 &VT, int option,
   const mesh_region &rg) {

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
      nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm1,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);
    contact_rigid_obstacle_nonlinear_term
      nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm2,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);
    contact_rigid_obstacle_nonlinear_term
      nterm3(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm3,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);
    contact_rigid_obstacle_nonlinear_term
      nterm4(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm4,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);

    getfem::generic_assembly assem;
    switch (option) {
    case 1: case 4: case 5: case 6:
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
  void asm_continuous_contact_rhs_Alart_Curnier
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

    contact_rigid_obstacle_nonlinear_term
      nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm1);
    contact_rigid_obstacle_nonlinear_term
      nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm2);

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
  void asm_continuous_contact_with_friction_rhs_Alart_Curnier
  (VECT1 &Ru, VECT1 &Rl, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda, const VECT1 &lambda,
   const getfem::mesh_fem &mf_obs, const VECT1 &obs, scalar_type r,
   scalar_type alpha, scalar_type gamma, const getfem::mesh_fem *mf_coeff,
   const VECT1 &f_coeff, const VECT1 &WT, const VECT1 &VT, int option,
   const mesh_region &rg) {

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
      nterm1(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm1,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);
    contact_rigid_obstacle_nonlinear_term
      nterm2(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm2,
             false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);

    getfem::generic_assembly assem;
    assem.set("V$1(#1)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1))(i,:,i); "
              "V$2(#2)+=comp(NonLin$2(#1,#1,#2,#3,#4).vBase(#2))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
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
    // option = 5 : De-Saxcé
    // option = 6 : New method based on De-Saxcé.

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
      GMM_ASSERT1(dl.size() >= 2 && dl.size() <= 7,
                  "Wrong number of data for continuous Coulomb friction "
                  << "brick, " << dl.size() << " should be between 2 and 7.");
      GMM_ASSERT1(matl.size() == size_type(3 + (option == 3)
                                           + (option == 2 && !contact_only)),
                  "Wrong number of terms for "
                  "continuous Coulomb friction brick");

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

      scalar_type gamma = 1;
      if (!contact_only && dl.size() >= 6) {
        gamma = md.real_variable(dl[5])[0];
        GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[5])) == 1,
                    "Parameter gamma should be a scalar");
      }

      const model_real_plain_vector &VT
        = (!contact_only && dl.size()>=7) ? md.real_variable(dl[6]) : voidvec;

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Continuous Coulomb friction tangent term");
        gmm::clear(matl[0]); gmm::clear(matl[1]); gmm::clear(matl[2]);
        if (matl.size() >= 4) gmm::clear(matl[3]);
        if (contact_only) {
          size_type fourthmat = (matl.size() >= 4) ? 3 : 1;
          asm_continuous_contact_tangent_matrix_Alart_Curnier
            (matl[0], matl[1], matl[2], matl[fourthmat], mim, mf_u, u,
             mf_lambda, lambda, mf_obstacle, obstacle, vr[0], option, rg);
        }
        else {
          size_type fourthmat = (matl.size() >= 4) ? 3 : 1;
          asm_continuous_contact_with_friction_tangent_matrix_Alart_Curnier
            (matl[0], matl[1], matl[2], matl[fourthmat], mim, mf_u, u,
             mf_lambda, lambda, mf_obstacle, obstacle, vr[0], alpha, gamma,
             mf_coeff, friction_coeff, WT, VT, option, rg);
        }
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]); gmm::clear(vecl[1]); gmm::clear(vecl[2]);
        if (matl.size() >= 4) gmm::clear(vecl[3]);

        if (contact_only)
          asm_continuous_contact_rhs_Alart_Curnier
            (vecl[0], vecl[2], mim, mf_u, u, mf_lambda, lambda,
             mf_obstacle, obstacle, vr[0], option, rg);
        else
          asm_continuous_contact_with_friction_rhs_Alart_Curnier
            (vecl[0], vecl[2], mim, mf_u, u, mf_lambda, lambda,
             mf_obstacle, obstacle, vr[0], alpha, gamma, mf_coeff,
             friction_coeff, WT, VT, option, rg);
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
   const std::string &dataname_alpha, const std::string &dataname_wt,
   const std::string &dataname_gamma, const std::string &dataname_vt) {

    pbrick pbr
      = new continuous_contact_rigid_obstacle_brick(false, option);

    model::termlist tl;

    switch (option) {
    case 1: case 4: case 5: case 6:
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
  void asm_penalized_contact_tangent_matrix
  (MAT &Kuu, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, int option, const mesh_region &rg) {

    contact_rigid_obstacle_nonlinear_term
      nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r,
            (option == 1) ? K_UU_V1 : K_UU_V2);

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
  void asm_penalized_contact_rhs
  (VECT1 &Ru, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U,  const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, int option, const mesh_region &rg) {

    contact_rigid_obstacle_nonlinear_term
      nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r,
            (option == 1) ? RHS_U_V5 : RHS_U_V2);

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
  void asm_penalized_contact_with_friction_tangent_matrix
  (MAT &Kuu, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U, const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, scalar_type alpha, scalar_type gamma,
   const getfem::mesh_fem *mf_coeff, const VECT1 &f_coeff, const VECT1 &WT,
   const VECT1 &VT, int option, const mesh_region &rg) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  K_UU_FRICT_V4; break;
    case 2 : subterm =  K_UU_FRICT_V3; break;
    case 3 : subterm =  K_UU_FRICT_V5; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm,
            false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);

    getfem::generic_assembly assem;
    assem.set
    ("M(#1,#1)+=comp(NonLin(#1,#1,#2,#3,#4).vBase(#1).vBase(#1))(i,j,:,i,:,j)");
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
  void asm_penalized_contact_with_friction_rhs
  (VECT1 &Ru, const mesh_im &mim, const getfem::mesh_fem &mf_u,
   const VECT1 &U,  const getfem::mesh_fem &mf_lambda,
   const VECT1 &lambda, const getfem::mesh_fem &mf_obs, const VECT1 &obs,
   scalar_type r, scalar_type alpha, scalar_type gamma,
   const getfem::mesh_fem *mf_coeff, const VECT1 &f_coeff, const VECT1 &WT,
   const VECT1 &VT, int option, const mesh_region &rg) {

    size_type subterm = 0;
    switch (option) {
    case 1 : subterm =  RHS_U_V7; break;
    case 2 : subterm =  RHS_U_V6; break;
    case 3 : subterm =  RHS_U_V8; break;
    }

    contact_rigid_obstacle_nonlinear_term
      nterm(mf_u, U, mf_lambda, lambda, mf_obs, obs, r, subterm,
            false, alpha, mf_coeff, &f_coeff, &WT, gamma, &VT);

    getfem::generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin$1(#1,#1,#2,#3,#4).vBase(#1))(i,:,i); ");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_lambda);
    assem.push_mf(mf_obs);
    assem.push_mf(mf_coeff ? *mf_coeff : mf_obs);
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

      // declare gamma and VT with implicit values to be consistent with the
      // augmented Lagrangian brick
      scalar_type gamma = 1;
      const model_real_plain_vector &VT = voidvec;

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        GMM_TRACE2("Penalized Coulomb friction tangent term");
        gmm::clear(matl[0]);
        if (contact_only)
          asm_penalized_contact_tangent_matrix
            (matl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
             vr[0], option, rg);
        else
          asm_penalized_contact_with_friction_tangent_matrix
            (matl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle,
             obstacle, vr[0], alpha, gamma, mf_coeff, friction_coeff, WT, VT,
             option, rg);
      }

      if (version & model::BUILD_RHS) {
        gmm::clear(vecl[0]);
        if (contact_only)
          asm_penalized_contact_rhs
            (vecl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle, obstacle,
             vr[0], option, rg);
        else
          asm_penalized_contact_with_friction_rhs
            (vecl[0], mim, mf_u, u, *mf_lambda, lambda, mf_obstacle,
             obstacle, vr[0], alpha, gamma, mf_coeff, friction_coeff, WT, VT,
             option, rg);
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


}  /* end of namespace getfem.                                             */
