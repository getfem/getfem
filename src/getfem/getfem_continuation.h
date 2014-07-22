/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2011-2014 Tomas Ligursky, Yves Renard
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/** @file getfem_continuation.h
    @author Tomas Ligursky <tomas.ligursky@gmail.com>
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date October 17, 2011.
    @brief inexact Moore-Penrose continuation method.
*/
#ifndef GETFEM_CONTINUATION_H__
#define GETFEM_CONTINUATION_H__

#include <getfem/getfem_model_solvers.h>

namespace getfem {


  //=========================================================================
  // Abstract Moore-Penrose continuation method
  //=========================================================================

  const double tau_bp_init = 1.e6;
  enum build_data { BUILD_F = 1, BUILD_F_x = 2, BUILD_ALL = 3 };

  /* Compute a unit tangent at (x, gamma) that is accute to the incoming
     vector. */
  template <typename CONT_S, typename VECT>
  void compute_tangent(CONT_S &S, const VECT &x, double gamma,
		       VECT &t_x, double &t_gamma) {
    double r;
    VECT g(x), y(x);
    S.F_gamma(x, gamma, g);
    S.solve_grad(x, gamma, y, g);
    t_gamma = 1. / (t_gamma - S.w_sp(t_x, y));
    S.scale(y, -t_gamma); S.copy(y, t_x);
    
    double no = S.w_norm(t_x, t_gamma);
    S.scale(t_x, 1./no); t_gamma /= no;

    S.mult_grad(x, gamma, t_x, y); S.scaled_add(y, g, t_gamma, y);
    r = S.norm(y);
    if (r > 1.e-10)
      GMM_WARNING1("Tangent computed with the residual " << r);
  }

  /* Calculate a tangent vector at (x, gamma) + h * (T_x, T_gamma) and test
     whether it is close to (T_x, T_gamma). Informatively, compare it with
     (t_x, t_gamma), as well. */
  template <typename CONT_S, typename VECT>
  bool test_tangent(CONT_S &S, const VECT &x, double gamma,
		    const VECT &T_x, double T_gamma,
		    const VECT &t_x, double t_gamma, double h) {
    bool res = false;
    double Gamma, T_Gamma = T_gamma, cang;
    VECT X(x), T_X(T_x);
    
    S.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
    S.set_build(BUILD_ALL);
    compute_tangent(S, X, Gamma, T_X, T_Gamma);
    
    cang = S.cosang(T_X, T_x, T_Gamma, T_gamma);
    if (S.noisy() > 1)
      cout << "cos of the angle with the tested tangent " << cang << endl;
    if (cang >= S.mincos()) res = true;
    else {    
      cang = S.cosang(T_X, t_x, T_Gamma, t_gamma);
      if (S.noisy() > 1)
	cout << "cos of the angle with the initial tangent " << cang << endl;
    }
    return res;
  }

  /* Simple tangent switch. */
  template <typename CONT_S, typename VECT>
  bool switch_tangent(CONT_S &S, const VECT &x, double gamma,
		      VECT &t_x, double &t_gamma, double &h) {    
    bool accepted;
    double t_gamma0 = t_gamma, T_gamma = t_gamma, Gamma;
    VECT t_x0(t_x), T_x(t_x), X(x);

    if (S.noisy() > 0) cout  << "trying simple tangent switch" << endl;
    if (S.noisy() > 0) cout << "starting computing a new tangent" << endl;
    h *= 1.5;
    S.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
    S.set_build(BUILD_ALL);
    compute_tangent(S, X, Gamma, T_x, T_gamma);
    // One can test the cosine of the angle between (T_x, T_gamma) and
    // (t_x, t_gamma), for sure, and increase h_min if it were greater or
    // equal to S.mincos(). However, this seems to be superfluous.
    
    if (S.noisy() > 0)
      cout << "starting testing the computed tangent" << endl;
    double h_test = (-0.9) * S.h_min();
    do {
      h_test = -h_test
	+ pow(10., floor(log10(-h_test / S.h_min()))) * S.h_min();
      accepted = test_tangent(S, x, gamma, T_x, T_gamma,
			      t_x, t_gamma, h_test);
      if (!accepted) {
	h_test *= -1.;
	accepted = test_tangent(S, x, gamma, T_x, T_gamma,
				t_x, t_gamma, h_test);
      }
    } while (!accepted && (h_test > -S.h_max()));
    
    if (accepted) {
      S.copy(T_x, t_x); t_gamma = T_gamma;
      if (h_test < 0) { 
	S.scale(t_x, -1.); t_gamma *= -1.; h_test *= -1.;
      }
      if (S.noisy() > 0)
	cout << "tangent direction switched, "
	     << "starting computing a suitable step size" << endl;
      bool h_adapted = false; h = S.h_init();
      while (!h_adapted && (h > h_test)) {
	h_adapted = test_tangent(S, x, gamma, t_x, t_gamma,
				 t_x0, t_gamma0, h);
	h *= S.h_dec();
      }
      h = (h_adapted) ? h / S.h_dec() : h_test;
    } else
      if (S.noisy() > 0) cout << "simple tangent switch has failed" << endl;
    
    return accepted;
  }

  /* Test for limit points (also called folds or turning points). */
  template <typename CONT_S>
  bool test_limit_point(CONT_S &S, double t_gamma) {
    double tau1 = S.get_tau_lp(), tau2 = t_gamma;
    S.set_tau_lp(tau2);
    return (tau2 * tau1 < 0);
  }


  /* Test function for bifurcation points for a given matrix. The first part
     of the solution of the augmented system is passed in
     (v_x, v_gamma). */
  template <typename CONT_S, typename MAT, typename VECT>
  double test_function_bp(CONT_S &S, const MAT &A, const VECT &g,
			  const VECT &t_x, double t_gamma,
			  VECT &v_x, double &v_gamma) {
    double q, r, tau;
    VECT y(g), z(g);

    S.solve(A, y, z, g, S.b_x());
    v_gamma = (S.b_gamma() - S.sp(t_x, z)) / (t_gamma - S.sp(t_x, y));
    S.scaled_add(z, y, -v_gamma, v_x);
    tau = 1. / (S.d() - S.sp(S.c_x(), v_x) - S.c_gamma() * v_gamma);
    S.scale(v_x, -tau); v_gamma *= -tau;

    // control of the norm of the residual
    S.mult(A, v_x, y);
    S.scaled_add(y, g, v_gamma, y); S.scaled_add(y, S.b_x(), tau, y); 
    r = S.sp(y, y);
    q = S.sp(t_x, v_x) + t_gamma * v_gamma + S.b_gamma() * tau; r += q * q;
    q = S.sp(S.c_x(), v_x) + S.c_gamma() * v_gamma + S.d() * tau - 1.;
    r += q * q; r = sqrt(r);
    if (r > 1.e-10)
      GMM_WARNING1("Test function evaluated with the residual " << r);

    return tau;
  }

  template <typename CONT_S, typename MAT, typename VECT>
  double test_function_bp(CONT_S &S, const MAT &A, const VECT &g,
			  const VECT &t_x, double t_gamma) {
    VECT v_x(g); double v_gamma;
    return test_function_bp(S, A, g, t_x, t_gamma, v_x, v_gamma);
  }

  /* Test function for bifurcation points for the gradient computed at
     (x, gamma). */
  template <typename CONT_S, typename VECT>
  double test_function_bp(CONT_S &S, const VECT &x, double gamma,
			  const VECT &t_x, double t_gamma,
			  VECT &v_x, double &v_gamma) {
    typename CONT_S::MAT A; S.F_x(x, gamma, A);
    VECT g(x); S.F_gamma(x, gamma, g);
    return test_function_bp(S, A, g, t_x, t_gamma, v_x, v_gamma);
  }

  template <typename CONT_S, typename VECT>
  double test_function_bp(CONT_S &S, const VECT &x, double gamma,
			  const VECT &t_x, double t_gamma) {
    VECT v_x(x); double v_gamma;
    return test_function_bp(S, x, gamma, t_x, t_gamma, v_x, v_gamma);
  }

  /* Test for smooth bifurcation points. */
  template <typename CONT_S, typename VECT>
  bool test_smooth_bifurcation(CONT_S &S, const VECT &x, double gamma,
			       const VECT &t_x, double t_gamma) {
    double tau0 = S.get_tau_bp_1(), tau1 = S.get_tau_bp_2(),
      tau2 = test_function_bp(S, x, gamma, t_x, t_gamma);
    S.set_tau_bp_1(tau1); S.set_tau_bp_2(tau2);
    return (tau2 * tau1 < 0) && (S.abs(tau1) < S.abs(tau0));
  }

  /* Test for non-smooth bifurcation points. */
  template <typename CONT_S, typename VECT>
  bool test_nonsmooth_bifurcation (CONT_S &S, const VECT &x1, double gamma1,
				   const VECT &t_x1, double t_gamma1,
				   const VECT &x2, double gamma2,
				   const VECT &t_x2, double t_gamma2) {
    unsigned long nb_changes = 0;
    double alpha = 0., delta = S.delta_min(),
      tau0 = tau_bp_init, tau1= S.get_tau_bp_2(), tau2, tau_var_ref, t_gamma;
    VECT g1(x1), g2(x1), g(x1), t_x(x1);

    // compute gradients at the two given points
    typename CONT_S::MAT A1, A2, A;
    S.F_x(x2, gamma2, A2); S.F_x(x2, gamma2, A); S.F_gamma(x2, gamma2, g2);
    S.F_x(x1, gamma1, A1); S.F_gamma(x1, gamma1, g1);
    S.init_tau_bp_graph();
    tau2 = test_function_bp(S, A2, g2, t_x2, t_gamma2);
    tau_var_ref = std::max(S.abs(tau2 - tau1),
			   (S.abs(tau1) + S.abs(tau2)) / 200);

    // monitor sign changes of the test function on the convex combination
    do {
      alpha = std::min(alpha + delta, 1.);
      S.scaled_add(A1, 1. - alpha, A2, alpha, A);
      S.scaled_add(g1, 1. - alpha, g2, alpha, g);
      S.scaled_add(t_x1, 1. - alpha, t_x2, alpha, t_x);
      t_gamma = (1. - alpha) * t_gamma1 + alpha * t_gamma2;
      
      tau2 = test_function_bp(S, A, g, t_x, t_gamma);
      if ((tau2 * tau1 < 0) && (S.abs(tau1) < S.abs(tau0))) ++nb_changes;
      S.insert_tau_bp_graph(alpha, tau2);

      if (S.abs(tau2 - tau1) < 0.5 * S.thrvar() * tau_var_ref)
	delta = std::min(2 * delta, S.delta_max());
      else if (S.abs(tau2 - tau1) > S.thrvar() * tau_var_ref) 
	delta = std::max(0.1 * delta, S.delta_min());
      tau0 = tau1; tau1 = tau2; 
    } while (alpha < 1.);
    
    S.set_tau_bp_1(tau_bp_init); S.set_tau_bp_2(tau2);
    return nb_changes % 2;
  }
  
  /* Newton-type corrections for the couple ((X, Gamma), (T_x, T_gamma)).
     The current direction of (T_x, T_gamma) is informatively compared with
     (t_x, t_gamma). */
  template <typename CONT_S, typename VECT>
  bool newton_corr(CONT_S &S, VECT &X, double &Gamma, VECT &T_x,
		   double &T_gamma, const VECT &t_x, double t_gamma,
		   unsigned long &it) {
    bool converged = false;
    double Delta_Gamma, no, res, diff;
    VECT F(X), g(X), Delta_X(X), y(X);

    if (S.noisy() > 0) cout << "starting correction " << endl;
    it = 0;
    S.F(X, Gamma, F);
    
    do {
      S.F_gamma(X, Gamma, g);
      S.solve_grad(X, Gamma, Delta_X, y, F, g);
      
      Delta_Gamma = S.sp(T_x, Delta_X) / (S.sp(T_x, y) - T_gamma);
      S.scaled_add(Delta_X, y, -Delta_Gamma, Delta_X);
      S.scaled_add(X, Delta_X, -1., X); Gamma -= Delta_Gamma;
      S.set_build(BUILD_ALL);
      
      T_gamma = 1. / (T_gamma - S.w_sp(T_x, y));
      S.scale(y, -T_gamma); S.copy(y, T_x);
      no = S.w_norm(T_x, T_gamma);
      S.scale(T_x, 1./no); T_gamma /= no;

      S.F(X, Gamma, F); res = S.norm(F); 
      diff = S.w_norm(Delta_X, Delta_Gamma);
      if (S.noisy() > 1)
	cout << " iter " << it << " residual " << res
	     << " difference " << diff 
	     << " cosang " << S.cosang(T_x, t_x, T_gamma, t_gamma) << endl;

      if (res <= S.maxres() && diff <= S.maxdiff()) {
	converged = true;
	// recalculate the final tangent, for sure
	compute_tangent(S, X, Gamma, T_x, T_gamma);
	break;
      }

      it++;      
    } while (it < S.maxit() && res < 1.e8);
    return converged;
  }
  
  template <typename CONT_S, typename VECT>
  bool newton_corr(CONT_S &S, VECT &X, double &Gamma, VECT &T_x,
		   double &T_gamma, const VECT &t_x, double t_gamma) {
    unsigned long it;
    return newton_corr(S, X, Gamma, T_x, T_gamma, t_x, t_gamma, it);
  }

  /* Try to perform one predictor-corrector step starting from the couple
     ((x, gamma), (t_x, t_gamma)). Return the resulting couple in the case of
     convergence. */
  template <typename CONT_S, typename VECT>
  bool test_predict_dir(CONT_S &S, VECT &x, double &gamma,
			VECT &t_x, double &t_gamma) {
    bool converged = false;
    double h =  S.h_init(), Gamma, T_gamma;
    VECT X(x), T_x(x);
    do { //step control
      
      // prediction
      if (S.noisy() > 0) cout << "prediction with h = " << h << endl;
      S.scaled_add(x, t_x, h, X); Gamma = gamma + h * t_gamma;
      S.set_build(BUILD_ALL);
      S.copy(t_x, T_x); T_gamma = t_gamma;

      //correction
      converged = newton_corr(S, X, Gamma, T_x, T_gamma, t_x, t_gamma);
      
      if (converged) {
	// check the direction of the tangent found
	S.scaled_add(X, x, -1., t_x); t_gamma = Gamma - gamma;
	if (S.sp(T_x, t_x, T_gamma, t_gamma) < 0)
	  { S.scale(T_x, -1.); T_gamma *= -1.; }
	S.copy(X, x); gamma = Gamma;
	S.copy(T_x, t_x); t_gamma = T_gamma;
      }
      else if (h > S.h_min())
	  h = (0.199 * S.h_dec() * h > S.h_min()) ?
	    0.199 * S.h_dec() * h : S.h_min();
      else break;

    } while(!converged);
    return converged;
  }

  /* A tool for approximating a smooth bifurcation point close to (x, gamma)
     and locating the two branches emanating from there. */
  template <typename CONT_S, typename VECT>
  void treat_smooth_bif_point(CONT_S &S, const VECT &x, double gamma,
			      const VECT &t_x, double t_gamma, double h) {
    unsigned long i = 0;
    double tau0 = S.get_tau_bp_1(), tau1 = S.get_tau_bp_2(),
      gamma0 = gamma, Gamma, t_gamma0 = t_gamma, T_gamma = t_gamma, v_gamma;
    VECT x0(x), X(x), t_x0(t_x), T_x(t_x), v_x(t_x);
    
    if (S.noisy() > 0)
      cout  << "starting locating the bifurcation point" << endl;

    // predictor-corrector steps with a secant-type step-length adaptation
    h *= tau1 / (tau0 - tau1);
    while ((S.abs(h) >= S.h_min()) && i < 10) {
      if (S.noisy() > 0) cout << "prediction with h = " << h << endl;
      S.scaled_add(x0, t_x0, h, X); Gamma = gamma0 + h * t_gamma0;
      S.set_build(BUILD_ALL);
      if (newton_corr(S, X, Gamma, T_x, T_gamma, t_x0, t_gamma0)) {
	S.copy(X, x0); gamma0 = Gamma;
	if (S.cosang(T_x, t_x0, T_gamma, t_gamma0) >= S.mincos())
	  { S.copy(T_x, t_x0); t_gamma0 = T_gamma; }
	tau0 = tau1;
	tau1 = test_function_bp(S, X, Gamma, t_x0, t_gamma0, v_x, v_gamma);
	h *= tau1 / (tau0 - tau1);
      }	else {
	S.scaled_add(x0, t_x0, h, x0); gamma0 += h * t_gamma0;
	test_function_bp(S, x0, gamma0, t_x0, t_gamma0, v_x, v_gamma);
	break;
      }
      ++i;
    }
    S.set_sing_point(x0, gamma0);
    S.insert_tangent_sing(t_x0, t_gamma0);

    if (S.noisy() > 0)
      cout  << "starting searching for the second branch" << endl;
    double no = S.w_norm(v_x, v_gamma);
    S.scale(v_x, 1./no); v_gamma /= no;
    if (test_predict_dir(S, x0, gamma0, v_x, v_gamma)
	&& S.insert_tangent_sing(v_x, v_gamma))
      { if (S.noisy() > 0) cout << "second branch found" << endl; }
    else if (S.noisy() > 0) cout << "Second branch not found!" << endl;
  }
  
  /* A tool for approximating a non-smooth point close to (x, gamma) and
     locating (preferably) all smooth one-sided solution branches emanating
     from there. It is supposed that (x, gamma) is a point on the previously
     traversed smooth solution branch within the distance of S.h_min() from
     the end point of this branch and (t_x, t_gamma) is the corresponding
     tangent that is directed towards the end point. */
  template <typename CONT_S, typename VECT>
  void treat_nonsmooth_point(CONT_S &S, const VECT &x, double gamma,
			     const VECT &t_x, double t_gamma, int version) {
    double gamma_end = gamma, Gamma, t_gamma0 = t_gamma, T_gamma = t_gamma,
      h = S.h_min(), cang, mcos = S.mincos();
    VECT x_end(x), X(x), t_x0(t_x), T_x(t_x);

    // approximate the non-smooth point by a bisection-like algorithm
    if (S.noisy() > 0)
      cout  << "starting locating a non-smooth point" << endl;
    S.scaled_add(x, t_x, h, X); Gamma = gamma + h * t_gamma;
    S.set_build(BUILD_ALL);
    if (newton_corr(S, X, Gamma, T_x, T_gamma, t_x0, t_gamma0)) {
      cang = S.cosang(T_x, t_x0, T_gamma, t_gamma0);
      if (cang >= mcos) mcos = (cang + 1.) / 2.;
    }

    S.copy(t_x0, T_x); T_gamma = t_gamma0;
    h /= 2.;
    for (unsigned long i = 0; i < 15; i++) {
      if (S.noisy() > 0) cout << "prediction with h = " << h << endl;
      S.scaled_add(x_end, t_x0, h, X); Gamma = gamma_end + h * t_gamma0;
      S.set_build(BUILD_ALL);
      if (newton_corr(S, X, Gamma, T_x, T_gamma, t_x0, t_gamma0)
	  && (S.cosang(T_x, t_x, T_gamma, t_gamma) >= mcos)) {
	S.copy(X, x_end); gamma_end = Gamma;
	S.copy(T_x, t_x0); t_gamma0 = T_gamma;
      } else {
	S.copy(t_x0, T_x); T_gamma = t_gamma0;
      }
      h /= 2.;
    }
    S.scaled_add(x_end, t_x0, h, x_end); gamma_end += h * t_gamma0;
    S.set_sing_point(x_end, gamma_end);

    // take two vectors to span a subspace of perturbations for the
    // non-smooth point
    if (S.noisy() > 0)
      cout  << "starting a thorough search for other branches" << endl;
    double t_gamma1 = t_gamma0, t_gamma2 = t_gamma0;
    VECT t_x1(t_x0), t_x2(t_x0);
    S.scale(t_x1, -1.); t_gamma1 *= -1.;
    S.insert_tangent_sing(t_x1, t_gamma1);

    h = S.h_min();
    S.scaled_add(x_end, t_x0, h, X); Gamma = gamma_end + h * t_gamma0;
    S.set_build(BUILD_ALL);
    compute_tangent(S, X, Gamma, t_x2, t_gamma2);

    // perturb the non-smooth point systematically to find new tangent
    // predictions
    unsigned long i1 = 0, i2 = 0, ncomb = 0;
    double a, a1, a2, no;
    S.clear(t_x0); t_gamma0 = 0.;

    do {
      for (unsigned long i = 0; i < S.nbdir(); i++) {
	a = (2 * M_PI * double(i)) / double(S.nbdir());
	a1 = h * sin(a); a2 = h * cos(a);
	S.scaled_add(x_end, t_x1, a1, X); Gamma = gamma_end + a1 * t_gamma1;
	S.scaled_add(X, t_x2, a2, X); Gamma += a2 * t_gamma2;
	S.set_build(BUILD_ALL);
	compute_tangent(S, X, Gamma, T_x, T_gamma);

	if (S.abs(S.cosang(T_x, t_x0, T_gamma, t_gamma0)) < S.mincos()) {
	  S.copy(T_x, t_x0); t_gamma0 = T_gamma;
	  if (S.insert_tangent_predict(T_x, T_gamma)) {
	    if (S.noisy() > 0)
	      cout << "new potential tangent vector found, "
		   << "trying one predictor-corrector step" << endl;
	    S.copy(x_end, X); Gamma = gamma_end;
	    
	    if (test_predict_dir(S, X, Gamma, T_x, T_gamma)) {
	      if (S.insert_tangent_sing(T_x, T_gamma)) {
		if ((a == 0) && (ncomb == 0)
		    && (S.abs(S.cosang(T_x, t_x0, T_gamma, t_gamma0))
			>= S.mincos())) { i2 = 1; ncomb = 1; }
		if (version) S.set_next_point(X, Gamma);
	      }
	      S.copy(x_end, X); Gamma = gamma_end;
	      S.copy(t_x0, T_x); T_gamma = t_gamma0;
	    }
	    
	    S.scale(T_x, -1.); T_gamma *= -1.;
	    if (test_predict_dir(S, X, Gamma, T_x, T_gamma)
		&& S.insert_tangent_sing(T_x, T_gamma) && version)
	      S.set_next_point(X, Gamma);
	  }
	}
      }
      
      // heuristics for varying the spanning vectors
      bool index_changed;
      if (i1 + 1 < i2) { ++i1; index_changed = true; }
      else if(i2 + 1 < S.nb_tangent_sing())
	{ ++i2; i1 = 0; index_changed = true; }
      else index_changed = false;
      if (index_changed) {
	S.copy(S.get_t_x_sing(i1), t_x1); t_gamma1 = S.get_t_gamma_sing(i1);
	S.copy(S.get_t_x_sing(i2), t_x2); t_gamma2 = S.get_t_gamma_sing(i2);
      } else {
	S.fill_random(T_x); T_gamma = S.random();
	no = S.w_norm(T_x, T_gamma);
	S.scaled_add(t_x2, T_x, 0.1/no, t_x2);
	t_gamma2 += 0.1/no * T_gamma;
	S.scaled_add(x_end, t_x2, h, X); Gamma = gamma_end + h * t_gamma2;
	S.set_build(BUILD_ALL);
	compute_tangent(S, X, Gamma, t_x2, t_gamma2);
      }
    } while (++ncomb < S.nbcomb());

    if (S.noisy() > 0)
      cout << "located branches " << S.nb_tangent_sing() << endl;
  }


  template <typename CONT_S, typename VECT>
  void init_test_functions(CONT_S &S, const VECT &x, double gamma,
			   const VECT &t_x, double t_gamma) {
    S.set_tau_lp(t_gamma);
    if (S.singularities() > 1) {
      if (S.noisy() > 0) cout << "starting computing an initial value of a "
			      << "test function for bifurcations" << endl;
      S.set_build(BUILD_ALL);
      double tau = test_function_bp(S, x, gamma, t_x, t_gamma);
      S.set_tau_bp_2(tau);
    }
  }

  template <typename CONT_S, typename VECT>
  void init_Moore_Penrose_continuation(CONT_S &S, const VECT &x,
				       double gamma, VECT &t_x,
				       double &t_gamma, double &h) {
    S.set_build(BUILD_ALL);
    S.clear(t_x); t_gamma = (t_gamma >= 0) ? 1. : -1.;
    if (S.noisy() > 0)
      cout << "starting computing an initial tangent" << endl;
    compute_tangent(S, x, gamma, t_x, t_gamma);
    h = S.h_init();
    if (S.singularities() > 0)
      init_test_functions(S, x, gamma, t_x, t_gamma);
  }

  
  /* Perform one step of the (non-smooth) Moore-Penrose continuation.
     NOTE: The new point need not to be saved in the model in the end! */
  template <typename CONT_S, typename VECT>
    void Moore_Penrose_continuation(CONT_S &S, VECT &x, double &gamma,
				    VECT &t_x, double &t_gamma, double &h) {
    bool converged, new_point = false, tangent_switched = false;
    unsigned long it, step_dec = 0;
    double t_gamma0 = t_gamma, Gamma, T_gamma;
    VECT t_x0(t_x), X(x), T_x(x);

    S.clear_tau_bp_currentstep(); S.clear_sing_data();

    do {
      // prediction
      if (S.noisy() > 0) cout << "prediction with h = " << h << endl;
      S.scaled_add(x, t_x, h, X); Gamma = gamma + h * t_gamma;
      S.set_build(BUILD_ALL);
      S.copy(t_x, T_x); T_gamma = t_gamma;
      
      // correction
      converged = newton_corr(S, X, Gamma, T_x, T_gamma, t_x, t_gamma, it);

      if (converged
	  && (S.cosang(T_x, t_x, T_gamma, t_gamma) >= S.mincos())) {
	new_point = true;
	if (S.singularities() > 0) {
	  if (test_limit_point(S, T_gamma)) {
	    S.set_sing_label("limit point");
	    if (S.noisy() > 0) cout << "Limit point detected!" << endl;
	  }
	  if (S.singularities() > 1) {
	    if (S.noisy() > 0)
	      cout << "new point found, starting computing a test function "
		   << "for bifurcations" << endl;
	    if (!tangent_switched) {
	      if(test_smooth_bifurcation(S, X, Gamma, T_x, T_gamma)) {
		S.set_sing_label("smooth bifurcation point");
		if (S.noisy() > 0)
		  cout << "Smooth bifurcation point detected!" << endl;
		treat_smooth_bif_point(S, X, Gamma, T_x, T_gamma, h);
	      }
	    } else if (test_nonsmooth_bifurcation(S, x, gamma, t_x0,
						  t_gamma0, X, Gamma, T_x,
						  T_gamma)) {
	      S.set_sing_label("non-smooth bifurcation point");
	      if (S.noisy() > 0)
		cout << "Non-smooth bifurcation point detected!" << endl;
	      treat_nonsmooth_point(S, x, gamma, t_x0, t_gamma0, 0);
	    }
	  }
	}
	
	if (step_dec == 0 && it < S.thrit())
	  h = (S.h_inc() * h < S.h_max()) ? S.h_inc() * h : S.h_max();
      } else if (h > S.h_min()) {
	h = (S.h_dec() * h > S.h_min()) ? S.h_dec() * h : S.h_min();
	step_dec++;
      } else if (S.non_smooth() && !tangent_switched) {
	if (S.noisy() > 0)
	  cout << "classical continuation has failed" << endl;
	if (switch_tangent(S, x, gamma, t_x, t_gamma, h)) {
	  tangent_switched = true;
	  step_dec = (h >= S.h_init()) ? 0 : 1;
	  if (S.noisy() > 0)
	    cout << "restarting the classical continuation" << endl;
	} else break;
      } else break;
    } while (!new_point);

    if (new_point) {
      S.copy(X, x); gamma = Gamma;
      S.copy(T_x, t_x); t_gamma = T_gamma;
    } else if (S.non_smooth()) {
      treat_nonsmooth_point(S, x, gamma, t_x0, t_gamma0, 1);
      if (S.next_point()) {
	if (S.singularities() > 0) {
	  if (test_limit_point(S, T_gamma)) {
	    S.set_sing_label("limit point");
	    if (S.noisy() > 0) cout << "Limit point detected!" << endl;
	  }
	  if (S.singularities() > 1) {
	    if (S.noisy() > 0)
	      cout << "starting computing a test function for bifurcations"
		   << endl;
	    S.set_build(BUILD_ALL);
	    bool bifurcation_detected = (S.nb_tangent_sing() > 2);
	    if (bifurcation_detected) {
	      // update the stored values of the test function only
	      S.set_tau_bp_1(tau_bp_init);
	      S.set_tau_bp_2(test_function_bp(S, S.get_x_next(),
					      S.get_gamma_next(),
					      S.get_t_x_sing(1),
					      S.get_t_gamma_sing(1)));
	    } else
	      bifurcation_detected
		= test_nonsmooth_bifurcation(S, x, gamma, t_x, t_gamma,
					     S.get_x_next(),
					     S.get_gamma_next(),
					     S.get_t_x_sing(1),
					     S.get_t_gamma_sing(1));
	    if (bifurcation_detected) {
	      S.set_sing_label("non-smooth bifurcation point");
	      if (S.noisy() > 0)
		cout << "Non-smooth bifurcation point detected!" << endl;
	    }
	  }
	}
	
	S.copy(S.get_x_next(), x); gamma = S.get_gamma_next();
	S.copy(S.get_t_x_sing(1), t_x); t_gamma = S.get_t_gamma_sing(1);
	h = S.h_init();
	new_point = true;
      }
    }
    
    if (!new_point) {
      cout << "Continuation has failed!" << endl;
      h = 0;
    }
  }
  

  //=========================================================================
  // Moore-Penrose continuation method for Getfem models
  //=========================================================================


#ifdef GETFEM_MODELS_H__
 
  struct cont_struct_getfem_model {
    
    typedef base_vector VECT;
    typedef model_real_sparse_matrix MAT;
    
  private:
    model *md;
    int singularities_;
    bool nonsmooth;
    std::string parameter_name_;
    bool with_parametrised_data;
    std::string initdata_name_, finaldata_name_, currentdata_name_;
    double scfac_;
    rmodel_plsolver_type lsolver;
    double h_init_, h_max_, h_min_, h_inc_, h_dec_;
    unsigned long maxit_, thrit_;
    double maxres_, maxdiff_, mincos_, maxres_solve_, delta_max_, delta_min_,
      thrvar_;
    unsigned long nbdir_, nbcomb_;
    int noisy_;
    VECT b_x_, c_x_;
    double b_gamma_, c_gamma_, d_;
    double tau_lp, tau_bp_1, tau_bp_2;
    VECT alpha_hist, tau_bp_hist;
    std::map<double, double> tau_bp_graph;
    std::string sing_label;
    VECT x_sing, x_next;
    double gamma_sing, gamma_next;
    std::vector<VECT> t_x_sing, t_x_predict;
    std::vector<double> t_gamma_sing, t_gamma_predict;
    build_data build;

  public:
    void init_border(void) {
      srand(unsigned(time(NULL)));
      unsigned long nbdof = md->nb_dof();
      gmm::resize(b_x_, nbdof); gmm::fill_random(b_x_);
      gmm::resize(c_x_, nbdof); gmm::fill_random(c_x_);
      b_gamma_ = gmm::random(1.); c_gamma_ = gmm::random(1.);
      d_ = gmm::random(1.);
    }

    cont_struct_getfem_model
    (model &m, const std::string &pn, double sfac, rmodel_plsolver_type ls,
     double hin = 1.e-2, double hmax = 1.e-1, double hmin = 1.e-5,
     double hinc = 1.3, double hdec = 0.5, unsigned long mit = 10,
     unsigned long tit = 4, double mres = 1.e-6, double mdiff = 1.e-6,
     double mcos = 0.9, double mress = 1.e-8, int noi = 0, int sing = 0,
     bool nonsm = false, double dmax = 0.005, double dmin = 0.00012,
     double tvar = 0.02, unsigned long ndir = 40, unsigned long ncomb = 1)
      : md(&m), singularities_(sing), nonsmooth(nonsm), parameter_name_(pn),
	with_parametrised_data(false), scfac_(sfac), lsolver(ls),
	h_init_(hin), h_max_(hmax), h_min_(hmin), h_inc_(hinc), h_dec_(hdec),
	maxit_(mit), thrit_(tit), maxres_(mres), maxdiff_(mdiff),
	mincos_(mcos), maxres_solve_(mress), delta_max_(dmax),
	delta_min_(dmin), thrvar_(tvar), nbdir_(ndir), nbcomb_(ncomb),
	noisy_(noi), tau_lp(0.), tau_bp_1(tau_bp_init),
	tau_bp_2(tau_bp_init), gamma_sing(0.), gamma_next(0.),
	build(BUILD_ALL)
    { GMM_ASSERT1(!md->is_complex(),
		  "Continuation has only a real version, sorry.");
      if (singularities_ > 1) init_border(); }
    
    cont_struct_getfem_model
    (model &m, const std::string &pn, const std::string &in,
     const std::string &fn, const std::string &cn, double sfac,
     rmodel_plsolver_type ls, double hin = 1.e-2, double hmax = 1.e-1,
     double hmin = 1.e-5, double hinc = 1.3, double hdec = 0.5,
     unsigned long mit = 10, unsigned long tit = 4, double mres = 1.e-6,
     double mdiff = 1.e-6, double mcos = 0.9, double mress = 1.e-8,
     int noi = 0, int sing = 0, bool nonsm = false, double dmax = 0.005,
     double dmin = 0.00012, double tvar = 0.02, unsigned long ndir = 40,
     unsigned long ncomb = 1)
      : md(&m), singularities_(sing), nonsmooth(nonsm), parameter_name_(pn),
	with_parametrised_data(true), initdata_name_(in),
	finaldata_name_(fn), currentdata_name_(cn), scfac_(sfac),
	lsolver(ls), h_init_(hin), h_max_(hmax), h_min_(hmin), h_inc_(hinc),
	h_dec_(hdec), maxit_(mit), thrit_(tit), maxres_(mres),
	maxdiff_(mdiff), mincos_(mcos), maxres_solve_(mress),
	delta_max_(dmax), delta_min_(dmin), thrvar_(tvar), nbdir_(ndir),
	nbcomb_(ncomb), noisy_(noi), tau_lp(0.), tau_bp_1(tau_bp_init),
	tau_bp_2(tau_bp_init), gamma_sing(0.), gamma_next(0.),
	build(BUILD_ALL)
    { GMM_ASSERT1(!md->is_complex(),
		  "Continuation has only a real version, sorry.");
      if (singularities_ > 1) init_border(); }

    cont_struct_getfem_model(void) {}
    

    // Linear algebra functions
    double abs(double a) { return gmm::abs(a); }
    void clear(VECT &v) { gmm::clear(v); }
    void copy(const VECT &v1, VECT &v) { gmm::copy(v1, v); }
    void scale(VECT &v, double a) { gmm::scale(v, a); }
    void scaled_add(const VECT &v1, const VECT &v2, double a, VECT &v)
    { gmm::add(v1, gmm::scaled(v2, a), v); }
    void scaled_add(const VECT &v1, double a1,
		    const VECT &v2, double a2, VECT &v)
    { gmm::add(gmm::scaled(v1, a1), gmm::scaled(v2, a2), v); }
    void scaled_add(const MAT &M1, double a1,
		    const MAT &M2, double a2, MAT &M)
    { gmm::add(gmm::scaled(M1, a1), gmm::scaled(M2, a2), M); }    
    void mult(const MAT &A, const VECT &v1, VECT &v)
    { gmm::mult(A, v1, v); }

    double sp(const VECT &v1, const VECT &v2)
    { return gmm::vect_sp(v1, v2); }
    double norm(const VECT &v)
    { return gmm::vect_norm2(v); }
    double w_sp(const VECT &v1, const VECT &v2)
    { return scfac_ * gmm::vect_sp(v1, v2); }
    double sp(const VECT &v1, const VECT &v2, double w1, double w2)
    { return sp(v1, v2) + w1 * w2; }
    double w_norm(const VECT &v, double w)
    { return sqrt(w_sp(v, v) + w * w); }
    double cosang(const VECT &v1, const VECT &v2, double w1, double w2) {
      double no = sqrt(sp(v1, v1, w1, w1) * sp(v2, v2, w2, w2)); 
      return ((no == 0) ? 0. : sp(v1, v2, w1, w2) / no);
    }
    
    double random(void) { return gmm::random(1.); }
    void fill_random(VECT &v) { gmm::fill_random(v); }

    void solve(const MAT &A, VECT &g, const VECT &L) { /* A * g = L */
      if (noisy_ > 2) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_, (noisy_ >= 2) ? noisy_ - 2 : 0,
			  40000);
      (*lsolver)(A, g, L, iter);
      if (noisy_ > 2) cout << "linear solver done" << endl;
    }

    void solve(const MAT &A, VECT &g1, VECT &g2,
	       const VECT &L1, const VECT &L2) { /* A * (g1|g2) = (L1|L2) */
      if (noisy_ > 2) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_, (noisy_ >= 2) ? noisy_ - 2 : 0,
			  40000);
      (*lsolver)(A, g1, L1, iter);
      iter.init(); (*lsolver)(A, g2, L2, iter); // (can be optimised)
      if (noisy_ > 2) cout << "linear solver done" << endl;
    }


    // Evaluation of  ...
    void set_variables(const VECT &x, double gamma) {
      md->set_real_variable(parameter_name_)[0] = gamma;
      if (with_parametrised_data) {
	gmm::add(gmm::scaled(md->real_variable(initdata_name_), 1. - gamma),
		 gmm::scaled(md->real_variable(finaldata_name_), gamma),
		 md->set_real_variable(currentdata_name_));
      }
      md->to_variables(x);
    }

    // F(x, gamma) --> f
    void F(const VECT &x, double gamma, VECT &f) {
      if (build == BUILD_ALL) set_variables(x, gamma);
      if (build & BUILD_F) {
	md->assembly(model::BUILD_RHS);
	build = build_data(build ^ BUILD_F);
      }
      gmm::copy(gmm::scaled(md->real_rhs(), -1.), f);
    }
    
    // (F(x, gamma + eps) - F(x, gamma)) / eps --> g
    void F_gamma(const VECT &x, double gamma, VECT &g) {
      const double eps = 1.e-8;
      VECT F0(x), F1(x);
      F(x, gamma, F0);
      build = BUILD_ALL; F(x, gamma + eps, F1); build = BUILD_ALL;
      gmm::add(F1, gmm::scaled(F0, -1.), g);
      gmm::scale(g, 1./eps);
    }

    void update_matrix(const VECT &x, double gamma) {
      if (build == BUILD_ALL) set_variables(x, gamma);
      if (build & BUILD_F_x) {
	if (noisy_ > 2) cout << "starting computing tangent matrix" << endl;
	md->assembly(model::BUILD_MATRIX);
	build = build_data(build ^ BUILD_F_x);
      }
    }
    // F_x(x, gamma) --> A
    void F_x(const VECT &x, double gamma, MAT &A) {
      update_matrix(x, gamma);
      unsigned long nbdof = md->nb_dof();
      gmm::resize(A, nbdof, nbdof);
      gmm::copy(md->real_tangent_matrix(), A);
    }

    // solve F_x(x, gamma) * g = L
    void solve_grad(const VECT &x, double gamma,
		    VECT &g, const VECT &L) {
      update_matrix(x, gamma);
      solve(md->real_tangent_matrix(), g, L);
    }

    // solve F_x(x, gamma) * (g1|g2) = (L1|L2)
    void solve_grad(const VECT &x, double gamma, VECT &g1,VECT &g2,
		    const VECT &L1, const VECT &L2) {
      update_matrix(x, gamma);
      solve(md->real_tangent_matrix(), g1, g2, L1, L2);
    }

    // F_x(x, gamma) * w --> y
    void mult_grad(const VECT &x, double gamma,
			const VECT &w, VECT &y) {
      update_matrix(x, gamma);
      mult(md->real_tangent_matrix(), w, y);
    }

    
    // Misc. for accessing private data
    model &linked_model(void) { return *md; }
    int singularities(void) { return singularities_; }
    bool non_smooth(void) { return nonsmooth; }
    std::string parameter_name(void) { return parameter_name_; }
    double scfac(void) { return scfac_; }
    double h_init(void) { return h_init_; }
    double h_min(void) { return h_min_; }
    double h_max(void) { return h_max_; }
    double h_dec(void) { return h_dec_; }
    double h_inc(void) { return h_inc_; }
    unsigned long maxit(void) { return maxit_; }
    unsigned long thrit(void) { return thrit_; }
    double maxres(void) { return maxres_; }
    double maxdiff(void) { return maxdiff_; }
    double mincos(void) { return mincos_; }
    double delta_max(void) { return delta_max_; }
    double delta_min(void) { return delta_min_; }
    double thrvar(void) { return thrvar_; }
    unsigned long nbdir(void) { return nbdir_; }
    unsigned long nbcomb(void) { return nbcomb_; }
    int noisy(void) { return noisy_; }
    VECT &b_x(void) { return b_x_; }
    VECT &c_x(void) { return c_x_; }
    double b_gamma(void) { return b_gamma_; }
    double c_gamma(void) { return c_gamma_; }
    double d(void) { return d_; }

    void set_tau_lp(double tau) { tau_lp = tau; }
    double get_tau_lp(void) { return tau_lp; }
    void set_tau_bp_1(double tau) { tau_bp_1 = tau; }
    double get_tau_bp_1(void) { return tau_bp_1; }
    void set_tau_bp_2(double tau) { tau_bp_2 = tau; }
    double get_tau_bp_2(void) { return tau_bp_2; }
    void clear_tau_bp_currentstep(void) {
      tau_bp_graph.clear();
      gmm::resize(alpha_hist, 0); gmm::resize(tau_bp_hist, 0);
    }
    void init_tau_bp_graph(void) { tau_bp_graph[0.] = tau_bp_2; }
    void insert_tau_bp_graph(double alpha, double tau) {
      tau_bp_graph[alpha] = tau;
    }
    VECT &get_alpha_hist(void) {
      unsigned long i = 0;
      gmm::resize(alpha_hist, tau_bp_graph.size());
      for (std::map<double, double>::iterator it = tau_bp_graph.begin();
	   it != tau_bp_graph.end(); it++) {
	alpha_hist[i] = (*it).first; i++;
      }	
      return alpha_hist;
    }
    VECT &get_tau_bp_hist(void) {
      unsigned long i = 0;
      gmm::resize(tau_bp_hist, tau_bp_graph.size());
      for (std::map<double, double>::iterator it = tau_bp_graph.begin();
	   it != tau_bp_graph.end(); it++) {
	tau_bp_hist[i] = (*it).second; i++; 
      }	
      return tau_bp_hist;
    }

    void clear_sing_data(void) {
      sing_label = "";
      gmm::resize(x_sing, 0); gmm::resize(x_next, 0);
      t_x_sing.clear(); t_gamma_sing.clear();
      t_x_predict.clear(); t_gamma_predict.clear();
    }
    void set_sing_label(std::string label) { sing_label = label; }
    std::string get_sing_label(void) { return sing_label; }
    void set_sing_point(const VECT &x, double gamma) {
      gmm::resize(x_sing, gmm::vect_size(x)); gmm::copy(x, x_sing);
      gamma_sing = gamma;
    }
    VECT &get_x_sing(void) { return x_sing; }
    double get_gamma_sing(void) { return gamma_sing; }
    unsigned long nb_tangent_sing(void) { return t_x_sing.size(); }
    bool insert_tangent_sing(const VECT &t_x, double t_gamma){
      bool is_included = false;
      unsigned long i = 0;
      double cang;
      while ((i < t_x_sing.size()) && (!is_included)){
	cang = cosang(t_x_sing[i], t_x, t_gamma_sing[i], t_gamma);
	is_included = (cang >= mincos_);
	++i;
      }
      if (!is_included) {
	t_x_sing.push_back(t_x); t_gamma_sing.push_back(t_gamma);
      }
      return !is_included;
    }
    VECT &get_t_x_sing(unsigned long i) { return t_x_sing[i]; }
    double get_t_gamma_sing(unsigned long i) { return t_gamma_sing[i]; }
    std::vector<VECT> &get_t_x_sing(void) { return t_x_sing; }
    std::vector<double> &get_t_gamma_sing(void) { return t_gamma_sing; }

    bool next_point(void) { return gmm::vect_size(x_next) > 0; }
    void set_next_point(const VECT &x, double gamma) {
      if (gmm::vect_size(x_next) == 0) {
	gmm::resize(x_next, gmm::vect_size(x)); gmm::copy(x, x_next);
	gamma_next = gamma;
      }
    }
    VECT &get_x_next(void) { return x_next; }
    double get_gamma_next(void) { return gamma_next; }

    bool insert_tangent_predict(const VECT &t_x, double t_gamma){
      bool is_included = false;
      unsigned long i = 0;
      double cang;
      while ((i < t_x_predict.size()) && (!is_included)){
	cang = gmm::abs(cosang(t_x_predict[i], t_x,
			       t_gamma_predict[i], t_gamma));
	is_included = (cang >= mincos_);
	++i;
      }
      if (!is_included) {
	t_x_predict.push_back(t_x); t_gamma_predict.push_back(t_gamma);
      }
      return !is_included;
    }

    void set_build(build_data build_) { build = build_; }
  };

#endif


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTINUATION_H__ */
