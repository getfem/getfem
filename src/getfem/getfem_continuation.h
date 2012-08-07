/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2011-2012 Tomas Ligursky, Yves Renard
 
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
    @brief (approximate) Moore-Penrose (also called Gauss-Newton) continuation method.

    NOTE: The bordered systems involved are solved by a block eliminiation
    although the bordered matrix may be ill-conditioned in some cases!
    Nevertheless, the algorithm seems to work well.
*/
#ifndef GETFEM_CONTINUATION_H__
#define GETFEM_CONTINUATION_H__

#include <getfem/getfem_model_solvers.h>

namespace getfem {


  //=========================================================================
  // Abstract Moore-Penrose continuation method
  //=========================================================================


  const double tau_init = 1.e6;
  enum build_data { BUILD_F = 1, BUILD_F_x = 2, BUILD_ALL = 3 };

  template <typename CONT_S, typename VECT> 
  double norm_(CONT_S &S, const VECT &x)
  { return sqrt(S.sp(x, x)); }
  
  template <typename CONT_S, typename VECT> 
  double w_sp_(CONT_S &S, const VECT &x1, const VECT &x2)
  { return S.scfac() * S.sp(x1, x2); }

  template <typename CONT_S, typename VECT> 
  double sp_(CONT_S &S, const VECT &x1, const VECT &x2,
	     double gamma1, double gamma2)
  { return w_sp_(S, x1, x2) + gamma1 * gamma2; }

  template <typename CONT_S, typename VECT> 
  double norm_(CONT_S &S, const VECT &x, double gamma)
  { return sqrt(sp_(S, x, x, gamma, gamma)); }


  template <typename CONT_S, typename VECT>
  void compute_tangent(CONT_S &S, const VECT &x, double gamma,
		       VECT &t_x, double &t_gamma) {
    VECT g(x), y(x);
    S.F_gamma(x, gamma, g);
    S.solve_grad(x, gamma, y, g);
    t_gamma = 1. / (t_gamma - w_sp_(S, t_x, y));
    S.scale(y, -t_gamma); S.copy(y, t_x);
    
    double no = norm_(S, t_x, t_gamma);
    S.scale(t_x, 1./no); t_gamma /= no;

    if (S.noisy() > 1) {
      S.mult_grad(x, gamma, t_x, y); S.scaled_add(y, g, t_gamma, y);
      cout << "tangent computed with the residual " 
	   << norm_(S, y) << endl;
    }
  }
  

  template <typename CONT_S, typename VECT>
  bool test_tangent(CONT_S &S, const VECT &x, double gamma,
		    const VECT &t_x, double t_gamma,
		    VECT &T_x, double T_gamma, double h) {
    bool res = false;
    double Gamma, T_Gamma = T_gamma, ang;
    VECT X(x), T_X(T_x);
    
    S.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
    S.set_build(BUILD_ALL);
    compute_tangent(S, X, Gamma, T_X, T_Gamma);
    
    ang = sp_(S, T_x, T_X, T_gamma, T_Gamma);
    if (S.noisy() > 1)
      cout << "the angle with the tested tangent " << ang << endl;
    if (ang >= S.minang()) res = true;
    else {
      ang = sp_(S, t_x, T_X, t_gamma, T_Gamma);
      if (S.noisy() > 1)
	cout << "the angle with the starting tangent " << ang << endl;
    }
    return res;
  }


  template <typename CONT_S, typename MAT, typename VECT>
  double test_function(CONT_S &S, const MAT &A, const VECT &g,
		       const VECT &t_x, double t_gamma) {
    double q, r, v_gamma, tau;
    VECT v_x(g), y(g), z(g);

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

  template <typename CONT_S, typename VECT>
  double test_function(CONT_S &S, const VECT &x, double gamma,
		       const VECT &t_x, double t_gamma) {
    VECT g(x); S.F_gamma(x, gamma, g);
    typename CONT_S::MAT A; S.F_x(x, gamma, A);
    return test_function(S, A, g, t_x, t_gamma);
  }

  template <typename CONT_S, typename VECT>
  bool test_smooth_bifurcation(CONT_S &S, const VECT &x, double gamma,
			       const VECT &t_x, double t_gamma) {
    double tau0 = S.get_tau1(), tau1 = S.get_tau2(),
      tau2 = test_function(S, x, gamma, t_x, t_gamma);
    S.set_tau1(tau1); S.set_tau2(tau2);
    return (tau2 * tau1 < 0) && (S.abs(tau1) < S.abs(tau0));
  }

  template <typename CONT_S, typename VECT>
  bool test_nonsmooth_bifurcation (CONT_S &S, const VECT &x1, double gamma1,
				   const VECT &t_x1, double t_gamma1,
				   const VECT &x2, double gamma2,
				   const VECT &t_x2, double t_gamma2) {
    unsigned long nb_changes = 0;
    double alpha = 0., delta = S.delta_min(),
      tau0 = tau_init, tau1= S.get_tau2(),tau2, tau_var_ref, t_gamma;
    VECT g1(x1), g2(x1), g(x1), t_x(x1);
    typename CONT_S::MAT A1, A2, A;
    S.F_gamma(x2, gamma2, g2);
    S.F_x(x1, gamma1, A1); S.F_gamma(x1, gamma1, g1);
    S.F_x(x2, gamma2, A2); S.F_x(x2, gamma2, A);
    S.init_tau_graph();
    tau2 = test_function(S, A2, g2, t_x2, t_gamma2);
    tau_var_ref = std::max(S.abs(tau2 - tau1),
			   (S.abs(tau1) + S.abs(tau2)) / 200);

    do {
      alpha = std::min(alpha + delta, 1.);
      S.scaled_add(A1, 1. - alpha, A2, alpha, A);
      S.scaled_add(g1, 1. - alpha, g2, alpha, g);
      S.scaled_add(t_x1, 1. - alpha, t_x2, alpha, t_x);
      t_gamma = (1. - alpha) * t_gamma1 + alpha * t_gamma2;
      
      tau2 = test_function(S, A, g, t_x, t_gamma);
      if ((tau2 * tau1 < 0) && (S.abs(tau1) < S.abs(tau0))) ++nb_changes;
      S.insert_tau_graph(alpha, tau2);

      if (S.abs(tau2 - tau1) < 0.5 * S.thrvar() * tau_var_ref)
	delta = std::min(2 * delta, S.delta_max());
      else if (S.abs(tau2 - tau1) > S.thrvar() * tau_var_ref) 
	delta = std::max(0.1 * delta, S.delta_min());
      tau0 = tau1; tau1 = tau2; 
    } while (alpha < 1.);
    
    S.set_tau1(tau_init); S.set_tau2(tau2);
    return nb_changes % 2;
  }

  template <typename CONT_S, typename VECT>
  void init_test_function(CONT_S &S, const VECT &x, double gamma,
			  const VECT &t_x, double t_gamma) {
    S.set_build(BUILD_ALL); S.init_border(x);
    if (S.noisy() > 0) cout << "starting computing initial value of the "
			    << "test function for bifurcations" << endl;
    double tau = test_function(S, x, gamma, t_x, t_gamma); S.set_tau2(tau);
  }

  template <typename CONT_S, typename VECT>
  void init_Moore_Penrose_continuation(CONT_S &S, const VECT &x,
				       double gamma, VECT &t_x,
				       double &t_gamma, double &h) {
    S.set_build(BUILD_ALL);
    S.clear(t_x); t_gamma = (t_gamma >= 0) ? 1. : -1.;
    if (S.noisy() > 0) cout << "starting computing initial tangent" << endl;
    compute_tangent(S, x, gamma, t_x, t_gamma);
    h = S.h_init();
    init_test_function(S, x, gamma, t_x, t_gamma);
  }

  
  /* Perform one step of the Moore-Penrose continuation. If a new point 
     (x, gamma) is found, it has to be saved in the model in the end! */
  template <typename CONT_S, typename VECT>
    void Moore_Penrose_continuation(CONT_S &S, VECT &x, double &gamma,
				    VECT &t_x, double &t_gamma, double &h) {
    bool converged, finished = false;
    int tangent_status = -1;
      /* -1: no manipulation with tangent direction (so far);
	  0: current direction not admitted;
	  1: current direction admitted; */
    unsigned long it, step_dec = 0;
    double t_gamma0 = t_gamma, Delta_Gamma, Gamma, T_gamma,
      r, no, res, diff, ang;
    VECT t_x0(t_x), F(x), g(x), Delta_X(x), X(x), T_x(x), y(x);

    do { // step control

      // prediction
      if (S.noisy() > 0) cout << "prediction with h = " << h << endl;
      S.scaled_add(x, t_x, h, X); Gamma = gamma + h * t_gamma;
      S.set_build(BUILD_ALL);
      S.copy(t_x, T_x); T_gamma = t_gamma;
      
      // correction
      if (S.noisy() > 0) cout << "starting correction " << endl;
      it = 0;
      S.F(X, Gamma, F);
      
      do { // Newton iterations
	S.F_gamma(X, Gamma, g);
	S.solve_grad(X, Gamma, Delta_X, y, F, g);
	r = w_sp_(S, T_x, y);

	Delta_Gamma = w_sp_(S, T_x, Delta_X) / (r - T_gamma);
	S.scaled_add(Delta_X, y, -Delta_Gamma, Delta_X);
	S.scaled_add(X, Delta_X, -1., X); Gamma -= Delta_Gamma;
	S.set_build(BUILD_ALL);
	
	T_gamma = 1. / (T_gamma - r);
	S.scale(y, -T_gamma); S.copy(y, T_x);
	no = norm_(S, T_x, T_gamma);
	S.scale(T_x, 1./no); T_gamma /= no;

	S.F(X, Gamma, F); res = norm_(S, F); 
	diff = norm_(S, Delta_X, Delta_Gamma);
	converged = (res <= S.maxres() && diff <= S.maxdiff());
	it++;

	if (S.noisy() > 0)
	  cout << "iter " << it << " residual " << res
	       << " difference " << diff
	       << " cos " << sp_(S, t_x, T_x, t_gamma, T_gamma) << endl;

      } while (!converged && it < S.maxit() && res < 1.e8);

      if (converged) {
// 	compute_tangent(S, X, Gamma, T_x, T_gamma);
	ang = sp_(S, t_x, T_x, t_gamma, T_gamma);
	if (S.noisy() > 1) {
	  S.F_gamma(X, Gamma, g);
	  S.update_matrix(X, Gamma); S.mult_grad(X, Gamma, T_x, y);
	  S.scaled_add(y, g, T_gamma, y);
	  cout << "the last tangent computed with the residual "
	       << norm_(S, y) << endl;
	}
	if (ang >= S.minang()) { // accept the new couple
	  S.set_message(""); S.clear_tau_graph();
	  if (S.noisy() > 0)
	    cout << "starting computing test function for bifurcations"
		 << endl;
	  bool bifurcation_detected = false;
	  if (tangent_status < 0)
	    bifurcation_detected =
	      test_smooth_bifurcation(S, X, Gamma, T_x, T_gamma);
	  else {
	    bifurcation_detected = test_nonsmooth_bifurcation
	      (S, x, gamma, t_x0, t_gamma0, X, Gamma, T_x, T_gamma);
	  }
	  if (bifurcation_detected) {
	    cout << "Bifurcation detected!" << endl;
	    S.set_message("Bifurcation detected!");
	  }
	  if (step_dec == 0 && it < S.thrit()) // elongate the step size
	    h = (S.h_inc() * h < S.h_max()) ? S.h_inc() * h : S.h_max();
	  finished = true;
	}
      }
      
      if (!finished) {
	if (h > S.h_min()) { // diminish the step size
	  h = (S.h_dec() * h > S.h_min()) ? S.h_dec() * h : S.h_min();
	  step_dec++;
	}
	else if (tangent_status < 0) {
	  if (S.noisy() > 0)
	    cout << "classical continuation has failed, "
		 << "starting seeking a new tangent" << endl;
	  S.copy(t_x, T_x); T_gamma = t_gamma;
	  S.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
	  S.set_build(BUILD_ALL);
	  compute_tangent(S, X, Gamma, T_x, T_gamma);

	  double h_test = (-0.9) * S.h_min();
	  do {
	    h_test = -h_test
	      + pow(10., floor(log10(- h_test / S.h_min()))) * S.h_min();
	    tangent_status = (int) test_tangent(S, x, gamma, t_x, t_gamma,
						T_x, T_gamma, h_test);
	    if (tangent_status == 0) {
	      h_test *= -1.;
	      tangent_status = (int) test_tangent(S, x, gamma, t_x, t_gamma,
						  T_x, T_gamma, h_test);
	    }
	  } while ((tangent_status == 0) && (h_test > -S.h_max()));
	  
	  if (tangent_status > 0) {
	    S.copy(T_x, t_x); t_gamma = T_gamma;
	    if (h_test < 0) { 
	      S.scale(t_x, -1.); t_gamma *= -1.; h_test *= -1.;
	    }
	    if (S.noisy() > 1)
	      cout << "starting computing a suitable stepsize" << endl;
	    bool h_adapted = false; h = S.h_init();
	    while ((!h_adapted) && (h > h_test)) {
	      h_adapted
		= test_tangent(S, x, gamma, t_x0, t_gamma0, t_x, t_gamma, h);
	      h *= S.h_dec();
	    }
	    h = (h_adapted) ? h / S.h_dec() : h_test;
	  } else break;
	} else break;
      }
    } while (!finished);

    if (finished) {
      S.copy(X, x); gamma = Gamma;
      S.copy(T_x, t_x); t_gamma = T_gamma;
    } else h = 0;
  }



  //=========================================================================
  // Moore-Penrose continuation method for Getfem models
  //=========================================================================


#ifdef GETFEM_MODELS_H__
 
  struct cont_struct_getfem_model {
    
    typedef base_vector VECT;
    typedef model_real_sparse_matrix MAT;
    
  private:
    model *md;  // for real models only
    std::string parameter_name_;
    bool with_parametrised_data;
    std::string initdata_name_, finaldata_name_, currentdata_name_;
    double scfac_;
    rmodel_plsolver_type lsolver;
    double h_init_, h_max_, h_min_, h_inc_, h_dec_;
    unsigned long maxit_, thrit_;
    double maxres_, maxdiff_, minang_, maxres_solve_, epsilon_, delta_max_,
      delta_min_, thrvar_;
    int noisy_;
    std::string message;
    VECT b_x_, c_x_;
    double b_gamma_, c_gamma_, d_;
    double tau1, tau2;
    VECT alpha_hist, tau_hist;
    std::map<double, double> tau_graph;
    build_data build;

  public:
    cont_struct_getfem_model
    (model &m, const std::string &pn, double sfac, rmodel_plsolver_type ls,
     double hin = 1.e-2, double hmax = 1.e-1, double hmin = 1.e-5,
     double hinc = 1.3, double hdec = 0.5, unsigned long mit = 10,
     unsigned long tit = 4, double mres = 1.e-6, double mdiff = 1.e-6,
     double mang = 0.9, double mress = 1.e-8, double eps = 1.e-8, 
     int noi = 0, double dmax = 0.005, double dmin = 0.00012,
     double tvar = 0.02)
      : md(&m), parameter_name_(pn), with_parametrised_data(false),
	scfac_(sfac), lsolver(ls), h_init_(hin), h_max_(hmax), h_min_(hmin),
	h_inc_(hinc), h_dec_(hdec),maxit_(mit), thrit_(tit), maxres_(mres),
	maxdiff_(mdiff), minang_(mang), maxres_solve_(mress),
	epsilon_(eps), delta_max_(dmax), delta_min_(dmin), thrvar_(tvar),
	noisy_(noi), message(""), tau1(tau_init), tau2(tau_init),
	alpha_hist(0), tau_hist(0), build(BUILD_ALL)
    {}
    
    cont_struct_getfem_model
    (model &m, const std::string &pn, const std::string &in,
     const std::string &fn, const std::string &cn, double sfac,
     rmodel_plsolver_type ls, double hin = 1.e-2, double hmax = 1.e-1,
     double hmin = 1.e-5, double hinc = 1.3, double hdec = 0.5,
     unsigned long mit = 10, unsigned long tit = 4, double mres = 1.e-6,
     double mdiff = 1.e-6, double mang = 0.9, double mress = 1.e-8,
     double eps = 1.e-8, int noi = 0, double dmax = 0.005,
     double dmin = 0.00012, double tvar = 0.02)
      : md(&m), parameter_name_(pn), with_parametrised_data(true),
	initdata_name_(in), finaldata_name_(fn), currentdata_name_(cn),
	scfac_(sfac), lsolver(ls), h_init_(hin), h_max_(hmax), h_min_(hmin),
	h_inc_(hinc), h_dec_(hdec), maxit_(mit), thrit_(tit), maxres_(mres),
	maxdiff_(mdiff), minang_(mang), maxres_solve_(mress), epsilon_(eps),
	delta_max_(dmax), delta_min_(dmin), thrvar_(tvar), noisy_(noi),
	message(""), tau1(tau_init), tau2(tau_init), alpha_hist(0),
	tau_hist(0), build(BUILD_ALL)
    {}

    cont_struct_getfem_model(void) {}
    

    // Linear algebra functions
    double abs(double a)
    { return gmm::abs(a); }
    void clear(VECT &v)
    { gmm::clear(v); }
    void copy(const VECT &v1, VECT &v)
    { gmm::copy(v1, v); }
    void scale(VECT &v, double a)
    { gmm::scale(v, a); }
    void scaled_add(const VECT &v1, const VECT &v2, double a, VECT &v)
    { gmm::add(v1, gmm::scaled(v2, a), v); }
    void scaled_add(const VECT &v1, double a1,
		    const VECT &v2, double a2, VECT &v)
    { gmm::add(gmm::scaled(v1, a1), gmm::scaled(v2, a2), v); }
    void scaled_add(const MAT &M1, double a1,
		    const MAT &M2, double a2, MAT &M)
    { gmm::add(gmm::scaled(M1, a1), gmm::scaled(M2, a2), M); }
    double sp(const VECT &v1, const VECT &v2)
    { return gmm::vect_sp(v1, v2); }
    void mult(const MAT &A, const VECT &v1, VECT &v)
    { gmm::mult(A, v1, v); }

    void solve(const MAT &A, VECT &g, const VECT &L) { /* A * g = L */
      if (noisy_ > 1) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_,
			  (noisy_ > 0) ? noisy_ - 1 : noisy_, 40000);
      (*lsolver)(A, g, L, iter);
      if (noisy_ > 1) cout << "linear solver done" << endl;
    }

    void solve(const MAT &A, VECT &g1, VECT &g2,
	       const VECT &L1, const VECT &L2) { /* A * (g1|g2) = (L1|L2) */
      if (noisy_ > 1) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_,
			  (noisy_ > 0) ? noisy_ - 1 : noisy_, 40000);
      (*lsolver)(A, g1, L1, iter);
      iter.init(); (*lsolver)(A, g2, L2, iter); // (can be optimised)
      if (noisy_ > 1) cout << "linear solver done" << endl;
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
    
    // (F(x, gamma + epsilon_) - F(x, gamma)) / epsilon_ --> g
    void F_gamma(const VECT &x, double gamma, VECT &g) {
      VECT F0(x), F1(x);
      F(x, gamma, F0);
      build = BUILD_ALL; F(x, gamma + epsilon_, F1); build = BUILD_ALL;
      gmm::add(F1, gmm::scaled(F0, -1.), g);
      gmm::scale(g, 1./epsilon_);
    }

    void update_matrix(const VECT &x, double gamma) {
      if (build == BUILD_ALL) set_variables(x, gamma);
      if (build & BUILD_F_x) {
	if (noisy_ > 1) cout << "starting computing tangent matrix" << endl;
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

    
    // Misc.
    model &linked_model(void) { return *md; }
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
    double minang(void) { return minang_; }
    double epsilon(void) { return epsilon_; }
    double delta_max(void) { return delta_max_; }
    double delta_min(void) { return delta_min_; }
    double thrvar(void) { return thrvar_; }
    int noisy(void) { return noisy_; }
    void set_message(std::string msg) { message = msg; }
    std::string get_message(void) { return message; }

    void init_border(const VECT &v) {
      srand(unsigned(time(NULL)));
      gmm::resize(b_x_, gmm::vect_size(v)); gmm::fill_random(b_x_);
      gmm::resize(c_x_, gmm::vect_size(v)); gmm::fill_random(c_x_);
      b_gamma_ = gmm::random(1.); c_gamma_ = gmm::random(1.);
      d_ = gmm::random(1.);
    }
    VECT &b_x(void) { return b_x_; }
    VECT &c_x(void) { return c_x_; }
    double b_gamma(void) { return b_gamma_; }
    double c_gamma(void) { return c_gamma_; }
    double d(void) { return d_; }
    void set_tau1(double tau) { tau1 = tau; }
    double get_tau1(void) { return tau1; }
    void set_tau2(double tau) { tau2 = tau; }
    double get_tau2(void) { return tau2; }

    void clear_tau_graph(void) {
      tau_graph.clear();
      gmm::resize(alpha_hist, 0); gmm::resize(tau_hist, 0);
    }
    void init_tau_graph(void) { tau_graph[0.] = tau2; }
    void insert_tau_graph(double alpha, double tau) {
      tau_graph[alpha] = tau;
    }
    VECT &get_alpha_hist(void) {
      gmm::resize(alpha_hist, tau_graph.size()); unsigned long i = 0;
      for (std::map<double, double>::iterator it = tau_graph.begin();
	   it != tau_graph.end(); it++) {
	alpha_hist[i] = (*it).first; i++;
      }	
      return alpha_hist;
    }
    VECT &get_tau_hist(void) {
      gmm::resize(tau_hist, tau_graph.size()); unsigned long i = 0;
      for (std::map<double, double>::iterator it = tau_graph.begin();
	   it != tau_graph.end(); it++) {
	tau_hist[i] = (*it).second; i++; 
      }	
      return tau_hist;
    }

    void set_build(build_data build_) { build = build_; }
  };

#endif


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTINUATION_H__ */
