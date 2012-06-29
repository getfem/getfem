/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2011-2012 Tomas Ligursky, Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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


  enum build_data { BUILD_F = 1, BUILD_F_x = 2, BUILD_ALL = 3 };

  template <typename CONT_S, typename VECT> 
  double norm_(CONT_S &s, const VECT &x) {
    double no = sqrt(s.sp(x, x));
    return no;
  }
  
  template <typename CONT_S, typename VECT> 
  double w_sp_(CONT_S &s, const VECT &x1, const VECT &x2) {
    double r = s.scfac() * s.sp(x1, x2);
    return r;
  }

  template <typename CONT_S, typename VECT> 
  double sp_(CONT_S &s, const VECT &x1, const VECT &x2,
	     double gamma1, double gamma2) {
    double r = w_sp_(s, x1, x2) + gamma1 * gamma2;
    return r;
  }

  template <typename CONT_S, typename VECT> 
  double norm_(CONT_S &s, const VECT &x, double gamma) {
    double no = sqrt(sp_(s, x, x, gamma, gamma));
    return no;
  }


  template <typename CONT_S, typename VECT>
  void compute_tangent(CONT_S &s, const VECT &x, double gamma,
		       VECT &t_x, double &t_gamma) {
    VECT g(x), y(x);
    s.F_gamma(x, gamma, g);
    s.solve_grad(x, gamma, g, y);
    t_gamma = 1. / (t_gamma - w_sp_(s, t_x, y));
    s.scale(y, -t_gamma); s.copy(y, t_x);
    
    double no = norm_(s, t_x, t_gamma);
    s.scale(t_x, 1./no); t_gamma /= no;

//     if (s.noisy() > 1) {
//       s.mult_grad(x, gamma, t_x, y); s.scaled_add(y, g, t_gamma, y);
//       cout << "new tangent computed with the residual " 
// 	   << norm_(s, y) << endl;
//     }
  }


  /* Compute the test function for bifurcations corresponding to the system
     where the augmented Jacobian is bordered by b = c = e_1, d = 0. */
  template <typename CONT_S, typename VECT>
  double test_function(CONT_S &s, const VECT &x, double gamma,
		       const VECT &t_x, double t_gamma) {
    double q, r, v_gamma, tau;
    VECT b_x(x), g(x), v_x(x), y(x), z(x);

    if (s.noisy() > 1) cout << "starting computing test function" << endl;
    s.clear(b_x); b_x[0] = 1.;
    s.F_gamma(x, gamma, g);
    s.solve_grad(x, gamma, g, b_x, y, z);
    v_gamma = s.sp(t_x, z) / (s.sp(t_x, y) - t_gamma);
    s.scaled_add(z, y, -v_gamma, v_x);
    tau = -1. / v_x[0]; s.scale(v_x, -tau); v_gamma *= -tau;

    s.mult_grad(x, gamma, v_x, y);
    s.scaled_add(y, g, v_gamma, y); y[0] += tau;
    r = s.sp(y, y);
    q = s.sp(t_x, v_x) + t_gamma * v_gamma; r += q * q;
    q = v_x[0] - 1.; r += q * q; r = sqrt(r);
    if (r > 1e-10)
      GMM_WARNING1("Test function evaluated with the residual " << r);

    return tau;
  }

  template <typename CONT_S, typename VECT>
  bool test_smooth_bifurcation(CONT_S &s, const VECT &x, double gamma,
			       const VECT &t_x, double t_gamma) {
    double tau1 = s.tau2(), tau2 = s.tau3(),
      tau3 = test_function(s, x, gamma, t_x, t_gamma);
    bool res = (tau3 * tau2 < 0) & (s.abs(tau2) < s.abs(tau1));
    s.set_tau1(tau1); s.set_tau2(tau2); s.set_tau3(tau3);
    return res;
  }
  

  template <typename CONT_S, typename VECT>
  int test_direction(CONT_S &s, const VECT &x, double gamma,
		     const VECT &t_x, double t_gamma,
		     VECT &T_x, double &T_gamma, double h) {
    int res = 1;
    double Gamma, T_Gamma = T_gamma, ang;
    VECT X(x), T_X(T_x);
    
    s.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
    s.set_build(BUILD_ALL);
    compute_tangent(s, X, Gamma, T_X, T_Gamma);
    
    ang = sp_(s, T_x, T_X, T_gamma, T_Gamma);
    if (s.noisy() > 1)
      cout << "the angle with the tested tangent " << ang << endl;
    if (ang >=  0.996) res = (h > 0) ? 3 : 4;
    else {
      ang = sp_(s, t_x, T_X, t_gamma, T_Gamma);
      if (s.noisy() > 1)
	cout << "the angle with the starting tangent " << ang << endl;
      if (ang < 0.86 && ang > -0.86) {
	res = 2;
	s.copy(T_X, T_x); T_gamma = T_Gamma; // try the new tangent next(?)
      }
    }
    return res;
  }


  template <typename CONT_S, typename VECT>
  void init_Moore_Penrose_continuation(CONT_S &s, const VECT &x,
				       double gamma, VECT &t_x,
				       double &t_gamma, double &h) {
    s.set_build(BUILD_ALL);
    s.clear(t_x); t_gamma = (t_gamma >= 0) ? 1. : -1.;
    if (s.noisy() > 0) cout << "computing initial tangent" << endl;
    compute_tangent(s, x, gamma, t_x, t_gamma);
    h = s.h_init();
    double tau = test_function(s, x, gamma, t_x, t_gamma); s.set_tau3(tau);
  }

  
  /* Perform one step of the Moore-Penrose continuation. If a new point 
     (x, gamma) is found, it has to be saved in the model in the end! */
  template <typename CONT_S, typename VECT>
    void Moore_Penrose_continuation(CONT_S &s, VECT &x, double &gamma,
				    VECT &t_x, double &t_gamma, double &h) {
    bool bifurcation = false, converged, finished = false;
    int tangent_status = 0;
      /* 0: no manipulation with tangent direction so far;
	 1: current direction neither admitted nor rejected;
	 2: direction rejected;
	 3: direction admitted with plus sign;
	 4: direction admitted with minus sign; */
    unsigned long it, step_dec = 0;
    double Delta_Gamma, Gamma, T_gamma, r, no, res, diff, ang;
    VECT F(x), g(x), Delta_X(x), X(x), T_x(x), y(x);

    do { // step control

      // prediction
      if (s.noisy() > 0) cout << "prediction with h = " << h << endl;
      s.scaled_add(x, t_x, h, X); Gamma = gamma + h * t_gamma;
      s.set_build(BUILD_ALL);
      s.copy(t_x, T_x); T_gamma = t_gamma;
      
      // correction
      if (s.noisy() > 0) cout << "starting correction " << endl;
      it = 0;
      s.F(X, Gamma, F);
      
      do { // Newton iterations
	s.F_gamma(X, Gamma, g);
	s.solve_grad(X, Gamma, F, g, Delta_X, y);
	r = w_sp_(s, T_x, y);

	Delta_Gamma = w_sp_(s, T_x, Delta_X) / (r - T_gamma);
	s.scaled_add(Delta_X, y, -Delta_Gamma, Delta_X);
	s.scaled_add(X, Delta_X, -1., X); Gamma -= Delta_Gamma;
	s.set_build(BUILD_ALL);
	
	T_gamma = 1. / (T_gamma - r);
	s.scale(y, -T_gamma); s.copy(y, T_x);
	no = norm_(s, T_x, T_gamma);
	s.scale(T_x, 1./no); T_gamma /= no;

	s.F(X, Gamma, F); res = norm_(s, F); 
	diff = norm_(s, Delta_X, Delta_Gamma);
	converged = (res <= s.maxres() && diff <= s.maxdiff());
	it++;

	if (s.noisy() > 1)
	  cout << "iter " << it << " residual " << res
	       << " difference " << diff
	       << " cos " << sp_(s, t_x, T_x, t_gamma, T_gamma) << endl;

      } while (!converged && it < s.maxit() && res < 1.e8);

      if (converged) {
	ang = sp_(s, t_x, T_x, t_gamma, T_gamma);
	if (s.noisy() > 0) cout << "cos " << ang << endl;
// 	if (s.noisy() > 1) {
// 	  s.F_gamma(X, Gamma, g);
// 	  s.update_matrix(X, Gamma); s.mult_grad(X, Gamma, T_x, y);
// 	  s.scaled_add(y, g, T_gamma, y);
// 	  cout << "final tangent computed with the residual "
// 	       << norm_(s, y) << endl;
// 	}
	if (ang >= s.minang()) { // accept the new couple
// 	  if (tangent_status == 0)
	    bifurcation = test_smooth_bifurcation(s, X, Gamma, T_x, T_gamma);
	  if (bifurcation) cout << "Bifurcation detected!" << endl;
	  if (step_dec == 0 && it < s.thrit()) // elongate the step size
	    h = (s.h_inc() * h < s.h_max()) ? s.h_inc() * h : s.h_max();
	  finished = true;
	}
      }
      
      if (!finished) {
	if (h > s.h_min()) { // diminish the step size
	  h = (s.h_dec() * h > s.h_min()) ? s.h_dec() * h : s.h_min();
	  step_dec++;
	}
	else if (tangent_status == 0) {
	  if (s.noisy() > 1)
	    cout << "Seeking a new tangent direction" << endl;
	  unsigned long tan = 0;
	  s.copy(t_x, T_x); T_gamma = t_gamma;
	  s.scaled_add(x, T_x, h, X); Gamma = gamma + h * T_gamma;
	  s.set_build(BUILD_ALL);
	  compute_tangent(s, X, Gamma, T_x, T_gamma);

	  do { // seek a new tangent
	    if (s.noisy() > 1)
	      cout << "Trying direction " << tan + 1 << endl;
	    h = s.h_min();

	    do { // test (T_x, T_gamma)
	      tangent_status =
		test_direction(s, x, gamma, t_x, t_gamma, T_x, T_gamma, h);
	      if (tangent_status == 1) {
		h *= -1.;
		tangent_status =
		  test_direction(s, x, gamma, t_x, t_gamma, T_x, T_gamma, h);
	        h *= -2.;
	      }
	    } while (tangent_status == 1 && h <= 1e5);

	    tan++;
	  } while (tangent_status <= 2 && tan < 1); // tan >= 1?
	  
	  if (tangent_status >= 3) {
	    if (s.noisy() > 1)
	      cout << "Direction " << tan << " accepted" << endl;
	    s.copy(T_x, t_x); t_gamma = T_gamma;
	    if (tangent_status == 4) { 
	      s.scale(t_x, -1.); t_gamma *= -1.; h /= 2.;
	    }
	    h = (h < s.h_init()) ? s.h_init() : h; step_dec = 0;
	  } else break;
	} else break;
      }
    } while (!finished);

    if (finished) {
      s.copy(X, x); gamma = Gamma;
      s.copy(T_x, t_x); t_gamma = T_gamma;
    } else h = 0;
  }



  //=========================================================================
  // Moore-Penrose continuation method for Getfem models
  //=========================================================================


#ifdef GETFEM_MODELS_H__
 
  struct cont_struct_getfem_model {

    model *md;  // for real models only
    std::string parameter_name;
    rmodel_plsolver_type lsolver;
    double scfac_;
    unsigned long maxit_, thrit_;
    double maxres_, maxdiff_, minang_, h_init_, h_max_, h_min_, h_inc_,
      h_dec_, epsilon_, maxres_solve_;
    int noisy_;
    bool with_parametrized_data;
    std::string initdata_name, finaldata_name, currentdata_name;
    build_data build;
    double tau1_, tau2_, tau3_;


    typedef base_vector VECT;

    cont_struct_getfem_model
    (model &m, const std::string &pn, rmodel_plsolver_type ls, double sfac,
     unsigned long mit = 10, unsigned long tit = 8, double mres = 1.e-6,
     double mdiff = 1.e-9, double mang = 0.9, double hin = 1.e-2,
     double hmax = 1.e-1, double hmin = 1.e-5, double hinc = 1.3,
     double hdec = 0.5, double eps = 1.e-8, double mress = 1.e-7,
     int noi = 0, double t1 = 1.e4, double t2 = 1.e4, double t3 = 1.e4)
      : md(&m), parameter_name(pn), lsolver(ls), scfac_(sfac), maxit_(mit),
	thrit_(tit), maxres_(mres), maxdiff_(mdiff), minang_(mang),
	h_init_(hin), h_max_(hmax), h_min_(hmin), h_inc_(hinc), h_dec_(hdec),
	epsilon_(eps), maxres_solve_(mress), noisy_(noi),
	with_parametrized_data(false), build(BUILD_ALL), tau1_(t1),
	tau2_(t2), tau3_(t3)
    {}
    
    cont_struct_getfem_model
    (model &m, const std::string &pn, const std::string &in,
     const std::string &fn, const std::string &cn, rmodel_plsolver_type ls,
     double sfac, unsigned long mit = 10, unsigned long tit = 8,
     double mres = 1.e-6, double mdiff = 1.e-9, double mang = 0.9,
     double hin = 1.e-2, double hmax = 1.e-1, double hmin = 1.e-5,
     double hinc = 1.3, double hdec = 0.5, double eps = 1.e-8,
     double mress = 1.e-7, int noi = 0, double t1 = 1.e4,
     double t2 = 1.e4, double t3 = 1.e4)
      : md(&m), parameter_name(pn), lsolver(ls), scfac_(sfac), maxit_(mit),
	thrit_(tit), maxres_(mres), maxdiff_(mdiff), minang_(mang),
	h_init_(hin), h_max_(hmax), h_min_(hmin), h_inc_(hinc), h_dec_(hdec),
	epsilon_(eps), maxres_solve_(mress), noisy_(noi),
	with_parametrized_data(true), initdata_name(in), finaldata_name(fn),
	currentdata_name(cn), build(BUILD_ALL), tau1_(t1),
	tau2_(t2), tau3_(t3)
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
    double sp(const VECT &v1, const VECT &v2)
    { return gmm::vect_sp(v1, v2); }


    // Evaluation of  ...
    void set_variables(const VECT &x, double gamma) {
      md->set_real_variable(parameter_name)[0] = gamma;
      if (with_parametrized_data) {
	gmm::add(gmm::scaled(md->real_variable(initdata_name), 1. - gamma),
		 gmm::scaled(md->real_variable(finaldata_name), gamma),
		 md->set_real_variable(currentdata_name));
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

    // solve F_x(x, gamma) * g = L
    void solve_grad(const VECT &x, double gamma,
		    const VECT &L, VECT &g) {
      update_matrix(x, gamma);
      if (noisy_ > 1) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_, noisy_, 40000);
      (*lsolver)(md->real_tangent_matrix(), g, L, iter);
      if (noisy_ > 1) cout << "linear solver done" << endl;
    }

    // solve F_x(x, gamma) * (g1|g2) = (L1|L2)
    // can be optimised!
    void solve_grad(const VECT &x, double gamma, const VECT &L1,
		    const VECT &L2, VECT &g1, VECT &g2) {
      update_matrix(x, gamma);
      if (noisy_ > 1) cout << "starting linear solver" << endl;
      gmm::iteration iter(maxres_solve_, noisy_, 40000);
      (*lsolver)(md->real_tangent_matrix(), g1, L1, iter);
      iter.init();
      (*lsolver)(md->real_tangent_matrix(), g2, L2, iter);
      if (noisy_ > 1) cout << "linear solver done" << endl;
    }

    // F_x(x, gamma) * w --> y
    void mult_grad(const VECT &x, double gamma,
			const VECT &w, VECT &y) {
      update_matrix(x, gamma);
      gmm::mult(md->real_tangent_matrix(), w, y);
    }

    
    // Misc.
    model &linked_model(void) { return *md; }
    double scfac(void) { return scfac_; }
    unsigned long thrit(void) { return thrit_; }
    unsigned long maxit(void) { return maxit_; }
    double epsilon(void) { return epsilon_; }
    double minang(void) { return minang_; }
    double maxres(void) { return maxres_; }
    double maxdiff(void) { return maxdiff_; }
    double h_init(void) { return h_init_; }
    double h_min(void) { return h_min_; }
    double h_max(void) { return h_max_; }
    double h_dec(void) { return h_dec_; }
    double h_inc(void) { return h_inc_; }
    int noisy(void) { return noisy_; }
    void set_build(build_data build_) { build = build_; }
    void set_tau1(double tau) { tau1_ = tau; }
    double tau1(void) { return tau1_; }
    void set_tau2(double tau) { tau2_ = tau; }
    double tau2(void) { return tau2_; }
    void set_tau3(double tau) { tau3_ = tau; }
    double tau3(void) { return tau3_; }

  };

#endif


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTINUATION_H__ */
