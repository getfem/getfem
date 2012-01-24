// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2011 Yves Renard
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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file gmm_solver_Newton.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @author  Michel Fournie <fournie@mip.ups-tlse.fr>
   @date January 24, 2006.
*/
#ifndef GMM_SOLVERS_NEWTON_H__
#define GMM_SOLVERS_NEWTON_H__ 

namespace gmm {

#include "gmm_kernel.h"

  /* ***************************************************************** */
  /*     Line search definition                                        */
  /* ***************************************************************** */

  struct abstract_newton_line_search {
    double conv_alpha, conv_r;
    size_t it, itmax, glob_it;
    //  size_t tot_it;
    virtual void init_search(double r, size_t git, double R0 = 0.0) = 0;
    virtual double next_try(void) = 0;
    virtual bool is_converged(double, double R1 = 0.0) = 0;
    virtual double converged_value(void) {
      // tot_it += it; cout << "tot_it = " << tot_it << endl; it = 0;
      return conv_alpha;
    };
    virtual double converged_residual(void) { return conv_r; };
    virtual ~abstract_newton_line_search() { }
  };


  struct quadratic_newton_line_search : public abstract_newton_line_search {
    double R0_, R1_;

    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    virtual void init_search(double r, size_t git, double R0 = 0.0) {
      GMM_ASSERT1(R0 != 0.0, "You have to specify R0");
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
      R0_ = R0;
    }
    virtual double next_try(void) {
      ++it;
      if (it == 1) return double(1);
      GMM_ASSERT1(R1_ != 0.0, "You have to specify R1");
      double a = R0_ / R1_;
      return (a < 0) ? (a*0.5 + sqrt(a*a*0.25-a)) : a*0.5;
    }
    virtual bool is_converged(double r, double R1 = 0.0) {
      conv_r = r;
      R1_ = R1; return ((gmm::abs(R1_) < gmm::abs(R0_*0.5)) || it >= itmax);
    }
    quadratic_newton_line_search(size_t imax = size_t(-1))  { itmax = imax; }
  };


  struct simplest_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r, double = 0.0) {
      conv_r = r;
      return ((it <= 1 && r < first_res)
	      || (r <= first_res * alpha_max_ratio)
	      || (conv_alpha <= alpha_min)
	      || it >= itmax);
    }
    simplest_newton_line_search
    (size_t imax = size_t(-1), double a_max_ratio = 6.0/5.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio), alpha_min(a_min)
      { itmax = imax; }
  };

  struct default_newton_line_search : public gmm::abstract_newton_line_search {
    // This line search try to detect where is the minimum, dividing the step
    // by a factor two each time.
    //    - it stops if the residual is less than the previous residual
    //      times alpha_min_ratio (= 0.9).
    //    - if the minimal step is reached with a residual greater than
    //      the previous residual times alpha_min_ratio then it decides
    //      between two options :
    //        - return with the step corresponding to the smallest residual
    //        - return with a greater residual corresponding to a residual
    //          less than the previous residual times alpha_max_ratio.
    //      the decision is taken regarding the previous iterations.
    //    - in order to shorten the line search, the process stops when
    //      the residual increases three times consecutively.
    // possible improvment : detect the curvature at the origin
    // (only one evaluation) and take it into account.
    // Fitted to some experiments in contrib/tests_newton

    double alpha, alpha_old, alpha_mult, first_res, alpha_max_ratio;
    double alpha_min_ratio, alpha_min;
    size_type count, count_pat;
    bool max_ratio_reached;
    double alpha_max_ratio_reached, r_max_ratio_reached;
    size_type it_max_ratio_reached;

    virtual void init_search(double r, size_t git, double = 0.0) {
      alpha_min_ratio = 0.9;
      alpha_min = 1e-8;
      alpha_max_ratio = 2.0;
      alpha_mult = 0.25;
      itmax = size_type(-1);
      glob_it = git; if (git <= 1) count_pat = 0;
      conv_alpha = alpha = alpha_old = 1.;
      conv_r = first_res = r; it = 0;
      count = 0;
      max_ratio_reached = false;
    }
    virtual double next_try(void) {
      alpha_old = alpha;
      if (alpha >= 0.4) alpha *= 0.5; else alpha *= alpha_mult; ++it;
      return alpha_old;
    }
    virtual bool is_converged(double r, double = 0.0) {
      // cout << "r = " << r << " alpha = " << alpha / alpha_mult << " count_pat = " << count_pat << endl;
      if (!max_ratio_reached && r < first_res * alpha_max_ratio) {
	alpha_max_ratio_reached = alpha_old; r_max_ratio_reached = r;
	it_max_ratio_reached = it; max_ratio_reached = true; 
      }
      if (max_ratio_reached && r < r_max_ratio_reached * 0.8
	  && r > first_res * 1.1 && it <= it_max_ratio_reached+1) {
	alpha_max_ratio_reached = alpha_old; r_max_ratio_reached = r;
	it_max_ratio_reached = it;
      }
      if (count == 0 || r < conv_r)
	{ conv_r = r; conv_alpha = alpha_old; count = 1; }
      if (conv_r < first_res) ++count;

      if (r < first_res *  alpha_min_ratio)
	{ count_pat = 0.; return true; }      
      if (count >= 5 || (alpha < alpha_min && max_ratio_reached)) {
	// double e = gmm::random() * 50.;
	if (conv_r < first_res * 0.999) count_pat = 0;
	if (/*e < -log(conv_alpha)-4.0 ||*/ count_pat >= 3)
	  { conv_r=r_max_ratio_reached; conv_alpha=alpha_max_ratio_reached; }
	if (conv_r >= first_res * 0.9999) count_pat++;
	return true;
      }
      return false;
    }
    default_newton_line_search(void) { count_pat = 0; }
  };

  /* the former default_newton_line_search */
  struct basic_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio;
    double alpha_min, prev_res, alpha_max_augment;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r, double = 0.0) {
      if (glob_it == 0 || (r < first_res / double(2))
	  || (conv_alpha <= alpha_min && r < first_res * alpha_max_augment)
	  || it >= itmax)
	{ conv_r = r; return true; }
      if (it > 1 && r > prev_res && prev_res < alpha_max_ratio * first_res)
	return true;
      prev_res = conv_r = r;
      return false;
    }
    basic_newton_line_search
    (size_t imax = size_t(-1),
     double a_max_ratio = 5.0/3.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0, double a_augm = 2.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio),
	alpha_min(a_min), alpha_max_augment(a_augm) { itmax = imax; }
  };


  struct systematic_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res;
    double alpha_min, prev_res;
    bool first;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0; first = true;
    }
    virtual double next_try(void)
    { double a = alpha; alpha *= alpha_mult; ++it; return a; }
    virtual bool is_converged(double r, double = 0.0) {
      // cout << "a = " << alpha / alpha_mult << " r = " << r << endl;
      if (r < conv_r || first)
	{ conv_r = r; conv_alpha = alpha / alpha_mult; first = false; }
      if ((alpha <= alpha_min*alpha_mult) || it >= itmax) return true;
      return false;
    }
    systematic_newton_line_search
    (size_t imax = size_t(-1),
     double a_min = 1.0/10000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_min(a_min)  { itmax = imax; }
  };

}


#endif
