// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_Newton.h : Some tools used for the Newton algorithms.
// Date    : January 19, 2006.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//


/**
   @file getfem_Newton.h
   @brief Tools for Newton algorithms
   @see getfem_mode_solvers.h
*/

#ifndef GETFEM_NEWTON_H__
#define GETFEM_NEWTON_H__

#include <getfem_config.h>

namespace getfem {


  struct abstract_newton_line_search {
    scalar_type conv_alpha, conv_r;
    virtual void init_search(scalar_type r) = 0;
    virtual scalar_type next_try(void) = 0;
    virtual bool is_converged(scalar_type) = 0;
    virtual scalar_type converged_value(void) { return conv_alpha; };
    virtual scalar_type converged_residual(void) { return conv_r; };
  };


  struct basic_newton_line_search : public abstract_newton_line_search {
    scalar_type alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    size_type k, itmax;
    virtual void init_search(scalar_type r)
    { conv_alpha = alpha = scalar_type(1); first_res = r; k = 0; }
    virtual scalar_type next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++k; return conv_alpha; }
    virtual bool is_converged(scalar_type r) {
      conv_r = r;
      return ((k <= 1 && r < first_res)
	      || (r <= first_res * alpha_max_ratio)
	      || (conv_alpha <= alpha_min)
	      || k >= itmax);
    }
    basic_newton_line_search
    (size_type imax = size_type(-1),
     scalar_type a_max_ratio = scalar_type(4)/scalar_type(3),
     scalar_type a_min = scalar_type(1)/scalar_type(1000),
     scalar_type a_mult = scalar_type(3)/scalar_type(5))
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio), alpha_min(a_min),
	itmax(imax) {}
  };

  struct basic_newton_line_search2 : public abstract_newton_line_search {
    scalar_type alpha, alpha_mult, first_res, alpha_max_ratio;
    scalar_type alpha_min, prev_res;
    size_type k, itmax;
    virtual void init_search(scalar_type r) {
      conv_alpha = alpha = (gmm::random() < 0.1)
	? scalar_type(1)/scalar_type(2) : scalar_type(1);
      first_res = r; k = 0;
    }
    virtual scalar_type next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++k; return conv_alpha; }
    virtual bool is_converged(scalar_type r) {
      if ((r < first_res / scalar_type(2)) || (conv_alpha <= alpha_min)
	  || k >= itmax)
	{ conv_r = r; return true; }
      if (k > 1 && r > prev_res && prev_res < alpha_max_ratio * first_res)
	return true;
      conv_r = r;
      return false;
    }
    basic_newton_line_search2
    (size_type imax = size_type(-1),
     scalar_type a_max_ratio = scalar_type(4)/scalar_type(3),
     scalar_type a_min = scalar_type(1)/scalar_type(1000),
     scalar_type a_mult = scalar_type(3)/scalar_type(5))
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio),
	alpha_min(a_min), itmax(imax)  {}
  };

  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NEWTON_H__  */
