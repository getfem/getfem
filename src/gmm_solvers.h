/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solvers.h : generic solvers.                             */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GMM++                                            */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#ifndef __GMM_SOLVERS_H
#define __GMM_SOLVERS_H

#include <gmm_iter.h>


namespace gmm {

  // Needed by ilut and choleskyt
  template<class T> struct _elt_rsvector_value_less {
    inline bool operator()(const _elt_rsvector<T>& a, 
			   const _elt_rsvector<T>& b) const
    { return (dal::abs(a.e) > dal::abs(b.e)); }
  };


  /** mixed method to find a zero of a real function G, a priori 
   * between a and b. If the zero is not between a and b, iterations
   * of secant are applied. When a convenient interval is found,
   * iterations of dichotomie and regula falsi are applied.
   */
  template <class FUNC, class T>
  T find_root(const FUNC &G, T a = 0.0, T b = a+1.0,
	      double tol = gmm::default_tol(T())) {
    
    T c, Ga = G(a), Gb = G(b), Gc, d;
    d = abs(b - a);
    for (int i = 0; i < 4; i++) { /* secant iterations.                   */
      if (d < tol) return (b + a) / 2.0;
      c = b - Gb * (b - a) / (Gb - Ga); Gc = G(c);
      a = b; b = c; Ga = Gb; Gb = Gc;
      d = abs(b - a);
    }
    while (Ga * Gb > 0.0) { /* secant iterations.                         */
      if (d < tol) return (b + a) / 2.0;
      c = b - Gb * (b - a) / (Gb - Ga); Gc = G(c);
      a = b; b = c; Ga = Gb; Gb = Gc;
      d = abs(b - a);
    }
    
    c = max(a, b); a = min(a, b); b = c;
    while (d > tol) {
      c = b - (b - a) * (Gb / (Gb - Ga)); Gc = G(c); /* regula falsi.     */
      if (c > b) c = b; if (c < a) c = a;
      if (Gc > 0) { b = c; Gb = Gc; } else { a = c; Ga = Gc; }
      c = (b + a) / 2.0 ; Gc = G(c); /* Dichotomie.                       */
      if (Gc > 0) { b = c; Gb = Gc; } else { a = c; Ga = Gc; }
      d = abs(b - a); c = (b + a) / 2.0; if ((c == a) || (c == b)) d = 0.0;
    }
    return (b + a) / 2.0;
  }
  
}

#include <gmm_precond_diagonal.h>
#include <gmm_precond_cholesky.h>
#include <gmm_precond_choleskyt.h>
#include <gmm_precond_mr_approx_inverse.h>
#include <gmm_precond_ilu.h>
#include <gmm_precond_ilut.h>



#include <gmm_solver_cg.h>
#include <gmm_solver_bicgstab.h>
#include <gmm_solver_qmr.h>
#include <gmm_solver_cheby.h>
#include <gmm_solver_constrained_cg.h>
#include <gmm_solver_Schwarz_additive.h>
#include <gmm_modified_gram_schmidt.h>
#include <gmm_tri_solve.h>
#include <gmm_solver_gmres.h>
// #include <gmm_solver_idgmres.h>
#include <gmm_superlu_interface.h>



#endif //  __GMM_SOLVERS_H
