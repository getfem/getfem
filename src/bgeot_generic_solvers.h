/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_generic_solver.h : generic algorithms on linear        */
/*                                      algebra                            */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __BGEOT_GENERIC_SOLVER_H
#define __BGEOT_GENERIC_SOLVER_H

#include <bgeot_abstract_linalg.h>

namespace bgeot {

  /* ******************************************************************** */
  /*		conjugate gradient  (unpreconditionned)        		  */
  /* ******************************************************************** */
  // Inspired from I.T.L. (http://www.osl.iu.edu/research/itl)

  template < class Matrix, class Vector>
  int cg(const Matrix& A, Vector& x, const Vector& b, int itemax, 
	 double residu, bool noisy = true) {
    typename linalg_traits<Vector>::base_type rho(0), rho_1(0), a(0), beta(0);
    Vector p(x.size()), q(x.size()), r(x.size());
    int iter = 0;
    mult(A, scaled(x, -1.0), b, r);

    rho = vect_sp(r,r);
    
    while (sqrt(rho) > residu) {

      if (iter == 0) copy(r, p);		  
      else { beta = rho / rho_1; add(r, scaled(p, beta), p); }

      mult(A, p, q);

      a = rho / vect_sp(p, q);
      
      add(scaled(p, a), x);
      add(scaled(q, -a), r);

      rho_1 = rho; rho = vect_sp(r, r);

      if (++iter >= itemax) return 1;
      if (noisy) cout << "iter " << iter << " residu " << sqrt(rho) << endl;
    }
    return 0;
  }
  
}


#endif //  __BGEOT_GENERIC_SOLVER_H
