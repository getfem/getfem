/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_cg.h : conjugated gradient.                       */
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


#ifndef __GMM_SOLVER_CG_H
#define __GMM_SOLVER_CG_H

namespace gmm {

  /* ******************************************************************** */
  /*		conjugate gradient                           		  */
  /* (preconditionned, with parametrable scalar product)        	  */
  /* ******************************************************************** */
  // Inspired from I.T.L. (http://www.osl.iu.edu/research/itl)

  template <class Matrix, class Matps, class Precond, 
            class Vector1, class Vector2>
  int cg(const Matrix& A, Vector1& x, const Vector2& b, const Matps& PS,
	 const Precond &P, int itemax, double residu, int noisy = 1) {

    typedef typename temporary_plain_vector<Vector1>::vector_type temp_vector;
    typedef typename linalg_traits<Vector1>::value_type value_type;

    value_type rho(0), rho_1(0), a(0), beta(0), norm_b = vect_norm2(b);
    temp_vector p(vect_size(x)), q(vect_size(x)), r(vect_size(x)),
      z(vect_size(x));
    int iter = 0;

    if (norm_b == value_type(0))
      clear(x);
    else {
      mult(A, scaled(x, -1.0), b, r);
      mult(P, r, z);
      rho = vect_sp(PS, r, z);
      
      while (sqrt(modulus(rho)) > residu * norm_b) {
	
	if (iter == 0) copy(r, p);		  
	else { beta = rho / rho_1; add(r, scaled(p, beta), p); }
	
	mult(A, p, q);
	
	a = rho / vect_sp(PS, p, q);
	
	add(scaled(p, a), x);
	add(scaled(q, -a), r);
	
	mult(P, r, z);
	rho_1 = rho;
	rho = vect_sp(PS, r, z);
	
	if (++iter >= itemax) return 1;
	if (noisy > 0)  cout << "iter " << iter << " residu "
			     << sqrt(modulus(rho)) / norm_b << endl;
      }
    }
    return 0;
  }

  template <class Matrix, class Matps, class Precond, 
            class Vector1, class Vector2>
  int cg(const Matrix& A, const Vector1& x, const Vector2& b, const Matps& PS,
	 const Precond &P, int itemax, double residu, int noisy = 1)
  { cg(A, linalg_cast(x), b, PS, P, itemax, residu, noisy); }

  template <class Matrix,  class Vector1, class Vector2> inline
  int cg(const Matrix& A, Vector1& x, const Vector2& b,
	 int itemax, double residu, int noisy = 1) {
    return cg(A, x, b, identity_matrix(), identity_matrix(), itemax,
	      residu, noisy);
  }
  
  template <class Matrix,  class Vector1, class Vector2> inline
  int cg(const Matrix& A, const Vector1& x, const Vector2& b,
	 int itemax, double residu, int noisy = 1) {
    return cg(A, linalg_cast(x), b, identity_matrix(), identity_matrix(),
	      itemax, residu, noisy);
  }
  
}


#endif //  __GMM_SOLVER_CG_H
