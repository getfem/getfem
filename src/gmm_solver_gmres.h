// -*- c++ -*-
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_gmres.h : from I.T.L.                                    */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
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

#ifndef ITL_KRYLOV_GMRES_H
#define ITL_KRYLOV_GMRES_H

#include <vector>
#include <algorithm>
#include "itl/itl.h"
#include "itl/givens_rotation.h"
#include "itl/number_traits.h"
#include <itl/modified_gram_schmidt.h>

namespace itl {

  // Generalized Minimum Residual
  //
  //   This solve the unsymmetric linear system Ax = b using restarted GMRES.
  //
  //   A return value of 0 indicates convergence within the
  //   maximum number of iterations (determined by the iter object).
  //   A return value of 1 indicates a failure to converge.
  //
  //   On instantiating Iteration object outer, the first parameter
  //   Vector w must be the precondtioned one, i.e., solve(M, b, w) where
  //   b is right side of linear system. See test_gmres.cc for example.
  //
  //   See: Y. Saad and M. Schulter. GMRES: A generalized minimum residual
  //   algorithm for solving nonsysmmetric linear systems, SIAM
  //   J. Sci. Statist. Comp.  7(1986), pp, 856-869
  //
  /* required operations: mult,copy,dot_conj,add,scaled,two_norm,tri_solve */
template < class Matrix, class Vector, class VectorB, class Preconditioner, 
	   class Iter, class Basis >
int 
gmres(const Matrix &A, Vector &x, const VectorB &b,
      const Preconditioner &M, int restart, Iter& outer, Basis& KS)
{
  typedef typename itl_traits<Vector>::value_type T;
  typedef typename itl_traits<Vector>::size_type size_type;

  typedef std::vector<T> TmpVec;
  typedef Vector internal_vector;

  //These must be arithmetically compatible with x
  internal_vector w(size(x)), r(size(x)), u(size(x));

  typedef typename number_traits<T>::magnitude_type Real;
  typedef typename internal_matrix_traits<T>::Matrix HMat;
  HMat H(restart+1, restart); //Elements in H must be real type

  TmpVec s(restart+1);
    
  std::vector< itl::givens_rotation<T> > rotations(restart+1);

  itl::mult(A, itl::scaled(x, -1.0), b, w);

  itl::solve(M, w, r);
  Real beta = std::abs(itl::two_norm(r));

  while (! outer.finished(beta)) {

    itl::copy(itl::scaled(r, 1./beta), KS[0]);
    std::fill(s.begin(), s.end(), 0.0);
    s[0] = beta;

    size_type i = 0;
    Iter inner(outer.normb(), restart, outer.tol(), outer.atol());

    do {
      size_type k;
      itl::mult(A, KS[i], u);
      itl::solve(M, u, KS[i+1]);

      itl::orthogonalize(KS, H[i], i);

      Real H_ip1_i = itl::two_norm(KS[i+1]);
      H(i+1, i) = H_ip1_i;
      itl::scale(KS[i+1], 1./H_ip1_i);

      for (k = 0; k < i; k++)
        rotations[k].scalar_apply(H(k,i), H(k+1,i));

      rotations[i] = itl::givens_rotation<T>(H(i,i), H(i+1,i));
      rotations[i].scalar_apply(H(i,i), H(i+1,i));
      rotations[i].scalar_apply(s[i], s[i+1]);
      
      ++inner, ++outer, ++i;

    } while (! inner.finished(std::abs(s[i]))); 

    itl::upper_tri_solve(H, s, i);

    itl::combine(KS, s, x, i);

    itl::mult(A, itl::scaled(x, -1.0), b, w);
    itl::solve(M, w, r);
    beta = std::abs(itl::two_norm(r));
  }
  
  return outer.error_code();
}


template < class Matrix, class Vector, class VectorB, class Preconditioner, 
	   class Iter >
int 
gmres(const Matrix &A, Vector &x, const VectorB &b,
      const Preconditioner &M, int restart, Iter& outer)
{
  itl::modified_gram_schmidt<Vector> orth(restart, size(x));
  return itl::gmres(A, x, b, M, restart, outer, orth); 
}

}

#endif
