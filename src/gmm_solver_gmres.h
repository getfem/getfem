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

#ifndef GMM_KRYLOV_GMRES_H
#define GMM_KRYLOV_GMRES_H

#include <gmm_solvers.h>

namespace gmm {

  // Generalized Minimum Residual
  //
  //   This solve the unsymmetric linear system Ax = b using restarted GMRES.
  //
  //   See: Y. Saad and M. Schulter. GMRES: A generalized minimum residual
  //   algorithm for solving nonsysmmetric linear systems, SIAM
  //   J. Sci. Statist. Comp.  7(1986), pp, 856-869
  //

  template < class Matrix, class Vector, class VectorB, class Preconditioner, 
	     class Basis >
  void gmres(const Matrix &A, Vector &x, const VectorB &b,
	     const Preconditioner &M, int restart,
	     iteration &outer, Basis& KS) {

    typedef typename linalg_traits<Vector>::value_type T;
    typedef size_t size_type;
    
    typedef std::vector<T> TmpVec;
    typedef typename temporary_vector<Vector>::vector_type internal_vector;
    
    internal_vector w(vect_size(x)), r(vect_size(x)), u(vect_size(x));
    
    typedef typename number_traits<T>::magnitude_type Real;
    typedef gmm::row_matrix<TmpVec> HMat;

    HMat H(restart+1, restart);
    TmpVec s(restart+1);
    outer.set_rhsnorm(gmm::vect_norm2(b));

    if (outerter.get_rhsnorm() == 0.0) { clear(x); return; }
    
    std::vector< gmm::givens_rotation<T> > rotations(restart+1);
    
    mult(A, scaled(x, -1.0), b, w);
    
    mult(M, w, r);
    Real beta = dal::abs(gmm::vect_norm2(r));
    
    while (! outer.finished(beta)) {
      
      gmm::copy(gmm::scaled(r, 1./beta), KS[0]);
      gmm::clear(s);
      s[0] = beta;
      
      size_type i = 0;
      iteration inner = outer;
      inner.reduce_noisy();
      inner.set_maxiter(restart);
      inner.init();
      
      do {
	size_type k;
	gmm::mult(A, KS[i], u);
	gmm::mult(M, u, KS[i+1]);
	
	orthogonalize(KS, mat_row(H, i), i);
	
	Real H_ip1_i = gmm::vect_norm2(KS[i+1]);
	H(i+1, i) = H_ip1_i;
	gmm::scale(KS[i+1], 1./H_ip1_i);
	
	for (k = 0; k < i; k++)
 	  rotations[k].scalar_apply(H(k,i), H(k+1,i));
	
	rotations[i] = givens_rotation<T>(H(i,i), H(i+1,i));
	rotations[i].scalar_apply(H(i,i), H(i+1,i));
	rotations[i].scalar_apply(s[i], s[i+1]);
	
	++inner, ++outer, ++i;
	
      } while (! inner.finished(dal::abs(s[i]))); 
      
      // i ou i+1 ?
      gmm::upper_tri_solve(gmm::sub_matrix(H, gmm::sub_interval(0, i+1),
					   gmm::sub_interval(0, i+1)), s);

      gmm::combine(KS, s, x, i);
      
      gmm::mult(A, gmm::scaled(x, -1.0), b, w);
      gmm::mult(M, w, r);
      beta = dal::abs(gmm::vect_norm2(r));
    }
  }


  template < class Matrix, class Vector, class VectorB, class Preconditioner >
  void gmres(const Matrix &A, Vector &x, const VectorB &b,
	     const Preconditioner &M, int restart, iteration& outer) {
    modified_gram_schmidt<Vector> orth(restart, vect_size(x));
    gmres(A, x, b, M, restart, outer, orth); 
  }

}

#endif
