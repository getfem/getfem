// -*- c++ -*-
//
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
/* File    :  gmm_solver_cheby.h : from I.T.L.                             */
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
#ifndef GMM_CHEBY_H
#define GMM_CHEBY_H

#include <gmm_solvers.h>


namespace gmm {

  //: Chebyshev Iteration
  //This solves the unsymmetric  linear system Ax = b using 
  //the Preconditioned Chebyshev Method.
  //<p>
  //A return value of 0 indicates convergence within the
  //maximum number of iterations (determined by the iter object).
  //A return value of 1 indicates a failure to converge.
  //<p>
  //See: T. Manteuffel, The Chebyshev iteration for nonsysmmetric 
  //     linear systems
  //Numer. Math. 28(1977), pp. 307-327
  //G. H. Golub and C. F. Van Loan, Matrix Computations, The Johns Hopkins 
  //University Press, Baltimore, Maryland, 1996 
  //
  //!category: itl,algorithms 
  //!component: function 
  //!definition: cheby.h 
  //!tparam: Matrix  - Matrix or multiplier for matrix free methods 
  //!tparam: Vector  - Vector  
  //!tparam: VectorB - Vector 
  //!tparam: Preconditioner -  Incomplete LU, Incomplete LU with threshold,
  //!SSOR or identity_preconditioner. 
  //!tparam: Iteration - Controls the stopping criteria 

template < class Matrix, class Vector, class VectorB, class Preconditioner>
void cheby(const Matrix &A, Vector &x, const VectorB &b,
	   const Preconditioner &M, iteration& iter,
	   typename Vector::value_type eigmin, 
	   typename Vector::value_type eigmax)
{
  typedef typename linalg_traits<Vector>::value_type Real;
  Real alpha, beta, c, d;
  typedef typename temporary_vector<Vector>::vector_type TmpVec;
  TmpVec p(vect_size(x)), q(vect_size(x)), z(vect_size(x)), r(vect_size(x));

  iter.set_rhsnorm(gmm::vect_norm2(b));
  if (iter.get_rhsnorm() == 0.0) { clear(x); return; }

  gmm::mult(A, gmm::scaled(x, -1.0), b, r);

  c = (eigmax - eigmin) / 2.0;
  d = (eigmax + eigmin) / 2.0;

  while ( ! iter.finished_vect(r) ) {
    gmm::mult(M, r, z);         

    if ( iter.first() ) {
      gmm::copy(z, p);          
      alpha = 2.0 / d;
    } else {
      beta = c * alpha / 2.0;    
      beta = beta * beta;
      alpha = 1.0 / (d - beta);  
      gmm::add(z, gmm::scaled(p, beta), p);
    }

    gmm::mult(A, p, q);
    gmm::add(x, gmm::scaled(p, alpha), x);
    gmm::add(r, gmm::scaled(q, -alpha), r);

    ++iter;
  }

}

}

#endif
