/* -*- c++ -*- (enables emacs c++ mode)                                    */
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
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
/* *********************************************************************** */
/*                                                                         */
/* Library : Generic Matrix Methods  (gmm)                                 */
/* File    : gmm_precond_mr_approx_inverse.h : modified version from I.T.L.*/
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
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
#ifndef GMM_PRECOND_MR_APPROX_INVERSE_H
#define GMM_PRECOND_MR_APPROX_INVERSE_H

/* Approximate inverse via MR iteration
 * see P301 of Saad book
*/

#include <gmm_solvers.h>
#include <gmm_tri_solve.h>
#include <gmm_interface.h>
#include <gmm_matrix.h>

namespace gmm {

  template <typename Matrix>
  struct mr_approx_inverse_precond {

    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef typename principal_orientation_type<typename
      linalg_traits<Matrix>::sub_orientation>::potype sub_orientation;
    typedef wsvector<value_type> VVector;
    typedef col_matrix<VVector> MMatrix;

    MMatrix M;

    mr_approx_inverse_precond(const Matrix& A, size_type nb_it,
			      value_type threshold);
  };

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const mr_approx_inverse_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { mult(P.M, v1, v2); }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const mr_approx_inverse_precond<Matrix>& P,
		       const V1 &v1,V2 &v2)
  { mult(gmm::transposed(P.M), v1, v2); }

  template <typename Matrix>
  mr_approx_inverse_precond<Matrix>::mr_approx_inverse_precond(const Matrix& A,
							 size_type nb_it,
							 value_type threshold)
    : M(mat_nrows(A), mat_ncols(A)) {

    VVector m(mat_ncols(A)),r(mat_ncols(A)),ei(mat_ncols(A)),Ar(mat_ncols(A)); 
    value_type alpha = gmm::mat_trace(A) / gmm::sqr(gmm::mat_norm2(A));
    
    for (size_type i = 0; i < mat_nrows(A); ++i) {
      gmm::clear(m); gmm::clear(ei); 
      m[i] = alpha;
      ei[i] = value_type(1);
      
      for (size_type j = 0; j < nb_it; ++j) {
	gmm::mult(A, gmm::scaled(m, -1.0), r);
	gmm::add(ei, r);
	gmm::mult(A, r, Ar);
	gmm::add(gmm::scaled(r, vect_sp(r, Ar) / vect_sp(Ar, Ar)), m);
	gmm::clean(m, threshold * gmm::vect_norm2(m));
      }
      gmm::copy(m, M.col(i));
    }
  }
}

#endif 

