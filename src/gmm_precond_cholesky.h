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
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_cholesky.h : modified version from I.T.L.        */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
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
#ifndef GMM_PRECOND_CHOLESKY_H
#define GMM_PRECOND_CHOLESKY_H

// Incomplete Level 0 Cholesky Preconditioner.
// For use with symmetric sparse matrices.
//
// Notes: The idea under a concrete Preconditioner such 
//        as Incomplete Cholesky is to create a Preconditioner
//        object to use in iterative methods. 
//

#include <gmm_solvers.h>
#include <gmm_tri_solve.h>
#include <gmm_interface.h>

namespace gmm {

  template <class Matrix>
  class cholesky_precond {

  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef csr_matrix_ref<value_type *, size_type *, size_type *, 0> tm_type;

    tm_type trimatrix;
  protected :
    std::vector<value_type> Tri_val;
    std::vector<size_type> Tri_ind, Tri_ptr;
 
    template<class M> void do_cholesky(const M& A, row_major);
    void do_cholesky(const Matrix& A, col_major);

  public:

    size_type nrows(void) const      { return mat_nrows(trimatrix); }
    size_type ncols(void) const { return mat_ncols(trimatrix); }
    
    cholesky_precond(const Matrix& A) : Tri_ptr(mat_nrows(A)+1) { 
      if (!is_sparse(A))
	DAL_THROW(failure_error,
		  "Matrix should be sparse for incomplete cholesky");
      do_cholesky(A, typename principal_orientation_type<typename
		  linalg_traits<Matrix>::sub_orientation>::potype());
    }
    cholesky_precond(void) {}
  };

  template <class Matrix> template<class M>
  void cholesky_precond<Matrix>::do_cholesky(const M& A, row_major) {
    size_type Tri_loc = 0, n = mat_nrows(A), d, g, h, i, j, k;
    value_type z;
    Tri_ptr[0] = 0;
    
    for (int count = 0; count < 2; ++count) {
      if (count) { Tri_val.resize(Tri_loc); Tri_ind.resize(Tri_loc); }
      for (Tri_loc = 0, i = 0; i < n; ++i) {
	typedef typename linalg_traits<M>::const_sub_row_type row_type;
	row_type row = mat_const_row(A, i);
        typename linalg_traits<row_type>::const_iterator
	  it = vect_const_begin(row), ite = vect_const_end(row);
	
	for (; it != ite; ++it)
	  if (it.index() >= i) {
	    if (count) { Tri_val[Tri_loc] = *it; Tri_ind[Tri_loc]=it.index(); }
	    ++Tri_loc;
	  }
	Tri_ptr[i+1] = Tri_loc;
      }
    }
    
    for (k = 0; k < n; k++) {
      d = Tri_ptr[k];
      if (Tri_val[d] < 0) { 
	DAL_WARNING(2, "negative value found");
	z = Tri_val[d] = 1.0;
      }
      else 
	z = Tri_val[d] = sqrt(Tri_val[d]);

      for (i = d + 1; i < Tri_ptr[k+1]; i++) Tri_val[i] /= z;
      for (i = d + 1; i < Tri_ptr[k+1]; i++) {
	z = Tri_val[i];
	h = Tri_ind[i];
	g = i;
	
	for (j = Tri_ptr[h] ; j < Tri_ptr[h+1]; j++)
	  for ( ; g < Tri_ptr[k+1] && Tri_ind[g+1] <= Tri_ind[j]; g++)
	    if (Tri_ind[g] == Tri_ind[j]) Tri_val[j] -= z * Tri_val[g];
      }
    }
    trimatrix = tm_type(&(Tri_val[0]), &(Tri_ind[0]), &(Tri_ptr[0]),
			n, mat_ncols(A));
  }
  
  template <class Matrix>
  void cholesky_precond<Matrix>::do_cholesky(const Matrix& A, col_major)
  { do_cholesky(gmm::transposed(A), row_major()); }

  template <class Matrix, class V1, class V2> inline
  void mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    gmm::lower_tri_solve(gmm::transposed(P.trimatrix), v2);
    gmm::upper_tri_solve(P.trimatrix, v2);
  }

  template <class Matrix, class V1, class V2> inline
  void transposed_mult(const cholesky_precond<Matrix>& P,const V1 &v1,V2 &v2)
  { mult(P, v1, v2); }

  template <class Matrix, class V1, class V2> inline
  void left_mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::transposed(P.trimatrix), v2); }

  template <class Matrix, class V1, class V2> inline
  void right_mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.trimatrix, v2);  }

  template <class Matrix, class V1, class V2> inline
  void transposed_left_mult(const cholesky_precond<Matrix>& P, const V1 &v1,
			    V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.trimatrix, v2); }

  template <class Matrix, class V1, class V2> inline
  void transposed_right_mult(const cholesky_precond<Matrix>& P, const V1 &v1,
			     V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::transposed(P.trimatrix), v2); }

}

#endif 

