// -*- c++ -*-
//
// Copyright 1997, 1998, 1999 University of Notre Dame.
// Authors: Andrew Lumsdaine, Jeremy G. Siek, Lie-Quan Lee
//
// You should have received a copy of the License Agreement for the
// Matrix Template Library along with the software;  see the
// file LICENSE.  If not, contact Office of Research, University of Notre
// Dame, Notre Dame, IN  46556.
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
//
//===========================================================================
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_lu.h : modified version from M.T.L.                */
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

#ifndef GMM_DENSE_LU_H
#define GMM_DENSE_LU_H

#include "gmm_solvers.h"

namespace gmm {

  // rank one update for dense matrices

  template <class Matrix, class VecX, class VecY>
  inline void rank_one_update(Matrix &A, const VecX& x,
			      const VecY& y, row_major) {
    size_type N = mat_nrows(A);
    if (N > vect_size(x) || mat_ncols(A) > vect_size(y))
      DAL_THROW(dimension_error,"dimensions mismatch");
    typename linalg_traits<VecX>::const_iterator itx = vect_const_begin(x);
    for (size_type i = 0; i < N; ++i, ++itx) {
      typedef typename linalg_traits<Matrix>::sub_row_type row_type;
      row_type row = mat_row(A, i);
      typename linalg_traits<row_type>::iterator
	it = vect_begin(row), ite = vect_end(row);
      typename linalg_traits<VecY>::const_iterator ity = vect_const_begin(y);
#   ifdef USING_BROKEN_GCC295
      typedef typename linalg_traits<Matrix>::value_type T;
      for (; it != ite; ++it, ++ity)
	const_cast<T &>(*it) += conj_product(*itx, *ity);
#   else
      for (; it != ite; ++it, ++ity) *it += conj_product(*itx, *ity);
#   endif
    }
  }

  template <class Matrix, class VecX, class VecY>
  inline void rank_one_update(Matrix &A, const VecX& x,
			      const VecY& y, col_major) {
    size_type M = mat_ncols(A);
    if (mat_nrows(A) > vect_size(x) || M > vect_size(y))
      DAL_THROW(dimension_error,"dimensions mismatch");
    typename linalg_traits<VecY>::const_iterator ity = vect_const_begin(y);
    for (size_type i = 0; i < M; ++i, ++ity) {
      typedef typename linalg_traits<Matrix>::sub_col_type col_type;
      col_type col = mat_col(A, i);
      typename linalg_traits<col_type>::iterator
	it = vect_begin(col), ite = vect_end(col);
      typename linalg_traits<VecX>::const_iterator itx = vect_const_begin(x);
#   ifdef USING_BROKEN_GCC295
      typedef typename linalg_traits<Matrix>::value_type T;
      for (; it != ite; ++it, ++itx)
	const_cast<T &>(*it) += conj_product(*itx, *ity);
#   else
      for (; it != ite; ++it, ++itx) *it += conj_product(*itx, *ity);
#   endif
    }
  }
  
  template <class Matrix, class VecX, class VecY>
  inline void rank_one_update(const Matrix &A_, const VecX& x,
			      const VecY& y) {
    Matrix& A = const_cast<Matrix&>(A_);
    if (is_sparse(A))
      DAL_THROW(failure_error,
		"Sorry, rank one update for sparse matrices does not exist");
    rank_one_update(A, x, y, typename principal_orientation_type<typename
		    linalg_traits<Matrix>::sub_orientation>::potype());
  }


  // LU Factorization of a general (dense) matrix
  //
  // This is the outer product (a level-2 operation) form of the LU
  // Factorization with pivoting algorithm . This is equivalent to
  // LAPACK's dgetf2. Also see "Matrix Computations" 3rd Ed.  by Golub
  // and Van Loan section 3.2.5 and especially page 115.
  // 
  // The pivot indices in ipvt are indexed starting from 1
  // so that this is compatible with LAPACK (Fortran).
  //
  template <class DenseMatrix, class Pvector>
  size_type lu_factor(DenseMatrix& A, Pvector& ipvt) {
    typedef typename linalg_traits<DenseMatrix>::value_type value_type;
    size_type info = 0, i, j, jp, M = A.nrows(), N = A.ncols();
    std::vector<value_type> c(M), r(N);
      
    if (M || N) {
      for (j = 0; j < std::min(M-1, N-1); ++j) {
	double max = dal::abs(A(j,j)); jp = j;
	for (i = j+1; i < M; ++i)		   /* find pivot.          */
	  if (dal::abs(A(i,j)) > max) { jp = i; max = dal::abs(A(i,j)); }
	ipvt[j] = jp + 1;
	
	if (A(jp, j) == value_type(0)) { info = j + 1; break; }
	if (jp != j) for (i = 0; i < N; ++i) std::swap(A(jp, i), A(j, i));
	
	for (i = j+1; i < M; ++i) { A(i, j) /= A(j,j); c[i-j-1] = -A(i, j); }
	for (i = j+1; i < N; ++i) r[i-j-1] = A(j, i);  // avoid the copy ?
	rank_one_update(sub_matrix(A, sub_interval(j+1, M-j-1),
				   sub_interval(j+1, N-j-1)), c, r);
      }
      ipvt[j] = j + 1;
    }
    
    return info;
  }
  
  //  LU Solve : Solve equation Ax=b, given an LU factored matrix.
  //  Thanks to Valient Gough for this routine!
  //
  template <class DenseMatrix, class VectorB, class VectorX, class Pvector>
  void lu_solve(const DenseMatrix &LU, const Pvector& pvector, 
		VectorX &x, const VectorB &b) {
    copy(b, x);
    /* use the permutation vector to modify the starting vector            */
    /*  to account for the permutations in LU.                             */
    for(size_type i = 0; i < pvector.size(); ++i) {
      size_type perm = pvector[i]-1;     // permutations stored in 1's offset
      if(i != perm) std::swap(x[i], x[perm]);
    }
    /* solve  Ax = b  ->  LUx = b  ->  Ux = L^-1 b.                        */
    lower_tri_solve(LU, x, true);
    upper_tri_solve(LU, x, false);
  }

  template <class DenseMatrix, class VectorB, class VectorX>
  void lu_solve(const DenseMatrix &A, VectorX &x, const VectorB &b) {
    DenseMatrix B(mat_nrows(A), mat_ncols(A));
    std::vector<size_type> ipvt(mat_nrows(A));
    gmm::copy(A, B);
    lu_factor(B, ipvt);
    lu_solve(B, ipvt, x, b);
  }
  
  template <class DenseMatrix, class VectorB, class VectorX, class Pvector>
  void lu_solve_transposed(const DenseMatrix &LU, const Pvector& pvector, 
			   VectorX &x, const VectorB &b) {
    copy(b, x);
    /* use the permutation vector to modify the starting vector            */
    /*  to account for the permutations in LU.                             */
    for(size_type i = 0; i < pvector.size(); ++i) {
      size_type perm = pvector[i]-1;       // permutations stored in 1's offset
      if(i != perm) std::swap(x[i], x[perm]);
    }
    /* solve  Ax = b  ->  LUx = b  ->  Ux = L^-1 b.                        */
    lower_tri_solve(transposed(LU), x, false);
    upper_tri_solve(transposed(LU), x, true);
  }


  // LU Inverse : Given an LU factored matrix, construct the inverse 
  //              of the matrix.
  //  Thanks to Valient Gough for this routine!
  template <class DenseMatrixLU, class DenseMatrix, class Pvector>
  void lu_inverse(const DenseMatrixLU& LU, const Pvector& pvector,
		  DenseMatrix& AInv, col_major) {
    typedef typename linalg_traits<DenseMatrixLU>::value_type value_type;
    std::vector<value_type> tmp(pvector.size(), value_type(0));
    std::vector<value_type> result(pvector.size());
    for(size_type i = 0; i < pvector.size(); ++i) {
      tmp[i] = 1.0;
      lu_solve(LU, pvector, result, tmp);
      copy(result, mat_col(AInv, i));
      tmp[i] = 0.0;
    }
  }

  template <class DenseMatrixLU, class DenseMatrix, class Pvector>
  void lu_inverse(const DenseMatrixLU& LU, const Pvector& pvector,
		  DenseMatrix& AInv, row_major) {
    typedef typename linalg_traits<DenseMatrixLU>::value_type value_type;
    std::vector<value_type> tmp(pvector.size(), value_type(0));
    std::vector<value_type> result(pvector.size());
    for(size_type i = 0; i < pvector.size(); ++i) {
      tmp[i] = 1.0;
      lu_solve_transposed(LU, pvector, result, tmp);
      copy(result, mat_row(AInv, i));
      tmp[i] = 0.0;
    }
  }
  
  template <class DenseMatrixLU, class DenseMatrix, class Pvector>
  void lu_inverse(const DenseMatrixLU& LU, const Pvector& pvector,
		  const DenseMatrix& AInv_) {
    DenseMatrix& AInv = const_cast<DenseMatrix&>(AInv_);
    lu_inverse(LU, pvector, AInv, typename principal_orientation_type<typename
	       linalg_traits<DenseMatrix>::sub_orientation>::potype());
  }

  template <class DenseMatrix>
  typename linalg_traits<DenseMatrix>::value_type
  lu_inverse(const DenseMatrix& A_) {
    DenseMatrix& A = const_cast<DenseMatrix&>(A_);
    DenseMatrix B(mat_nrows(A), mat_ncols(A));
    std::vector<size_type> ipvt(mat_nrows(A));
    gmm::copy(A, B);
    if (lu_factor(B, ipvt)) DAL_THROW(failure_error, "Non invertible matrix");
    lu_inverse(B, ipvt, A);
    return lu_det(B, ipvt);
  }

  template <class DenseMatrixLU, class Pvector>
  typename linalg_traits<DenseMatrixLU>::value_type
  lu_det(const DenseMatrixLU& LU, const Pvector&) {
    typename linalg_traits<DenseMatrixLU>::value_type det(1);
    for (size_type j = 0; j < std::min(mat_nrows(LU), mat_ncols(LU)); ++j)
      det *= LU(j,j);
    return det;
  }

  template <class DenseMatrix>
  typename linalg_traits<DenseMatrix>::value_type
  lu_det(const DenseMatrix& A) {
    DenseMatrix B(mat_nrows(A), mat_ncols(A));
    std::vector<size_type> ipvt(mat_nrows(A));
    gmm::copy(A, B);
    lu_factor(B, ipvt);
    return lu_det(B, ipvt);
  }

}

#endif

