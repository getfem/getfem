/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_Householder.h : Householder calculus.              */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Caroline Lecalvez, Caroline.Lecalvez@gmm.insa-tlse.fr          */
/*          Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
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

#ifndef GMM_DENSE_HOUSEHOLDER_H
#define GMM_DENSE_HOUSEHOLDER_H

#include "gmm_solvers.h"

namespace gmm {

  /* ********************************************************************* */
  /*    Rank one update  (complex and real version)                        */
  /* ********************************************************************* */

  template <class Matrix, class VecX, class VecY>
  inline void rank_one_update(Matrix &A, const VecX& x,
			      const VecY& y, row_major) {
    typedef typename linalg_traits<Matrix>::value_type value_type;
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
      value_type tx = *itx;
      for (; it != ite; ++it, ++ity) *it += conj_product(tx, *ity);
#   endif
    }
  }

  template <class Matrix, class VecX, class VecY>
  inline void rank_one_update(Matrix &A, const VecX& x,
			      const VecY& y, col_major) {
    typedef typename linalg_traits<Matrix>::value_type value_type;
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
      value_type ty = *ity;
      for (; it != ite; ++it, ++itx) *it += conj_product(*itx, ty); 
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

  /* ********************************************************************* */
  /*    Householder vector computation (complex and real version)          */
  /* ********************************************************************* */

  template <class VECT> void house_vector(const VECT &VV) {
    VECT &V = const_cast<VECT &>(VV);
    typedef typename linalg_traits<VECT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    
    magnitude_type mu = vect_norm2(V);
    value_type beta;
    if (mu != magnitude_type(0)) {
      if (dal::abs(V[0]) != magnitude_type(0))
	beta = V[0] + mu * V[0] / dal::abs(V[0]);
      else
	beta = mu;
      gmm::scale(V, value_type(1) / beta);
    }
    V[0] = value_type(1);
  }
  
  /* ********************************************************************* */
  /*    Householder updates  (complex and real version)                    */
  /* ********************************************************************* */

  template <class MAT, class VECT1, class VECT2> inline
    void row_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW); MAT &A = const_cast<MAT &>(AA);

    gmm::mult(conjugated(transposed(A)), scaled(V, -2.0/vect_norm2_sqr(V)), W);
    rank_one_update(A, V, W);
  }

  template <class MAT, class VECT1, class VECT2> inline
    void col_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW); MAT &A = const_cast<MAT &>(AA);
    
    gmm::mult(A, scaled(V, -2.0 / vect_norm2_sqr(V)), W);
    rank_one_update(A, W, V);
  }

  /* ********************************************************************* */
  /*    Hessemberg reduction with Householder.                             */
  /* ********************************************************************* */

  template <class MAT1, class MAT2>
    void Hessenberg_reduction(const MAT1& AA, const MAT2 &QQ, bool compute_Q){
    MAT1& A = const_cast<MAT1&>(AA); MAT1& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(A);
    std::vector<value_type> v(n), w(n);
    if (compute_Q) gmm::copy(identity_matrix(), Q);
    for (size_type k = 1; k < n-1; ++k) {
      sub_interval SUBI(k, n-k), SUBJ(k-1,n-k+1), SUBK(0,n);
      v.resize(n-k);
      for (size_type j = k; j < n; ++j) v[j-k] = A(j, k-1);
      house_vector(v);
      row_house_update(sub_matrix(A, SUBI, SUBJ), v, sub_vector(w, SUBJ));
      col_house_update(sub_matrix(A, SUBK, SUBI), v, w);
      if (compute_Q) col_house_update(sub_matrix(Q, SUBK, SUBI), v, w);
    }
  }

  /* ********************************************************************* */
  /*    Givens rotations                                                   */
  /* ********************************************************************* */

  template <class T> void Givens_rotation(T a, T b, T &c, T &s) {
    if (b == T(0)) { c = T(1); s = T(0); return; }
    if (dal::abs(b) > dal::abs(a))
      { T tau = -a/b; s = T(1) / sqrt(1+tau*tau); c = s*tau; }
    else
      { T tau = -b/a; c = T(1) / sqrt(1+tau*tau); s = c*tau; }
  }

  template <class T> inline void Apply_Givens_rotation(T &x, T &y, T c, T s)
  { T t1=x, t2=y; x = c*t1 - s*t2; y = c*t2 + s*t1; }

  template <class T>
  void Givens_rotation(std::complex<T> a, std::complex<T> b,
		       std::complex<T> &c, std::complex<T> &s) {
    T aa = dal::abs(a);
    if (aa == T(0)) { c = std::complex<T>(0); s = std::complex<T>(1); return; }
    T norm = ::sqrt(dal::abs_sqr(a) + dal::abs_sqr(b));
    c = std::complex<T>(aa / norm);
    s = a/aa * std::conj(b)/norm;
  }

  template <class T> inline
  void Apply_Givens_rotation(std::complex<T> &x, std::complex<T> &y,
			     std::complex<T>  c, std::complex<T>  s)
  { std::complex<T> t1=x, t2=y; x = c*t1 - s*t2; y = c*t2 + std::conj(s)*t1; }

  template <class MAT, class T>
  void row_rot(const MAT &AA, T c, T s, size_type i, size_type k) {
    MAT &A = const_cast<MAT &>(AA); // can be specialized for row matrices
    for (size_type j = 0; j < mat_ncols(A); ++j)
      Apply_Givens_rotation(A(i,j), A(k,j), c, s);
  }

  template <class MAT, class T>
  void col_rot(const MAT &AA, T c, T s, size_type i, size_type k) {
    MAT &A = const_cast<MAT &>(AA); // can be specialized for column matrices
    for (size_type j = 0; j < mat_nrows(A); ++j)
      Apply_Givens_rotation(A(j,i), A(j,k), c, s);
  }

  /* ********************************************************************* */
  /*    Householder tridiagonalisation for symmetric matrices              */
  /* ********************************************************************* */

  template <class MAT1, class MAT2> 
  void Householder_tridiagonalisation(const MAT1 &AA, const MAT2 &QQ,
				     bool compute_q) {
    MAT1 &A = const_cast<MAT1 &>(AA); MAT2 &Q = const_cast<MAT2 &>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    size_type n = mat_nrows(A); 
    std::vector<value_type> v(n), p(n), w(n);
    for (size_type k = 1; k < n-1; ++k) { // A pas mal optimiser ...
      sub_interval SUBI(k, n-k);
      v.resize(n-k); p.resize(n-k); w.resize(n-k); 
      for (size_type l = k; l < n; ++l) 
	{ v[l-k] = A(l, k-1); A(l, k-1)=0.0; A(k-1, l)=0.0; }
      magnitude_type norm_x = -vect_norm2(v) * dal::sgn(v[0]);
      house_vector(v);
      magnitude_type norm = vect_norm2_sqr(v);
      A(k, k-1) = A(k-1, k) = norm_x;
      gmm::mult(sub_matrix(A, SUBI), v, p);
      gmm::scale(p, value_type(2) / norm);
      gmm::add(gmm::scaled(p, value_type(-1)),
	       gmm::scaled(v, vect_sp(p, v) / norm), w);
      rank_one_update(sub_matrix(A, SUBI), v, w);
      rank_one_update(sub_matrix(A, SUBI), w, v);
      // + act de Q ...
      if (compute_q)
	cout << "WARNING : Q not computed in Householder_tridiagonalisation\n";
    }
  }

}

#endif

