/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_ilutp.h : ILUTP preconditionner for sparse       */
/*                                  matrices                               */
/*     									   */
/* Date : October 14, 2004.                                                */
/* Author : Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
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
#ifndef GMM_PRECOND_ILUTP_H
#define GMM_PRECOND_ILUTP_H

//: ILUTP:  Incomplete LU with threshold and K fill-in Preconditioner and
//          partial pivoting (See Yousef Saad, Iterative Methods for
//          sparse linear systems, PWS Publishing Company, section 10.4.4

#include <gmm_precond_ilut.h>

namespace gmm {

  template <typename Matrix>
  class ilutp_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;
    typedef row_matrix<svector> LU_Matrix;

    bool invert;
    LU_Matrix L, U;
    std::vector<size_type> ipvt;
    std::vector<size_type> ipvtinv;
    gmm::unsorted_sub_index indperm;
    gmm::unsorted_sub_index indperminv;    
    mutable std::vector<value_type> temporary;

  protected:
    int K;
    double eps;    

    template<typename M> void do_ilutp(const M&, row_major);
    void do_ilutp(const Matrix&, col_major);

  public:
    void build_with(const Matrix& A) {
      invert = false;
      gmm::resize(L, mat_nrows(A), mat_ncols(A));
      gmm::resize(U, mat_nrows(A), mat_ncols(A));
      do_ilutp(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    ilutp_precond(const Matrix& A, int k_, double eps_) 
      : L(mat_nrows(A), mat_ncols(A)), U(mat_nrows(A), mat_ncols(A)),
	K(k_), eps(eps_) { build_with(A); }
    ilutp_precond(int k_, double eps_) :  K(k_), eps(eps_) {}
    ilutp_precond(void) { K = 10; eps = 1E-7; }
    size_type memsize() const { 
      return sizeof(*this) + (nnz(U)+nnz(L))*sizeof(value_type);
    }
  };


  template<typename V, typename R>
  size_type find_max_after_i(const V &v, size_type i, R val, abstract_dense) {
    typename linalg_traits<V>::const_iterator it = vect_const_begin(v) + i,
      ite = vect_const_end(v);
    size_type j = i;
    for ( ; it != ite; ++it, ++j)
      if (gmm::abs(*it) > val && j > i)
	{ i = j; val = gmm::abs(*it); }
    return i;
  }

  template<typename V, typename R>
  size_type find_max_after_i(const V &v, size_type i, R val, abstract_sparse) {
    typename linalg_traits<V>::const_iterator it = vect_const_begin(v),
      ite = vect_const_end(v);
    for ( ; it != ite; ++it)
      if (gmm::abs(*it) > val && it.index() > i)
	{ i = it.index(); val = gmm::abs(*it); }
    return i;
  }

  template<typename V, typename R>
  size_type find_max_after_i(const V &v, size_type i, R val, abstract_skyline)
  { return find_max_after_i(v, i, val, abstract_sparse()); }

  template<typename V, typename R>
  size_type find_max_after_i(const V &v, size_type i, R val) {
    return find_max_after_i(v, i, val,
			    typename linalg_traits<V>::storage_type());
  }

  template<typename Matrix> template<typename M> 
  void ilutp_precond<Matrix>::do_ilutp(const M& A, row_major) {
    typedef value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    size_type n = mat_nrows(A);
    if (n == 0) return;
    std::vector<T> indiag(n);
    ipvt.resize(n); ipvtinv.resize(n); temporary.resize(n);
    for (size_type i = 0; i < n; ++i) ipvt[i] = i;
    svector w(mat_ncols(A));
    T tmp;
    gmm::clear(U); gmm::clear(L);
    R prec = default_tol(R()); 
    R max_pivot = gmm::abs(A(0,0)) * prec;

    for (size_type i = 0; i < n; ++i) {

      // To be optimized, computation of sub_index is made twice ...
      //     the reverse index could be updated at each iteration
      size_type ip= find_max_after_i(gmm::sub_vector(mat_const_row(A, i),
						     unsorted_sub_index(ipvt)),
				     i, gmm::abs(A(i,i)));
      if (ip != i) std::swap(ipvt[i], ipvt[ip]);
      
      gmm::copy(gmm::sub_vector(mat_const_row(A, i),
				unsorted_sub_index(ipvt)), w);

      double norm_row = gmm::vect_norm2(w);

      size_type nL = 0, nU = 0;
      if (is_sparse(A)) {
	typename linalg_traits<svector>::iterator it = vect_begin(w),
	  ite = vect_end(w);
	for (; it != ite; ++it) if (i > it.index()) nL++;
	nU = w.nb_stored() - nL - 1;
      }

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	tmp = (wk->e) * indiag[k];
	if (gmm::abs(tmp) < eps * norm_row) { w.sup(k); --krow; } 
	else { wk->e += tmp; gmm::add(scaled(mat_row(U, k), -tmp), w); }
      }
      tmp = w[i];

      if (gmm::abs(tmp) <= max_pivot)
	{ DAL_WARNING(2, "pivot " << i << " is too small"); tmp = T(1); }

      max_pivot = std::max(max_pivot, std::min(gmm::abs(tmp) * prec, R(1)));
      indiag[i] = T(1) / tmp;
      U(i,i) = tmp; gmm::clean(w, eps * norm_row); w[i] = T(0);
      std::sort(w.begin(), w.end(), elt_rsvector_value_less_<T>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      size_type nnl = 0, nnu = 0;
      for (; wit != wite; ++wit) // copy to be optimized ...
	if (wit->c < i) { if (nnl < nL+K) L(i, wit->c) = wit->e; ++nnl; }
	else            { if (nnu < nU+K) U(i, wit->c) = wit->e; ++nnu; }
    }
    indperm = unsorted_sub_index(ipvt);
    cout << "pivots = " << ipvt << endl;
    for (size_type i = 0; i < n; ++i) ipvtinv[ipvt[i]] = i;
    indperminv = unsorted_sub_index(ipvtinv);
  }

  template<typename Matrix> 
  void ilutp_precond<Matrix>::do_ilutp(const Matrix& A, col_major) {
    do_ilutp(gmm::transposed(A), row_major());
    invert = true;
  }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      gmm::copy(gmm::sub_vector(v1, P.indperminv), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      gmm::copy(v1, P.temporary);
      gmm::lower_tri_solve(P.L, P.temporary, true);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperm), v2);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const ilutp_precond<Matrix>& P,const V1 &v1,V2 &v2) {
    if (P.invert) {
      gmm::copy(v1, P.temporary);
      gmm::lower_tri_solve(P.L, P.temporary, true);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperm), v2);
    }
    else {
      gmm::copy(gmm::sub_vector(v1, P.indperminv), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      gmm::copy(gmm::sub_vector(v1, P.indperminv), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    }
    else {
      copy(v1, v2);
      gmm::lower_tri_solve(P.L, v2, true);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      copy(v1, v2);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      copy(v1, P.temporary);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperm), v2);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const ilutp_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    if (P.invert) {
      copy(v1, P.temporary);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperm), v2);
    }
    else {
      copy(v1, v2);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }
  
  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const ilutp_precond<Matrix>& P, const V1 &v1,
			     V2 &v2) {
    if (P.invert) {
      copy(v1, v2);
      gmm::lower_tri_solve(P.L, v2, true);
    }
    else {
      gmm::copy(gmm::sub_vector(v1, P.indperminv), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    }
  }

}

#endif 

