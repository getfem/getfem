/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_qr.h : QR decomposition and QR method for dense    */
/*                             matrices with Householder method.           */
/*                                                                         */
/* ref :  G.H. Golub, C.F. Van Loan, Matrix Computations, second edition   */
/*        The Johns Hopkins University Press, 1989.                        */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
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

#ifndef GMM_DENSE_QR_H
#define GMM_DENSE_QR_H

#include "gmm_dense_lu.h"

namespace gmm {

  template <class VECT1, class VECT2>
    void house_vector(const VECT1 &X, const VECT2 &VV) {
    VECT2 &V = const_cast<VECT2 &>(VV);
    typedef typename linalg_traits<VECT1>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    size_type n = vect_size(X);
    if (vect_size(V) != n) DAL_THROW(dimension_error, "dimensions mismatch");
    magnitude_type mu = gmm::vect_norm2(X);
    gmm::copy(X, V);
    if (mu != magnitude_type(0)) {
      value_type beta = X[0] + dal::sgn(X[0]) * mu; // with complexes ??
      gmm::scale(V, value_type(1)/beta);
    }
    V[0] = value_type(1);
  }
  
  template <class MAT, class VECT1, class VECT2>
    void row_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW);
    MAT &A = const_cast<MAT &>(AA);
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    magnitude_type beta = magnitude_type(2) / vect_norm2_sqr(V);
    gmm::mult(gmm::transposed(A), gmm::scaled(V, -beta), W);
    rank_one_update(A, V, W);
  }

  template <class MAT, class VECT1, class VECT2>
    void col_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW);
    MAT &A = const_cast<MAT &>(AA);
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    magnitude_type beta = magnitude_type(2) / vect_norm2_sqr(V);
    gmm::mult(A, gmm::scaled(V, -beta), W);
    rank_one_update(A, W, V);
  }
  
  template <class MAT1, class MAT2, class MAT3>
    void qr_factor(const MAT1 &A, MAT2 &Q, MAT3 &R) { 
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef std::vector<value_type> temp_vector;
    // To be optimized ..

    gmm::copy(A, R);
    size_type m = mat_nrows(A), n = mat_ncols(A);
    
    if (m < n) DAL_THROW(dimension_error, "dimensions mismatch");
    std::vector<value_type> W(m);
    dense_matrix<value_type> VV(m, n);
    

    for (size_type j = 0; j < n; ++j) {
      sub_interval SUBI(j, m-j), SUBJ(j, n-j);

      for (size_type i = j; i < m; ++i) W[i] = R(i, j);
      house_vector(sub_vector(W, SUBI), sub_vector(mat_col(VV,j), SUBI));

      row_house_update(sub_matrix(R, SUBI, SUBJ),
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBJ));
  
      for (size_type i = j+1; i < m; ++i) R(i, j) = 0.0;
    }

    gmm::copy(identity_matrix(), Q);
    for (size_type j = n-1; j != size_type(-1); --j) {
      sub_interval SUBI(j, m-j);
      row_house_update(sub_matrix(Q, SUBI, SUBI), 
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBI));
    }

  }

  template <class MAT1, class VECT, class MAT2>
    void qr_method(const MAT1 &a, VECT &eigval, MAT2 &eigvect,
		   bool vect = true) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    
    size_type n = mat_nrows(a);
    MAT1 q(n, n), r(n,n), a1(n,n); 
    gmm::copy(a, a1);
    if (vect) gmm::copy(identity_matrix(), eigvect);
    
    for (size_type ite = 0; ite < n*100; ++ite) { // trop d'itérations ?
      qr_factor(a1, q, r);
      if (vect) { gmm::mult(eigvect, q, a1); gmm::copy(a1, eigvect); }
      gmm::mult(r, q, a1);
      
      magnitude_type vmax(0);
      for (size_type i = 0; i < n; ++i)
	vmax = std::max(vmax, dal::abs(a1(i,i)));
      
      bool ok = true;
      for (size_type i = 0; i < n && ok; ++i) { // to be optimized
        for (size_type j = 0; j < i; ++j)
	  if (dal::abs(a1(i,j)) > vmax * 1E-12) // critère à revoir ?
	    { ok = false; break; }
      }
      if (ok) {
	for (size_type i = 0; i < n; ++i) eigval[i] = a1(i,i);
	return;
      }
    }
    DAL_THROW(failure_error, "QR method failed");
  }

  template <class MAT1, class VECT> void qr_method(const MAT1 &a, VECT &eigval)
  {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    
    dense_matrix<value_type> m(0,0);
    qr_method(a, eigval, m, false); 
  }


}

#endif

