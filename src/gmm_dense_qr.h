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
/* Authors : Caroline Lecalvez, Caroline.Lecalvez@gmm.insa-tlse.fr         */
/*           Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
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
    temp_vector V(m), W(m);

    for (size_type j = 0; j < n; ++j) {
      sub_interval SUBI(j, m-j), SUBJ(j, n-j);

      for (size_type i = j; i < m; ++i) W[i] = R(i, j);
      house_vector(sub_vector(W, SUBI), sub_vector(V, SUBI));

      row_house_update(sub_matrix(R, SUBI, SUBJ),
		       sub_vector(V, SUBI), sub_vector(W, SUBJ));
  
      for (size_type i = j+1; i < m; ++i) R(i, j) = V[i];
    }

    cout << "R prov = " << R << endl;

    gmm::copy(identity_matrix(), Q);
    for (size_type j = n-2; j != size_type(-1); --j) {
      sub_interval SUBI(j, m-j);

      V[j] = value_type(1);
      for (size_type i = j+1; i < m; ++i)
	{ V[i] = R(i, j); R(i, j) = value_type(0); }

      row_house_update(sub_matrix(Q, SUBI, SUBI), 
		       sub_vector(V, SUBI), sub_vector(W, SUBI));
    }

  }

}

#endif

