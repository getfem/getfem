/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_qr.h : QR decomposition and QR method for dense    */
/*                             matrices with Householder method.           */
/*                                                                         */
/* ref :  G.H. Golub, C.F. Van Loan, Matrix Computations, second edition   */
/*        The Johns Hopkins University Press, 1989.                        */
/*     									   */
/* Date : September 11, 2003.                                              */
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
      value_type beta = X[0] + dal::sgn(X[0]) * mu;
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

  // QR method for real square matrices based on QR factorisation.
  // eigval has to be a complex vector.
  template <class MAT1, class VECT, class MAT2>
    void qr_algorithm(const MAT1 &a, VECT &eigval, MAT2 &eigvect,
		      double tol = 1E-12, bool compvect = true) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    typedef typename linalg_traits<VECT>::value_type vect_value_type;

    size_type n = mat_nrows(a);
    MAT1 q(n, n), r(n,n), a1(n,n); 
    gmm::copy(a, a1);
    if (compvect) gmm::copy(identity_matrix(), eigvect);
    
    for (size_type ite = 0; ite < n*1000; ++ite) {
      qr_factor(a1, q, r);
      if (compvect) { gmm::mult(eigvect, q, a1); gmm::copy(a1, eigvect); }
      gmm::mult(r, q, a1);
      
      gmm::clean(a1, 1E-10);
      cout << "iteration " << ite << " a1 = " << a1 << endl;

      magnitude_type vmax(0);
      for (size_type i = 0; i < n; ++i)
	vmax = std::max(vmax, dal::abs(a1(i,i)));
      
      bool ok = true;
      for (size_type i = 1; i < n && ok; ++i) { // to be optimized
        for (size_type j = 0; j < i-1; ++j)
	  if (dal::abs(a1(i,j)) > vmax * 1E-12) // validity of criterion ?
	    { ok = false; break; }
      }
      if (ok) {
	for (size_type i = 0; i < n; ++i)
	  if ((i == n-1) || dal::abs(a1(i+1,i)) <= vmax * tol)
	    eigval[i] = vect_value_type(a1(i,i));
	  else {
	    // ne marche pas pour les matrices complexes
	    value_type tr = a1(i,i) + a1(i+1, i+1);
	    value_type det = a1(i,i)*a1(i+1, i+1) - a1(i,i+1)*a1(i+1, i);
	    value_type delta = tr*tr - 4 * det;
	    if (delta < value_type(0)) {
	      eigval[i] = vect_value_type(tr / 2, sqrt(-delta) / 2);
	      eigval[i+1] = vect_value_type(tr / 2, -sqrt(-delta) / 2);
	    }
	    else {
	      eigval[i]   = vect_value_type(a1(i,i));
	      eigval[i+1] = vect_value_type(a1(i+1, i+1));
	    }
	    ++i;
	  }
	return;
      }
    }
    DAL_THROW(failure_error, "QR algorithm failed");
  }

  template <class MAT1, class VECT>
    void qr_algorithm(const MAT1 &a, VECT &eigval, double tol = 1E-12) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    
    dense_matrix<value_type> m(0,0);
    qr_algorithm(a, eigval, m, tol, false); 
  }

  template <class MAT1, class MAT2>
    void Francis_qr_step(const MAT1& HH, const MAT2 &QQ, bool compute_Q) {
    MAT1& H = const_cast<MAT1&>(HH);
    MAT1& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(H), m = n - 1;
    std::vector<value_type> v(3), w(3);
    if (compute_Q) gmm::copy(identity_matrix(), Q);

    value_type s = H(m-1, m-1) + H(n-1, n-1);
    value_type t = H(m-1,m-1)*H(n-1,n-1) - H(m-1,n-1)*H(n-1,m-1);
    value_type x = H(0, 0) + H(1, 1) + H(0,1) * H(1, 0) - s * H(0,0) + t;
    value_type y = H(1,0) * (H(0,0) + H(1,1) - s);
    value_type z = H(1, 0) * H(2, 1);
    for (size_type k = 0; k < n - 2; ++k) {
      w[0] = x; w[1] = y; w[2] = z;
      house_vector(w, v);
      sub_interval SUBI(k, 3);
      row_house_update(sub_matrix(H, SUBI, sub_interval(k, n-k)),  v, w);
      size_type r = std::min(k+4, n);
      col_house_update(sub_matrix(H, sub_interval(0,r), SUBI),  v, w);
      if (compute_Q)
	col_house_update(sub_matrix(Q, sub_interval(0,r), SUBI),  v, w);
      x = H(k+1, k);
      y = H(k+2, k);
      if (k < n-3) z = H(k+3, k);
    }
    std::vector<value_type> vv(2), ww(2);
    ww[0] = x; ww[1] = y;
    house_vector(ww, vv);
    row_house_update(sub_matrix(H,sub_interval(n-2,2),sub_interval(n-3,3)),
		     vv, ww);
    col_house_update(sub_matrix(H,sub_interval(0,n)  ,sub_interval(n-2,2)),
		     vv, ww);
  }

  template <class MAT1, class MAT2>
    void Hessenberg_reduction(const MAT1& AA, const MAT2 &QQ, bool compute_Q) {
    MAT1& A = const_cast<MAT1&>(AA);
    MAT1& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(A);
    std::vector<value_type> v(n), w(n);
    if (compute_Q) gmm::copy(identity_matrix(), Q);
    for (size_type k = 1; k < n-1; ++k) {
      for (size_type j = k; j < n; ++j) w[j] = A(j, k-1);
      sub_interval SUBI(k, n-k), SUBJ(k-1,n-k+1), SUBK(0,n);
      house_vector(sub_vector(w, SUBI), sub_vector(v, SUBI));
       row_house_update(sub_matrix(A, SUBI, SUBJ), sub_vector(v, SUBI),
 		       sub_vector(w, SUBI));
       col_house_update(sub_matrix(A, SUBK, SUBI), sub_vector(v, SUBI),
 		       sub_vector(w, SUBI));
       if (compute_Q)
 	col_house_update(sub_matrix(Q, SUBK, SUBI), sub_vector(v, SUBI),
			 sub_vector(w, SUBI));
       for (size_type j = k+1; j < n; ++j) A(j, k-1) = v[j];
    }
  }


  // implicit QR method for real square matrices based on QR factorisation.
  // eigval has to be a complex vector.
  template <class MAT1, class VECT, class MAT2>
    void implicit_qr_algorithm(const MAT1 &a, VECT &eigval, MAT2 &eigvect,
			       double tol = 1E-12, bool compvect = true) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    typedef typename linalg_traits<VECT>::value_type vect_value_type;

    size_type n = mat_nrows(a), q = 0, p;
    dense_matrix<value_type> z(n,n), aux(n,n), h(n,n);
    gmm::copy(a, h);
    Hessenberg_reduction(h, eigvect, compvect);
    
    while (q != n-1) {
      
       for (size_type i = 1; i < n; ++i)
 	if (dal::abs(h(i,i-1)) < (dal::abs(h(i,i)) + dal::abs(h(i+1,i+1)))*tol)
 	  h(i,i-1) = value_type(0);

       cout << "h = " << h << endl;

       p = 0;
       while (p < n-1 && h(p+1,p) == value_type(0)) ++p;
       q = 0;
       while (q < n-1 && h(n-1-q, n-2-q) == value_type(0)) ++q;

       if (q < n-1) {
 	sub_interval SUBI(p, n-p-q), SUBJ(0,n);
 	Francis_qr_step(sub_matrix(h, SUBI, SUBI),
 			sub_matrix(z, SUBI, SUBI), compvect);
 	if (compvect) {
 	  gmm::mult(sub_matrix(eigvect, SUBJ, SUBI),
 		    sub_matrix(z, SUBI, SUBI),
 		    sub_matrix(aux, SUBJ, SUBI));
 	  gmm::copy(sub_matrix(aux, SUBJ, SUBI),
 		    sub_matrix(eigvect, SUBJ, SUBI));
 	}
       } 
    } 
  }


  template <class MAT1, class VECT>
    void implicit_qr_algorithm(const MAT1 &a, VECT &eigval, 
			       double tol = 1E-12) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    
    dense_matrix<value_type> m(0,0);
    implicit_qr_algorithm(a, eigval, m, tol, false); 
  }


}

#endif

