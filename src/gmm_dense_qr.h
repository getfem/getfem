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
/* Date : September 12, 2003.                                              */
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


  /* ********************************************************************* */
  /*    Householder vector computation                                     */
  /* ********************************************************************* */

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
  
  /* ********************************************************************* */
  /*    Householder updates                                                */
  /* ********************************************************************* */

  template <class MAT, class VECT1, class VECT2>
    void row_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW); MAT &A = const_cast<MAT &>(AA);
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    magnitude_type beta = magnitude_type(2) / vect_norm2_sqr(V);
    gmm::mult(gmm::transposed(A), gmm::scaled(V, -beta), W);
    rank_one_update(A, V, W);
  }

  template <class MAT, class VECT1, class VECT2>
    void col_house_update(const MAT &AA, const VECT1 &V, const VECT2 &WW) {
    VECT2 &W = const_cast<VECT2 &>(WW); MAT &A = const_cast<MAT &>(AA);
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    magnitude_type beta = magnitude_type(2) / vect_norm2_sqr(V);
    gmm::mult(A, gmm::scaled(V, -beta), W);
    rank_one_update(A, W, V);
  }

  /* ********************************************************************* */
  /*    QR factorization using Householder method.                         */
  /* ********************************************************************* */
  
  template <class MAT1, class MAT2, class MAT3>
    void qr_factor(const MAT1 &A, MAT2 &Q, MAT3 &R) { 
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef std::vector<value_type> temp_vector;

    size_type m = mat_nrows(A), n = mat_ncols(A);
    gmm::copy(A, R);
    
    if (m < n) DAL_THROW(dimension_error, "dimensions mismatch");
    std::vector<value_type> W(m);
    dense_matrix<value_type> VV(m, n);

    for (size_type j = 0; j < n; ++j) {
      sub_interval SUBI(j, m-j), SUBJ(j, n-j);

      for (size_type i = j; i < m; ++i) W[i] = R(i, j);
      house_vector(sub_vector(W, SUBI), sub_vector(mat_col(VV,j), SUBI));
      row_house_update(sub_matrix(R, SUBI, SUBJ),
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBJ));
      for (size_type i = j+1; i < m; ++i) R(i, j) = 0.0; // à garder ?
    }

    gmm::copy(identity_matrix(), Q);
    for (size_type j = n-1; j != size_type(-1); --j) {
      sub_interval SUBI(j, m-j);
      row_house_update(sub_matrix(Q, SUBI, SUBI), 
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBI));
    }
  }

  /* ********************************************************************* */
  /*    Hessember reduction with Householder.                              */
  /* ********************************************************************* */

  template <class MAT1, class MAT2>
    void Hessenberg_reduction(const MAT1& AA, const MAT2 &QQ, bool compute_Q){
    MAT1& A = const_cast<MAT1&>(AA); MAT1& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(A);
    std::vector<value_type> v(n), w(n);
    if (compute_Q) gmm::copy(identity_matrix(), Q);
    for (size_type k = 1; k < n-1; ++k) {
      for (size_type j = k; j < n; ++j) w[j] = A(j, k-1);
      sub_interval SUBI(k, n-k), SUBJ(k-1,n-k+1), SUBK(0,n);
      house_vector(sub_vector(w, SUBI), sub_vector(v, SUBI));
      row_house_update(sub_matrix(A, SUBI, SUBJ), sub_vector(v, SUBI),
 		       sub_vector(w, SUBJ));
      col_house_update(sub_matrix(A, SUBK, SUBI), sub_vector(v, SUBI),
 		       sub_vector(w, SUBK));
       if (compute_Q)
       col_house_update(sub_matrix(Q, SUBK, SUBI), sub_vector(v, SUBI),
                       sub_vector(w, SUBK));
       // for (size_type j = k+1; j < n; ++j) A(j, k-1) = v[j];
    }
  }

  /* ********************************************************************* */
  /*    Compute eigenvalue vector.                                         */
  /* ********************************************************************* */

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, VECT &V, double tol, TA, TV) {
    tol *= 2.0;
    for (size_type i = 0; i < mat_nrows(A); ++i)
      if ((i == n-1) ||
	  dal::abs(A(i+1,i)) < (dal::abs(A(i,i))+dal::abs(A(i+1,i+1)))*tol)
	V[i] = TV(A(i,i));
      else {
	TA tr = A(i,i) + A(i+1, i+1);
	TA det = A(i,i)*A(i+1, i+1) - A(i,i+1)*A(i+1, i);
	TA delta = tr*tr - 4 * det;
	if (delta < TA(0)) {
	  DAL_WARNING(2, "A Ccomplex eigenvalue has been detected");
	  V[i+1] = V[i] =std::abs(std::complex<TV>(tr/2.0, sqrt(-delta)/2.0));
	  
	}
	else {
	  V[i]   = TV(A(i,i));
	  V[i+1] = TV(A(i+1, i+1));
	}
	++i;
      }
  }

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, VECT &V, double tol,
		   TA, std::complex<TV>) {
    size_type n = mat_nrows(A);
    tol *= 2.0;
    for (size_type i = 0; i < n; ++i)
      if ((i == n-1) ||
	  dal::abs(A(i+1,i)) < (dal::abs(A(i,i))+dal::abs(A(i+1,i+1)))*tol)
	V[i] = std::complex<TV>(A(i,i));
      else {
	TA tr = A(i,i) + A(i+1, i+1);
	TA det = A(i,i)*A(i+1, i+1) - A(i,i+1)*A(i+1, i);
	TA delta = tr*tr - 4 * det;
	if (delta < TA(0)) {
	  V[i] = std::complex<TV>(tr / 2.0, sqrt(-delta) / 2.0);
	  V[i+1] = std::complex<TV>(tr / 2.0, -sqrt(-delta) / 2.0);
	}
	else {
	  V[i]   = std::complex<TV>(A(i,i));
	  V[i+1] = std::complex<TV>(A(i+1, i+1));
	}
	++i;
      }
  }

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, const VECT &VV, double tol,
		   std::complex<TA>, TV)
  { for (size_type i = 0; i < mat_nrows(A); ++i) V[i] = TV(A1(i,i)); }

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, const VECT &VV, double tol,
		   std::complex<TA>, std::complex<TV>) {
    for (size_type i = 0; i < mat_nrows(A); ++i)
      V[i] = std::complex<TV>(A1(i,i));
  }

  template <class MAT, class VECT> inline
  void extract_eig(const MAT &A, const VECT &V, double tol) {
    extract_eig(A, const_cast<VECT&>(V), tol,
		typename linalg_traits<MAT>::value_type(),
		typename linalg_traits<VECT>::value_type());
  }

  /* ********************************************************************* */
  /*    Stop criterion for QR algorithms                                   */
  /* ********************************************************************* */

  template <class MAT, class T>
  void stop_criterion(MAT &A, size_type &p, size_type &q,
		      double tol, T) {

    size_type n = mat_nrows(A);
    for (size_type i = 1; i < n; ++i)
      if (dal::abs(A(i,i-1)) < (dal::abs(A(i,i))+ dal::abs(A(i-1,i-1)))*tol)
	A(i,i-1) = T(0);
       
    q = 0;
    while ((q < n-1 && A(n-1-q, n-2-q) == T(0)) ||
	   (q < n-2 && A(n-2-q, n-3-q) == T(0))) ++q;
    if (q >= n-2) q = n;
    p = n-q; if (p) --p; if (p) --p;
    while (p > 0 && A(p,p-1) != T(0)) --p;
  }

  template <class MAT, class T> // complex version, To be verified
  void stop_criterion(MAT &A, size_type &p, size_type &q,
		      double tol, std::complex<T>) {

    for (size_type i = 1; i < n; ++i)
      if (dal::abs(A(i,i-1)) < (dal::abs(A(i,i))+ dal::abs(A(i-1,i-1)))*tol)
	A(i,i-1) = T(0);
       
    q = 0;
    while (q < n-1 && A(n-1-q, n-2-q) == T(0)) ++q;
    if (q >= n-1) q = n;
    p = n-q; if (p) --p; if (p) --p;
    while (p > 0 && A(p,p-1) != T(0)) --p;
  }

  template <class MAT> inline
  void stop_criterion(const MAT &A, size_type &p, size_type &q, double tol) {
    stop_criterion(const_cast<MAT&>(A), p, q, tol,
		   typename linalg_traits<MAT>::value_type());
  }
  

  /* ********************************************************************* */
  /*    Power qr algorihm.                                                 */
  /* ********************************************************************* */

  // QR method for real square matrices based on QR factorisation.
  // eigval has to be a complex vector.
  template <class MAT1, class VECT, class MAT2>
    void power_qr_algorithm(const MAT1 &A, VECT &eigval, MAT2 &eigvect,
		      double tol = 1E-12, bool compvect = true) {
    typedef typename linalg_traits<MAT1>::value_type value_type;

    size_type n = mat_nrows(A), p, q, ite = 0;
    MAT1 Q(n, n), R(n,n), A1(n,n); 
    gmm::copy(A, A1);

    Hessenberg_reduction(A1, eigvect, compvect);
    stop_criterion(A1, p, q, tol);

    while (q < n) {
      qr_factor(A1, Q, R);
      gmm::mult(R, Q, A1);
      if (compvect) { gmm::mult(eigvect, Q, R); gmm::copy(R, eigvect); }
      
      stop_criterion(A1, p, q, tol);
      if (++ite > n*1000) DAL_THROW(failure_error, "QR algorithm failed");
    }
    extract_eig(A1, eigval, tol); 
  }

  template <class MAT1, class VECT>
    void power_qr_algorithm(const MAT1 &a, VECT &eigval, double tol = 1E-12) {
    typedef typename linalg_traits<MAT1>::value_type value_type;
    
    dense_matrix<value_type> m(0,0);
    power_qr_algorithm(a, eigval, m, tol, false); 
  }

  /* ********************************************************************* */
  /*    Francis QR step.                                                   */
  /* ********************************************************************* */

  template <class MAT1, class MAT2>
    void Francis_qr_step(const MAT1& HH, const MAT2 &QQ, bool compute_Q) {
    MAT1& H = const_cast<MAT1&>(HH); MAT1& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(H); 
    
    std::vector<value_type> v(3), w(n);
    if (compute_Q) gmm::copy(identity_matrix(), Q);

    value_type s = H(n-2, n-2) + H(n-1, n-1);
    value_type t = H(n-2, n-2) * H(n-1, n-1) - H(n-2, n-1) * H(n-1, n-2);
    value_type x = H(0, 0)*H(0, 0) + H(0,1) * H(1, 0) - s * H(0,0) + t;
    value_type y = H(1,0) * (H(0,0) + H(1,1) - s);
    value_type z = H(1, 0) * H(2, 1);
    for (size_type k = 0; k < n - 2; ++k) {
      w[0] = x; w[1] = y; w[2] = z;
      house_vector(sub_vector(w, sub_interval(0, 3)), v);
      size_type r = std::min(k+4, n), q = (k==0) ? 0 : k-1;
      sub_interval SUBI(k, 3), SUBJ(0, r), SUBK(q, n-q);
      
      row_house_update(sub_matrix(H, SUBI, SUBK),  v, sub_vector(w, SUBK));
      col_house_update(sub_matrix(H, SUBJ, SUBI),  v, sub_vector(w, SUBJ));
      
      if (compute_Q)
       	col_house_update(sub_matrix(Q, SUBJ, SUBI),  v, sub_vector(w, SUBJ));
      x = H(k+1, k);
      y = H(k+2, k);
      if (k < n-3) z = H(k+3, k);
    }
    std::vector<value_type> vv(2);
    w[0] = x; w[1] = y;
    house_vector(sub_vector(w, sub_interval(0, 2)), vv);
    row_house_update(sub_matrix(H,sub_interval(n-2,2), sub_interval(n-3,3)),
		     vv, sub_vector(w, sub_interval(0, 3)));
    col_house_update(sub_matrix(H,sub_interval(0,n), sub_interval(n-2,2)),
		     vv, w);
    if (compute_Q)
      col_house_update(sub_matrix(Q,sub_interval(0,n), sub_interval(n-2,2)),
		       vv, w);
  }

  /* ********************************************************************* */
  /*    Implicit QR algorithm.                                             */
  /* ********************************************************************* */

  // implicit QR method for real square matrices based on QR factorisation.
  // eigval has to be a complex vector.
  template <class MAT1, class VECT, class MAT2>
    void implicit_qr_algorithm(const MAT1 &A, VECT &eigval, MAT2 &eigvect,
			       double tol = 1E-12, bool compvect = true) {
    typedef typename linalg_traits<MAT1>::value_type value_type;

    size_type n = mat_nrows(A), q = 0, p;
    dense_matrix<value_type> Z(n,n), B(n,n), H(n,n);
    gmm::copy(A, H);
    Hessenberg_reduction(H, eigvect, compvect);
    stop_criterion(H, p, q, tol);

    while (q < n) {
      sub_interval SUBI(p, n-p-q), SUBJ(0,n);
      Francis_qr_step(sub_matrix(H, SUBI, SUBI),
		      sub_matrix(Z, SUBI, SUBI), compvect);
      if (compvect) {
	gmm::mult(sub_matrix(eigvect, SUBJ, SUBI),
		  sub_matrix(Z, SUBI, SUBI), sub_matrix(B, SUBJ, SUBI));
	gmm::copy(sub_matrix(B, SUBJ, SUBI), sub_matrix(eigvect, SUBJ, SUBI));
      }
      stop_criterion(H, p, q, tol);
    }
    extract_eig(H, eigval, tol);
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

