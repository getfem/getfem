/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_qr.h : QR decomposition and QR method for dense    */
/*                             matrices.                                   */
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

#ifndef GMM_DENSE_QR_H
#define GMM_DENSE_QR_H

namespace gmm {


  /* ********************************************************************* */
  /* Default tolerance.                                                    */
  /* ********************************************************************* */
  
  template<class T> inline double default_tol(T)
  { int i=sizeof(T)/4; double tol(2); while(i-- > 0) tol*=1E-8; return tol; }
  template<class T> inline double default_tol(std::complex<T>)
  { return default_tol(T()); }

  /* ********************************************************************* */
  /* QR factorization using Householder method (complex and real version). */
  /* ********************************************************************* */
  
  template <class MAT1, class MAT2, class MAT3>
    void qr_factor(const MAT1 &A, const MAT2 &QQ, const MAT3 &RR) { 
    MAT2 &Q = const_cast<MAT2 &>(QQ); MAT3 &R = const_cast<MAT3 &>(RR); 
    typedef typename linalg_traits<MAT1>::value_type value_type;
    typedef std::vector<value_type> temp_vector;

    size_type m = mat_nrows(A), n = mat_ncols(A);
    if (m < n) DAL_THROW(dimension_error, "dimensions mismatch");
    dense_matrix<value_type> *MWORK;
    std::auto_ptr< dense_matrix<value_type> > auto_MWORK;
    if (m == n)
      { gmm::copy(A, R); MWORK = &R; }
    else
      { MWORK = new dense_matrix<value_type>(m,n); auto_MWORK.reset(MWORK); gmm::copy(A, *MWORK); }
    gmm::copy(identity_matrix(), Q);
    
    std::vector<value_type> W(m);
    dense_matrix<value_type> VV(m, n);

    for (size_type j = 0; j < n; ++j) {
      sub_interval SUBI(j, m-j), SUBJ(j, n-j) /* , SUBK(0, m) */;

      for (size_type i = j; i < m; ++i) VV(i,j) = (*MWORK)(i, j);
      house_vector(sub_vector(mat_col(VV,j), SUBI));

      row_house_update(sub_matrix(*MWORK, SUBI, SUBJ),
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBJ));
      //      col_house_update(sub_matrix(Q, SUBK, SUBI), 
      //	       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBK));
    }

    for (size_type j = n-1; j != size_type(-1); --j) {
      sub_interval SUBI(j, m-j), SUBJ(j, n-j);
      row_house_update(sub_matrix(Q, SUBI, SUBJ), 
		       sub_vector(mat_col(VV,j), SUBI), sub_vector(W, SUBJ));
    }
    
    if (m != n) {
      gmm::copy(sub_matrix(*MWORK, sub_interval(0, n), sub_interval(0, n)), R);
    }

  }

  /* ********************************************************************* */
  /*    Compute eigenvalue vector.                                         */
  /* ********************************************************************* */

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, VECT &V, double tol, TA, TV) {
    size_type n = mat_nrows(A);
    tol *= 2.0;
    for (size_type i = 0; i < n; ++i) {
      if ((i < n-1) &&
	  dal::abs(A(i+1,i)) >= (dal::abs(A(i,i))+dal::abs(A(i+1,i+1)))*tol) {
	TA tr = A(i,i) + A(i+1, i+1);
	TA det = A(i,i)*A(i+1, i+1) - A(i,i+1)*A(i+1, i);
	TA delta = tr*tr - TA(4.0) * det;
	if (delta < TA(0)) {
	  DAL_WARNING(2, "A complex eigenvalue has been detected");
	  V[i] = V[i+1] = tr / TA(2.0);
	}
	else {
	  V[i  ] = TA(tr + sqrt(delta))/ TA(2.0);
	  V[i+1] = TA(tr -  sqrt(delta))/ TA(2.0);
	}
	++i;
      }
      else
	V[i] = TV(A(i,i));
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
	TA delta = tr*tr - TA(4.0) * det;
	if (delta < TA(0)) {
	  V[i] = std::complex<TV>(tr / TA(2.0), sqrt(-delta) / TA(2.0));
	  V[i+1] = std::complex<TV>(tr / TA(2.0), -sqrt(-delta) / TA(2.0));
	}
	else {
	  V[i  ] = TA(tr + sqrt(delta)) / TA(2.0);
	  V[i+1] = TA(tr -  sqrt(delta)) / TA(2.0);
	}
	++i;
      }
  }

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, const VECT &VV, double tol,
		   std::complex<TA>, TV)
  { DAL_THROW(failure_error, "Sorry, not allowed"); }

  template <class TA, class TV, class MAT, class VECT>
  void extract_eig(const MAT &A, VECT &V, double tol,
		   std::complex<TA>, std::complex<TV>) {
    size_type n = mat_nrows(A);
    tol *= 2.0;
    for (size_type i = 0; i < n; ++i)
      if ((i == n-1) ||
	  dal::abs(A(i+1,i)) < (dal::abs(A(i,i))+dal::abs(A(i+1,i+1)))*tol)
	V[i] = std::complex<TV>(A(i,i));
      else {
	std::complex<TA> tr = A(i,i) + A(i+1, i+1);
	std::complex<TA> det = A(i,i)*A(i+1, i+1) - A(i,i+1)*A(i+1, i);
	std::complex<TA> delta = tr*tr - TA(4.0) * det;
	V[i] = (tr + sqrt(delta)) / TA(2.0);
	V[i+1] = (tr - sqrt(delta)) / TA(2.0);
	++i;
      }
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

  template <class MAT>
  void stop_criterion(MAT &A, size_type &p, size_type &q, double tol) {
    typedef typename linalg_traits<MAT>::value_type value_type;
    size_type n = mat_nrows(A);
    for (size_type i = 1; i < n-q; ++i)
      if (dal::abs(A(i,i-1)) < (dal::abs(A(i,i))+ dal::abs(A(i-1,i-1)))*tol)
	A(i,i-1) = value_type(0);

    while ((q < n-1 && A(n-1-q, n-2-q) == value_type(0)) ||
	   (q < n-2 && A(n-2-q, n-3-q) == value_type(0))) ++q;
    if (q >= n-2) q = n;
    p = n-q; if (p) --p; if (p) --p;
    while (p > 0 && A(p,p-1) != value_type(0)) --p;
  }
  
  template <class MAT> inline
  void symmetric_stop_criterion(const MAT &AA, size_type &p, size_type &q,
				double tol) {
    typedef typename linalg_traits<MAT>::value_type value_type;
    MAT& A = const_cast<MAT&>(AA);
    size_type n = mat_nrows(A);
    for (size_type i = 1; i < n-q; ++i)
      if (dal::abs(A(i,i-1)) < (dal::abs(A(i,i))+ dal::abs(A(i-1,i-1)))*tol)
	A(i,i-1) = value_type(0);

    while (q < n-1 && A(n-1-q, n-2-q) == value_type(0)) ++q;
    if (q >= n-1) q = n;
    p = n-q; if (p) --p; if (p) --p;
    while (p > 0 && A(p,p-1) != value_type(0)) --p;
  }


  /* ********************************************************************* */
  /*    Basic qr algorithm.                                                */
  /* ********************************************************************* */

  // QR method for real or complex square matrices based on QR factorisation.
  // eigval has to be a complex vector if A has complex eigeinvalues.
  // Very slow method. Use implicit_qr_method instead.
  template <class MAT1, class VECT, class MAT2>
    void rudimentary_qr_algorithm(const MAT1 &A, const VECT &eigval_,
				  const MAT2 &eigvect_, double tol =
		 gmm::default_tol(typename linalg_traits<MAT1>::value_type()),
				  bool compvect = true) {
    VECT &eigval = const_cast<VECT &>(eigval_);
    MAT2 &eigvect = const_cast<MAT2 &>(eigvect_);

    typedef typename linalg_traits<MAT1>::value_type value_type;

    size_type n = mat_nrows(A), p, q = 0, ite = 0;
    dense_matrix<value_type> Q(n, n), R(n,n), A1(n,n); 
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
    void rudimentary_qr_algorithm(const MAT1 &a, VECT &eigval,
				  double tol = gmm::default_tol(typename
				  linalg_traits<MAT1>::value_type())) {
    dense_matrix<typename linalg_traits<MAT1>::value_type> m(0,0);
    rudimentary_qr_algorithm(a, eigval, m, tol, false); 
  }

  /* ********************************************************************* */
  /*    Francis QR step.                                                   */
  /* ********************************************************************* */

  template <class MAT1, class MAT2>
    void Francis_qr_step(const MAT1& HH, const MAT2 &QQ, bool compute_Q) {
    MAT1& H = const_cast<MAT1&>(HH); MAT2& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    size_type n = mat_nrows(H), nq = mat_nrows(Q); 
    
    std::vector<value_type> v(3), w(std::max(n, nq));

    value_type s = H(n-2, n-2) + H(n-1, n-1);
    value_type t = H(n-2, n-2) * H(n-1, n-1) - H(n-2, n-1) * H(n-1, n-2);
    value_type x = H(0, 0) * H(0, 0) + H(0,1) * H(1, 0) - s * H(0,0) + t;
    value_type y = H(1, 0) * (H(0,0) + H(1,1) - s);
    value_type z = H(1, 0) * H(2, 1);

    sub_interval SUBQ(0, nq);

    for (size_type k = 0; k < n - 2; ++k) {
      v[0] = x; v[1] = y; v[2] = z;
      house_vector(v);
      size_type r = std::min(k+4, n), q = (k==0) ? 0 : k-1;
      sub_interval SUBI(k, 3), SUBJ(0, r), SUBK(q, n-q);
      
      row_house_update(sub_matrix(H, SUBI, SUBK),  v, sub_vector(w, SUBK));
      col_house_update(sub_matrix(H, SUBJ, SUBI),  v, sub_vector(w, SUBJ));
      
      if (compute_Q)
       	col_house_update(sub_matrix(Q, SUBQ, SUBI),  v, sub_vector(w, SUBQ));

      x = H(k+1, k); y = H(k+2, k);
      if (k < n-3) z = H(k+3, k);
    }
    sub_interval SUBI(n-2,2), SUBJ(0, n), SUBK(n-3,3), SUBL(0, 3);
    v.resize(2);
    v[0] = x; v[1] = y;
    house_vector(v);
    row_house_update(sub_matrix(H, SUBI, SUBK), v, sub_vector(w, SUBL));
    col_house_update(sub_matrix(H, SUBJ, SUBI), v, sub_vector(w, SUBJ));
    if (compute_Q)
      col_house_update(sub_matrix(Q, SUBQ, SUBI), v, sub_vector(w, SUBQ));
  }

  /* ********************************************************************* */
  /*    Wilkinson Double shift QR step (from Lapack).                      */
  /* ********************************************************************* */

  template <class MAT1, class MAT2>
  void Wilkinson_double_shift_qr_step(const MAT1& HH, const MAT2 &QQ,
				      double tol, bool exc, bool compute_Q) {
    MAT1& H = const_cast<MAT1&>(HH); MAT2& Q = const_cast<MAT2&>(QQ);
    typedef typename linalg_traits<MAT1>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type n = mat_nrows(H), nq = mat_nrows(Q), m;
    std::vector<T> v(3), w(std::max(n, nq));
    const R dat1 = 0.75, dat2 = -0.4375;
    T h33, h44, h43h34, v1 = 0, v2 = 0, v3 = 0;
    
    if (exc) {                    /* Exceptional shift.                    */
      R s = dal::abs(H(n-1, n-2)) + dal::abs(H(n-2, n-3));
      h33 = h44 = dat1 * s;
      h43h34 = dat2*s*s;
    }
    else {                        /* Wilkinson double shift.               */
      h44 = H(n-1,n-1); h33 = H(n-2, n-2);
      h43h34 = H(n-1, n-2) * H(n-2, n-1);
    }

    /* Look for two consecutive small subdiagonal elements.                */
    /* Determine the effect of starting the double-shift QR iteration at   */
    /* row m, and see if this would make H(m-1, m-2) negligible.           */
    for (m = n-2; m != 0; --m) {
      T h11  = H(m-1, m-1), h22  = H(m, m);
      T h21  = H(m, m-1),   h12  = H(m-1, m);
      T h44s = h44 - h11,   h33s = h33 - h11;
      v1 = (h33s*h44s-h43h34) / h21 + h12;
      v2 = h22 - h11 - h33s - h44s;
      v3 = H(m+1, m);
      R s = dal::abs(v1) + dal::abs(v2) + dal::abs(v3);
      v1 /= s; v2 /= s; v3 /= s;
      if (m == 1) break;
      T h00 = H(m-2, m-2);
      T h10 = H(m-1, m-2);
      R tst1 = dal::abs(v1)*(dal::abs(h00)+dal::abs(h11)+dal::abs(h22));
      if (dal::abs(h10)*(dal::abs(v2)+dal::abs(v3)) <= tol * tst1) break;
    }

    /* Double shift QR step.                                               */
    sub_interval SUBQ(0, nq);
    for (size_type k = m-1; k < n-2; ++k) {
      v[0] = v1; v[1] = v2; v[2] = v3;
      house_vector(v);
      size_type r = std::min(k+4, n), q = (k==0) ? 0 : k-1;
      sub_interval SUBI(k, 3), SUBJ(0, r), SUBK(q, n-q);
      
      row_house_update(sub_matrix(H, SUBI, SUBK),  v, sub_vector(w, SUBK));
      col_house_update(sub_matrix(H, SUBJ, SUBI),  v, sub_vector(w, SUBJ));
      if (k > m-1) { H(k+1, k-1) = T(0); if (k < n-3) H(k+2, k-1) = T(0); }
      
      if (compute_Q)
       	col_house_update(sub_matrix(Q, SUBQ, SUBI),  v, sub_vector(w, SUBQ));

      v1 = H(k+1, k); v2 = H(k+2, k);
      if (k < n-3) v3 = H(k+3, k);
    }
    sub_interval SUBI(n-2,2), SUBJ(0, n), SUBK(n-3,3), SUBL(0, 3);
    v.resize(2); v[0] = v1; v[1] = v2;
    house_vector(v);
    row_house_update(sub_matrix(H, SUBI, SUBK), v, sub_vector(w, SUBL));
    col_house_update(sub_matrix(H, SUBJ, SUBI), v, sub_vector(w, SUBJ));
    if (compute_Q)
      col_house_update(sub_matrix(Q, SUBQ, SUBI), v, sub_vector(w, SUBQ));
  }
    


  /* ********************************************************************* */
  /*    Implicit QR algorithm.                                             */
  /* ********************************************************************* */

  // QR method for real or complex square matrices based on an
  // implicit QR factorisation. eigval has to be a complex vector
  // if A has complex eigeinvalues. complexity about 10n^3, 25n^3 if
  // eigenvectors are computed
  template <class MAT1, class VECT, class MAT2>
    void implicit_qr_algorithm(const MAT1 &A, const VECT &eigval_,
			       const MAT2 &eigvect_,
			       double tol = gmm::default_tol(typename
		   linalg_traits<MAT1>::value_type()), bool compvect = true) {
    VECT &eigval = const_cast<VECT &>(eigval_);
    MAT2 &eigvect = const_cast<MAT2 &>(eigvect_);
    typedef typename linalg_traits<MAT1>::value_type value_type;

    size_type n = mat_nrows(A), q = 0, q_old, p, ite = 0, its = 0;
    dense_matrix<value_type> H(n,n);

    gmm::copy(A, H);
    Hessenberg_reduction(H, eigvect, compvect);
    stop_criterion(H, p, q, tol);

    while (q < n) {
      sub_interval SUBI(p, n-p-q), SUBJ(0,n);
//       Francis_qr_step(sub_matrix(H, SUBI),
// 		      sub_matrix(eigvect, SUBJ, SUBI), compvect);
      Wilkinson_double_shift_qr_step(sub_matrix(H, SUBI), 
				     sub_matrix(eigvect, SUBJ, SUBI),
				     tol, (its == 10 || its == 20), compvect);

      q_old = q;
      stop_criterion(H, p, q, tol);
      if (q != q_old) its = 0;
      ++its;
      if (++ite > n*100) DAL_THROW(failure_error, "QR algorithm failed");
    }
    extract_eig(H, eigval, tol);
  }


  template <class MAT1, class VECT>
    void implicit_qr_algorithm(const MAT1 &a, VECT &eigval, 
			       double tol = gmm::default_tol(typename 
			       linalg_traits<MAT1>::value_type())) {
    dense_matrix<typename linalg_traits<MAT1>::value_type> m(0,0);
    implicit_qr_algorithm(a, eigval, m, tol, false); 
  }

  /* ********************************************************************* */
  /*    Implicit symmetric QR step with Wilkinson Shift.                   */
  /* ********************************************************************* */

  template <class MAT1, class MAT2> 
    void symmetric_Wilkinson_qr_step(const MAT1& TT, const MAT2 &ZZ,
				     bool compute_z) {
    MAT1& T = const_cast<MAT1&>(TT); MAT2& Z = const_cast<MAT2&>(ZZ);
    typedef typename linalg_traits<MAT1>::value_type value_type;
    // to be optimized (use a real tridiag matrix).
    
    size_type n = mat_nrows(T);
    value_type d = (T(n-2, n-2) - T(n-1, n-1)) / value_type(2);
    value_type e = dal::sqr(T(n-1, n-2));
    value_type mu = T(n-1, n-1) - e / (d + dal::sgn(d) * sqrt(d*d + e));
    value_type x = T(0,0) - mu, z = T(1, 0), c, s;

    for (size_type k = 1; k < n; ++k) {
      Givens_rotation(x, z, c, s);

      if (k > 1) Apply_Givens_rotation_left(T(k-1,k-2), T(k,k-2), c, s);
      Apply_Givens_rotation_left(T(k-1,k-1), T(k,k-1), c, s);
      Apply_Givens_rotation_left(T(k-1,k  ), T(k,k  ), c, s);
      if (k < n-1) Apply_Givens_rotation_left(T(k-1,k+1), T(k,k+1), c, s);
      
      if (k > 1) Apply_Givens_rotation_right(T(k-2,k-1), T(k-2,k), c, s);
      Apply_Givens_rotation_right(T(k-1,k-1), T(k-1,k), c, s);
      Apply_Givens_rotation_right(T(k  ,k-1), T(k,k)  , c, s);
      if (k < n-1) Apply_Givens_rotation_right(T(k+1,k-1), T(k+1,k), c, s);
      
      if (compute_z) col_rot(Z, c, s, k-1, k);
      if (k < n-1) { x = T(k, k-1); z = T(k+1, k-1); }
    }

  }

  /* ********************************************************************* */
  /*    Implicit QR algorithm for symmetric matrices.                      */
  /* ********************************************************************* */

  // implicit QR method for real square symmetric matrices.
  // eigval has to be a complex vector if A has complex eigeinvalues.
  // complexity about 4n^3/3, 9n^3 if eigenvectors are computed
  template <class MAT1, class VECT, class MAT2>
  void symmetric_qr_algorithm(const MAT1 &A, const VECT &eigval_,
			      const MAT2 &eigvect_, double tol =
		 gmm::default_tol(typename linalg_traits<MAT1>::value_type()),
			      bool compvect = true) {
    VECT &eigval = const_cast<VECT &>(eigval_);
    MAT2 &eigvect = const_cast<MAT2 &>(eigvect_);
    typedef typename linalg_traits<MAT1>::value_type value_type;

    size_type n = mat_nrows(A), q = 0, p, ite = 0;
    dense_matrix<value_type> T(n,n);
    gmm::copy(A, T);
    Householder_tridiagonalization(T, eigvect, compvect);

//     dense_matrix<value_type> aux1(n,n), aux2(n,n);
//     gmm::mult(eigvect, T, aux1);
//     gmm::mult(aux1, conjugated(eigvect), aux2);
//     gmm::add(scaled(A, -1), aux2);
//     cout << "Ca donne : " << mat_norm2(aux2) << endl;

    symmetric_stop_criterion(T, p, q, tol);

    while (q < n) {
      sub_interval SUBI(p, n-p-q), SUBJ(0,n);
      symmetric_Wilkinson_qr_step(sub_matrix(T, SUBI), 
				  sub_matrix(eigvect, SUBJ, SUBI), compvect);
      symmetric_stop_criterion(T, p, q, tol);
      if (++ite > n*100) DAL_THROW(failure_error, "QR algorithm failed");
    }
    extract_eig(T, eigval, tol);
  }


  template <class MAT1, class VECT>
    void symmetric_qr_algorithm(const MAT1 &a, VECT &eigval, 
				double tol = gmm::default_tol(typename
				linalg_traits<MAT1>::value_type())) {
    dense_matrix<typename linalg_traits<MAT1>::value_type> m(0,0);
    symmetric_qr_algorithm(a, eigval, m, tol, false);
  }


}

#endif

