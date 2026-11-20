/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2003-2026 Yves Renard, Julien Pommier

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program. If not, see https://www.gnu.org/licenses/.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file gmm_condition_number.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date August 27, 2003.
   @brief computation of the condition number of dense matrices.
*/
#ifndef GMM_CONDITION_NUMBER_H__
#define GMM_CONDITION_NUMBER_H__

#include "gmm_dense_qr.h"

namespace gmm {

  /* Scalar type used internally for the eigenvalue computation of
     condition_number. Single precision is promoted to double precision:
     with the default tolerance of the single precision type, the
     symmetric QR algorithm may fail to converge on moderately large
     matrices ("QR algorithm failed."). */
  template <typename T> struct condition_number_work_type
  { typedef T type; };
  template <> struct condition_number_work_type<float>
  { typedef double type; };
  template <> struct condition_number_work_type<std::complex<float> >
  { typedef std::complex<double> type; };

  /** computation of the condition number of matrices (ratio of the
      extreme singular values).

      Uses symmetric_qr_algorithm on a dense copy of the matrix. Single
      precision matrices are handled in double precision.

      @param M a matrix.
      @param emin smallest (in magnitude) eigenvalue
      @param emax largest eigenvalue.
   */
  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condition_number(const MAT& M, 
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emin,
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emax) {
    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    typedef typename condition_number_work_type<T>::type TW;
    typedef typename number_traits<TW>::magnitude_type RW;

    // Added because of errors in complex with zero det
    if (sizeof(T) != sizeof(R) && gmm::abs(gmm::lu_det(M)) == R(0))
      return  gmm::default_max(R());
      
    size_type m = mat_nrows(M), n = mat_ncols(M);
    emax = emin = R(0);
    std::vector<RW> eig(m+n);

    if (m+n == 0) return R(0);
    if (is_hermitian(M)) {
      dense_matrix<TW> A(m, m);
      gmm::copy(M, A);
      eig.resize(m);
      gmm::symmetric_qr_algorithm(A, eig);
    }
    else {
      dense_matrix<TW> B(m+n, m+n); // not very efficient ??
      gmm::copy(conjugated(M), sub_matrix(B, sub_interval(m, n), sub_interval(0, m)));
      gmm::copy(M, sub_matrix(B, sub_interval(0, m),
					  sub_interval(m, n)));
      gmm::symmetric_qr_algorithm(B, eig);
    }
    RW eminw = gmm::abs(eig[0]), emaxw = eminw;
    for (size_type i = 1; i < eig.size(); ++i) {
      RW e = gmm::abs(eig[i]); 
      eminw = std::min(eminw, e);
      emaxw = std::max(emaxw, e);
    }
    emin = R(eminw); emax = R(emaxw);
    if (eminw == RW(0)) return gmm::default_max(R());
    RW c = emaxw / eminw;
    // avoid an overflow when converting back to single precision
    if (c >= RW(gmm::default_max(R()))) return gmm::default_max(R());
    return R(c);
  }

  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condition_number(const MAT& M) { 
    typename number_traits<typename
      linalg_traits<MAT>::value_type>::magnitude_type emax, emin;
    return condition_number(M, emin, emax);
  }

  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  Frobenius_condition_number_sqr(const MAT& M) { 
    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    size_type m = mat_nrows(M), n = mat_ncols(M);
    dense_matrix<T> B(std::min(m,n), std::min(m,n));
    if (m < n) mult(M,gmm::conjugated(M),B);
    else       mult(gmm::conjugated(M),M,B);
    R trB = abs(mat_trace(B));
    lu_inverse(B);
    return trB*abs(mat_trace(B));
  }

  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  Frobenius_condition_number(const MAT& M)
  { return sqrt(Frobenius_condition_number_sqr(M)); }

  /** estimation of the condition number (TO BE DONE...)
   */
  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condest(const MAT& M, 
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emin,
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emax) {
    return condition_number(M, emin, emax);
  }
  
  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condest(const MAT& M) { 
    typename number_traits<typename
      linalg_traits<MAT>::value_type>::magnitude_type emax, emin;
    return condest(M, emin, emax);
  }
}

#endif
