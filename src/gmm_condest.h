/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_condest.h : condition number estimation                  */
/*     				                		           */
/* Date : August 27, 2003.                                                 */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr,                     */
/*          Julien Pommier, Julien.Pommier@gmm.insa-tlse.fr.               */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard, Julien Pommier.                        */
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

#ifndef __GMM_CONDEST_H
#define __GMM_CONDEST_H

#include "gmm_dense_qr.h"

namespace gmm {

//  /** estimation of the magnitude of the largest eigenvalue 
//   *  works also with non-square matrices
//   */
//   template <typename MAT> 
//   typename number_traits<typename
//   linalg_traits<MAT>::value_type>::magnitude_type
//   norm_lin2_est(const MAT& M) {
//     typedef typename number_traits<typename
//       linalg_traits<MAT>::value_type>::magnitude_type magnitude_type;
//     typedef typename linalg_traits<MAT>::value_type T;
//     typedef typename temporary_dense_vector<MAT>::vector_type vector_type;
    
//     int d = mat_nrows(M) - mat_ncols(M);
//     int vsz = std::min(mat_nrows(M), mat_ncols(M));
//     int vsz2 = std::max(mat_nrows(M), mat_ncols(M));
//     vector_type v(vsz), tmp(vsz), tmp2(vsz2); fill_random(v);
//     magnitude_type e = vect_norm2(v), e0 = 0, pert = 1E-1;
//     while (gmm::abs(e-e0) > 1e-6 * e0) {
//       e0 = e;
//       if (d == 0)
// 	mult(B,v,tmp);
//       else if (d > 0)
// 	{ mult(B,v,tmp2); mult(transposed(B),tmp2,tmp); }
//       else 
// 	{ mult(transposed(B),v,tmp2); mult(B,tmp2,tmp); }
      
//       e = vect_norm2(tmp);
//       scale(tmp, 1.0 / e);
//       copy(tmp,v);
//       clear(tmp); fill_random(tmp, 0.05);
//       add(scaled(tmp, pert), v);
//       pert /= 10.0;
//     }
//     if (d == 0)
//       return e;
//     else return sqrt(e);
//  }

//  /** estimation of the condition number 
//   * (using lu_inverse => dense matrix only)
//   */
//   template <typename MAT> 
//   typename number_traits<typename
//   linalg_traits<MAT>::value_type>::magnitude_type
//   condest(const MAT& M) {
//     typedef typename linalg_traits<MAT>::value_type T;
    
//     int d = mat_nrows(M) - mat_ncols(M);
//     int vsz = std::min(mat_nrows(M), mat_ncols(M));
//     dense_matrix<T> B(vsz, vsz);
//     if (d == 0)
//       copy(M, B);
//     else if (d > 0)
//       mult(transposed(M), M, B);
//     else 
//       mult(M,transposed(M), B);

//     gmm::lu_inverse(B);
//     return  norm_lin2_est(M) 
//       * ((d == 0) ? norm_lin2_est(B) : sqrt(norm_lin2_est(B)));
//   }


  /** estimation of the condition number 
   * (using symmetric_qr_algorithm => dense matrix only)
   */

  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condest(const MAT& M, 
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emin,
	  typename number_traits<typename
	  linalg_traits<MAT>::value_type>::magnitude_type& emax) {
    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename number_traits<T>::magnitude_type magnitude_type;

    int d = mat_nrows(M) - mat_ncols(M);
    int vsz = std::min(mat_nrows(M), mat_ncols(M));
    
    dense_matrix<T> B(vsz, vsz);
    std::vector<magnitude_type> eig(vsz);
    if (d >= 0) mult(transposed(M), M, B); else mult(M, transposed(M), B);
    
    gmm::symmetric_qr_algorithm(B, eig);
    
    emin = emax = gmm::abs(eig[0]);
    for (int i = 1; i < vsz; ++i) {
      emin = std::min(emin, gmm::abs(eig[i]));
      emax = std::max(emax, gmm::abs(eig[i]));
    }
    return sqrt(emax / emin);
  }
  
  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condest(const MAT& M) { 
    typename number_traits<typename
      linalg_traits<MAT>::value_type>::magnitude_type emax, emin;
    return condest(M, emax, emin);
  }
}

#endif
