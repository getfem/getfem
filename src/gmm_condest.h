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

#ifndef __GMM_CONDEST_H
#define __GMM_CONDEST_H

namespace gmm {

  // pas vraiment des estimations en l'état ...

  /** estimation of the magnitude of the largest eigenvalue 
   *  works also with non-square matrices
   */
  template <typename MAT> 
  typename number_traits<typename
  linalg_traits<MAT>::value_type>::magnitude_type
  norm_lin2_est(const MAT& M) {
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    
    int d = mat_nrows(M) - mat_ncols(M);
    int vsz = std::min(mat_nrows(M), mat_ncols(M));
    
    dense_matrix<value_type> B(vsz, vsz);
    std::vector<magnitude_type> eig(vsz);
    if (d >= 0) mult(transposed(M), M, B); else mult(M,transposed(M), B);
    
    gmm::symmetric_qr_algorithm(B, eig);
    
    magnitude_type emax = dal::abs(eig[0]);
    for (int i = 1; i < vsz; ++i) emax = std::max(emax, eig[i]); 
    
    return sqrt(emax);
  }

  /** estimation of the condition number 
   * (using impicit_qr_algorithm => dense matrix only)
   */

  template <typename MAT> 
  typename number_traits<typename 
  linalg_traits<MAT>::value_type>::magnitude_type
  condest(const MAT& M) {
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;

    int d = mat_nrows(M) - mat_ncols(M);
    int vsz = std::min(mat_nrows(M), mat_ncols(M));
    
    dense_matrix<value_type> B(vsz, vsz);
    std::vector<magnitude_type> eig(vsz);
    if (d >= 0) mult(transposed(M), M, B); else mult(M,transposed(M), B);
    
    gmm::symmetric_qr_algorithm(B, eig);
    
    magnitude_type emin = dal::abs(eig[0]), emax = dal::abs(eig[0]);
    for (int i = 1; i < vsz; ++i)
      { emin = std::min(emin, eig[i]); emax = std::max(emax, eig[i]); } 
    
    return sqrt(emax / emin);
  }
  

}

#endif
