/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_diagonal.h : diagonal preconditoner.             */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
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
#ifndef GMM_PRECOND_DIAGONAL_H
#define GMM_PRECOND_DIAGONAL_H

#include <gmm_solvers.h>

namespace gmm {

  template<typename Matrix> struct diagonal_precond {
    typedef typename linalg_traits<Matrix>::value_type value_type;
    std::vector<value_type> diag;

    diagonal_precond(const Matrix &M) : diag(mat_nrows(M))
    { for (size_type i = 0; i < mat_nrows(M); ++i) diag[i] = 1.0 / M(i, i); }

    diagonal_precond(void) {}
  };

  template <typename Matrix, typename V2> inline
  void mult_diag_p(const diagonal_precond<Matrix>& P, V2 &v2, abstract_sparse){
    typename linalg_traits<V2>::iterator it = vect_begin(v2),
      ite = vect_end(v2);
    for (; it != ite; ++it) *it *= P.diag[it.index()];
  }

  template <typename Matrix, typename V2> inline
  void mult_diag_p(const diagonal_precond<Matrix>& P,V2 &v2, abstract_skyline)
    { mult_diag_p(P, v2, abstract_sparse()); }

  template <typename Matrix, typename V2> inline
  void mult_diag_p(const diagonal_precond<Matrix>& P, V2 &v2, abstract_dense){
    for (size_type i = 0; i < P.diag.size(); ++i) v2[i] *= P.diag[i];
  }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const diagonal_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.diag.size() != vect_size(v2))
      DAL_THROW(dimension_error, "dimensions mismatch");
    copy(v1, v2);
    mult_diag_p(P, v2, typename linalg_traits<V2>::storage_type());
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const diagonal_precond<Matrix>& P,const V1 &v1,V2 &v2) {
    mult(P, v1, v2);
  }
  
  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const diagonal_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.diag.size() != vect_size(v2))
      DAL_THROW(dimension_error, "dimensions mismatch");
    copy(v1, v2);
    for (size_type i = 0; i < P.diag.size(); ++i)
      v2[i] *= sqrt(dal::abs(P.diag[i]));
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const diagonal_precond<Matrix>& P,
			    const V1 &v1, V2 &v2)
    { left_mult(P, v1, v2); }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const diagonal_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.diag.size() != vect_size(v2))
      DAL_THROW(dimension_error, "dimensions mismatch");
    copy(v1, v2);
    for (size_type i = 0; i < P.diag.size(); ++i)
      if (P.diag[i] < 0.0) 
	v2[i] *= -sqrt(dal::abs(P.diag[i]));
      else
	v2[i] *= sqrt(dal::abs(P.diag[i]));
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const diagonal_precond<Matrix>& P,
			    const V1 &v1, V2 &v2)
    { right_mult(P, v1, v2); }


}

#endif 

