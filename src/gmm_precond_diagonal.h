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
#ifndef GMM_PRECOND_DIAGONAL_H
#define GMM_PRECOND_DIAGONAL_H

#include <gmm_solvers.h>

namespace gmm {

  template<class T> struct diagonal_precond {
    std::vector<T> diag;

    template<class Matrix> diagonal_precond(const Matrix &M)
      : diag(mat_nrows(M)) {
      for (size_type i = 0; i < mat_nrows(M); ++i)
	diag[i] = 1.0 / A(i, i);
    }
  };

  template<class T>
  const &diagonal_precond<T> &gmm::transposed(const &diagonal_precond<T> &P)
  { return P; }

  template <class T, class V1, class V2> inline
  void mult(const diagonal_precond<T>& P, const V1 &v1, V2 &v2) {
    if (P.diag.size() != vect_size(v2))
      DAL_THROW(dimension_error, "dimensions mismatch");
    copy(v1, v2); 
    for (size_type i = 0; i < P.diag.size(); ++i)
      v2[i] *= p.diag[i];
  }
}

#endif 

