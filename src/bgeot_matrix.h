/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_matrix.h : small matrices.                             */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1995-2002  Yves Renard.                                   */
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


#ifndef __BGEOT_MATRIX_H
#define __BGEOT_MATRIX_H

#include <bgeot_vector.h>

namespace bgeot
{

  // Some further functions on gmm::dense_matrix<T>

  template<class T> vsvector<T>& operator *=(vsvector<T>& v,
					     const gmm::dense_matrix<T>& m)
  { gmm::mult(m, v, v); return v; }

  template<class T>
    vsvector<T> operator *(const gmm::dense_matrix<T>& m, const vsvector<T>& v)
  { vsvector<T> res(m.nrows()); gmm::mult(m, v, res); return res; }

  template<class T>
  gmm::dense_matrix<T> operator *(const gmm::dense_matrix<T>& m,
				  const gmm::dense_matrix<T>& n) {
    gmm::dense_matrix<T> res(m.nrows(), n.ncols()); 
    gmm::mult(m, n, res);
    return res;
  }

  typedef gmm::dense_matrix<scalar_type> base_matrix;


}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_MATRIX_H */
