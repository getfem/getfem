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
  
  template <class T>  T mat_det(gmm::dense_matrix<T> &M) {
    T *p = &(M(0,0));
    switch (mat_nrows(M)) {
    case 1 : return (*p);
    case 2 : return (*p) * (*(p+3)) - (*(p+1)) * (*(p+2));
    case 3 : return (*p) * ((*(p+4)) * (*(p+8)) - (*(p+5)) * (*(p+7)))
	       - (*(p+1)) * ((*(p+3)) * (*(p+8)) - (*(p+5)) * (*(p+6)))
	       + (*(p+2)) * ((*(p+3)) * (*(p+7)) - (*(p+4)) * (*(p+6)));
    default : return gmm::lu_det(M);
    }
    return T(0);
  }

  template <class T> T mat_inverse(gmm::dense_matrix<T> &M) {
    size_type N = mat_nrows(M);
    T *p = &(M(0,0)), det(0);
    if (N <= 3) {
      switch (N) {
      case 1 : det = *p; *p = T(1) / det; break;
      case 2 : det = (*p) * (*(p+3)) - (*(p+1)) * (*(p+2));
	std::swap(*p, *(p+3));
	*p++ /= det; *p++ /= -det; *p++ /= -det; *p++ /= det; break;
      case 3 :
	{
	  T a, b, c, d, e, f, g, h, i;
	  a =   (*(p+4)) * (*(p+8)) - (*(p+5)) * (*(p+7));
	  b = - (*(p+1)) * (*(p+8)) + (*(p+2)) * (*(p+7));
	  c =   (*(p+1)) * (*(p+5)) - (*(p+2)) * (*(p+4));
	  d = - (*(p+3)) * (*(p+8)) + (*(p+5)) * (*(p+6));
	  e =   (*(p+0)) * (*(p+8)) - (*(p+2)) * (*(p+6));
	  f = - (*(p+0)) * (*(p+5)) + (*(p+2)) * (*(p+3));
	  g =   (*(p+3)) * (*(p+7)) - (*(p+4)) * (*(p+6));
	  h = - (*(p+0)) * (*(p+7)) + (*(p+1)) * (*(p+6));
	  i =   (*(p+0)) * (*(p+4)) - (*(p+1)) * (*(p+3));
	  det = (*p) * a + (*(p+1)) * d + (*(p+2)) * g;
	  *p++ = a / det; *p++ = b / det; *p++ = c / det; 
	  *p++ = d / det; *p++ = e / det; *p++ = f / det; 
	  *p++ = g / det; *p++ = h / det; *p++ = i / det; 
	}
      }
    }
    else det = gmm::lu_inverse(M);
    return det;
  }

  typedef gmm::dense_matrix<scalar_type> base_matrix;


}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_MATRIX_H */
