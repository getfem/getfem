/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_mtl.h : Interface with some mtl types.               */
/*                                                                         */
/*     									   */
/* Date : February 6, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __GENSOLV_MTL_H
#define __GENSOLV_MTL_H

#include <mtl/mtl.h>
#include <mtl/matrix.h>
#include <mtl/utils.h>
#include <gensolv_secant.h>

namespace gensolv
{ 

  template<class T> struct mtl_sym_sparse_matrix
  {
    typedef typename mtl::matrix<T,
             mtl::symmetric<mtl::upper>, 
             mtl::array< mtl::compressed<> >,
             mtl::row_major >::type type;
  };

  template<class T> struct mtl_gen_sparse_matrix
  {
    typedef typename mtl::matrix<T,
             mtl::rectangle<>, 
             mtl::array< mtl::compressed<> >,
             mtl::row_major >::type type;
  };

  template<class T> struct mtl_dense_vector
  {
    typedef mtl::dense1D<T> type;
  };

  template<class T> T dot(const mtl::dense1D<T> &x,
			  const mtl::dense1D<T> &y)
    { return mtl::dot(x, y); }
  // template<class T> void copy(const mtl::dense1D<T> &x,
  //                             const mtl::dense1D<T> &y)
  // { mtl::copy(x, y); }
  template<class T> void add_mult(mtl::dense1D<T> &x, T a,
				  const mtl::dense1D<T> &y)
  { mtl::add(x, mtl::scaled(y, a), x); }
  template<class T> void mult(mtl::dense1D<T> &x, T a)
  { mtl::scale(x, a); }
  template<class T> void add(mtl::dense1D<T> &x, const mtl::dense1D<T> &y)
  { mtl::add(x, y, x); }

}  /* end of namespace gensolv.                                           */


#endif /* __GENSOLV_MTL_H                                                 */
