/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot) version 1.0                    */
/* File    :  bgeot_config.h : basic configuration.                        */
/*     									   */
/* Date : December 20, 1999.                                               */
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


#ifndef __BGEOT_CONFIG_H
#define __BGEOT_CONFIG_H

#include <bgeot_tensor.h>
#include <bgeot_poly.h>

namespace bgeot
{
  static const size_t ST_NIL = size_t(-1);
  typedef dal::uint8_type dim_type;
  typedef dal::uint16_type short_type;
  typedef size_t size_type;
  typedef double scalar_type;

  typedef vsvector<scalar_type> base_vector;
  typedef vsmatrix<scalar_type> base_matrix;
  typedef tensor<scalar_type> base_tensor;
  typedef polynomial<scalar_type> base_poly;
  typedef PT<base_vector> base_node;


  /* usual constant polynomials  */

  inline base_poly null_poly(short_type n)
    { return base_poly(n, 0); }
  inline base_poly one_poly(short_type n)
    { base_poly res=base_poly(n, 0); res.one(); return res;  }
  inline base_poly one_var_poly(short_type n, short_type k)
    { return base_poly(n, 1, k); }


}  /* end of namespace bgeot.                                             */


#endif /* __BGEOT_CONFIG_H                                                */
