/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_def.h : basic definitions of G.M.M.                      */
/*     									   */
/* Date : October 13, 2002.                                                */
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

#ifndef __GMM_DEF_H
#define __GMM_DEF_H

#include <dal_ref.h>
#include <complex>

namespace gmm {


  typedef size_t size_type;
  using dal::dimension_error;
  using dal::file_not_found_error;
  using dal::internal_error;
  using dal::not_linear_error;
  using dal::to_be_done_error;
  using dal::failure_error;

  /* ******************************************************************** */
  /*		Specifier types                             		  */
  /* ******************************************************************** */

  

  struct abstract_null_type {}; // specify an information lake.

  struct linalg_true {};
  struct linalg_false {};

  struct abstract_vector {};
  struct abstract_matrix {};
  
  struct abstract_sparse {};    // sparse matrix or vector
  struct abstract_plain {};     // plain matrix or vector
  struct abstract_indirect {};  // matrix given by the product with a vector

  struct row_major {};          // matrix with a row access.
  struct col_major {};       // matrix with a column access
  struct row_and_col {};     // both accesses but row preference
  struct col_and_row {};     // both accesses but column preference

  template <class T> struct transposed_type;
  template<> struct transposed_type<row_major>   {typedef col_major   t_type;};
  template<> struct transposed_type<col_major>   {typedef row_major   t_type;};
  template<> struct transposed_type<row_and_col> {typedef col_and_row t_type;};
  template<> struct transposed_type<col_and_row> {typedef row_and_col t_type;};

  template <class T> struct principal_orientation_type;
  template<> struct principal_orientation_type<row_major>
  { typedef row_major potype; };
  template<> struct principal_orientation_type<col_major>
  { typedef col_major potype; };
  template<> struct principal_orientation_type<row_and_col>
  { typedef row_major potype; };
  template<> struct principal_orientation_type<col_and_row>
  { typedef col_major potype; };

  template <class V> struct linalg_traits;

  template <class PT, class V> struct vect_ref_type;
  template <class P, class V> struct vect_ref_type<P *, V> {
    typedef typename linalg_traits<V>::reference_type access_type;
    typedef typename linalg_traits<V>::iterator iterator;
  };
  template <class P, class V> struct vect_ref_type<const P *, V> {
    typedef typename linalg_traits<V>::value_type access_type;
    typedef typename linalg_traits<V>::const_iterator iterator;
  };
  template <class PT, class V> struct mat_ref_type;
  template <class P, class V> struct mat_ref_type<P *, V> {
    typedef typename linalg_traits<V>::reference_type access_type;
  };
  template <class P, class V> struct mat_ref_type<const P *, V> {
    typedef typename linalg_traits<V>::value_type access_type;
  };
  
  template <class PT> struct const_pointer;
  template <class P> struct const_pointer<P *>
  { typedef const P* pointer; };
  template <class P> struct const_pointer<const P *>
  { typedef const P* pointer; };
  
  /* ******************************************************************** */
  /*		Operations on scalars                         		  */
  /* ******************************************************************** */

  template <class T> struct number_traits
  { typedef T magnitude_type; };
 
  template <class T> struct number_traits<std::complex<T> >
  { typedef T magnitude_type; };

  template <class T> T conj_product(T a, T b) { return a * b; }
  template <class T> std::complex<T> conj_product(const std::complex<T> &a,
						  const std::complex<T> &b)
  { return std::conj(a) * b; } // to be optimized
  
  template <class T> T modulus(T a)
  { return dal::abs(a); }
  template <class T> T modulus(const std::complex<T> &a)
  { return std::abs(a); }

  /* ******************************************************************** */
  /*   Selects a temporary vector type                                    */
  /*   V if V is a valid vector type,                                     */
  /*   svector if V is a reference on a sparse vector,                    */
  /*   std::vector if V is a reference on a plain vector.                 */
  /* ******************************************************************** */

  template <class T> class wsvector;
  template <class R, class S, class L, class V> struct _temporary_vector {};
  template <class V, class L>
  struct _temporary_vector<linalg_true, abstract_sparse, L, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V, class L>
  struct _temporary_vector<linalg_true, abstract_plain, L, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class S, class V>
  struct _temporary_vector<linalg_false, S, abstract_vector, V>
  { typedef V vector_type; };
  template <class V>
  struct _temporary_vector<linalg_false, abstract_plain, abstract_matrix, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_vector<linalg_false, abstract_sparse, abstract_matrix, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };

  template <class V> struct temporary_vector {
    typedef typename _temporary_vector<typename linalg_traits<V>::is_reference,
				       typename linalg_traits<V>::storage_type,
				       typename linalg_traits<V>::linalg_type,
				       V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary plain vector type                              */
  /*   V if V is a valid plain vector type,                               */
  /*   std::vector if V is a reference or a sparse vector                 */
  /* ******************************************************************** */

  template <class R, class S, class V> struct _temporary_plain_vector;
  template <class S, class V> struct _temporary_plain_vector<linalg_true, S, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_plain_vector<linalg_false, abstract_sparse, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_plain_vector<linalg_false, abstract_plain, V>
  { typedef V vector_type; };

  template <class V> struct temporary_plain_vector {
    typedef typename _temporary_plain_vector<typename
    linalg_traits<V>::is_reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary sparse vector type                             */
  /*   V if V is a valid sparse vector type,                              */
  /*   std::vector if V is a reference or a plain vector                  */
  /* ******************************************************************** */

  template <class R, class S, class V> struct _temporary_sparse_vector;
  template <class S, class V>
  struct _temporary_sparse_vector<linalg_true, S, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };
  template <class V>
  struct _temporary_sparse_vector<linalg_false, abstract_sparse, V>
  { typedef V vector_type; };
  template <class V>
  struct _temporary_sparse_vector<linalg_false, abstract_plain, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };

  template <class V> struct temporary_sparse_vector {
    typedef typename _temporary_sparse_vector<typename
    linalg_traits<V>::is_reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };


}

#endif //  __GMM_DEF_H
