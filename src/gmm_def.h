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
/* Copyright (C) 2002  Yves Renard.                                        */
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
  /* not perfectly null, required by aCC 3.33                             */
  struct abstract_null_type { 
    abstract_null_type(int=0) {}
    template <typename A,typename B,typename C> void operator()(A,B,C) {}
  }; // specify an information lake.

  struct linalg_true {};
  struct linalg_false {};

  struct linalg_const {};       // A reference is either linalg_const,
  struct linalg_modifiable {};  //  linalg_modifiable or linalg_false.

  struct abstract_vector {};    // The object is a vector
  struct abstract_matrix {};    // The object is a matrix
  
  struct abstract_sparse {};    // sparse matrix or vector
  struct abstract_skyline {};   // 'sky-line' matrix or vector
  struct abstract_dense {};     // dense matrix or vector
  struct abstract_indirect {};  // matrix given by the product with a vector

  struct row_major {};          // matrix with a row access.
  struct col_major {};          // matrix with a column access
  struct row_and_col {};        // both accesses but row preference
  struct col_and_row {};        // both accesses but column preference

  template <typename T> struct transposed_type;
  template<> struct transposed_type<row_major>   {typedef col_major   t_type;};
  template<> struct transposed_type<col_major>   {typedef row_major   t_type;};
  template<> struct transposed_type<row_and_col> {typedef col_and_row t_type;};
  template<> struct transposed_type<col_and_row> {typedef row_and_col t_type;};

  template <typename T> struct principal_orientation_type
  { typedef abstract_null_type potype; };
  template<> struct principal_orientation_type<row_major>
  { typedef row_major potype; };
  template<> struct principal_orientation_type<col_major>
  { typedef col_major potype; };
  template<> struct principal_orientation_type<row_and_col>
  { typedef row_major potype; };
  template<> struct principal_orientation_type<col_and_row>
  { typedef col_major potype; };

  //  template <typename V> struct linalg_traits;
  template <typename V> struct linalg_traits {    
    typedef abstract_null_type this_type;
    typedef abstract_null_type linalg_type;
    typedef abstract_null_type value_type;
    typedef abstract_null_type is_reference;
    typedef abstract_null_type& reference;
    typedef abstract_null_type* iterator;
    typedef const abstract_null_type* const_iterator;
    typedef abstract_null_type storage_type;
    typedef abstract_null_type clear_type;
    typedef abstract_null_type access_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_row_iterator;
    typedef abstract_null_type row_iterator;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_col_iterator;
    typedef abstract_null_type col_iterator;
    typedef abstract_null_type sub_orientation;    
  };

  template <typename PT, typename V> struct vect_ref_type;
  template <typename P, typename V> struct vect_ref_type<P *, V> {
    typedef typename linalg_traits<V>::reference access_type;
    typedef typename linalg_traits<V>::iterator iterator;
  };
  template <typename P, typename V> struct vect_ref_type<const P *, V> {
    typedef typename linalg_traits<V>::value_type access_type;
    typedef typename linalg_traits<V>::const_iterator iterator;
  };
  
  template <typename PT> struct const_pointer;
  template <typename P> struct const_pointer<P *>
  { typedef const P* pointer; };
  template <typename P> struct const_pointer<const P *>
  { typedef const P* pointer; };

  template <typename PT> struct modifiable_pointer;
  template <typename P> struct modifiable_pointer<P *>
  { typedef P* pointer; };
  template <typename P> struct modifiable_pointer<const P *>
  { typedef P* pointer; };

  inline bool _is_sparse(abstract_sparse)   { return true;  }
  inline bool _is_sparse(abstract_dense)    { return false; }
  inline bool _is_sparse(abstract_skyline)  { return true;  }
  inline bool _is_sparse(abstract_indirect) { return false; }

  template <typename L> inline bool is_sparse(const L &) 
  { return _is_sparse(typename linalg_traits<L>::storage_type()); }

  inline bool _is_row_matrix(row_major)     { return true;  }
  inline bool _is_row_matrix(col_major)     { return false; }
  inline bool _is_row_matrix(row_and_col)   { return true;  }
  inline bool _is_row_matrix(col_and_row)   { return true;  }

  template <typename L> inline bool is_row_matrix(const L &) 
  { return _is_row_matrix(typename linalg_traits<L>::sub_orientation()); }

  inline bool _is_col_matrix(row_major)     { return false; }
  inline bool _is_col_matrix(col_major)     { return true;  }
  inline bool _is_col_matrix(row_and_col)   { return true;  }
  inline bool _is_col_matrix(col_and_row)   { return true;  }

  template <typename L> inline bool is_col_matrix(const L &) 
  { return _is_col_matrix(typename linalg_traits<L>::sub_orientation()); }

  inline bool is_col_matrix(row_major) { return false; }
  inline bool is_col_matrix(col_major) { return true; }
  inline bool is_row_matrix(row_major) { return true; }
  inline bool is_row_matrix(col_major) { return false; }

  template <typename L> inline bool is_const_reference(L) { return false; }
  inline bool is_const_reference(linalg_const) { return true; }

  /* ******************************************************************** */
  /*  types to deal with const object representing a modifiable reference */
  /* ******************************************************************** */
  
  template <typename PT, typename R> struct _mref_type 
  { typedef abstract_null_type return_type; };
  template <typename L, typename R> struct _mref_type<L *, R>
  { typedef L & return_type; };
  template <typename L, typename R> struct _mref_type<const L *, R>
  { typedef const L & return_type; };
  template <typename L> struct _mref_type<L *, linalg_const>
  { typedef const L & return_type; };
  template <typename L> struct _mref_type<const L *, linalg_const>
  { typedef const L & return_type; };
  template <typename L> struct _mref_type<const L *, linalg_modifiable>
  { typedef L & return_type; };
  template <typename L> struct _mref_type<L *, linalg_modifiable>
  { typedef L & return_type; };

  template <typename PT> struct mref_type {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _mref_type<PT, 
      typename linalg_traits<L>::is_reference>::return_type return_type;
  };

  template <typename L> typename mref_type<const L *>::return_type 
  linalg_cast(const L &l)
  { return const_cast<typename mref_type<const L *>::return_type>(l); }

  template <typename L> typename mref_type<L *>::return_type linalg_cast(L &l)
  { return const_cast<typename mref_type<L *>::return_type>(l); }

  template <typename L, typename R> struct _cref_type
  { typedef abstract_null_type return_type; };
  template <typename L> struct _cref_type<L, linalg_modifiable>
  { typedef L & return_type; };
  template <typename L> struct cref_type {
    typedef typename _cref_type<L, 
      typename linalg_traits<L>::is_reference>::return_type return_type;
  };

  template <typename L> typename cref_type<L>::return_type 
  linalg_const_cast(const L &l)
  { return const_cast<typename cref_type<L>::return_type>(l); }


  // To be used to select between a reference or a const refercence for
  // the return type of a function
  // select_return<C1, C2, L *> return C1 if L is a const reference,
  //                                   C2 otherwise.
  // select_return<C1, C2, const L *> return C2 if L is a modifiable reference
  //                                         C1 otherwise. 
  template <typename C1, typename C2, typename REF> struct _select_return {
    typedef abstract_null_type return_type;
  };
  template <typename C1, typename C2, typename L>
  struct _select_return<C1, C2, const L &> { typedef C1 return_type; };
  template <typename C1, typename C2, typename L>
  struct _select_return<C1, C2, L &> { typedef C2 return_type; };
  template <typename C1, typename C2, typename PT> struct select_return {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _select_return<C1, C2, 
      typename mref_type<PT>::return_type>::return_type return_type;
  };

  
  // To be used to select between a reference or a const refercence inside
  // a structure or a linagl_traits
  // select_ref<C1, C2, L *> return C1 if L is a const reference,
  //                                C2 otherwise.
  // select_ref<C1, C2, const L *> return C2 in any case. 
  template <typename C1, typename C2, typename REF> struct _select_ref {
    typedef abstract_null_type ref_type;
  };
  template <typename C1, typename C2, typename L>
  struct _select_ref<C1, C2, const L &> { typedef C1 ref_type; };
  template <typename C1, typename C2, typename L>
  struct _select_ref<C1, C2, L &> { typedef C2 ref_type; };
  template <typename C1, typename C2, typename PT> struct select_ref {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _select_ref<C1, C2, 
      typename mref_type<PT>::return_type>::ref_type ref_type;
  };
  template <typename C1, typename C2, typename L> struct select_ref<C1, C2, const L *> {
    typedef C1 ref_type;
  };


  template<typename R> struct _is_a_reference
  { typedef linalg_true reference; };
  template<> struct _is_a_reference<linalg_false>
  { typedef linalg_false reference; };

  template<typename L> struct is_a_reference {
    typedef typename _is_a_reference<typename linalg_traits<L>::is_reference>
      ::reference reference;
  };

  template <typename PT> struct which_reference {
    typedef abstract_null_type is_reference;
  };
  template <typename PT> struct which_reference<PT *>
  { typedef linalg_modifiable is_reference; };
  template <typename PT> struct which_reference<const PT *>
  { typedef linalg_const is_reference; };


  template <typename C1, typename C2, typename R> struct _select_orientation {
    typedef abstract_null_type return_type;
  };
  template <typename C1, typename C2> struct _select_orientation<C1, C2, row_major>
  { typedef C1 return_type; };
   template <typename C1, typename C2> struct _select_orientation<C1, C2, col_major>
   { typedef C2 return_type; };
  template <typename C1, typename C2, typename L> struct select_orientation {
    typedef typename _select_orientation<C1, C2,
      typename principal_orientation_type<typename
      linalg_traits<L>::sub_orientation>::potype>::return_type return_type;
  };
  
  /* ******************************************************************** */
  /*		Operations on scalars                         		  */
  /* ******************************************************************** */

  template <typename T> struct number_traits
  { typedef T magnitude_type; };
 
  template <typename T> struct number_traits<std::complex<T> >
  { typedef T magnitude_type; };

  template <typename T> inline T conj_product(T a, T b) { return a * b; }
  template <typename T> inline std::complex<T> conj_product(std::complex<T> a,
							 std::complex<T> b)
  { return std::conj(a) * b; } // to be optimized ?

  template <typename T> inline bool is_complex(T a) { return false; }
  template <typename T> inline bool is_complex(std::complex<T> a) { return true; }
  

  /* ******************************************************************** */
  /*		Basic vectors used                         		  */
  /* ******************************************************************** */
  
  template<typename T> struct dense_vector_type 
  { typedef std::vector<T> vector_type; };

  template <typename T> class wsvector;
  template<typename T> struct sparse_vector_type 
  { typedef wsvector<T> vector_type; };

  template <typename T> class slvector;
  template <typename T> class dense_matrix;
  template <typename VECT> class row_matrix;
  

  /* ******************************************************************** */
  /*   Selects a temporary vector type                                    */
  /*   V if V is a valid vector type,                                     */
  /*   wsvector if V is a reference on a sparse vector,                   */
  /*   std::vector if V is a reference on a dense vector.                 */
  /* ******************************************************************** */

  
  template <typename R, typename S, typename L, typename V> struct _temporary_vector {
    typedef abstract_null_type vector_type;
  };
  template <typename V, typename L>
  struct _temporary_vector<linalg_true, abstract_sparse, L, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V, typename L>
  struct _temporary_vector<linalg_true, abstract_skyline, L, V>
  { typedef slvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V, typename L>
  struct _temporary_vector<linalg_true, abstract_dense, L, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename S, typename V>
  struct _temporary_vector<linalg_false, S, abstract_vector, V>
  { typedef V vector_type; };
  template <typename V>
  struct _temporary_vector<linalg_false, abstract_dense, abstract_matrix, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_vector<linalg_false, abstract_sparse, abstract_matrix, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };

  template <typename V> struct temporary_vector {
    typedef typename _temporary_vector<typename is_a_reference<V>::reference,
				       typename linalg_traits<V>::storage_type,
				       typename linalg_traits<V>::linalg_type,
				       V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary matrix type                                    */
  /*   M if M is a valid matrix type,                                     */
  /*   row_matrix<wsvector> if M is a reference on a sparse matrix,       */
  /*   dense_matrix if M is a reference on a dense matrix.                */
  /* ******************************************************************** */

  
  template <typename R, typename S, typename L, typename V> struct _temporary_matrix {
    typedef abstract_null_type matrix_type;
  };
  template <typename V, typename L>
  struct _temporary_matrix<linalg_true, abstract_sparse, L, V>
  { typedef row_matrix<wsvector<typename linalg_traits<V>::value_type> > matrix_type; };
  template <typename V, typename L>
  struct _temporary_matrix<linalg_true, abstract_skyline, L, V>
  { typedef row_matrix<slvector<typename linalg_traits<V>::value_type> > matrix_type; };
  template <typename V, typename L>
  struct _temporary_matrix<linalg_true, abstract_dense, L, V>
  { typedef dense_matrix<typename linalg_traits<V>::value_type> matrix_type; };
  template <typename S, typename V>
  struct _temporary_matrix<linalg_false, S, abstract_matrix, V>
  { typedef V matrix_type; };

  template <typename V> struct temporary_matrix {
    typedef typename _temporary_matrix<typename is_a_reference<V>::reference,
				       typename linalg_traits<V>::storage_type,
				       typename linalg_traits<V>::linalg_type,
				       V>::matrix_type matrix_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary dense vector type                              */
  /*   V if V is a valid dense vector type,                               */
  /*   std::vector if V is a reference or another type of vector          */
  /* ******************************************************************** */

  template <typename R, typename S, typename V> struct _temporary_dense_vector {
    typedef abstract_null_type vector_type;
  };
  template <typename S, typename V> struct _temporary_dense_vector<linalg_true, S, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_dense_vector<linalg_false, abstract_sparse, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_dense_vector<linalg_false, abstract_skyline, V>
  { typedef std::vector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_dense_vector<linalg_false, abstract_dense, V>
  { typedef V vector_type; };

  template <typename V> struct temporary_dense_vector {
    typedef typename _temporary_dense_vector<typename
    is_a_reference<V>::reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary sparse vector type                             */
  /*   V if V is a valid sparse vector type,                              */
  /*   wsvector if V is a reference or another type of vector             */
  /* ******************************************************************** */

  template <typename R, typename S, typename V> struct _temporary_sparse_vector {
    typedef abstract_null_type vector_type;
  };
  template <typename S, typename V>
  struct _temporary_sparse_vector<linalg_true, S, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_sparse_vector<linalg_false, abstract_sparse, V>
  { typedef V vector_type; };
  template <typename V>
  struct _temporary_sparse_vector<linalg_false, abstract_dense, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_sparse_vector<linalg_false, abstract_skyline, V>
  { typedef wsvector<typename linalg_traits<V>::value_type> vector_type; };

  template <typename V> struct temporary_sparse_vector {
    typedef typename _temporary_sparse_vector<typename
    is_a_reference<V>::reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary sky-line vector type                           */
  /*   V if V is a valid sky-line vector type,                            */
  /*   slvector if V is a reference or another type of vector             */
  /* ******************************************************************** */

  template <typename R, typename S, typename V> struct _temporary_skyline_vector {
    typedef abstract_null_type vector_type;
  };
  template <typename S, typename V>
  struct _temporary_skyline_vector<linalg_true, S, V>
  { typedef slvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_skyline_vector<linalg_false, abstract_skyline, V>
  { typedef V vector_type; };
  template <typename V>
  struct _temporary_skyline_vector<linalg_false, abstract_dense, V>
  { typedef slvector<typename linalg_traits<V>::value_type> vector_type; };
  template <typename V>
  struct _temporary_skyline_vector<linalg_false, abstract_sparse, V>
  { typedef slvector<typename linalg_traits<V>::value_type> vector_type; };

  template <typename V> struct temporary_skylines_vector {
    typedef typename _temporary_skyline_vector<typename
    is_a_reference<V>::reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };


  /* ********************************************************************* */
  /*		Standard access and clear objects             		   */
  /* ********************************************************************* */

  template <typename IT, typename CIT> struct dense_access {
    
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::reference reference;
    typedef typename std::iterator_traits<CIT>::pointer const_pointer;

    typedef typename select_return<const value_type &, reference,
			 const_pointer>::return_type const_reference;
  
    reference operator()(const void *, const IT &_begin,
			 const IT &, size_type i)
    { return _begin[i]; }
    const_reference operator()(const void *, const CIT &_begin,
			       const CIT &, size_type i)
    { return _begin[i]; }
  };

  template <typename IT> struct dense_access<IT, IT> {
    typedef typename std::iterator_traits<IT>::value_type value_type;
  
    value_type operator()(const void *, const IT &_begin,
			 const IT &, size_type i)
    { return _begin[i]; }
  };


  template <typename IT> struct dense_clear {
    typedef typename std::iterator_traits<IT>::value_type value_type;
    inline void operator()(const void *,const IT &_begin,const IT &_end);
  };
  
  template <typename IT> inline void dense_clear<IT>::operator()(const void *,
			  const IT &_begin, const IT &_end)
  { std::fill(_begin, _end, value_type(0)); }

  /* ********************************************************************* */
  /* Set to begin end set to end for iterators on non-const sparse vectors.*/
  /* ********************************************************************* */

  template <typename IT, typename VECT> inline
  void set_to_begin(IT &it, const void *o, VECT *)
  { it = vect_begin(*(const_cast<VECT *>((const VECT *)(o)))); }

  template <typename IT, typename VECT> inline
  void set_to_begin(IT &it, const void *o, const VECT *) 
  { it = vect_const_begin(*((const VECT *)(o))); }

  template <typename IT, typename VECT> inline
  void set_to_end(IT &it, const void *o, VECT *)
  { it = vect_end(*(const_cast<VECT *>((const VECT *)(o)))); }
  
  template <typename IT, typename VECT> inline
  void set_to_end(IT &it, const void *o, const VECT *)
  { it = vect_const_end(*((const VECT *)(o))); }

}

#endif //  __GMM_DEF_H
