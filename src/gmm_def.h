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
  struct linalg_const {};
  struct linalg_modifiable {};

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
    typedef typename linalg_traits<V>::reference access_type;
    typedef typename linalg_traits<V>::iterator iterator;
  };
  template <class P, class V> struct vect_ref_type<const P *, V> {
    typedef typename linalg_traits<V>::value_type access_type;
    typedef typename linalg_traits<V>::const_iterator iterator;
  };
  
  template <class PT> struct const_pointer;
  template <class P> struct const_pointer<P *>
  { typedef const P* pointer; };
  template <class P> struct const_pointer<const P *>
  { typedef const P* pointer; };

  template <class PT> struct modifiable_pointer;
  template <class P> struct modifiable_pointer<P *>
  { typedef P* pointer; };
  template <class P> struct modifiable_pointer<const P *>
  { typedef P* pointer; };

  inline bool _is_sparse(abstract_sparse)  { return true;  }
  inline bool _is_sparse(abstract_plain)   { return false; }
  inline bool _is_sparse(abstract_indirect)  { return false; }

  template <class L> inline bool is_sparse(const L &) 
  { return _is_sparse(typename linalg_traits<L>::storage_type()); }


  /* ******************************************************************** */
  /*  types to deal with const object representing a modifiable reference */
  /* ******************************************************************** */
  
  template <class PT, class R> struct _mref_type;

  template <class L, class R> struct _mref_type<L *, R>
  { typedef L & return_type; };
  template <class L, class R> struct _mref_type<const L *, R>
  { typedef const L & return_type; };
  template <class L> struct _mref_type<L *, linalg_const>
  { typedef const L & return_type; };
  template <class L> struct _mref_type<const L *, linalg_const>
  { typedef const L & return_type; };
  template <class L> struct _mref_type<const L *, linalg_modifiable>
  { typedef L & return_type; };
  template <class L> struct _mref_type<L *, linalg_modifiable>
  { typedef L & return_type; };

  template <class PT> struct mref_type {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _mref_type<PT, 
      typename linalg_traits<L>::is_reference>::return_type return_type;
  };

  template <class L> typename mref_type<const L *>::return_type 
  linalg_cast(const L &l)
  { return const_cast<typename mref_type<const L *>::return_type>(l); }
  template <class L> typename mref_type<L *>::return_type linalg_cast(L &l)
  { return const_cast<typename mref_type<L *>::return_type>(l); }



  template <class L, class R> struct _cref_type;
  template <class L> struct _cref_type<L, linalg_modifiable>
  { typedef L & return_type; };

  template <class L> struct cref_type {
    typedef typename _cref_type<L, 
      typename linalg_traits<L>::is_reference>::return_type return_type;
  };

  template <class L> typename cref_type<L>::return_type 
  linalg_const_cast(const L &l)
  { return const_cast<typename cref_type<L>::return_type>(l); }



  template <class C1, class C2, class REF> struct _select_return;
  template <class C1, class C2, class L>
  struct _select_return<C1, C2, const L &> { typedef C1 return_type; };
  template <class C1, class C2, class L>
  struct _select_return<C1, C2, L &> { typedef C2 return_type; };
  template <class C1, class C2, class PT> struct select_return {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _select_return<C1, C2, 
      typename mref_type<PT>::return_type>::return_type return_type;
  };

  template<class R> struct _is_a_reference
  { typedef linalg_true reference; };
  template<> struct _is_a_reference<linalg_false>
  { typedef linalg_false reference; };

  template<class L> struct is_a_reference {
    typedef typename _is_a_reference<typename linalg_traits<L>::is_reference>
      ::reference reference;
  };

  template <class PT> struct which_reference;
  template <class PT> struct which_reference<PT *>
  { typedef linalg_modifiable is_reference; };
  template <class PT> struct which_reference<const PT *>
  { typedef linalg_const is_reference; };


  template <class C1, class C2, class R> struct _select_orientation;
  template <class C1, class C2> struct _select_orientation<C1, C2, row_major>
  { typedef C1 return_type; };
   template <class C1, class C2> struct _select_orientation<C1, C2, col_major>
   { typedef C2 return_type; };
  template <class C1, class C2, class L> struct select_orientation {
    typedef typename _select_orientation<C1, C2,
      typename principal_orientation_type<typename
      linalg_traits<L>::sub_orientation>::potype>::return_type return_type;
  };
  
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

  /* ******************************************************************** */
  /*		Basic vectors used                         		  */
  /* ******************************************************************** */
  
  template<class T> struct plain_vector_type 
  { typedef std::vector<T> vector_type; };

  template <class T> class wsvector;
  template<class T> struct sparse_vector_type 
  { typedef wsvector<T> vector_type; };

  /* ******************************************************************** */
  /*   Selects a temporary vector type                                    */
  /*   V if V is a valid vector type,                                     */
  /*   wsvector if V is a reference on a sparse vector,                   */
  /*   std::vector if V is a reference on a plain vector.                 */
  /* ******************************************************************** */

  
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
    typedef typename _temporary_vector<typename is_a_reference<V>::reference,
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
    is_a_reference<V>::reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ******************************************************************** */
  /*   Selects a temporary sparse vector type                             */
  /*   V if V is a valid sparse vector type,                              */
  /*   wsvector if V is a reference or a plain vector                     */
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
    is_a_reference<V>::reference,
    typename linalg_traits<V>::storage_type, V>::vector_type vector_type;
  };

  /* ********************************************************************* */
  /*		Standard access and clear objects             		   */
  /* ********************************************************************* */

  template <class V> struct plain_access {
    typedef typename linalg_traits<V>::value_type value_type;
    typedef value_type &reference;
    typedef typename linalg_traits<V>::iterator iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    reference operator()(const void *, const iterator &_begin,
			 const iterator &, size_type i)
    { return _begin[i]; }
    value_type operator()(const void *, const const_iterator &_begin,
			 const const_iterator &, size_type i)
    { return _begin[i]; }
  };

  template <class V> struct plain_clear {
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::iterator iterator;
    
    void operator()(const void *, const iterator &_begin, const iterator &_end)
    { std::fill(_begin, _end, value_type(0)); }
  };

  /* ******************************************************************** */
  /*		sub indexes                               		  */
  /* ******************************************************************** */

  struct sub_index {

    typedef std::vector<size_t> base_type;
    typedef base_type::const_iterator const_iterator;

    std::vector<size_t> ind;
    std::vector<size_t> rind;

    size_type size(void) const { return ind.size(); }
    size_type index(size_type i) const { return ind[i]; }
    size_type rindex(size_type i) const { return rind[i]; }
   
    const_iterator  begin(void) const { return  ind.begin(); }
    const_iterator    end(void) const { return  ind.end();   }
    const_iterator rbegin(void) const { return rind.begin(); }
    const_iterator   rend(void) const { return rind.end();   }

    sub_index(void) {}
    template <class IT> void init(IT, IT, size_type, abstract_sparse);
    template <class IT> void init(IT, IT, size_type, abstract_plain);
    template <class IT, class L> sub_index(IT it, IT ite,
					       size_type n, const L &)
    { init(it, ite, n, typename linalg_traits<L>::storage_type()); }
    template <class CONT, class L> sub_index(const CONT &c,
					     size_type n, const L &)
    { init(c.begin(), c.end(), n, typename linalg_traits<L>::storage_type()); }
    template <class CONT> sub_index(const CONT &c, size_type n)
    { init(c.begin(), c.end(), n, abstract_sparse()); }

  };

  template <class IT>
  void sub_index::init(IT it, IT ite, size_type n, abstract_sparse) {
      rind.resize(n); std::fill(rind.begin(), rind.end(), size_t(-1));
      IT itaux = it;
      for (size_type i = 0; it != ite; ++it, ++i) rind[*it] = i;
      //      cout << "i = " << i << endl;
      ind.resize(ite - itaux); std::copy(itaux, ite, ind.begin());
  }

  template <class IT>
  void sub_index::init(IT it, IT ite, size_type, abstract_plain)
  { ind.resize(ite - it); std::copy(it, ite, ind.begin()); }


  struct sub_interval {
    size_type min, max; 

    size_type size(void) const { return max - min + 1; }
    size_type index(size_type i) const { return min + i; }
    size_type rindex(size_type i) const
    { if (i >= min && i <= max) return i - min; return size_type(-1); }
    sub_interval(size_type mi, size_type l) : min(mi), max(mi+l-1) {}
    sub_interval() {}
  };

  struct sub_slice {
    size_type min, max, N; 

    size_type size(void) const { return (max - min) / N + 1; }
    size_type index(size_type i) const { return min + N * i; }
    size_type rindex(size_type i) const { 
      if (i >= min && i <= max)
	{ size_type j = (i - min); if (j % N == 0) return j / N; }
      return size_type(-1);
    }
    sub_slice(size_type mi, size_type l, size_type n)
      : min(mi), max(mi+(l-1)*n), N(n) {}
    sub_slice(void) {}
  };

}

#endif //  __GMM_DEF_H
