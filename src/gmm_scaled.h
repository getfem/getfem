/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_transposed.h : generic transposed matrices.              */
/*     									   */
/* Date : November 10, 2002.                                               */
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

//   . faire scale et scaled sur les matrices aussi
//

#ifndef __GMM_SCALED_H
#define __GMM_SCALED_H

namespace gmm {

  /* ********************************************************************* */
  /*		Scaled references on vectors            		   */
  /* ********************************************************************* */

  template <class IT> struct scaled_const_iterator {
    typedef typename std::iterator_traits<IT>::value_type      value_type;
    typedef typename std::iterator_traits<IT>::pointer         pointer;
    typedef typename std::iterator_traits<IT>::reference       reference;
    typedef typename std::iterator_traits<IT>::difference_type difference_type;
    typedef typename std::iterator_traits<IT>::iterator_category
    iterator_category;

    IT it;
    value_type r;
    
    scaled_const_iterator(void) {}
    scaled_const_iterator(const IT &i, value_type x) : it(i), r(x) {}
    
    inline size_type index(void) const { return it.index(); }
    scaled_const_iterator operator ++(int)
    { scaled_const_iterator tmp = *this; ++it; return tmp; }
    scaled_const_iterator operator --(int) 
    { scaled_const_iterator tmp = *this; --it; return tmp; }
    scaled_const_iterator &operator ++() { ++it; return *this; }
    scaled_const_iterator &operator --() { --it; return *this; }
    scaled_const_iterator &operator +=(difference_type i)
      { it += i; return *this; }
    scaled_const_iterator &operator -=(difference_type i)
      { it -= i; return *this; }
    scaled_const_iterator operator +(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb += i); }
    scaled_const_iterator operator -(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb -= i); }
    difference_type operator -(const scaled_const_iterator &i) const
      { return difference_type(it - i.it); }
    
    value_type operator  *() const { return (*it) * r; }
    value_type operator [](size_type ii) const { return it[ii] * r; }
    
    bool operator ==(const scaled_const_iterator &i) const
      { return (i.it == it); }
    bool operator !=(const scaled_const_iterator &i) const
      { return (i.it != it); }
    bool operator < (const scaled_const_iterator &i) const
      { return (it < i.it); }
  };

  template <class V> struct scaled_vector_const_ref {
    typedef scaled_vector_const_ref<V> this_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::const_iterator iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<V>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    size_type _size;
    value_type r;

    scaled_vector_const_ref(const V &v, value_type rr)
      : _begin(vect_const_begin(v)), _end(vect_const_end(v)), origin(linalg_origin(v)),
	_size(vect_size(v)), r(rr) {}

    reference operator[](size_type i) const
    { return r * access_type()(origin, _begin, _end, i); }
  };

  template <class V> struct scaled_vector_const_access {
    typedef scaled_vector_const_ref<V> this_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<this_type>::iterator iterator;
    
    value_type operator()(const void *o, const iterator &_begin,
			  const iterator &_end, size_type i) {
      return _begin.r 
	* typename linalg_traits<V>::access_type()(o, _begin.it, _end.it, i);
    }
  };

  template <class V> struct linalg_traits<scaled_vector_const_ref<V> > {
    typedef scaled_vector_const_ref<V> this_type;
    typedef linalg_const is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef value_type reference;
    typedef scaled_const_iterator<typename linalg_traits<V>::const_iterator>
            iterator;
    typedef iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef scaled_vector_const_access<V> access_type;
    typedef abstract_null_type clear_type;
    static size_type size(const this_type &v) { return v._size; }
    static iterator begin(this_type &v)
    { return iterator(v._begin, v.r); }
    static const_iterator begin(const this_type &v)
    { return const_iterator(v._begin, v.r); }
    static iterator end(this_type &v)
    { return iterator(v._end, v.r); }
    static const_iterator end(const this_type &v)
    { return const_iterator(v._end, v.r); }
    static const void* origin(const this_type &v) { return v.origin; }
  };

  // for GCC 2.95
  template <class V> struct linalg_traits<const scaled_vector_const_ref<V> > 
    : public linalg_traits<scaled_vector_const_ref<V> > {};


  /* ********************************************************************* */
  /*		Scaled references on matrices            		   */
  /* ********************************************************************* */

  template <class M> struct scaled_row_const_iterator {
    typedef scaled_row_const_iterator<M> iterator;
    typedef typename linalg_traits<M>::const_row_iterator ITER;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;

    ITER it;
    value_type r;

    iterator operator ++(int) { iterator tmp = *this; it++; return tmp; }
    iterator operator --(int) { iterator tmp = *this; it--; return tmp; }
    iterator &operator ++()   { it++; return *this; }
    iterator &operator --()   { it--; return *this; }
    iterator &operator +=(difference_type i) { it += i; return *this; }
    iterator &operator -=(difference_type i) { it -= i; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return it - i.it; }

    ITER operator *() const { return it; }
    ITER operator [](int i) { return it + i; }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    scaled_row_const_iterator(void) {}
    scaled_row_const_iterator(const ITER &i, value_type rr)
      : it(i), r(rr) { }

  };

  template <class M> struct  scaled_row_matrix_const_ref {
    
    typedef scaled_row_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<M>::const_row_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    value_type r;

    scaled_row_matrix_const_ref(const M &m, value_type rr)
      : _begin(mat_row_begin(m)), 
      _end(mat_row_end(m)), origin(linalg_origin(m)), r(rr) {}

    value_type operator()(size_type i, size_type j) const
    { return r * access_type()(_begin+i, j); }
  };

  template <class M> struct scaled_row_matrix_access {
    typedef scaled_row_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<this_type>::row_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::access_type access_type;
    
    value_type operator()(const iterator &itrow, size_type i)
    { return itrow.r * access_type(itrow.it, i); }
  };

  template <class M> struct linalg_traits<scaled_row_matrix_const_ref<M> > {
    typedef scaled_row_matrix_const_ref<M> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef value_type reference;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef typename linalg_traits<M>::const_sub_row_type vector_type;
    typedef scaled_vector_const_ref<vector_type> sub_row_type;
    typedef scaled_vector_const_ref<vector_type> const_sub_row_type;
    typedef scaled_row_const_iterator<M> row_iterator;
    typedef scaled_row_const_iterator<M> const_row_iterator;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_col_iterator;
    typedef abstract_null_type col_iterator;
    typedef row_major sub_orientation;
    typedef scaled_row_matrix_access<M> access_type;
    static size_type nrows(const this_type &m) { return m._end - m._begin; }
    static size_type ncols(const this_type &m)
    { return (nrows(m) == 0) ? 0 : vect_size(mat_row(m, 0)); }
    static const_sub_row_type row(const const_row_iterator &it)
    { return scaled(linalg_traits<M>::row(it.it), it.r); }
    static const_row_iterator row_begin(const this_type &m)
    { return const_row_iterator(m._begin, m.r); }
    static const_row_iterator row_end(const this_type &m)
    { return const_row_iterator(m._end, m.r); }
    static const void* origin(const this_type &m) { return m.origin; }
  };

  // for GCC 2.95
  template <class M>
  struct linalg_traits<const scaled_row_matrix_const_ref<M> > 
    : public linalg_traits<scaled_row_matrix_const_ref<M> > {};


  template <class M> struct scaled_col_const_iterator {
    typedef scaled_col_const_iterator<M> iterator;
    typedef typename linalg_traits<M>::const_col_iterator ITER;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;

    ITER it;
    value_type r;

    iterator operator ++(int) { iterator tmp = *this; it++; return tmp; }
    iterator operator --(int) { iterator tmp = *this; it--; return tmp; }
    iterator &operator ++()   { it++; return *this; }
    iterator &operator --()   { it--; return *this; }
    iterator &operator +=(difference_type i) { it += i; return *this; }
    iterator &operator -=(difference_type i) { it -= i; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return it - i.it; }

    ITER operator *() const { return it; }
    ITER operator [](int i) { return it + i; }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    scaled_col_const_iterator(void) {}
    scaled_col_const_iterator(const ITER &i, value_type rr)
      : it(i), r(rr) { }

  };

  template <class M> struct  scaled_col_matrix_const_ref {
    
    typedef scaled_col_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<M>::const_col_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    value_type r;

    scaled_col_matrix_const_ref(const M &m, value_type rr)
      : _begin(mat_col_begin(m)), 
      _end(mat_col_end(m)), origin(linalg_origin(m)), r(rr) {}

    value_type operator()(size_type i, size_type j) const
    { return r * access_type()(_begin+i, j); }
  };

  template <class M> struct scaled_col_matrix_access {
    typedef scaled_col_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<this_type>::col_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::access_type access_type;
    
    value_type operator()(const iterator &itcol, size_type i)
    { return itcol.r * access_type(itcol.it, i); }
  };

  template <class M> struct linalg_traits<scaled_col_matrix_const_ref<M> > {
    typedef scaled_col_matrix_const_ref<M> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef value_type reference;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef typename linalg_traits<M>::const_sub_col_type vector_type;
    typedef scaled_vector_const_ref<vector_type> sub_col_type;
    typedef scaled_vector_const_ref<vector_type> const_sub_col_type;
    typedef scaled_col_const_iterator<M> col_iterator;
    typedef scaled_col_const_iterator<M> const_col_iterator;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_row_iterator;
    typedef abstract_null_type row_iterator;
    typedef col_major sub_orientation;
    typedef scaled_col_matrix_access<M> access_type;
    static size_type ncols(const this_type &m) { return m._end - m._begin; }
    static size_type nrows(const this_type &m)
    { return (ncols(m) == 0) ? 0 : vect_size(mat_col(m, 0)); }
    static const_sub_col_type col(const const_col_iterator &it)
    { return scaled(linalg_traits<M>::col(it.it), it.r); }
    static const_col_iterator col_begin(const this_type &m)
    { return const_col_iterator(m._begin, m.r); }
    static const_col_iterator col_end(const this_type &m)
    { return const_col_iterator(m._end, m.r); }
    static const void* origin(const this_type &m) { return m.origin; }
  };

  // for GCC 2.95
  template <class M>
  struct linalg_traits<const scaled_col_matrix_const_ref<M> > 
    : public linalg_traits<scaled_col_matrix_const_ref<M> > {};

  template <class L, class R> struct __scaled_return;
  template <class L> struct __scaled_return<L, row_major> 
  { typedef scaled_row_matrix_const_ref<L> return_type; };
  template <class L> struct __scaled_return<L, col_major> 
  { typedef scaled_col_matrix_const_ref<L> return_type; };
  

  template <class L, class LT> struct _scaled_return;
  template <class L> struct _scaled_return<L, abstract_vector> 
  { typedef scaled_vector_const_ref<L> return_type; };
  template <class L> struct _scaled_return<L, abstract_matrix> {
    typedef typename __scaled_return<L, 
      typename principal_orientation_type<typename
      linalg_traits<L>::sub_orientation>::potype>::return_type return_type;
  };

  template <class L> struct scaled_return {
    typedef typename _scaled_return<L, typename
      linalg_traits<L>::linalg_type>::return_type return_type;
  };

  template <class L> inline
  typename scaled_return<L>::return_type
  scaled(const L &v, typename linalg_traits<L>::value_type x)
  { return scaled(v, x, typename linalg_traits<L>::linalg_type()); }

  template <class V> inline
  typename scaled_return<V>::return_type
  scaled(const V &v, typename linalg_traits<V>::value_type x, abstract_vector)
  { return scaled_vector_const_ref<V>(v, x); }

  template <class M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x,abstract_matrix) {
    return scaled(m, x,  typename principal_orientation_type<typename
		  linalg_traits<M>::sub_orientation>::potype());
  }

  template <class M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x, row_major) {
    return scaled_row_matrix_const_ref<M>(m, x);
  }

  template <class M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x, col_major) {
    return scaled_col_matrix_const_ref<M>(m, x);
  }

  
  /* ******************************************************************** */
  /*	matrix or vector scale                                	          */
  /* ******************************************************************** */

  template <class L> inline
  void scale(L& l, typename linalg_traits<L>::value_type a)
  { scale(l, a, typename linalg_traits<L>::linalg_type()); }

  template <class L> inline
  void scale(const L& l, typename linalg_traits<L>::value_type a)
  { scale(linalg_const_cast(l), a); }

  template <class L> inline
  void scale(L& l, typename linalg_traits<L>::value_type a, abstract_vector) {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for ( ; it != ite; ++it) *it *= a;
  }

  template <class L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, abstract_matrix) {
    scale(l, a, typename principal_orientation_type<typename
	  linalg_traits<L>::sub_orientation>::potype());
  }

  template <class L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, row_major) {
    typename linalg_traits<L>::row_iterator it = mat_row_begin(l),
      ite = mat_row_end(l);
    for ( ; it != ite; ++it) scale(linalg_traits<L>::row(*it), a);
  }

  template <class L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, col_major) {
    typename linalg_traits<L>::col_iterator it = mat_col_begin(l),
      ite = mat_col_end(l);
    for ( ; it != ite; ++it) scale(linalg_traits<L>::col(*it), a);
  }

}

#endif //  __GMM_SCALED_H
