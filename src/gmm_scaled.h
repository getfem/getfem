/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_scaled.h : generic scaled vectors and matrices.          */
/*     									   */
/* Date : November 10, 2002.                                               */
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

#ifndef __GMM_SCALED_H
#define __GMM_SCALED_H

#include <gmm_def.h>

namespace gmm {

  /* ********************************************************************* */
  /*		Scaled references on vectors            		   */
  /* ********************************************************************* */

  template <typename IT> struct scaled_const_iterator {
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
    inline scaled_const_iterator operator ++(int)
    { scaled_const_iterator tmp = *this; ++it; return tmp; }
    inline scaled_const_iterator operator --(int) 
    { scaled_const_iterator tmp = *this; --it; return tmp; }
    inline scaled_const_iterator &operator ++() { ++it; return *this; }
    inline scaled_const_iterator &operator --() { --it; return *this; }
    inline scaled_const_iterator &operator +=(difference_type i)
      { it += i; return *this; }
    inline scaled_const_iterator &operator -=(difference_type i)
      { it -= i; return *this; }
    inline scaled_const_iterator operator +(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb += i); }
    inline scaled_const_iterator operator -(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb -= i); }
    inline difference_type operator -(const scaled_const_iterator &i) const
      { return difference_type(it - i.it); }
    
    inline value_type operator  *() const { return (*it) * r; }
    inline value_type operator [](size_type ii) const { return it[ii] * r; }
    
    inline bool operator ==(const scaled_const_iterator &i) const
      { return (i.it == it); }
    inline bool operator !=(const scaled_const_iterator &i) const
      { return (i.it != it); }
    inline bool operator < (const scaled_const_iterator &i) const
      { return (it < i.it); }
  };

  template <typename V> struct scaled_vector_const_ref {
    typedef scaled_vector_const_ref<V> this_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::const_iterator iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::origin_type origin_type;

    iterator _begin, _end;
    const origin_type *origin;
    size_type _size;
    value_type r;

    scaled_vector_const_ref(const V &v, value_type rr)
      : _begin(vect_const_begin(v)), _end(vect_const_end(v)),
	origin(linalg_origin(v)), _size(vect_size(v)), r(rr) {}

    reference operator[](size_type i) const
    { return r * linalg_traits<V>::access(origin, _begin, _end, i); }
  };

  template <typename V> struct linalg_traits<scaled_vector_const_ref<V> > {
    typedef scaled_vector_const_ref<V> this_type;
    typedef linalg_const is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::origin_type origin_type;
    typedef value_type reference;
    typedef abstract_null_type iterator;
    typedef scaled_const_iterator<typename linalg_traits<V>::const_iterator>
            const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    static size_type size(const this_type &v) { return v._size; }
    static iterator begin(this_type &v)
    { return iterator(v._begin, v.r); }
    static const_iterator begin(const this_type &v)
    { return const_iterator(v._begin, v.r); }
    static iterator end(this_type &v)
    { return iterator(v._end, v.r); }
    static const_iterator end(const this_type &v)
    { return const_iterator(v._end, v.r); }
    static const origin_type* origin(const this_type &v) { return v.origin; }
    static value_type access(const origin_type *o, const const_iterator &it,
			     const const_iterator &ite, size_type i)
    { return it.r * (linalg_traits<V>::access(o, it.it, ite.it, i)); }

  };

   template<typename V> std::ostream &operator <<
  (std::ostream &o, const scaled_vector_const_ref<V>& m)
  { gmm::write(o,m); return o; }

#ifdef USING_BROKEN_GCC295
  template <typename V> struct linalg_traits<const scaled_vector_const_ref<V> > 
    : public linalg_traits<scaled_vector_const_ref<V> > {};
#endif


  /* ********************************************************************* */
  /*		Scaled references on matrices            		   */
  /* ********************************************************************* */

  template <typename M> struct scaled_row_const_iterator {
    typedef scaled_row_const_iterator<M> iterator;
    typedef typename linalg_traits<M>::const_row_iterator ITER;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;

    ITER it;
    value_type r;

    inline iterator operator ++(int) { iterator tmp=*this; it++; return tmp; }
    inline iterator operator --(int) { iterator tmp=*this; it--; return tmp; }
    inline iterator &operator ++()   { it++; return *this; }
    inline iterator &operator --()   { it--; return *this; }
    iterator &operator +=(difference_type i) { it += i; return *this; }
    iterator &operator -=(difference_type i) { it -= i; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return it - i.it; }

    inline ITER operator *() const { return it; }
    inline ITER operator [](int i) { return it + i; }

    inline bool operator ==(const iterator &i) const { return (it == i.it); }
    inline bool operator !=(const iterator &i) const { return !(i == *this); }
    inline bool operator < (const iterator &i) const { return (it < i.it); }

    scaled_row_const_iterator(void) {}
    scaled_row_const_iterator(const ITER &i, value_type rr)
      : it(i), r(rr) { }

  };

  template <typename M> struct  scaled_row_matrix_const_ref {
    
    typedef scaled_row_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<M>::const_row_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<this_type>::origin_type origin_type;

    iterator _begin, _end;
    const origin_type *origin;
    value_type r;
    size_type nr, nc;

    scaled_row_matrix_const_ref(const M &m, value_type rr)
      : _begin(mat_row_begin(m)), _end(mat_row_end(m)),
	origin(linalg_origin(m)), r(rr), nr(mat_ncols(m)), nc(mat_nrows(m)) {}

    value_type operator()(size_type i, size_type j) const
    { return r * linalg_traits<M>::access(_begin+i, j); }
  };

  template <typename M> struct linalg_traits<scaled_row_matrix_const_ref<M> > {
    typedef scaled_row_matrix_const_ref<M> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::origin_type origin_type;
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
    static size_type nrows(const this_type &m)
    { return m.nr; }
    static size_type ncols(const this_type &m)
    { return m.nc; }
    static const_sub_row_type row(const const_row_iterator &it)
    { return scaled(linalg_traits<M>::row(it.it), it.r); }
    static const_row_iterator row_begin(const this_type &m)
    { return const_row_iterator(m._begin, m.r); }
    static const_row_iterator row_end(const this_type &m)
    { return const_row_iterator(m._end, m.r); }
    static const origin_type* origin(const this_type &m) { return m.origin; }
    static value_type access(const const_row_iterator &it, size_type i)
    { return it.r * (linalg_traits<M>::access(it.it, i)); }
  };

  template<typename M> std::ostream &operator <<
  (std::ostream &o, const scaled_row_matrix_const_ref<M>& m)
  { gmm::write(o,m); return o; }

#ifdef USING_BROKEN_GCC295
  template <typename M>
  struct linalg_traits<const scaled_row_matrix_const_ref<M> > 
    : public linalg_traits<scaled_row_matrix_const_ref<M> > {};
#endif


  template <typename M> struct scaled_col_const_iterator {
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

  template <typename M> struct  scaled_col_matrix_const_ref {
    
    typedef scaled_col_matrix_const_ref<M> this_type;
    typedef typename linalg_traits<M>::const_col_iterator iterator;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<this_type>::origin_type origin_type;

    iterator _begin, _end;
    const origin_type *origin;
    value_type r;
    size_type nr, nc;

    scaled_col_matrix_const_ref(const M &m, value_type rr)
      : _begin(mat_col_begin(m)), _end(mat_col_end(m)),
	origin(linalg_origin(m)), r(rr), nr(mat_nrows(m)), nc(mat_ncols(m)) {}

    value_type operator()(size_type i, size_type j) const
    { return r * linalg_traits<M>::access(_begin+j, i); }
  };

  template <typename M> struct linalg_traits<scaled_col_matrix_const_ref<M> > {
    typedef scaled_col_matrix_const_ref<M> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::origin_type origin_type;
    typedef value_type reference;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef typename linalg_traits<M>::const_sub_col_type vector_type;
    typedef abstract_null_type sub_col_type;
    typedef scaled_vector_const_ref<vector_type> const_sub_col_type;
    typedef abstract_null_type  col_iterator;
    typedef scaled_col_const_iterator<M> const_col_iterator;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_row_iterator;
    typedef abstract_null_type row_iterator;
    typedef col_major sub_orientation;
    static size_type ncols(const this_type &m)
    { return m.nc; }
    static size_type nrows(const this_type &m)
    { return m.nr; }
    static const_sub_col_type col(const const_col_iterator &it)
    { return scaled(linalg_traits<M>::col(it.it), it.r); }
    static const_col_iterator col_begin(const this_type &m)
    { return const_col_iterator(m._begin, m.r); }
    static const_col_iterator col_end(const this_type &m)
    { return const_col_iterator(m._end, m.r); }
    static const origin_type* origin(const this_type &m) { return m.origin; }
    static value_type access(const const_col_iterator &it, size_type i)
    { return it.r * (linalg_traits<M>::access(it.it, i)); }
  };

  template<typename M> std::ostream &operator <<
  (std::ostream &o, const scaled_col_matrix_const_ref<M>& m)
  { gmm::write(o,m); return o; }


#ifdef USING_BROKEN_GCC295
  template <typename M>
  struct linalg_traits<const scaled_col_matrix_const_ref<M> > 
    : public linalg_traits<scaled_col_matrix_const_ref<M> > {};
#endif

  template <typename L, typename R> struct __scaled_return {
    typedef abstract_null_type return_type;
  };
  template <typename L> struct __scaled_return<L, row_major> 
  { typedef scaled_row_matrix_const_ref<L> return_type; };
  template <typename L> struct __scaled_return<L, col_major> 
  { typedef scaled_col_matrix_const_ref<L> return_type; };
  

  template <typename L, typename LT> struct _scaled_return {
    typedef abstract_null_type return_type;
  };
  template <typename L> struct _scaled_return<L, abstract_vector> 
  { typedef scaled_vector_const_ref<L> return_type; };
  template <typename L> struct _scaled_return<L, abstract_matrix> {
    typedef typename __scaled_return<L, 
      typename principal_orientation_type<typename
      linalg_traits<L>::sub_orientation>::potype>::return_type return_type;
  };

  template <typename L> struct scaled_return {
    typedef typename _scaled_return<L, typename
      linalg_traits<L>::linalg_type>::return_type return_type;
  };

  template <typename L> inline
  typename scaled_return<L>::return_type
  scaled(const L &v, typename linalg_traits<L>::value_type x)
  { return scaled(v, x, typename linalg_traits<L>::linalg_type()); }

  template <typename V> inline
  typename scaled_return<V>::return_type
  scaled(const V &v, typename linalg_traits<V>::value_type x, abstract_vector)
  { return scaled_vector_const_ref<V>(v, x); }

  template <typename M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x,abstract_matrix) {
    return scaled(m, x,  typename principal_orientation_type<typename
		  linalg_traits<M>::sub_orientation>::potype());
  }

  template <typename M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x, row_major) {
    return scaled_row_matrix_const_ref<M>(m, x);
  }

  template <typename M> inline
  typename scaled_return<M>::return_type
  scaled(const M &m, typename linalg_traits<M>::value_type x, col_major) {
    return scaled_col_matrix_const_ref<M>(m, x);
  }

  
  /* ******************************************************************** */
  /*	matrix or vector scale                                	          */
  /* ******************************************************************** */

  template <typename L> inline
  void scale(L& l, typename linalg_traits<L>::value_type a)
  { scale(l, a, typename linalg_traits<L>::linalg_type()); }

  template <typename L> inline
  void scale(const L& l, typename linalg_traits<L>::value_type a)
  { scale(linalg_const_cast(l), a); }

  template <typename L> inline
  void scale(L& l, typename linalg_traits<L>::value_type a, abstract_vector) {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for ( ; it != ite; ++it) *it *= a;
  }

  template <typename L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, abstract_matrix) {
    scale(l, a, typename principal_orientation_type<typename
	  linalg_traits<L>::sub_orientation>::potype());
  }

  template <typename L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, row_major) {
    typename linalg_traits<L>::row_iterator it = mat_row_begin(l),
      ite = mat_row_end(l);
    for ( ; it != ite; ++it) scale(linalg_traits<L>::row(it), a);
  }

  template <typename L> 
  void scale(L& l, typename linalg_traits<L>::value_type a, col_major) {
    typename linalg_traits<L>::col_iterator it = mat_col_begin(l),
      ite = mat_col_end(l);
    for ( ; it != ite; ++it) scale(linalg_traits<L>::col(it), a);
  }

}

#endif //  __GMM_SCALED_H
