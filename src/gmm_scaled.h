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
    value_type operator [](size_type ii) const { return it[in+ii] * r; }
    
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
      : _begin(vect_begin(v)), _end(vect_end(v)), origin(linalg_origin(v)),
	_size(vect_size(v)), r(rr) {}

    reference operator[](size_type i) const
    { return access_type()(origin, _begin, _end, i); }
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
    size_type size(const this_type &v) { return v._size; }
    iterator begin(this_type &v)
    { return iterator(v._begin, v.r); }
    const_iterator begin(const this_type &v)
    { return const_iterator(v._begin, v.r); }
    iterator end(this_type &v)
    { return iterator(v._end, v.r); }
    const_iterator end(const this_type &v)
    { return const_iterator(v._end, v.r); }
    const void* origin(const this_type &v) { return v.origin; }
  };

  // for GCC 2.95
  template <class V> struct linalg_traits<const scaled_vector_const_ref<V> > 
    : public linalg_traits<scaled_vector_const_ref<V> > {};

  template <class V> inline
  scaled_vector_const_ref<V> scaled(const V &v,
				    typename linalg_traits<V>::value_type x)
  { return scaled_vector_const_ref<V>(v, x); }
  
  /* ******************************************************************** */
  /*		scale                                    	          */
  /* ******************************************************************** */

  template <class L> void scale(L& l, typename linalg_traits<L>::value_type a)
  {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for ( ; it != ite; ++it) *it *= a;
  }

  template <class L> inline
  void scale(const L& l, typename linalg_traits<L>::value_type a)
  { scale(linalg_const_cast(l), a); }

}

#endif //  __GMM_SCALED_H
