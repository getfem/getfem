/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_interface_bgeot.h : interface to bgeot vsvectors.        */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2003  Yves Renard.                                   */
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

#ifndef __GMM_INTERFACE_BGEOT_H
#define __GMM_INTERFACE_BGEOT_H


namespace gmm {

  /* ********************************************************************* */
  /*		                                         	 	   */
  /*		Traits for bgeot objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class T> struct linalg_traits<bgeot::vsvector<T> > {
    typedef bgeot::vsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_dense storage_type;
    typedef dense_access<iterator,const_iterator> access_type;
    typedef dense_clear<iterator> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return v.begin(); }
    static const_iterator begin(const this_type &v) { return v.begin(); }
    static iterator end(this_type &v) { return v.end(); }
    static const_iterator end(const this_type &v) { return v.end(); }
    static const void* origin(const this_type &v) { return &v; }
    static void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

#ifdef USING_BROKEN_GCC295
  template <class T> struct linalg_traits<const bgeot::vsvector<T> > 
    : public linalg_traits<bgeot::vsvector<T> > {};
#endif

  template <class VECT> struct linalg_traits<bgeot::PT<VECT> > {
    typedef bgeot::PT<VECT> this_type;
    typedef linalg_false is_reference;
    typedef typename linalg_traits<VECT>::value_type value_type;
    typedef typename linalg_traits<VECT>::reference reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<VECT>::iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    typedef dense_access<iterator,const_iterator> access_type;
    typedef dense_clear<iterator> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return v.begin(); }
    static const_iterator begin(const this_type &v) { return v.begin(); }
    static iterator end(this_type &v) { return v.end(); }
    static const_iterator end(const this_type &v) { return v.end(); }
    static const void* origin(const this_type &v) { return &v; }
    static void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

#ifdef USING_BROKEN_GCC295
  template <class VECT> struct linalg_traits<const bgeot::PT<VECT> >
  : public linalg_traits<bgeot::PT<VECT> > {};
#endif

  template<class ITER, class MIT> struct dense_compressed_iterator
  {
    typedef ITER value_type;
    typedef ITER *pointer;
    typedef ITER &reference;
    typedef ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t size_type;
    typedef dense_compressed_iterator<ITER, MIT> iterator;

    ITER it;
    size_type N, nrows, ncols;
    const void *origin;
    
    iterator operator ++(int) { iterator tmp = *this; it += N; return tmp; }
    iterator operator --(int) { iterator tmp = *this; it -= N; return tmp; }
    iterator &operator ++()   { it += N; return *this; }
    iterator &operator --()   { it -= N; return *this; }
    iterator &operator +=(difference_type i) { it += i * N; return *this; }
    iterator &operator -=(difference_type i) { it -= i * N; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return (it - i.it) / N; }

    ITER operator *() const { return it; }
    ITER operator [](int ii) const { return it + ii * N; }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    dense_compressed_iterator(void) {}
    dense_compressed_iterator(const dense_compressed_iterator<MIT, MIT> &ii)
    : it(ii.it), N(ii.N), nrows(ii.nrows),ncols(ii.ncols),origin(ii.origin) {}
    dense_compressed_iterator(const ITER &iter, size_type n, size_type r,
			      size_type c, const void *o)
      : it(iter), N(n), nrows(r), ncols(c), origin(o) { }
    
  };

}


#endif //  __GMM_INTERFACE_BGEOT_H
