/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_vector.h : vectors.                                      */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002 Yves Renard.                                         */
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

#ifndef __GMM_VECTOR_H
#define __GMM_VECTOR_H

#include <map>

namespace gmm
{

  /*************************************************************************/
  /*                                                                       */
  /* Class ref_elt_vector: reference on a vector component.                */
  /*                                                                       */
  /*************************************************************************/


  template<class T, class V> class ref_elt_vector {

    V *pm;
    size_type l;
    
    public :

    operator T() const { return pm->r(l); }
    ref_elt_vector(V *p, size_type ll) : pm(p), l(ll) {}
    inline ref_elt_vector operator =(T v)
      { (*pm).w(l,v); return *this; }
    inline bool operator ==(T v) const
      { return ((*pm).r(l) == v); }
    inline bool operator !=(T v) const
      { return ((*pm).r(l) != v); }
    inline ref_elt_vector operator +=(T v)
      { (*pm).w(l,(*pm).r(l) + v); return *this; }
    inline ref_elt_vector operator -=(T v)
      { (*pm).w(l,(*pm).r(l) - v); return *this; }
    inline ref_elt_vector operator /=(T v)
      { (*pm).w(l,(*pm).r(l) / v); return *this; }
    inline ref_elt_vector operator *=(T v)
      { (*pm).w(l,(*pm).r(l) * v); return *this; }
    inline ref_elt_vector operator =(const ref_elt_vector &re)
      { *this = T(re); return *this; }
  };  
  
  template<class T, class V> T operator +(const ref_elt_vector<T, V> &re)
    { return T(re); }
  template<class T, class V> T operator -(const ref_elt_vector<T, V> &re)
    { return -T(re); }
  template<class T, class V> T operator +(const ref_elt_vector<T, V> &re, T v)
    { return T(re)+ v; }
  template<class T, class V> T operator +(T v, const ref_elt_vector<T, V> &re)
    { return v+ T(re); }
  template<class T, class V> T operator -(const ref_elt_vector<T, V> &re, T v)
    { return T(re)- v; }
  template<class T, class V> T operator -(T v, const ref_elt_vector<T, V> &re)
    { return v- T(re); }
  template<class T, class V> T operator *(const ref_elt_vector<T, V> &re, T v)
    { return T(re)* v; }
  template<class T, class V> T operator *(T v, const ref_elt_vector<T, V> &re)
    { return v* T(re); }
  template<class T, class V> T operator /(const ref_elt_vector<T, V> &re, T v)
    { return T(re)/ v; }
  template<class T, class V> T operator /(T v, const ref_elt_vector<T, V> &re)
    { return v/ T(re); }

  /*************************************************************************/
  /*                                                                       */
  /* Class wsvector: sparse vector optimized for random write operations.  */
  /*                                                                       */
  /*************************************************************************/
  
  template<class T> struct wsvector_iterator
    : public std::map<size_type, T>::iterator {
    typedef typename std::map<size_type, T>::iterator base_it_type;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    // typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    reference operator *() const { return (base_it_type::operator*()).second; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).first; }

    wsvector_iterator(void) {}
    wsvector_iterator(const base_it_type &it) : base_it_type(it) {}
  };

  template<class T> struct wsvector_const_iterator
    : public std::map<size_type, T>::const_iterator {
    typedef typename std::map<size_type, T>::const_iterator base_it_type;
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    // typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    reference operator *() const { return (base_it_type::operator*()).second; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).first; }

    wsvector_const_iterator(void) {}
    wsvector_const_iterator(const wsvector_iterator<T> &it)
      : base_it_type(it) {}
    wsvector_const_iterator(const base_it_type &it) : base_it_type(it) {}
  };


  template<class T> class wsvector : public std::map<size_type, T> {
  public:
    
    typedef std::map<size_type, T> base_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    //  typedef typename base_type::size_type size_type;

  protected:
    size_type nbl;
    
  public:
    void out_of_range_error(void) const;
    void clean(double eps);
    void resize(size_type l) { nbl = l;  /* + suprimer les elements en trop */}
    
    inline ref_elt_vector<T, wsvector<T> > operator [](size_type c)
    { return ref_elt_vector<T, wsvector<T> >(this, c); }

    inline void w(size_type c, const T &e) {
#   ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#   endif
      if (e == T(0)) erase(c);
      else base_type::operator [](c) = e;
    }

    inline T r(size_type c) const {
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      const_iterator it = lower_bound(c);
      if (it != end() && c == it->first) return it->second;
      else return T(0);
    }

    inline T operator [](size_type c) const { return r(c); }
    
    size_type nb_stored(void) const { return base_type::size(); }
    size_type size(void) const { return nbl; }

    /* Constructeurs */
    void init(size_type l) { nbl = l; this->clear(); }
    explicit wsvector(size_type l){ init(l); }
    wsvector(void) { init(0); }
  };

  template<class T>  void wsvector<T>::clean(double eps) {
    iterator it = begin(), itf = ++it, ite = end();
    for ( ; it != ite; ++itf)
      { if (dal::abs((*it).e) <= eps) { erase(it); it = itf; } else ++it; }
  }

  template<class T>  void wsvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  template <class T> struct wsvector_access {
    typedef wsvector<T> V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::iterator iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    reference operator()(const void *o, const iterator &,
			 const iterator &, size_type i)
    { return (*(const_cast<V *>((const V *)(o))))[i]; }

    value_type operator()(const void *o, const const_iterator &,
			 const const_iterator &, size_type i)
    { return (*(const V *)(o))[i]; }

  };

  template <class T> struct wsvector_clear {
    typedef wsvector<T> V;
    typedef typename linalg_traits<V>::iterator iterator;
    
    void operator()(const void *o, const iterator &, const iterator &)
    { (const_cast<V *>((const V *)(o)))->clear(); }
  };

  template <class T> struct linalg_traits<wsvector<T> > {
    typedef wsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_vector<T, wsvector<T> > reference;
    typedef wsvector_iterator<T>  iterator;
    typedef wsvector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef wsvector_access<T> access_type;
    typedef wsvector_clear<T> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return v.begin(); }
    static const_iterator begin(const this_type &v) { return v.begin(); }
    static iterator end(this_type &v) { return v.end(); }
    static const_iterator end(const this_type &v) { return v.end(); }
    static const void* origin(const this_type &v) { return &v; }
    static void do_clear(this_type &v)
    { clear_type()(origin(v), begin(v), end(v)); }
  };

  template<class T> std::ostream &operator <<
  (std::ostream &o, const wsvector<T>& v) { gmm::write(o,v); return o; }

#ifdef USING_BROKEN_GCC295
  template <class T> struct linalg_traits<const wsvector<T> >
    : public linalg_traits<wsvector<T> > {};
#endif

  /******* Optimized BLAS for wsvector<T> **********************************/

  template <class T> inline void copy(const wsvector<T> &v1, wsvector<T> &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    v2 = v1;
  }
  template <class T> inline
  void copy(const wsvector<T> &v1, const simple_vector_ref<wsvector<T> *> &v2){
    simple_vector_ref<wsvector<T> *>
      *svr = const_cast<simple_vector_ref<wsvector<T> *> *>(&v2);
    wsvector<T>
      *pv = const_cast<wsvector<T> *>((const wsvector<T> *)(v2.origin));
    if (vect_size(v1) != vect_size(v2))
	DAL_THROW(dimension_error,"dimensions mismatch");
    *pv = v1; svr->_begin = vect_begin(*pv); svr->_end = vect_end(*pv);
  }
  template <class T> inline
  void copy(const simple_vector_ref<const wsvector<T> *> &v1,
	    wsvector<T> &v2)
  { copy(*(const wsvector<T> *)(v1.origin), v2); }
  template <class T> inline
  void copy(const simple_vector_ref<wsvector<T> *> &v1, wsvector<T> &v2)
  { copy(*(const wsvector<T> *)(v1.origin), v2); }

  template <class T> inline void clean(wsvector<T> &v, double eps) {
    typename wsvector<T>::iterator it = v.begin(), ite = v.end(), itc;
    while (it != ite) 
      if (dal::abs((*it).second) <= eps)
	{ itc=it; ++it; v.erase(itc); } else ++it; 
  }

  template <class T>
  inline void clean(const simple_vector_ref<wsvector<T> *> &l, double eps) {
    simple_vector_ref<wsvector<T> *>
      *svr = const_cast<simple_vector_ref<wsvector<T> *> *>(&l);
    wsvector<T>
      *pv = const_cast<wsvector<T> *>((const wsvector<T> *)(l.origin));
    clean(*pv, eps);
    svr->_begin = vect_begin(*pv); svr->_end = vect_end(*pv);
  }

  /*************************************************************************/
  /*                                                                       */
  /*    rsvector: sparse vector optimized for linear algebra operations.   */
  /*                                                                       */
  /*************************************************************************/

  template<class T> struct _elt_rsvector {
    size_type c; T e;
    _elt_rsvector(void) {  }
    _elt_rsvector(size_type cc) { c = cc; }
    _elt_rsvector(size_type cc, const T &ee) { c = cc; e = ee; }
    bool operator < (const _elt_rsvector &a) const { return c < a.c; }
    bool operator == (const _elt_rsvector &a) const { return c == a.c; }
    bool operator != (const _elt_rsvector &a) const { return c != a.c; }
  };

  template<class T> struct rsvector_iterator {
    typedef typename std::vector<_elt_rsvector<T> >::iterator IT;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef rsvector_iterator<T> iterator;

    IT it;

    reference operator *() const { return it->e; }
    pointer operator->() const { return &(operator*()); }

    iterator &operator ++() { ++it; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator &operator --() { --it; return *this; }
    iterator operator --(int) { iterator tmp = *this; --(*this); return tmp; }

    bool operator ==(const iterator &i) const { return it == i.it; }
    bool operator !=(const iterator &i) const { return !(i == *this); }

    size_type index(void) const { return it->c; }
    rsvector_iterator(void) {}
    rsvector_iterator(const IT &i) : it(i) {}
  };

  template<class T> struct rsvector_const_iterator {
    typedef typename std::vector<_elt_rsvector<T> >::const_iterator IT;
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef rsvector_const_iterator<T> iterator;

    IT it;

    reference operator *() const { return it->e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return it->c; }

    iterator &operator ++() { ++it; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator &operator --() { --it; return *this; }
    iterator operator --(int) { iterator tmp = *this; --(*this); return tmp; }

    bool operator ==(const iterator &i) const { return it == i.it; }
    bool operator !=(const iterator &i) const { return !(i == *this); }

    rsvector_const_iterator(void) {}
    rsvector_const_iterator(const rsvector_iterator<T> &i) : it(i.it) {}
    rsvector_const_iterator(const IT &i) : it(i) {}
  };

  template<class T> class rsvector : public std::vector<_elt_rsvector<T> > {
  public:
    
    typedef std::vector<_elt_rsvector<T> > _base_type;
    typedef typename _base_type::iterator iterator;
    typedef typename _base_type::const_iterator const_iterator;
    typedef typename _base_type::size_type size_type;
    typedef T value_type;

  protected:
    size_type nbl;    	/* Nombre d'elements max.	        	  */
    
  public:

    void sup(size_type j);
    void out_of_range_error(void) const;
    void base_resize(size_type n) { _base_type::resize(n); }
    void resize(size_type l) { nbl = l; /* + suprimer les elements en trop*/ }
    
    ref_elt_vector<T, rsvector<T> > operator [](size_type c)
    { return ref_elt_vector<T, rsvector<T> >(this, c); }

    void w(size_type c, const T &e);
    T r(size_type c) const;

    inline T operator [](size_type c) const { return r(c); }
    
    size_type nb_stored(void) const { return _base_type::size(); }
    size_type size(void) const { return nbl; }
    void clear(void) { _base_type::resize(0); }
    /* Constructeurs */
    explicit rsvector(size_type l) : nbl(l) { }
    rsvector(void) : nbl(0) { }
  };

  template <class T> void rsvector<T>::sup(size_type j) {
    if (nb_stored() != 0) {
      _elt_rsvector<T> ev(j);
      iterator it = std::lower_bound(this->begin(), this->end(), ev);
      if (it != this->end() && it->c == j) {
	for (iterator ite = this->end() - 1; it != ite; ++it) *it = *(it+1);
	_base_type::resize(nb_stored()-1);
      }
    }
  }

  template <class T> void rsvector<T>::w(size_type c, const T &e) {
#   ifdef __GETFEM_VERIFY
    if (c >= nbl) out_of_range_error();
#   endif
    if (e == T(0)) sup(c);
    else {
      _elt_rsvector<T> ev(c, e);
      if (nb_stored() == 0) {
	_base_type::resize(1);
	*(this->begin()) = ev;
      }
      else {
	iterator it = std::lower_bound(this->begin(), this->end(), ev);
	if (it != this->end() && it->c == c) it->e = e;
	else {
	  size_type ind = it - this->begin();
	  _base_type::resize(nb_stored()+1);
	  it = this->begin() + ind;
	  for (iterator ite = this->end() - 1; ite != it; --ite) *ite = *(ite-1);
	  *it = ev;  // à verifier
	}
      }
    }
  }
  
  template <class T> T rsvector<T>::r(size_type c) const {
#   ifdef __GETFEM_VERIFY
    if (c >= nbl) out_of_range_error();
#   endif
    if (nb_stored() != 0) {
      _elt_rsvector<T> ev(c);
      const_iterator it = std::lower_bound(this->begin(), this->end(), ev);
      if (it != this->end() && it->c == c) return it->e;
    }
    return T(0);
  }

  template<class T>  void rsvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  template <class T> struct rsvector_access {
    typedef rsvector<T> V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::iterator iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    reference operator()(const void *o, const iterator &,
			 const iterator &, size_type i)
    { return (*(const_cast<V *>((const V *)(o))))[i]; }

    value_type operator()(const void *o, const const_iterator &,
			 const const_iterator &, size_type i)
    { return (*(const V *)(o))[i]; }

  };

  template <class T> struct rsvector_clear {
    typedef rsvector<T> V;
    typedef typename linalg_traits<V>::iterator iterator;
    
    void operator()(const void *o, const iterator &, const iterator &)
    { (*(const_cast<V *>((const V *)(o)))).clear(); }
  };

  template <class T> struct linalg_traits<rsvector<T> > {
    typedef rsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_vector<T, rsvector<T> > reference;
    typedef rsvector_iterator<T>  iterator;
    typedef rsvector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef rsvector_access<T> access_type;
    typedef rsvector_clear<T> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return iterator(v.begin()); }
    static const_iterator begin(const this_type &v)
    { return const_iterator(v.begin()); }
    static iterator end(this_type &v) { return iterator(v.end()); }
    static const_iterator end(const this_type &v)
      { return const_iterator(v.end()); }
    static const void* origin(const this_type &v) { return &v; }
    static void do_clear(this_type &v)
      { clear_type()(origin(v), begin(v), end(v)); }
  };

  template<class T> std::ostream &operator <<
  (std::ostream &o, const rsvector<T>& v) { gmm::write(o,v); return o; }

#ifdef USING_BROKEN_GCC295
  template <class T> struct linalg_traits<const rsvector<T> >
    : public linalg_traits<rsvector<T> > {};
#endif

  /******* Optimized BLAS for rsvector<T> **********************************/

  template <class T> inline void copy(const rsvector<T> &v1, rsvector<T> &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    v2 = v1;
  }
  template <class T> inline
  void copy(const rsvector<T> &v1, const simple_vector_ref<rsvector<T> *> &v2){
    simple_vector_ref<rsvector<T> *>
      *svr = const_cast<simple_vector_ref<rsvector<T> *> *>(&v2);
    rsvector<T>
      *pv = const_cast<rsvector<T> *>((const rsvector<T> *)(v2.origin));
    if (vect_size(v1) != vect_size(v2))
	DAL_THROW(dimension_error,"dimensions mismatch");
    *pv = v1; svr->_begin = vect_begin(*pv); svr->_end = vect_end(*pv);
  }
  template <class T> inline
  void copy(const simple_vector_ref<const rsvector<T> *> &v1,
	    rsvector<T> &v2)
  { copy(*(const rsvector<T> *)(v1.origin), v2); }
  template <class T> inline
  void copy(const simple_vector_ref<rsvector<T> *> &v1, rsvector<T> &v2)
  { copy(*(const rsvector<T> *)(v1.origin), v2); }

  template <class V, class T> inline void add(const V &v1, rsvector<T> &v2) {
    if ((const void *)(&v1) != (const void *)(&v2)) {
      if (vect_size(v1) != vect_size(v2))
	DAL_THROW(dimension_error,"dimensions mismatch");
	add_rsvector(v1, v2, typename linalg_traits<V>::storage_type());
    }
  }

  template <class V, class T> 
  inline void add_rsvector(const V &v1, rsvector<T> &v2, abstract_plain)
  { add(v1, v2, abstract_plain(), abstract_sparse()); }

  template <class V, class T> 
  inline void add_rsvector(const V &v1, rsvector<T> &v2, abstract_skyline)
  { add(v1, v2, abstract_skyline(), abstract_sparse()); }

  template <class V, class T> 
  void add_rsvector(const V &v1, rsvector<T> &v2, abstract_sparse) {
    typename linalg_traits<V>::const_iterator it1 = vect_const_begin(v1),
      ite1 = vect_const_end(v1);
    typename rsvector<T>::iterator it2 = v2.begin(), ite2 = v2.end(), it3;
    size_type nbc = 0, old_nbc = v2.nb_stored();
    for (; it1 != ite1 && it2 != ite2 ; ++nbc)
      if (it1.index() == it2->c) { ++it1; ++it2; }
      else if (it1.index() < it2->c) ++it1; else ++it2;
    for (; it1 != ite1; ++it1) ++nbc;
    for (; it2 != ite2; ++it2) ++nbc;

    v2.base_resize(nbc);
    it3 = v2.begin() + old_nbc; --it3;
    it2 = v2.end(); --it2; ite2 = v2.begin(); --ite2;
    it1 = vect_end(v1); --it1; ite1 = vect_const_begin(v1); --ite1;

    for (; it1 != ite1 && it3 != ite2; --it2) {
      if (it3->c > it1.index()) { *it2 = *it3; --it3; }
      else if (it3->c == it1.index()) { *it2=*it3; it2->e+=*it1; --it3; --it1;}
      else { it2->c = it1.index(); it2->e = *it1; --it1; }
    }
    for (; it1 != ite1; --it2) { it2->c = it1.index(); it2->e = *it1; --it1; }
  }

  template <class V, class T> void copy(const V &v1, rsvector<T> &v2) {
    if ((const void *)(&v1) != (const void *)(&v2)) {
      if (vect_size(v1) != vect_size(v2))
	DAL_THROW(dimension_error,"dimensions mismatch");
#       ifdef __GETFEM_VERIFY
        if (linalg_origin(v1) == linalg_origin(v2))
	  DAL_WARNING(2, "a conflict is possible in vector copy\n");
#       endif
	copy_rsvector(v1, v2, typename linalg_traits<V>::storage_type());
    }
  }

  template <class V, class T> 
  void copy_rsvector(const V &v1, rsvector<T> &v2, abstract_plain)
  { copy_vect(v1, v2, abstract_plain(), abstract_sparse()); }

  template <class V, class T> 
  void copy_rsvector(const V &v1, rsvector<T> &v2, abstract_skyline)
  { copy_vect(v1, v2, abstract_skyline(), abstract_sparse()); }

  template <class V, class T> // à refaire
  void copy_rsvector(const V &v1, rsvector<T> &v2, abstract_sparse) {
     typename linalg_traits<V>::const_iterator it = vect_const_begin(v1),
      ite = vect_const_end(v1);
    std::vector<size_type> tab(100);
    size_type i = 0;
    for (; it != ite; ++it)
      if ((*it) != typename linalg_traits<V>::value_type(0)) {
	tab[i++] = it.index();
	if (i >= tab.size()) tab.resize(i + 100);
      }
    v2.base_resize(i);
    if (i > 0) {
      typename rsvector<T>::iterator it2 = v2.begin(), ite2 = v2.end();
      for (i = 0; it2 != ite2; ++it2, ++i)
	{ it2->c = tab[i]; it2->e = v1[tab[i]]; }
    }
  }
  
  template <class T> inline void clean(rsvector<T> &v, double eps) {
    typename rsvector<T>::iterator it = v.begin(), ite = v.end();
    for (; it != ite; ++it) if (dal::abs((*it).e) <= eps) break;
    if (it != ite) {
      typename rsvector<T>::iterator itc = it;
      size_type erased = 1;
      for (++it; it != ite; ++it)
	{ *itc = *it; if (dal::abs((*it).e) <= eps) ++erased; else ++itc; }
      v.base_resize(v.nb_stored() - erased);
    }
  }

  template <class T>
  inline void clean(const simple_vector_ref<rsvector<T> *> &l, double eps) {
    simple_vector_ref<rsvector<T> *>
      *svr = const_cast<simple_vector_ref<rsvector<T> *> *>(&l);
    rsvector<T>
      *pv = const_cast<rsvector<T> *>((const rsvector<T> *)(l.origin));
    clean(*pv, eps);
    svr->_begin = vect_begin(*pv); svr->_end = vect_end(*pv);
  }
  

  /*************************************************************************/
  /*                                                                       */
  /* Class slvector: 'sky-line' vector.                                    */
  /*                                                                       */
  /*************************************************************************/

  template<class T> struct slvector_iterator
  {
    typedef T value_type;
    typedef T *pointer;
    typedef T &reference;
    typedef ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t size_type;
    typedef slvector_iterator<T> iterator;
    typedef typename std::vector<T>::iterator base_iterator;

    base_iterator it;
    size_type shift;
    
   
    iterator &operator ++()
    { ++it; ++shift; return *this; }
    iterator &operator --()
    { --it; --shift; return *this; }
    iterator operator ++(int)
    { iterator tmp = *this; ++(*(this)); return tmp; }
    iterator operator --(int)
    { iterator tmp = *this; --(*(this)); return tmp; }
    iterator &operator +=(difference_type i)
    { it += i; shift += i; return *this; }
    iterator &operator -=(difference_type i)
    { it -= i; shift -= i; return *this; }
    iterator operator +(difference_type i) const
    { iterator tmp = *this; return (tmp += i); }
    iterator operator -(difference_type i) const
    { iterator tmp = *this; return (tmp -= i); }
    difference_type operator -(const iterator &i) const
    { return it - i.it; }
	
    reference operator *() const
    { return *it; }
    reference operator [](int ii)
    { return *(it + ii); }
    
    bool operator ==(const iterator &i) const
    { return it == i.it; }
    bool operator !=(const iterator &i) const
    { return !(i == *this); }
    bool operator < (const iterator &i) const
    { return it < i.it; }
    size_type index(void) const { return shift; }

    slvector_iterator(void) {}
    slvector_iterator(const base_iterator &iter, size_type s)
      : it(iter), shift(s) {}
  };

  template<class T> struct slvector_const_iterator
  {
    typedef T value_type;
    typedef const T *pointer;
    typedef const T &reference;
    typedef ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t size_type;
    typedef slvector_const_iterator<T> iterator;
    typedef typename std::vector<T>::const_iterator base_iterator;

    base_iterator it;
    size_type shift;
    
   
    iterator &operator ++()
    { ++it; ++shift; return *this; }
    iterator &operator --()
    { --it; --shift; return *this; }
    iterator operator ++(int)
    { iterator tmp = *this; ++(*(this)); return tmp; }
    iterator operator --(int)
    { iterator tmp = *this; --(*(this)); return tmp; }
    iterator &operator +=(difference_type i)
    { it += i; shift += i; return *this; }
    iterator &operator -=(difference_type i)
    { it -= i; shift -= i; return *this; }
    iterator operator +(difference_type i) const
    { iterator tmp = *this; return (tmp += i); }
    iterator operator -(difference_type i) const
    { iterator tmp = *this; return (tmp -= i); }
    difference_type operator -(const iterator &i) const
    { return it - i.it; }
	
    value_type operator *() const
    { return *it; }
    value_type operator [](int ii)
    { return *(it + ii); }
    
    bool operator ==(const iterator &i) const
    { return it == i.it; }
    bool operator !=(const iterator &i) const
    { return !(i == *this); }
    bool operator < (const iterator &i) const
    { return it < i.it; }
    size_type index(void) const { return shift; }

    slvector_const_iterator(void) {}
    slvector_const_iterator(const slvector_iterator<T>& iter)
      : it(iter.it), shift(iter.shift) {}
    slvector_const_iterator(const base_iterator &iter, size_type s)
      : it(iter), shift(s) {}
  };


  template <class T> class slvector {
    
  public :
    typedef slvector_iterator<T> iterators;
    typedef slvector_const_iterator<T> const_iterators;
    typedef typename std::vector<T>::size_type size_type;
    typedef T value_type;

  protected :
    std::vector<T> data;
    size_type shift;
    size_type _size;


  public :

    void out_of_range_error(void) const;
    size_type size(void) const { return _size; }
    size_type first(void) const { return shift; }
    size_type last(void) const { return shift + data.size(); }
    ref_elt_vector<T, slvector<T> > operator [](size_type c)
    { return ref_elt_vector<T, slvector<T> >(this, c); }

    typename std::vector<T>::iterator data_begin(void) { return data.begin(); }
    typename std::vector<T>::iterator data_end(void) { return data.end(); }
    typename std::vector<T>::const_iterator data_begin(void) const
      { return data.begin(); }
    typename std::vector<T>::const_iterator data_end(void) const
      { return data.end(); }

    void w(size_type c, const T &e);
    T r(size_type c) const {
#   ifdef __GETFEM_VERIFY
      if (c >= _size) out_of_range_error();
#   endif
      if (c < shift || c >= shift + data.size()) return T(0);
      return data[c - shift];
    }

    inline T operator [](size_type c) const { return r(c); }

    void clear(void) { data.resize(0); shift = 0; }

    slvector(void) : data(0), shift(0), _size(0) {}
    explicit slvector(size_type l) : data(0), shift(0), _size(l) {}
    slvector(size_type l, size_type d, size_type s)
      : data(d), shift(d), _size(l) {}

  };

  template<class T>  void slvector<T>::w(size_type c, const T &e) {
    // cout << "vecteur avant : " << *this << " ajout à l'indice " << c << " de " << e << endl;
#   ifdef __GETFEM_VERIFY
      if (c >= _size) out_of_range_error();
#   endif
      size_type s = data.size();
      if (c < shift) { // à verifier
	data.resize(s + shift - c); 
	typename std::vector<T>::iterator it = data.begin(), ite = data.end();
	typename std::vector<T>::iterator it2 = ite - 1,
	  it3 = ite - shift + c - 1;
	for (; it3 != it; --it3, --it2) *it2 = *it3;
	std::fill(it, it + shift - c, T(0));
	shift = c;
      }
      else if (c >= shift + s && s) {
	data.resize(c - shift + 1);
	std::fill(data.begin() + s, data.end(), T(0));
      }
      else { data.resize(1); shift = c; }
      data[c - shift] = e;
      // cout << "résultat : " << *this << endl;
    }

  template<class T>  void slvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  template <class T> struct slvector_access {
    typedef slvector<T> V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::iterator iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    reference operator()(const void *o, const iterator &,
			 const iterator &, size_type i)
    { return (*(const_cast<V *>((const V *)(o))))[i]; }

    value_type operator()(const void *o, const const_iterator &,
			 const const_iterator &, size_type i)
    { return (*(const V *)(o))[i]; }

  };

  template <class T> struct slvector_clear {
    typedef slvector<T> V;
    typedef typename linalg_traits<V>::iterator iterator;
    
    void operator()(const void *o, const iterator &, const iterator &)
    { (*(const_cast<V *>((const V *)(o)))).clear(); }
  };

  template <class T> struct linalg_traits<slvector<T> > {
    typedef slvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_vector<T, slvector<T> > reference;
    typedef slvector_iterator<T>  iterator;
    typedef slvector_const_iterator<T> const_iterator;
    typedef abstract_skyline storage_type;
    typedef slvector_access<T> access_type;
    typedef slvector_clear<T> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v)
      { return iterator(v.data_begin(), v.first()); }
    static const_iterator begin(const this_type &v)
      { return const_iterator(v.data_begin(), v.first()); }
    static iterator end(this_type &v)
      { return iterator(v.data_end(), v.last()); }
    static const_iterator end(const this_type &v)
      { return const_iterator(v.data_end(), v.last()); }
    static const void* origin(const this_type &v) { return &v; }
    static void do_clear(this_type &v)
      { clear_type()(origin(v), begin(v), end(v)); }
  };

  template<class T> std::ostream &operator <<
  (std::ostream &o, const slvector<T>& v) { gmm::write(o,v); return o; }

#ifdef USING_BROKEN_GCC295
  template <class T> struct linalg_traits<const slvector<T> >
    : public linalg_traits<slvector<T> > {};
#endif


}

#endif /* __GMM_VECTOR_H */
