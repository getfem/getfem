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

#include <dal_tree_sorted.h>

namespace gmm
{
  /************************************************************************/
  /* Class wsvector: sparse vector optimized for random write operations. */
  /************************************************************************/

  template<class T> struct _elt_wsvector {
    size_type c; T e;
    _elt_wsvector(void) { c = size_type(-1); }
    _elt_wsvector(size_type cc) { c = cc; }
    _elt_wsvector(size_type cc, const T &ee) { c = cc; e = ee; }
  };
  
  template<class T> struct comp_elt_wsvector
    : public std::binary_function<_elt_wsvector<T>, _elt_wsvector<T>, int> {
    int operator()(const _elt_wsvector<T>& m, const _elt_wsvector<T>& n) const
    { int d = m.c-n.c; if (d<0) return -1; if (d>0) return 1; return 0; }
  };
  
  template<class T> class ref_elt_wsvector;

  template<class T> struct wsvector_iterator
    : public dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::tas_iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::tas_iterator base_it_type;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::forward_iterator_tag iterator_category;

    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    
    size_type index(void) const { return (base_it_type::operator*()).c; }
    wsvector_iterator(void) {}
    wsvector_iterator(const base_it_type &it) : base_it_type(it) {}
  };

  template<class T> struct wsvector_const_iterator
    : public dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::const_tas_iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::const_tas_iterator base_it_type;
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::forward_iterator_tag iterator_category;
    
    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).c; }
    wsvector_const_iterator(void) {}
    wsvector_const_iterator(const wsvector_iterator<T> &it)
      : base_it_type(it) {}
    wsvector_const_iterator(const base_it_type &it) : base_it_type(it) {}
  };


  template<class T> class wsvector
    : public dal::dynamic_tree_sorted<_elt_wsvector<T>, comp_elt_wsvector<T>,3>
  {
  public:
    
    typedef dal::dynamic_tree_sorted<_elt_wsvector<T>, comp_elt_wsvector<T>,3>
      _base_type;
    typedef typename _base_type::tas_iterator tas_iterator;
    typedef typename _base_type::const_tas_iterator const_tas_iterator;
    typedef typename _base_type::sorted_iterator sorted_iterator;
    typedef typename _base_type::const_sorted_iterator const_sorted_iterator;
    typedef typename _base_type::size_type size_type;
    typedef T value_type;

  protected:
    size_type nbl;    	/* Nombre d'elements max.	        	  */
    
  public:

    void out_of_range_error(void) const;
    void clean(double eps) {
      tas_iterator it = tas_begin(), ite = tas_end();
      for ( ; it != ite; ++it)
	if (dal::abs((*it).e) <= eps) sup(it.index());
    }
    
    void resize(size_type l) { nbl = l;  /* + suprimer les elements en trop */}
    
    inline ref_elt_wsvector<T> operator [](size_type c)
    { return ref_elt_wsvector<T>(this, c); }

    inline void w(size_type c, const T &e) {
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      _elt_wsvector<T> ev(c, e);
      if (e == T(0)) {
	size_type i = search(ev);
	if (i != size_type(-1)) sup(i);
      }
      else
	add_norepeat(ev, true);
    }

    inline T r(size_type c) const {  
      _elt_wsvector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = search(ev);
      if (i == size_type(-1)) return T(0);
      return (_base_type::operator[](i)).e;
    }

    inline T operator [](size_type c) const { return r(c); }
    
    size_type nb_stored(void) const { return card(); }
    size_type size(void) const { return nbl; }

    /* Constructeurs */
    void init(size_type l) { nbl = l; clear(); }
    wsvector(size_type l){ init(l); }
    wsvector(void) { init(0); }

  };

  template<class T>  void wsvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  /*********  intermediary structure for r/w operation.	*******************/
  
  /* revoir cette stratégie pour l'interfacage générique de abstract_linalg */

  template<class T> class ref_elt_wsvector { /* ref. on an matrix element.  */

    wsvector<T> *pm;
    size_type l;
    
    public :

      operator T() const { return pm->r(l); }
    ref_elt_wsvector(wsvector<T> *p, size_type ll)
      { pm = p; l = ll; }
    inline ref_elt_wsvector operator =(T v)
      { (*pm).w(l,v); return *this; }
    inline bool operator ==(T v) const
      { return ((*pm).r(l) == v); }
    inline bool operator !=(T v) const
      { return ((*pm).r(l) != v); }
    inline ref_elt_wsvector operator +=(T v)
      { (*pm).w(l,(*pm).r(l) + v); return *this; }
    inline ref_elt_wsvector operator -=(T v)
      { (*pm).w(l,(*pm).r(l) - v); return *this; }
    inline ref_elt_wsvector operator /=(T v)
      { (*pm).w(l,(*pm).r(l) / v); return *this; }
    inline ref_elt_wsvector operator *=(T v)
      { (*pm).w(l,(*pm).r(l) * v); return *this; }
    inline ref_elt_wsvector operator =(const ref_elt_wsvector &re)
      { *this = T(re); return *this; }
  };  
  
  template<class T> T operator +(const ref_elt_wsvector<T> &re)
    { return T(re); }
  template<class T> T operator -(const ref_elt_wsvector<T> &re)
    { return -T(re); }
  template<class T> T operator +(const ref_elt_wsvector<T> &re, T v)
    { return T(re)+ v; }
  template<class T> T operator +(T v, const ref_elt_wsvector<T> &re)
    { return v+ T(re); }
  template<class T> T operator -(const ref_elt_wsvector<T> &re, T v)
    { return T(re)- v; }
  template<class T> T operator -(T v, const ref_elt_wsvector<T> &re)
    { return v- T(re); }
  template<class T> T operator *(const ref_elt_wsvector<T> &re, T v)
    { return T(re)* v; }
  template<class T> T operator *(T v, const ref_elt_wsvector<T> &re)
    { return v* T(re); }
  template<class T> T operator /(const ref_elt_wsvector<T> &re, T v)
    { return T(re)/ v; }
  template<class T> T operator /(T v, const ref_elt_wsvector<T> &re)
    { return v/ T(re); }

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
    typedef ref_elt_wsvector<T> reference;
    typedef wsvector_iterator<T>  iterator;
    typedef wsvector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef wsvector_access<T> access_type;
    typedef wsvector_clear<T> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.tas_begin(); }
    const_iterator begin(const this_type &v) { return v.tas_begin(); }
    iterator end(this_type &v) { return v.tas_end(); }
    const_iterator end(const this_type &v) { return v.tas_end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

  template<class T> std::ostream &operator <<
  (std::ostream &o, const wsvector<T>& v) { gmm::write(o,v); return o; }

  // for GCC 2.95
  template <class T> struct linalg_traits<const wsvector<T> >
    : public linalg_traits<wsvector<T> > {};

  /************************************************************************/
  /*    rsvector: sparse vector optimized for linear algebra operations.  */
  /************************************************************************/

  template<class T> struct _elt_rsvector {
    size_type c; T e;
    _elt_rsvector(void) {  }
    _elt_rsvector(size_type cc) { c = cc; }
    _elt_rsvector(size_type cc, const T &ee) { c = cc; e = ee; }
    bool operator < (const _elt_rsvector &a) const { return c < a.c; }
    bool operator == (const _elt_rsvector &a) const { return c == a.c; }
    bool operator != (const _elt_rsvector &a) const { return c != a.c; }
  };
  
  template<class T> class ref_elt_rsvector;

  template<class T> struct rsvector_iterator {
    typedef typename std::vector<_elt_rsvector<T> >::iterator IT;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef rsvector_iterator<T> iterator;

    IT it;

    reference operator *() const { return it->e; }
    pointer operator->() const { return &(operator*()); }

    iterator &operator ++() { ++it; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }

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
    void clean(double eps) {
      for (size_type i = 0; i < nb_stored(); ++i)
	if (dal::abs(_base_type::operator[](i)) <= eps) sup(i);
    }
    void base_resize(size_type n) { _base_type::resize(n); }
    void resize(size_type l) { nbl = l; /* + suprimer les elements en trop*/ }
    
    ref_elt_rsvector<T> operator [](size_type c)
    { return ref_elt_rsvector<T>(this, c); }

    void w(size_type c, const T &e);
    T r(size_type c) const;

    inline T operator [](size_type c) const { return r(c); }
    
    size_type nb_stored(void) const { return _base_type::size(); }
    size_type size(void) const { return nbl; }
    void clear(void) { _base_type::resize(0); }
    /* Constructeurs */
    rsvector(size_type l) : nbl(l) { }
    rsvector(void) : nbl(0) { }
  };

  template <class T> void rsvector<T>::sup(size_type j) {
    if (nb_stored() != 0) {
      _elt_rsvector<T> ev(j);
      iterator it = std::lower_bound(begin(), end(), ev);
      if (it != end() && it->c == j) {
	for (iterator ite = end() - 1; it != ite; ++it) *it = *(it+1);
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
	*(begin()) = ev;
      }
      else {
	iterator it = std::lower_bound(begin(), end(), ev);
	if (it != end() && it->c == c) it->e = e;
	else {
	  size_type ind = it - begin();
	  _base_type::resize(nb_stored()+1);
	  it = begin() + ind;
	  for (iterator ite = end() - 1; ite != it; --ite) *ite = *(ite-1);
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
      const_iterator it = std::lower_bound(begin(), end(), ev);
      if (it != end() && it->c == c) return it->e;
    }
    return T(0);
  }


  template<class T>  void rsvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  /*********  intermediary structure for r/w operation.	*******************/
  
  /* revoir cette stratégie pour l'interfacage générique de abstract_linalg */

  template<class T> class ref_elt_rsvector { /* ref. on an matrix element.  */

    rsvector<T> *pm;
    size_type l;
    
    public :

      operator T() const { return pm->r(l); }
    ref_elt_rsvector(rsvector<T> *p, size_type ll)
      { pm = p; l = ll; }
    inline ref_elt_rsvector operator =(T v)
      { (*pm).w(l,v); return *this; }
    inline bool operator ==(T v) const
      { return ((*pm).r(l) == v); }
    inline bool operator !=(T v) const
      { return ((*pm).r(l) != v); }
    inline ref_elt_rsvector operator +=(T v)
      { (*pm).w(l,(*pm).r(l) + v); return *this; }
    inline ref_elt_rsvector operator -=(T v)
      { (*pm).w(l,(*pm).r(l) - v); return *this; }
    inline ref_elt_rsvector operator /=(T v)
      { (*pm).w(l,(*pm).r(l) / v); return *this; }
    inline ref_elt_rsvector operator *=(T v)
      { (*pm).w(l,(*pm).r(l) * v); return *this; }
    inline ref_elt_rsvector operator =(const ref_elt_rsvector &re)
      { *this = T(re); return *this; }
  };  
  
  template<class T> T operator +(const ref_elt_rsvector<T> &re)
    { return T(re); }
  template<class T> T operator -(const ref_elt_rsvector<T> &re)
    { return -T(re); }
  template<class T> T operator +(const ref_elt_rsvector<T> &re, T v)
    { return T(re)+ v; }
  template<class T> T operator +(T v, const ref_elt_rsvector<T> &re)
    { return v+ T(re); }
  template<class T> T operator -(const ref_elt_rsvector<T> &re, T v)
    { return T(re)- v; }
  template<class T> T operator -(T v, const ref_elt_rsvector<T> &re)
    { return v- T(re); }
  template<class T> T operator *(const ref_elt_rsvector<T> &re, T v)
    { return T(re)* v; }
  template<class T> T operator *(T v, const ref_elt_rsvector<T> &re)
    { return v* T(re); }
  template<class T> T operator /(const ref_elt_rsvector<T> &re, T v)
    { return T(re)/ v; }
  template<class T> T operator /(T v, const ref_elt_rsvector<T> &re)
    { return v/ T(re); }

  template <class T> struct rsvector_access {
    typedef rsvector<T> V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::iterator iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    reference operator()(const void *o, const iterator &,
			 const iterator &, size_type i)
    { return (*(V *)(o))[i]; }

    value_type operator()(const void *o, const const_iterator &,
			 const const_iterator &, size_type i)
    { return (*(const V *)(o))[i]; }

  };

  template <class T> struct rsvector_clear {
    typedef rsvector<T> V;
    typedef typename linalg_traits<V>::iterator iterator;
    
    void operator()(const void *o, const iterator &, const iterator &)
    { (*(V *)(o)).clear(); }
  };

  template <class T> struct linalg_traits<rsvector<T> > {
    typedef rsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_rsvector<T> reference;
    typedef rsvector_iterator<T>  iterator;
    typedef rsvector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef rsvector_access<T> access_type;
    typedef rsvector_clear<T> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return iterator(v.begin()); }
    const_iterator begin(const this_type &v)
    { return const_iterator(v.begin()); }
    iterator end(this_type &v) { return iterator(v.end()); }
    const_iterator end(const this_type &v) { return const_iterator(v.end()); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

  template<class T> std::ostream &operator <<
  (std::ostream &o, const rsvector<T>& v) { gmm::write(o,v); return o; }

  // for GCC 2.95
  template <class T> struct linalg_traits<const rsvector<T> >
    : public linalg_traits<rsvector<T> > {};

  template <class V, class T> void copy(const V &v1, rsvector<T> &v2) {
    if ((const void *)(&v1) != (const void *)(&v2)) {
      if (vect_size(v1) != vect_size(v2))
	DAL_THROW(dimension_error,"dimensions mismatch");
#       ifdef __GETFEM_VERIFY
        if (linalg_origin(v1) == linalg_origin(v2))
	  cerr << "Warning : a conflict is possible in vector copy\n";
#       endif
	copy_rsvector(v1, v2, typename linalg_traits<V>::storage_type());
    }
  }

  template <class V, class T> 
  void copy_rsvector(const V &v1, rsvector<T> &v2, abstract_plain) {
    cout << "routine à verifier\n";
    typename linalg_traits<V>::const_iterator it = vect_begin(v1),
      ite = vect_end(v1);
    std::vector<size_type> tab(100);
    size_type i = 0, j = 0;
    for (; it != _end; ++it, ++j)
      if ((*it) != typename linalg_traits<V>::value_type(0)) {
	tab[i++] = j;
	if (i >= tab.size()) tab.resize(i + 100);
      }
    v2.base_resize(i);
    if (i > 0) {
      typename rsvector<T>::iterator it2 = v2.begin(), ite2 = v2.end();
      for (i = 0; it2 != ite2; ++it2, ++i) { it2->c = tab[i]; it2->e = v1[tab[i]]; }
    }
  }

  template <class V, class T> 
  void copy_rsvector(const V &v1, rsvector<T> &v2, abstract_sparse) {
     typename linalg_traits<V>::const_iterator it = vect_begin(v1),
      ite = vect_end(v1);
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
      for (i = 0; it2 != ite2; ++it2, ++i) { it2->c = tab[i]; it2->e = v1[tab[i]]; }
      std::sort(v2.begin(), v2.end());
    }
  }
  
  template <class T> 
  void copy_rsvector(const wsvector<T> &v1, rsvector<T> &v2, abstract_sparse) {
     v2.base_resize(v1.nb_stored());
     std::copy(v1.sorted_begin(), v1.sorted_end(), v2.begin());
  }
  
}

#endif /* __GMM_VECTOR_H */
