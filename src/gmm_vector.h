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

#ifndef __GMM_VECTOR_H
#define __GMM_VECTOR_H

#include <dal_tree_sorted.h>

namespace gmm
{
  /************************************************************************/
  /*		Class wsvector: sparse vector.				  */
  /************************************************************************/

  template<class T> struct _elt_wsvector {
    size_type c; T e;
    _elt_wsvector(void) { c = size_type(-1); e = T(0); }
    _elt_wsvector(size_type cc) { c = cc; e = T(0); }
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
      comp_elt_wsvector<T>,3>::iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::iterator base_it_type;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;

    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).c; }
    wsvector_iterator(void) {}
    wsvector_iterator(const base_it_type &it) : base_it_type(it) {}
  };

  template<class T> struct wsvector_const_iterator
    : public dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::const_iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_wsvector<T>,
      comp_elt_wsvector<T>,3>::const_iterator base_it_type;
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).c; }
    wsvector_const_iterator(void) {}
    wsvector_const_iterator(const wsvector_iterator<T> &it) : base_it_type(it) {}
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
    typedef typename _base_type::size_type size_type;
    typedef T value_type;

  protected:
    size_type nbl;    	/* Nombre d'elements max.	        	  */
    
    void init(size_type l) { nbl = l; clear(); }
    
  public:

    void out_of_range_error(void) const;
    void clean(double eps) {
      tas_iterator it = tas_begin(), ite = tas_end();
      for ( ; it != ite; ++it)
	if (dal::abs((*it).e) <= eps)
	  { it->e = T(0); sup(it.index()); it->c = size_type(-1); }
    }
    
    void resize(size_type l) { nbl = l; }
    
    inline ref_elt_wsvector<T> operator [](size_type c)
    { return ref_elt_wsvector<T>(this, c); }

    inline void w(size_type c, const T &e) { 
      _elt_wsvector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = add_norepeat(ev);
      _base_type::operator[](i).e = e;
      if (e == T(0)) { sup(i); _base_type::operator[](i).c = size_type(-1); } 
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

    inline T &ref(size_type c) {  
      _elt_wsvector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = add_norepeat(ev);
      return (_base_type::operator[](i)).e;
    }

    inline T operator [](size_type c) const { return r(c); }
    
    bool stored(size_type c) {
      _elt_wsvector<T> ev(c); 
      size_type i=search(ev); 
      return (i!=size_type(-1));
    }
    size_type nb_stored(void) { return card(); }
      
    /* Operations algebriques sur les vecteurs. */

    size_type size(void) const { return nbl; }

    /* Constructeurs */

    wsvector(size_type l){ init(l); }
    wsvector(void) { init(1); }

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

  template <class T> struct linalg_traits<wsvector<T> > {
    typedef wsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_wsvector<T> reference_type;
    typedef wsvector_iterator<T>  iterator;
    typedef wsvector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };
}

#endif /* __GMM_VECTOR_H */
