/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_svector.h : Sparse vectors.                            */
/*     									   */
/*                                                                         */
/* Date : February 01, 1998.                                               */
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

#ifndef __BGEOT_SVECTOR_H
#define __BGEOT_SVECTOR_H

#include <dal_tree_sorted.h>
#include <bgeot_matrix.h>

namespace bgeot
{
  /************************************************************************/
  /*		Class svector: sparse vector.				  */
  /************************************************************************/

  template<class T> struct _elt_svector {
    size_type c; T e;
    _elt_svector(void) { c = size_type(-1); e = T(0); }
    _elt_svector(size_type cc) { c = cc; e = T(0); }
    _elt_svector(size_type cc, const T &ee) { c = cc; e = ee; }
  };
  
  template<class T> struct comp_elt_svector
    : public std::binary_function<_elt_svector<T>, _elt_svector<T>, int> {
    int operator()(const _elt_svector<T>& m, const _elt_svector<T>& n) const
    { int d = m.c-n.c; if (d<0) return -1; if (d>0) return 1; return 0; }
  };
  
  template<class T> class ref_elt_svector;

  template<class T> struct svector_iterator
    : public dal::dynamic_tree_sorted<_elt_svector<T>,
      comp_elt_svector<T>,3>::iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_svector<T>,
      comp_elt_svector<T>,3>::iterator base_it_type;
    typedef T                   value_type;
    typedef value_type*         pointer;
    typedef value_type&         reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;

    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).c; }
    svector_iterator(void) {}
    svector_iterator(const base_it_type &it) : base_it_type(it) {}
  };

  template<class T> struct svector_const_iterator
    : public dal::dynamic_tree_sorted<_elt_svector<T>,
      comp_elt_svector<T>,3>::const_iterator {
    typedef typename dal::dynamic_tree_sorted<_elt_svector<T>,
      comp_elt_svector<T>,3>::const_iterator base_it_type;
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    reference operator *() const { return (base_it_type::operator*()).e; }
    pointer operator->() const { return &(operator*()); }
    size_type index(void) const { return (base_it_type::operator*()).c; }
    svector_const_iterator(void) {}
    svector_const_iterator(const svector_iterator<T> &it) : base_it_type(it) {}
    svector_const_iterator(const base_it_type &it) : base_it_type(it) {}
  };


  template<class T> class svector
    : public dal::dynamic_tree_sorted<_elt_svector<T>, comp_elt_svector<T>,3>
  {
  public:
    
    typedef dal::dynamic_tree_sorted<_elt_svector<T>, comp_elt_svector<T>,3>
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
    
    inline ref_elt_svector<T> operator [](size_type c)
    { return ref_elt_svector<T>(this, c); }

    inline void w(size_type c, const T &e) { 
      _elt_svector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = add_norepeat(ev);
      _base_type::operator[](i).e = e;
      if (e == T(0)) { sup(i); _base_type::operator[](i).c = size_type(-1); } 
    }

    inline T r(size_type c) const {  
      _elt_svector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = search(ev);
      if (i == size_type(-1)) return T(0);
      return (_base_type::operator[](i)).e;
    }

    inline T &ref(size_type c) {  
      _elt_svector<T> ev(c);
#ifdef __GETFEM_VERIFY
      if (c >= nbl) out_of_range_error();
#endif
      size_type i = add_norepeat(ev);
      return (_base_type::operator[](i)).e;
    }

    inline T operator [](size_type c) const { return r(c); }
    
    bool stored(size_type c) {
      _elt_svector<T> ev(c); 
      size_type i=search(ev); 
      return (i!=size_type(-1));
    }
    size_type nb_stored(void) { return card(); }
      
    /* Operations algebriques sur les vecteurs. */

    size_type nbline(void) const { return nbl; }
    size_type size(void) const { return nbl; } 
    void addmul(T, const svector<T>&);
    void fill(T);

    vsvector<T> full() const { 
      vsvector<T> v(size()); v.fill(0);
      const_tas_iterator it = tas_begin(), end = tas_end();
      for ( ; it != end; ++it) v[(*it).c] = (*it).e;
      return v;
    }
    
    /* Constructeurs */

    svector(size_type l){ init(l); }
    svector(void) { init(1); }

  };

  template<class T>  void svector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }



  /************************  Operations arithmetiques  **********************/

  template<class T>  void svector<T>::addmul(T a, const svector<T> &v) {
    if(v.size() != size()) DAL_THROW(dimension_error, "dimensions mismatch");
    const_tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*this)[(*it).c] += a * (*it).e;
  }

  template<class T>  void svector<T>::fill(T xx) {
    if (xx == T(0))
      elt.clear();
    else
      for (size_type i = 0; i < nbl; i++) (*this)[i] = xx;
  }

  template<class T>  svector<T>& operator *=(svector<T> &v, T xx) {
    typename svector<T>::tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*it).e *= xx;
    return v;
  }

  template<class T>  svector<T>& operator /=(svector<T> &v, T xx) {
    typename svector<T>::tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*it).e /= xx;
    return v;
  }

  template<class T>  svector<T>& operator +=(svector<T> &v, 
					     const svector<T>& w) {
    if (v.size()!= w.size()) DAL_THROW(dimension_error, "dimensions mismatch");

    typename svector<T>::const_tas_iterator it = w.tas_begin(),
      end = w.tas_end();
    for ( ; it != end; ++it) v[(*it).c] += (*it).e;
    return v;
  }

  template<class T>  svector<T>& operator -=(svector<T> &v,
					     const svector<T>& w) {
    if (v.size()!= w.size()) DAL_THROW(dimension_error, "dimensions mismatch");

    typename svector<T>::const_tas_iterator it = w.tas_begin(),
      end = w.tas_end();
    for ( ; it != end; ++it) v[(*it).c] -= (*it).e;
    return v;                      
  }


  template<class T>  bool operator ==(const svector<T> &v,
				      const svector<T> &w) {
    if (v.size() != w.size()) return false;
    typename svector<T>::const_tas_iterator it = w.tas_begin(),
      end = w.tas_end();
    for ( ; it != end; ++it) if (v[(*it).c] != (*it).e) return false; 
    it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) if (w[(*it).c] != (*it).e) return false;
    return true;
  }

  template<class T>  bool operator !=(const svector<T> &v, const svector<T> &w)
  { return ( !(v == w)); }

  template<class T> inline svector<T> operator *(const svector<T>& m, T x)
  { svector<T> p = m; p *= x; return p; }

  template<class T> inline svector<T> operator *(T x, const svector<T>& m)
  { svector<T> p = m; p *= x; return p; }

  template<class T> inline svector<T> operator /(const svector<T>& m, T x)
  { svector<T> p = m; p /= x; return p; }

  template<class T> inline svector<T> operator +(const svector<T>& m,
						 const svector<T>& n)
  { svector<T> p = m; p += n; return p; }

  template<class T> inline svector<T> operator -(const svector<T>& m,
						 const svector<T>& n)
  { svector<T> p = m; p -= n; return p; }

  template<class T> inline svector<T> operator -(const svector<T>& m)
  { svector<T> p = m; p *= -1; return p; }

  template<class T> inline svector<T> operator +(const svector<T>& p)
  { return p; }

  template<class T> T vect_norm2(const svector<T>&v) {
    /* pas bon sur les complexes.                                           */
    register T res = 0;
    typename svector<T>::const_tas_iterator it = v.tas_begin(),
      end = v.tas_end();
    for ( ; it != end; ++it) res += dal::sqr((*it).e);
    return sqrt(res);
  }

  template<class T, class VECT> T vect_sp(const svector<T>&v, const VECT &w) {
    /* pas bon sur les complexes.                                           */
    register T res = 0;
    typename svector<T>::const_tas_iterator it = v.tas_begin(),
      end = v.tas_end();
    for ( ; it != end; ++it) res += w[(*it).c] * (*it).e;
    return res;
  }

  /*********  intermediary structure for r/w operation.	*******************/
  
  /* revoir cette stratégie pour l'interfacage générique de abstrct_linalg */

  template<class T> class ref_elt_svector { /* ref. on an matrix element.  */

    svector<T> *pm;
    size_type l;
    
    public :

      operator T() const { return pm->r(l); }
    ref_elt_svector(svector<T> *p, size_type ll)
      { pm = p; l = ll; }
    inline ref_elt_svector operator =(T v)
      { (*pm).w(l,v); return *this; }
    inline bool operator ==(T v) const
      { return ((*pm).r(l,v) == v); }
    inline bool operator !=(T v) const
      { return ((*pm).r(l,v) != v); }
    inline ref_elt_svector operator +=(T v)
      { (*pm).w(l,(*pm).r(l) + v); return *this; }
    inline ref_elt_svector operator -=(T v)
      { (*pm).w(l,(*pm).r(l) - v); return *this; }
    inline ref_elt_svector operator /=(T v)
      { (*pm).w(l,(*pm).r(l) / v); return *this; }
    inline ref_elt_svector operator *=(T v)
      { (*pm).w(l,(*pm).r(l) * v); return *this; }
    inline ref_elt_svector operator =(const ref_elt_svector &re)
      { *this = T(re); return *this; }
  };  
  
  template<class T> T operator +(const ref_elt_svector<T> &re)
    { return T(re); }
  template<class T> T operator -(const ref_elt_svector<T> &re)
    { return -T(re); }
  template<class T> T operator +(const ref_elt_svector<T> &re, T v)
    { return T(re)+ v; }
  template<class T> T operator +(T v, const ref_elt_svector<T> &re)
    { return v+ T(re); }
  template<class T> T operator -(const ref_elt_svector<T> &re, T v)
    { return T(re)- v; }
  template<class T> T operator -(T v, const ref_elt_svector<T> &re)
    { return v- T(re); }
  template<class T> T operator *(const ref_elt_svector<T> &re, T v)
    { return T(re)* v; }
  template<class T> T operator *(T v, const ref_elt_svector<T> &re)
    { return v* T(re); }
  template<class T> T operator /(const ref_elt_svector<T> &re, T v)
    { return T(re)/ v; }
  template<class T> T operator /(T v, const ref_elt_svector<T> &re)
    { return v/ T(re); }
  

 
}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_SVECTOR_H */
