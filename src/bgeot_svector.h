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
  /* classe a tester ... */

  /************************************************************************/
  /*		Class svector: sparse vector.				  */
  /************************************************************************/

  template<class T> struct _elt_svector
  {
    size_type c; T e;
    _elt_svector(void) { c = size_type(-1); e = T(0); }
    _elt_svector(size_type cc) { c = cc; e = T(0); }
    _elt_svector(size_type cc, const T &ee) { c = cc; e = ee; }
  };

  template<class T> struct comp_elt_svector
    : public std::binary_function<_elt_svector<T>, _elt_svector<T>, int>
  {
    int operator()(const _elt_svector<T>& m, const _elt_svector<T>& n) const
    { int d = m.c-n.c; if (d<0) return -1; if (d>0) return 1; return 0; }
  };
  

  template<class T> class svector
    : public dal::dynamic_tree_sorted<_elt_svector<T>, comp_elt_svector<T>,3>
  {
  public:

    typedef dal::dynamic_tree_sorted<_elt_svector<T>, comp_elt_svector<T>,3> _base_type;
    typedef typename _base_type::tas_iterator tas_iterator;
    typedef typename _base_type::const_tas_iterator const_tas_iterator;
    typedef typename _base_type::size_type size_type;

  protected:
    size_type nbl;    	/* Nombre d'elements max.	        	  */

    void __init_vp(size_type l) { nbl = l; clear(); }

  public:

    void clean(double eps)
    {
      tas_iterator it = tas_begin();
      tas_iterator end = tas_end();
      for ( ; it != end; ++it) if (dal::abs((*it).e) <= eps) sup(it.index());
    }

    void resize(size_type l) {
      nbl = l;
    }
      
    typedef T value_type;

    /* operations de lecture et ecriture. pas la bonne methode,         */
    /* faire une structure intermediaire.                               */
    inline T& operator [](size_type c)
    { 
      _elt_svector<T> ev(c);
#ifdef __GETFEM_VERIFY
      assert(c < nbl);
#endif
      size_type i = add_norepeat(ev, false);
      return (_base_type::operator[](i)).e;
    }
    inline const T& operator [](size_type c) const
    {  
      _elt_svector<T> ev(c);
      static T zero = T(0);
#ifdef __GETFEM_VERIFY
      assert(c < nbl);
#endif
      size_type i = search(ev);
      if (i == size_type(-1)) return zero;
      return (_base_type::operator[](i)).e;
    }
     
    bool stored(size_type c)
    { _elt_svector<T> ev(c); size_type i=search(ev); return (i!=size_type(-1)); }
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

    svector(size_type l);
    svector(void) { __init_vp(1); };

  };


  /************************  Constructeurs.	*************************/

  template<class T>  svector<T>::svector(size_type l)
  {
#ifdef __GETFEM_VERIFY
    assert(l > 0);
#endif
    __init_vp(l);
  }

  /************************  Operations arithmetiques  **********************/

  template<class T>  void svector<T>::addmul(T a, const svector<T> &v)
  {
#ifdef __GETFEM_VERIFY
    assert(v.size() == this->size());
#endif
    const_tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*this)[(*it).c] += a * (*it).e;
  }

  template<class T>  void svector<T>::fill(T xx)
  {
    if (xx == T(0))
      elt.clear();
    else
      for (size_type i = 0; i < nbl; i++) (*this)[i] = xx;
  }

  template<class T>  svector<T>& operator *=(svector<T> &v, T xx)
  {
    typename svector<T>::tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*it).e *= xx;
    return v;
  }

  template<class T>  svector<T>& operator /=(svector<T> &v, T xx)
  {
    typename svector<T>::tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) (*it).e /= xx;
    return v;
  }

  template<class T>  svector<T>& operator +=(svector<T> &v, const svector<T>& w)
  {
#ifdef __GETFEM_VERIFY
    assert(v.size() == w.size());
#endif

    typename svector<T>::const_tas_iterator it = w.tas_begin(), end = w.tas_end();
    for ( ; it != end; ++it) v[(*it).c] += (*it).e;
    return v;
  }

  template<class T>  svector<T>& operator -=(svector<T> &v, const svector<T>& w)
  {
#ifdef __GETFEM_VERIFY
    assert(v.size() == w.size());
#endif

    typename svector<T>::const_tas_iterator it = w.tas_begin(), end = w.tas_end();
    for ( ; it != end; ++it) v[(*it).c] -= (*it).e;
    return v;                      
  }


  template<class T>  bool operator ==(const svector<T> &v, const svector<T> &w)
  {
    if (v.size() != w.size()) return false;
    typename svector<T>::const_tas_iterator it = w.tas_begin(), end = w.tas_end();
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

  template<class T> T vect_norm2(const svector<T>&v)
  { /* pas bon sur les complexes.                                           */
    register T res = 0;
    typename svector<T>::const_tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) res += dal::sqr((*it).e);
    return sqrt(res);
  }

  template<class T, class VECT> T vect_sp(const svector<T>&v, const VECT &w)
  { /* pas bon sur les complexes.                                           */
    register T res = 0;
    typename svector<T>::const_tas_iterator it = v.tas_begin(), end = v.tas_end();
    for ( ; it != end; ++it) res += w[(*it).c] * (*it).e;
    return res;
  }
 
}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_SVECTOR_H */
