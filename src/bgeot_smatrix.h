/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  bgeot_smatrix.h : a example of non optimized sparse matrix   */
/*                              type.                                      */
/*            Please use a specialized matrix package instead.             */
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

#ifndef __BGEOT_SMATRIX_H
#define __BGEOT_SMATRIX_H

#include <bgeot_svector.h>

namespace bgeot
{

  template<class T> class ref_elt_smatrix;
  
  template<class T> class smatrix  /* General Sparse Matrix.		   */
  {
  public :

    typedef size_t size_type;
    typedef T value_type;
    typedef vsvector<T> vector_type;
    typedef svector<T> _line_type;

    std::vector<_line_type> li; /* array of lines.                         */
    
    /* Initialisation and copy.                                   	   */
    
    void init(size_type l, size_type c)
    { li.resize(l); std::fill(li.begin(), li.end(), svector<T>(c)); }
    
    void resize(size_type l, size_type c)
    { li.resize(l); for (size_type i=0; i < l; ++i) li[i].resize(c); }

    /* read and write operations.                                          */
    void out_of_range_error(void) const;
    inline ref_elt_smatrix<T> operator ()(size_type l, size_type c)
    { return ref_elt_smatrix<T>(this, l, c); }
    
    inline void w(size_type l, size_type c, const T &e) {
      #ifdef __GETFEM_VERIFY
        if (l >= li.size()) out_of_range_error();
      #endif
      li[l].w(c,e);
    }
    inline T r(size_type l, size_type c) const {
      #ifdef __GETFEM_VERIFY
        if (l >= li.size()) out_of_range_error();
      #endif
      return li[l].r(c);
    }
    inline T operator ()(size_type l, size_type c) const { return r(l,c);}

    /* Manipulation of lines and columns functions.               	*/

    inline void clear_line(size_type i) { li[i].clear(); }
    void clear() { for (size_type i=0; i < nrows(); ++i) clear_line(i); }

    inline _line_type& row(size_type i) { return li[i]; }
    inline const _line_type& row(size_type i) const { return li[i]; }

    /* Operations algebriques sur les matrices. */

    void l_mul(size_type l, T x) { li[l] *= x; }

    inline size_type nrows(void) const { return li.size(); }
    inline size_type ncols(void) const
    { return (nrows() == 0) ? 0 : li[0].size(); }
    smatrix<T>& operator *=(T x);
    smatrix<T>& operator /=(T x);
    inline void inv(void)
      { mat_gauss_inverse(*this); }
    inline smatrix<T>& operator /=(const smatrix<T>&m)
      { smatrix<T> n = m; n.inv(); return ((*this) *= n); }
    smatrix<T>& operator *=(const smatrix<T> &m);
    smatrix<T>& operator +=(const smatrix<T> &m);
    smatrix<T>& operator -=(const smatrix<T> &m);
    bool operator ==(const smatrix<T> &m) const;
    inline bool operator !=(const smatrix<T> &m) const
      { return ( !(m == *this) ); }

    /* Constructors / Destructor.                                         */
    smatrix(size_type l, size_type c = size_type(-6))
      { if (c == size_type(-6)) c = l; init(l, c); }
    smatrix(void) { init(1, 1); }

    /* Fonctions diverses. */
    template<class R> friend std::ostream& operator << (std::ostream& o,
							const smatrix<R>& m) { 
      mat_out(o,m);
      size_type i,j;
      for (i = 0,  j = 0 ; i < m.nrows(); i++) j += m.li[i].card();
      cout << "total : " << j << " elements\n";
      return o;
    }
  };

  template<class T>  void smatrix<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }


  template<class T> vsvector<T>& operator *=(vsvector<T>& v,
					     const smatrix<T>& m) {
    /* a optimiser.							*/
    if ( (v.size() != m.nrows()) || (v.size() != m.ncols()) )
      DAL_THROW(dimension_error, "dimensions mismatch");
    vsvector<T> tmp = v;
    for (size_type l = 0; l < m.nrows(); l++) v(l) = lv_product(m, l, tmp);
    return v;
  }
  template<class R, class VECT> VECT operator *(const smatrix<R>& m,
						const VECT& v) {
    /* optimiser l' acces a la matrice.				*/
    if (v.size() != m.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");
    VECT res(m.nrows());
    for (size_type l = 0; l < m.nrows(); l++) res[l] = lv_product(m, l, v);
    return res;
  }
  
  template<class T> T lc_product(const smatrix<T>& m1, size_type i,
				 const smatrix<T>& m2, size_type j) {
    typename smatrix<T>::_line_type::iterator it = m1.row(i).begin(),
      ite = m1.row(i).end();
    T res = T(0);
    for (; it != ite; ++it) res += (it->e) * m2.r(it->c, j);
    return res;
  }
  
  template<class T> T partial_lc_product(const smatrix<T>& m1, size_type i,
					 const smatrix<T>& m2, size_type j,
					 size_type k) {
    typename smatrix<T>::_line_type::iterator it = m1.row(i).begin(),
      ite = m1.row(i).end();
    T res = T(0);
    for (; it != ite; ++it) if (it->c <= k) res += (it->e) * m2.r(it->c, j);
    return res;
  }
  
  template<class T> T ll_product(const smatrix<T>& m1, size_type i,
				 const smatrix<T>& m2, size_type j)
  { return vect_sp(m1.row(i), m2.row(j)); }
  
  
  template<class T> T partial_ll_product(const smatrix<T>& m1, size_type i,
					 const smatrix<T>& m2, size_type j,
					 size_type k) {
    typename smatrix<T>::_line_type::iterator it = m1.row(i).begin(),
      ite = m1.row(i).end();
    T res = T(0);
    for (; it != ite; ++it) if (it->c <= k) res += (it->e) * m2(j, it->c);
    return res;
  }
  
  template<class T, class VECT> inline T lv_product(const smatrix<T> &m1,
						    size_type i, const VECT& v)
  { return vect_sp(m1.row(i), v); }
  
  template<class T> smatrix<T>& smatrix<T>::operator *=(T x) {
    for (size_type l = 0; l < nrows(); l++) l_mul(l, x);
    return *this;
  }

  template<class T> smatrix<T>& smatrix<T>::operator /=(T x)
  { (*this) *= 1 / x; return *this; }

  template<class T> smatrix<T>& smatrix<T>::operator *=(const smatrix<T> &m) {
   			/* a optimiser.					*/
    size_type l, c, i;	/* traiter le cas ou les deux matrices sont egales. */
    typename smatrix<T>::_line_type t;
    
    if ( (ncols() != m.nrows()) || (ncols() != m.ncols()) )
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    for (l = 0; l < nrows(); l++) {
      t = li[l];
      li[l].clear();
      
      for (c = 0; c < ncols(); c++) {
	T aux = 0.0;
	typename smatrix<T>::_line_type::iterator it = t.begin(),
	  ite = t.end();
	for ( ; it != ite; ++it) aux += (it->e) * m.r(t[i].c, c);
	this->w(l,c,aux);
      }
    }
    return *this;
  }
  
  template<class T> smatrix<T>& smatrix<T>::operator +=(const smatrix<T> &m) {
    if ( (ncols() != m.ncols()) || (nrows() != m.nrows()) )
      DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type l = 0; l < m.nrows(); l++) row(i) += m.row(i);
    return *this;
  }
  
  template<class T> smatrix<T>& smatrix<T>::operator -=(const smatrix<T> &m) {
    if ( (ncols() != m.ncols()) || (nrows() != m.nrows()) )
      DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type l = 0; l < m.nrows(); l++) row(i) -= m.row(i);
    return *this;
  }
  
  template<class T> bool smatrix<T>::operator ==(const smatrix<T> &m) const {
    if ( (ncols() != m.ncols()) || (nrows() != m.nrows()) ) return false;
    for (size_type l = 0; l < m.nrows(); l++)
      if (row(i) != m.row(i)) return false;
    return true;
  }
  
  template<class T> smatrix<T> operator +(const smatrix<T> &m)
  { return m; }
  template<class T> smatrix<T> operator -(const smatrix<T> &m)
  { smatrix<T> p = m; p *= -1.0; return p; }
  template<class T> smatrix<T> operator +(const smatrix<T> &m,
					  const smatrix<T> &n)
  { smatrix<T> p = m; p += n; return p; }
  template<class T> smatrix<T> operator -(const smatrix<T> &m,
					  const smatrix<T> &n)
  { smatrix<T> p = m; p -= n; return p; }
  template<class T> smatrix<T> operator *(const smatrix<T> &m,
					  const smatrix<T> &n) {
    smatrix<T> p(m.nrows(), n.ncols());
    if (m.ncols() != n.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");
      
    for (size_type l = 0; l < m.nrows(); l++)
      for (size_type c = 0; c < n.ncols(); c++) {
	T e = lc_product(m, l, n, c);
	if (e != T(0)) p.w(l,c, e); 
      }
    return p;
  }
  template<class T> smatrix<T> operator *(const smatrix<T> &m, T x)
  { smatrix<T> p = m; p *= x; return p; }
  template<class T> smatrix<T> operator *(T x, const smatrix<T> &m)
  { smatrix<T> p = m; p *= x; return p; }
  template<class T> smatrix<T> operator /(const smatrix<T>& m,
					  const smatrix<T>& n)
  { smatrix<T> p = n; p.inv(); return m * p; }
  template<class T> smatrix<T> operator /(const smatrix<T>& m, T x)
  { smatrix<T> p = m; p /= x; return p; }
  template<class T> smatrix<T> operator /(T x, const smatrix<T>& m)
  { smatrix<T> p = m; p.inv(); return p * x; }
  

  /*********  intermediary structure for r/w operation.	*******************/
  
  template<class T> class ref_elt_smatrix { /* ref. on an matrix element.  */

    smatrix<T> *pm;
    size_type l,c;
    
    public :
      
      operator T() const { return (*pm).r(l,c); }
    ref_elt_smatrix(smatrix<T> *p, size_type ll, size_type cc)
      { pm = p; l = ll; c = cc; }
    inline ref_elt_smatrix operator =(T v)
      { (*pm).w(l,c,v); return *this; }
    inline ref_elt_smatrix operator +=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) + v); return *this; }
    inline ref_elt_smatrix operator -=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) - v); return *this; }
    inline ref_elt_smatrix operator /=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) / v); return *this; }
    inline ref_elt_smatrix operator *=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) * v); return *this; }
    inline ref_elt_smatrix operator =(const ref_elt_smatrix &re)
      { *this = T(re); return *this; }
  };  
  
  template<class T> T operator +(const ref_elt_smatrix<T> &re)
    { return T(re); }
  template<class T> T operator -(const ref_elt_smatrix<T> &re)
    { return -T(re); }
  template<class T> T operator +(const ref_elt_smatrix<T> &re, T v)
    { return T(re)+ v; }
  template<class T> T operator +(T v, const ref_elt_smatrix<T> &re)
    { return v+ T(re); }
  template<class T> T operator -(const ref_elt_smatrix<T> &re, T v)
    { return T(re)- v; }
  template<class T> T operator -(T v, const ref_elt_smatrix<T> &re)
    { return v- T(re); }
  template<class T> T operator *(const ref_elt_smatrix<T> &re, T v)
    { return T(re)* v; }
  template<class T> T operator *(T v, const ref_elt_smatrix<T> &re)
    { return v* T(re); }
  template<class T> T operator /(const ref_elt_smatrix<T> &re, T v)
    { return T(re)/ v; }
  template<class T> T operator /(T v, const ref_elt_smatrix<T> &re)
    { return v/ T(re); }
  

  /*********  generic solver.	*******************/


  template < class Matrix, class Vector>
    int cg(const Matrix& A, Vector& x, const Vector& b, int itemax, 
	   double residu, bool noisy = true)
  {
    typedef typename Vector::value_type value_type;
    value_type rho(0), rho_1(0), alpha(0), beta(0);
    Vector p(x.size()), q(x.size()), r(x.size());
    int iter = 0;
    r = b - A * x;
    rho = vect_sp(r, r);

    while (::sqrt(rho) > residu) {
      if (iter == 0) p = r;		  
      else { beta = rho / rho_1; p = r + p * beta; }
      
      q = A * p;
      alpha = rho / vect_sp(p, q);
      
      x += alpha * p;
      r -= alpha * q;
      
      rho_1 = rho;
      
      ++iter;
      rho = vect_sp(r, r);
      if (noisy) cout << "iter " << iter << " residu " << ::sqrt(rho) << endl;
      if (iter >= itemax) return 1;
    }
    return 0;
  }

}

#endif /* __BGEOT_SMATRIX_H */
