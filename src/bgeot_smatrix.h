/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  bgeot_smatrix.h : a example of non optimized sparse matrix   */
/*                              type.                                      */
/*            Please use a specialized matrix modulus instead.             */
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

#include <dal_tree_sorted.h>
#include <bgeot_matrix.h>

namespace bgeot
{

  template<class T> class ref_elt_smatrix;
  
  template<class T> class smatrix  /* General Sparse Matrix.		*/
  {
  public :

    typedef size_t size_type;
    typedef T value_type;
    typedef vsvector<T> vector_type;

    
    struct elt_m	       	/* basic element of the matrix.		*/
    {
      T val;	       	        /* value of the element.	        */
      size_type rang;	        /* column number.			*/
      elt_m(T v, size_type r) { val = v; rang = r; }
      elt_m(size_type r) { rang = r; }
      elt_m(void) {}
    };
    
    struct _elt_m_comp
    {
      int operator() (const elt_m &m, const elt_m &n) const
	{ return (int(m.rang) - int(n.rang)); }
    };
    
    typedef dal::dynamic_tree_sorted<elt_m, _elt_m_comp, 3> _line_m;

  protected:
    
    size_type nbc, nbl;      /* Number of columns and lines.           	*/
    std::vector<_line_m> li; /* array of lines.                         */
    
    /* Initialisation and copy.                                   	*/
    
    void init(size_type l, size_type c) { nbl = l; nbc = c; li.resize(l); }
    
  public:
    
    /* read and write operations.                                       */
    
    inline _line_m & operator[](size_type i) { return li[i]; }
    inline const _line_m & operator[](size_type i) const { return li[i]; }

    inline ref_elt_smatrix<T> operator ()(size_type l, size_type c)
    { return ref_elt_smatrix<T>(this, l, c); }
    
    inline void w(size_type l, size_type c, T e) { 
      elt_m el(e,c); _line_m *pl = &(li[l]); size_type num = pl->search(el);
      if (num != size_type(-1)) (*pl)[num] = el; else pl->add(el);
    }
    
    inline T r(size_type l, size_type c) const { 
      elt_m el(c); size_type num; const _line_m *pl = &(li[l]);
      if ((num = pl->search(el)) != size_type(-1)) return (*pl)[num].val;
      else return T(0);
    }

    inline ref_elt_smatrix<T> operator ()(size_type l, size_type c) const
    { r(l,c); }


    /* Manipulation of lines and columns functions.               	*/

    inline void clear_line(size_type i) { li[i].clear(); }

    /* Operations algebriques sur les matrices. */

    template<class R> friend R lc_product(const smatrix<R>& m1, size_type i,
			const smatrix<R>& m2, size_type j);
    template<class R> friend R partial_lc_product(const smatrix<R>& m1, size_type i,
				 const smatrix<R>& m2, size_type j, size_type k);
    template<class R> friend R ll_product(const smatrix<R>& m1, size_type i,
			 const smatrix<R>& m2, size_type j);
    template<class R> friend R partial_ll_product(const smatrix<R>& m1, size_type i,
				 const smatrix<R>& m2, size_type j, size_type k);
    template<class R> friend R lv_product(const smatrix<R> &m1, size_type i,
				   const typename smatrix<R>::vector_type &v);

    void l_mul(size_type l, T x);

    inline size_type nrows(void) const
      { return nbl; }
    inline size_type ncols(void) const
      { return nbc; }
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
    template<class R> friend smatrix<R> operator +(const smatrix<R> &m);
    template<class R> friend smatrix<R> operator -(const smatrix<R> &m);
    template<class R> friend smatrix<R> operator +(const smatrix<R> &m,
						   const smatrix<R> &n);
    template<class R> friend smatrix<R> operator -(const smatrix<R> &m,
						   const smatrix<R> &n);
    template<class R> friend smatrix<R> operator *(const smatrix<R> &m,
						   const smatrix<R> &n);
    template<class R> friend smatrix<R> operator *(const smatrix<R> &m, R x);
    template<class R> friend smatrix<R> operator *(R x, const smatrix<T> &m);
    template<class R> friend smatrix<R> operator /(const smatrix<R>& m,
						   const smatrix<R>& n);
    template<class R> friend smatrix<R> operator /(const smatrix<R>& m, R x);
    template<class R> friend smatrix<R> operator /(R x, const smatrix<R>& m);

    template<class R> friend vsvector<R>& operator *=(vsvector<R>& v, 
						     const smatrix<R>& m);

    template<class R> friend vsvector<R> operator *(const smatrix<R>& m,
						   vsvector<R>& v);
   
    /* Constructors / Destructor.                                         */
    smatrix(size_type l, size_type c = -6) { if (c == -6) c = l; init(l, c); }
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

  template<class T> vsvector<T>& operator *=(vsvector<T>& v,const smatrix<T>& m)
  {	/* a optimiser.							*/
    if ( (v.size() != m.nbl) || (v.size() != m.nbc) )
      throw dimension_error
	("vsvector<T>& smatrix<T>::operator *= dimensions mismatch");
    vsvector<T> tmp = v;
    for (size_type l = 0; l < m.nbl; l++) v(l) = lv_product(m, l, tmp);
    return v;
  }
  
  template<class T> vsvector<T> operator *(const smatrix<T>& m, vsvector<T>& v)
  {	/* optimiser l' acces a la matrice.				*/
    if (v.size() != m.ncols())
      throw dimension_error
	("vsvector<T>& smatrix<T>::operator *= dimensions mismatch");
    vsvector<T> res(m.nrows());
    for (size_type l = 0; l < m.nbl; l++) res[l] = lv_product(m, l, v);
    return res;
  }
  
  template<class T> T lc_product(const smatrix<T>& m1, size_type i,
				 const smatrix<T>& m2, size_type j) {
    smatrix<T>::_line_m *t = &(m1.li[i]); smatrix<T>::elt_m *e;
    T res = T(0);  
    for (size_type k = t->card()-1; k != size_type(-1); k--)
      { e = &((*t)[k]); res += e->val*m2.r(e->rang, j); }
    return res;
  }
  
  template<class T> T partial_lc_product(const smatrix<T>& m1, size_type i,
					 const smatrix<T>& m2, size_type j, size_type k)
  { /* Somme pour l <= k de m1(i,l) * m2(l,j).		*/
    smatrix<T>::_line_m *t = &(m1.li[i]); smatrix<T>::elt_m *e;
    T res = 0;
    for (size_type l = t->card()-1; l != size_type(-1); l--)
      { e =&((*t)[l]); if (e->rang <= k) res += e->val * m2.r(e->rang, j);}
    return res;
  }
  
  template<class T> T ll_product(const smatrix<T>& m1, size_type i,
				 const smatrix<T>& m2, size_type j)
  { /* Produit de la ligne i de m1 par la ligne j de m2.	*/
    smatrix<T>::_line_m *t = &(m1.li[i]); smatrix<T>::elt_m *e;
    T res = 0;
    for (size_type k = t->card()-1; k != size_type(-1); k--)
      {  e = &((*t)[k]); res += e->val * m2.r(j, e->rang); }
    return res;
  }
  
  template<class T> T partial_ll_product(const smatrix<T>& m1, size_type i,
					 const smatrix<T>& m2, size_type j, size_type k)
  { /* Somme pour l <= k de m1(i,l) * m2(j,l).		*/
    const typename smatrix<T>::_line_m *t = &(m1.li[i]);
    const typename smatrix<T>::elt_m *e;
    T res = 0;
    for (size_type l = t->card()-1; l != size_type(-1); l--)
      { e =&((*t)[l]); if (e->rang <= k) res += e->val * m2.r(j, e->rang);}
    return res;
  }
  
  template<class T> T lv_product(const smatrix<T> &m1, size_type i,
				 const typename smatrix<T>::vector_type &v)
  { 
    const typename smatrix<T>::_line_m *t = &(m1.li[i]); T res = T(0);  
    for (size_type k = t->card()-1; k != size_type(-1); k--)
    {
      const typename smatrix<T>::elt_m *e = &((*t)[k]);
      res += e->val * v[e->rang];
    }
    return res;
  }
  
  template<class T> void smatrix<T>::l_mul(size_type l, T x)
  { for (size_type j = li[l].card() - 1; j != size_type(-1); j--) li[l][j].val *= x; }
  
 
  
  template<class T> smatrix<T>& smatrix<T>::operator *=(T x)
  {
    for (size_type l = 0; l < nbl; l++) l_mul(l, x);
    return *this;
  }

  template<class T> smatrix<T>& smatrix<T>::operator /=(T x)
  {
    for (size_type l = 0; l < nbl; l++)
      for (size_type j = li[l].card() - 1; j != size_type(-1); j--) li[l][j].val /= x;
    return *this;
  }

  template<class T> smatrix<T>& smatrix<T>::operator *=(const smatrix<T> &m)
  {			/* a optimiser.					*/
    size_type l, c, i;	/* traiter le cas ou les deux matrices sont egales. */
    smatrix<T>::_line_m t;
    
    if ( (nbc != m.nbl) || (nbc != m.nbc) )
      throw dimension_error
	("smatrix<T>& smatrix<T>::operator *= dimensions mismatch");
      
    
    for (l = 0; l < nbl; l++)
    {
      t = li[l];
      li[l].clear();
      
      for (c = 0; c < nbc; c++)
	{
	  T aux = 0.0;
	  for (i = 0; i < t.card(); i++) aux += t[i].val*m.r(t[i].rang, c);
	  this->w(l,c,aux);
	}
    }
    return *this;
  }
  
  template<class T> smatrix<T>& smatrix<T>::operator +=(const smatrix<T> &m)
  { /* a optimiser. */
    if ( (nbc != m.nbc) || (nbl != m.nbl) )
      throw dimension_error
	("smatrix<T>& smatrix<T>::operator += dimensions mismatch");
      
    
    for (size_type l = 0; l < m.nbl; l++)
      for (size_type j = m.li[l].card() - 1; j != size_type(-1); j--)
      {
	elt_m *d = &(m.li[l][j]);
	this->w(l, d->rang, this->r(l, d->rang) + d->val);
      }
    return *this;
  }

  template<class T> smatrix<T>& smatrix<T>::operator -=(const smatrix<T> &m) {
    if ( (nbc != m.nbc) || (nbl != m.nbl) )
      throw dimension_error
	("smatrix<T>& smatrix<T>::operator -= dimensions mismatch");
     
    
    for (size_type l = 0; l < m.nbl; l++)
      for (size_type j = m.li[l].card() - 1; j != size_type(-1); j--)
	{
	  elt_m *d = &(m.li[l][j]);
	  this->w(l, d->rang, this->r(l, d->rang) - d->val);
	} 
    return *this;
  }

  template<class T> bool smatrix<T>::operator ==(const smatrix<T> &m) const { 
    size_type l, j; elt_m *d;
    
    if ( ( nbc != m.nbc ) || ( nbl != m.nbl ) ) return false;
    
    for (l = 0; l < nbl; l++)
      {
	for (j = li[l].card() - 1; j != size_type(-1); j--)
	  { d = &(li[l][j]); if ( d->val != m.r(l,d->rang)) return false; }
	for (j = m.li[l].card() - 1; j != size_type(-1); j--)
	  { d = &(m.li[l][j]); if (d->val != r(l,d->rang)) return false; }
      }
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
    smatrix<T> p(m.nbl, n.nbc);
    
    if (m.nbc != n.nbl)
      throw dimension_error
	("smatrix<T> smatrix<T>::operator * dimensions mismatch");
      
    
    for (size_type l = 0; l < m.nbl; l++)
      for (size_type c = 0; c < n.nbc; c++)
	p.w(l,c, lc_product(m, l, n, c) );
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

  template<class T> class ref_elt_smatrix /* ref. on an matrix element.  */
  { 
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
      { *this = T(re); }
    template<class R> friend R operator +(const ref_elt_smatrix<R> &re);
    template<class R> friend R operator -(const ref_elt_smatrix<R> &re);
    template<class R> friend R operator +(const ref_elt_smatrix<R> &re, R v);
    template<class R> friend R operator +(R v, const ref_elt_smatrix<R> &re);
    template<class R> friend R operator -(const ref_elt_smatrix<R> &re, R v);
    template<class R> friend R operator -(R v, const ref_elt_smatrix<R> &re);
    template<class R> friend R operator *(const ref_elt_smatrix<R> &re, R v);
    template<class R> friend R operator *(R v, const ref_elt_smatrix<R> &re);
    template<class R> friend R operator /(const ref_elt_smatrix<R> &re, R v);
    template<class R> friend R operator /(R v, const ref_elt_smatrix<R> &re);
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
  




  template < class Matrix, class Vector>
    int cg(const Matrix& A, Vector& x, const Vector& b, int itemax, 
	   double residu)
  {
    typedef typename Vector::value_type value_type;
    value_type rho(0), rho_1(0), alpha(0), beta(0);
    Vector p(x.size()), q(x.size()), r(x.size());
    int iter = 0;
    r = b - A * x;
    rho = vect_sp(r, r);

    while (::sqrt(rho) > residu) {

      if (iter == 0)
	p = r;		  
      else {
	beta = rho / rho_1;
	p = r + p * beta;
      }
      
      q = A * p;
      alpha = rho / vect_sp(p, q);
      
      x += alpha * p;
      r -= alpha * q;
      
      rho_1 = rho;
      
      ++iter;
      rho = vect_sp(r, r);
      cout << "iter " << iter << " residu " << ::sqrt(rho) << endl;
      if (iter >= itemax) return 1;
    }
    
    return 0;
  }

}

#endif /* __BGEOT_SMATRIX_H */
