/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  bgeot_smatrix.h : a sparse matrix type.                      */
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


template<class T> class ref_elt_matrix_gsp;

template<class T> class matrix_gsp  /* General Sparse Matrix.		*/
{
  protected:

    typedef short int IND;

    struct elt_m	       	/* basic element of the matrix.		*/
    {
          T val;	       	/* value of the element.	        */
      int rang;	        	/* column number.			*/
      elt_m(T v, int r) { val = v; rang = r; }
      elt_m(int r) { rang = r; }
      elt_m(void) {}
    };

    struct _elt_m_comp
    {
      static int compare(const matrix_gsp<T>::elt_m &m,
			 const matrix_gsp<T>::elt_m &n)
      { return (m.rang - n.rang); }
    };

    /* le tri peut se faire par insertion, le nombre d'elements est petit. */
    typedef dal_tree_sorted<elt_m, IND, _elt_m_comp> _line_m;

    int nbc, nbl;    	/* Number of columns and lines.           	*/
    _line_m *li;        /* array of lines.                              */

    /* Initialisation and copy.                                   	*/

    friend void __init_creuse(matrix_gsp<T> &m, int l, int c);
    friend void __copy_creuse(matrix_gsp<T> &d, const matrix_gsp<T> &s);

  public:

    typedef T base_type;
    typedef vectorp<T> vector_type;

    /* read and write operations.                                       */

    inline ref_elt_matrix_gsp<T> operator ()(int l, int c)
      { return ref_elt_matrix_gsp<T>(this, l, c); }
    
    inline void w(int l, int c, T e)
      { 
	elt_m el(e,c); _line_m *pl = &(li[l]); IND num = pl->search(el);
	// if (e == T(0)) { if (num != IND(-1)) pl->sup(num); }
	// else
	  if (num != IND(-1)) (*pl)[num] = el; else pl->add(el);
      } /* suppression pose probleme ... a voir ... */
   
    inline T r(int l, int c) const
      { 
	elt_m el(c); int num; _line_m *pl = &(li[l]);
	if ((num = pl->search(el)) != -1) return (*pl)[num].val;
	                             else return T(0);
      }

    /* Manipulation of lines and columns functions.               	*/

    inline void clear_line(int i)
      { li[i].clear(); }

    void line_exchange(int i, int j);
    void column_exchange(int i, int j);
    void transpose( void );
    void del_line(int l);
    void del_column(int co);

    /* Operations algebriques sur les matrices. */

    friend T lc_product(const matrix_gsp<T>& m1, int i,
			const matrix_gsp<T>& m2, int j);
    friend  T partial_lc_product(const matrix_gsp<T>& m1, int i,
				 const matrix_gsp<T>& m2, int j, int k);
    friend  T ll_product(const matrix_gsp<T>& m1, int i,
			 const matrix_gsp<T>& m2, int j);
    friend  T partial_ll_product(const matrix_gsp<T>& m1, int i,
				 const matrix_gsp<T>& m2, int j, int k);
    friend T lv_product(const matrix_gsp<T> &m1, int i,
			const matrix_gsp<T>::vector_type &v);

    void l_mul(int l, T x);

    matrix_gsp<T>& operator =(const matrix_gsp<T> &m);
    matrix_gsp<T>& operator =(T v);
    inline int nbline(void) const
      { return nbl; }
    inline int nbcolumn(void) const
      { return nbc; }
    matrix_gsp<T>& operator *=(T x);
    matrix_gsp<T>& operator /=(T x);
    inline void inv(void)
      { mat_gauss_inverse(*this); }
    inline matrix_gsp<T>& operator /=(const matrix_gsp<T>&m)
      { matrix_gsp<T> n = m; n.inv(); return ((*this) *= n); }
    matrix_gsp<T>& operator *=(const matrix_gsp<T> &m);
    matrix_gsp<T>& operator +=(const matrix_gsp<T> &m);
    matrix_gsp<T>& operator -=(const matrix_gsp<T> &m);
    bool operator ==(const matrix_gsp<T> &m) const;
    inline bool operator !=(const matrix_gsp<T> &m) const
      { return ( !(m == *this) ); }
    friend matrix_gsp<T> operator +(const matrix_gsp<T> &m);
    friend matrix_gsp<T> operator -(const matrix_gsp<T> &m);
    friend matrix_gsp<T> operator +(const matrix_gsp<T> &m,
				    const matrix_gsp<T> &n);
    friend matrix_gsp<T> operator -(const matrix_gsp<T> &m,
				    const matrix_gsp<T> &n);
    friend matrix_gsp<T> operator *(const matrix_gsp<T> &m,
				    const matrix_gsp<T> &n);
    friend matrix_gsp<T> operator *(const matrix_gsp<T> &m, T x);
    friend matrix_gsp<T> operator *(T x, const matrix_gsp<T> &m);
    friend matrix_gsp<T> operator /(const matrix_gsp<T>& m,
				    const matrix_gsp<T>& n);
    friend matrix_gsp<T> operator /(const matrix_gsp<T>& m, T x);
    friend matrix_gsp<T> operator /(T x, const matrix_gsp<T>& m);

    friend vectorp<T>& operator *=(vectorp<T>& v, const matrix_gsp<T>& m);

    friend vectorp<T> operator *(const matrix_gsp<T>& m, vectorp<T>& v);
   
    /* Constructors / Destructor.                                         */
    matrix_gsp(const matrix_gsp<T>&m)
      { __init_creuse(*this, m.nbl, m.nbc); __copy_creuse(*this, m); }
    matrix_gsp(int l, int c = -6)
      { if (c == -6) c = l; __init_creuse(*this, l, c); }
    matrix_gsp(void)
      { __init_creuse(*this, 1, 1); }
    ~matrix_gsp(void)
      { delete[] li; li = NULL; }

    /* Fonctions diverses. */
    friend ostream& operator <<(ostream& o, const matrix_gsp<T>& m)
      { 
	mat_out(o,m);
	int i,j;
	for (i = 0,  j = 0 ; i < m.nbline(); i++) j += m.li[i].nb_elt();
	cout << "total : " << j << " elements\n";
	return o;
      }
};

template<class T> void __init_creuse(matrix_gsp<T> &m, int l, int c)
{
  m.nbl = l; m.nbc = c;
  m.li = new matrix_gsp<T>::_line_m[l](3,3,true,true);
}

template<class T> void __copy_creuse(matrix_gsp<T> &d, const matrix_gsp<T> &s)
{

  if (d.nbl != s.nbl) { delete[] d.li; __init_creuse(d, s.nbl, s.nbc); }


  d.nbc = s.nbc; for (int i = 0; i < d.nbl; i++) d.li[i] = s.li[i];

}

template<class T>
  matrix_gsp<T>& matrix_gsp<T>::operator =(const matrix_gsp<T> &m)
{ 
  __copy_creuse(*this, m);
  return *this;
}

template<class T> void matrix_gsp<T>::line_exchange(int i, int j)
{
  char *p = (char *)(&(li[j])), *q = (char *)(&(li[i]));
  for (int i = 0; i < sizeof(*li); i++, p++, q++) exx(*p, *q);
}

template<class T> void matrix_gsp<T>::transpose( void )
{ /* a optimiser pour matrices pleines non carrees : inutile de copier */
  _line_m *t = new _line_m[nbc](3,3,true,true);
  
  for (int l = 0; l < nbl; l++)
    for (int j = li[l].nb_elt() - 1; j >= 0; j--)
      {
	matrix_gsp<T>::elt_m *d = &(li[l][j]);
	t[ d->rang ].add( matrix_gsp<T>::elt_m(d->val, l) );
      }
  exx(nbl, nbc); delete[] li; li = t;
}

template<class T> void matrix_gsp<T>::del_line(int l)
{ nbl--; for (int i=l; i < nbl-1; i++) li[i]=li[i+1]; li[nbl].clear(); }

template<class T> void matrix_gsp<T>::del_column(int co)
{
  for (int l = 0; l < nbl; l++)
    {
      int num = -1;
      for (int j = li[l].nb_elt() - 1; j >= 0; j--)
	{
	  int r = (li[l][j]).rang;
	  if (r > co) (li[l][j]).rang--; else if (r == co) num = j;
	}
      if (num != -1) li[l].sup(num);
    }
  nbc--;
} 

template<class T> void matrix_gsp<T>::column_exchange(int i, int j)
/* hum !! */
{ for (int l=0; l<nbl;l++) { T a=r(l,i); w(l,i,this->r(l,j)); w(l,j,a);}}

template<class T> vectorp<T>& operator *=(vectorp<T>& v,const matrix_gsp<T>& m)
{	/* a optimiser.							*/
  if ( (v.size() != m.nbl) || (v.size() != m.nbc) )
    err("MATRIX: *= with incompatible dimensions.");
  vectorp<T> tmp = v;
  for (int l = 0; l < m.nbl; l++) v(l) = lv_product(m, l, tmp);
  return v;
}

template<class T> vectorp<T> operator *(const matrix_gsp<T>& m, vectorp<T>& v)
{	/* optimiser l' acces a la matrice.				*/
  if (v.size() != m.nbcolumn())
    err("MATRIX: * with incompatible dimensions.");
  vectorp<T> res(m.nbline());
  for (int l = 0; l < m.nbl; l++) res(l) = lv_product(m, l, v);
  return res;
}

template<class T> T lc_product(const matrix_gsp<T>& m1, int i,
			const matrix_gsp<T>& m2, int j)
{
  matrix_gsp<T>::_line_m *t = &(m1.li[i]); matrix_gsp<T>::elt_m *e;
  T res = T(0);  
  for (int k = t->nb_elt()-1; k >= 0; k--)
    { e = &((*t)[k]); res += e->val*m2.r(e->rang, j); }
  return res;
}

template<class T> T partial_lc_product(const matrix_gsp<T>& m1, int i,
				       const matrix_gsp<T>& m2, int j, int k)
{ /* Somme pour l <= k de m1(i,l) * m2(l,j).		*/
  matrix_gsp<T>::_line_m *t = &(m1.li[i]); matrix_gsp<T>::elt_m *e;
  T res = 0;
  for (int l = t->nb_elt()-1; l >= 0; l--)
    { e =&((*t)[l]); if (e->rang <= k) res += e->val * m2.r(e->rang, j);}
  return res;
}

template<class T> T ll_product(const matrix_gsp<T>& m1, int i,
			       const matrix_gsp<T>& m2, int j)
{ /* Produit de la ligne i de m1 par la ligne j de m2.	*/
  matrix_gsp<T>::_line_m *t = &(m1.li[i]); matrix_gsp<T>::elt_m *e;
  T res = 0;
  for (int k = t->nb_elt()-1; k >= 0; k--)
    {  e = &((*t)[k]); res += e->val * m2.r(j, e->rang); }
  return res;
}

template<class T> T partial_ll_product(const matrix_gsp<T>& m1, int i,
				       const matrix_gsp<T>& m2, int j, int k)
{ /* Somme pour l <= k de m1(i,l) * m2(j,l).		*/
  matrix_gsp<T>::_line_m *t = &(m1.li[i]); matrix_gsp<T>::elt_m *e;
  T res = 0;
  for (int l = t->nb_elt()-1; l >= 0; l--)
    { e =&((*t)[l]); if (e->rang <= k) res += e->val * m2.r(j, e->rang);}
  return res;
}

template<class T> T lv_product(const matrix_gsp<T> &m1, int i,
	     const matrix_gsp<T>::vector_type &v)
{ 
  matrix_gsp<T>::_line_m *t = &(m1.li[i]); T res = T(0);  
  for (int k = t->nb_elt()-1; k >= 0; k--)
    { matrix_gsp<T>::elt_m *e = &((*t)[k]); res += e->val * v(e->rang); }
  return res;
}

template<class T> void matrix_gsp<T>::l_mul(int l, T x)
{ for (int j = li[l].nb_elt() - 1; j >= 0; j--) li[l][j].val *= x; }

template<class T> matrix_gsp<T>& matrix_gsp<T>::operator =(T v)
{ 
  int m = min(nbl, nbc);
  for(int i = 0; i < nbl; i++)
    { li[i].clear(); if (i < m) this->w(i,i,v); }
  return *this;
}

template<class T> matrix_gsp<T>& matrix_gsp<T>::operator *=(T x)
{
  for (int l = 0; l < nbl; l++) l_mul(l, x);
  return *this;
}
template<class T> matrix_gsp<T>& matrix_gsp<T>::operator /=(T x)
{
  for (int l = 0; l < nbl; l++)
    for (int j = li[l].nb_elt() - 1; j >= 0; j--) li[l][j].val /= x;
  return *this;
}
template<class T> matrix_gsp<T>& matrix_gsp<T>::operator *=(const matrix_gsp<T> &m)
{			/* a optimiser.					*/
  int l, c, i;	/* traiter le cas ou les deux matrices sont egales. */
  matrix_gsp<T>::_line_m t;
  
  if ( (nbc != m.nbl) || (nbc != m.nbc) )
    err("MATRIX: *= with incompatible dimensions.");
  
  for (l = 0; l < nbl; l++)
    {
      t = li[l];
      li[l].clear();
      
      for (c = 0; c < nbc; c++)
	{
	  T aux = 0.0;
	  for (i = 0; i < t.nb_elt(); i++) aux += t[i].val*m.r(t[i].rang, c);
	  this->w(l,c,aux);
	}
    }
  return *this;
}

template<class T> matrix_gsp<T>& matrix_gsp<T>::operator +=(const matrix_gsp<T> &m)
{ /* a optimiser. */
  if ( (nbc != m.nbc) || (nbl != m.nbl) )
    err("MATRIX: operation += with incompatible dimensions.");
  
  for (int l = 0; l < m.nbl; l++)
    for (int j = m.li[l].nb_elt() - 1; j >= 0; j--)
      {
	elt_m *d = &(m.li[l][j]);
	this->w(l, d->rang, this->r(l, d->rang) + d->val);
      }
  return *this;
}
template<class T> matrix_gsp<T>& matrix_gsp<T>::operator -=(const matrix_gsp<T> &m)
{
  if ( (nbc != m.nbc) || (nbl != m.nbl) )
    err("MATRIX: operation -= with incompatibles dimensions.");
  
  for (int l = 0; l < m.nbl; l++)
    for (int j = m.li[l].nb_elt() - 1; j >= 0; j--)
      {
	elt_m *d = &(m.li[l][j]);
	this->w(l, d->rang, this->r(l, d->rang) - d->val);
      } 
  return *this;
}
template<class T> bool matrix_gsp<T>::operator ==(const matrix_gsp<T> &m) const
{ 
  int l, j; elt_m *d;
  
  if ( ( nbc != m.nbc ) || ( nbl != m.nbl ) ) return false;
  
  for (l = 0; l < nbl; l++)
    {
      for (j = li[l].nb_elt() - 1; j >= 0; j--)
	{ d = &(li[l][j]); if ( d->val != m.r(l,d->rang)) return false; }
      for (j = m.li[l].nb_elt() - 1; j >= 0; j--)
	{ d = &(m.li[l][j]); if (d->val != r(l,d->rang)) return false; }
    }
  return true;
}

template<class T> matrix_gsp<T> operator +(const matrix_gsp<T> &m)
{ return m; }
template<class T> matrix_gsp<T> operator -(const matrix_gsp<T> &m)
{ matrix_gsp<T> p = m; p *= -1; return p; }
template<class T> matrix_gsp<T> operator +(const matrix_gsp<T> &m,
						  const matrix_gsp<T> &n)
{ matrix_gsp<T> p = m; p += n; return p; }
template<class T> matrix_gsp<T> operator -(const matrix_gsp<T> &m,
						  const matrix_gsp<T> &n)
{ matrix_gsp<T> p = m; p -= n; return p; }
template<class T> matrix_gsp<T> operator *(const matrix_gsp<T> &m,
						  const matrix_gsp<T> &n)
{
  matrix_gsp<T> p(m.nbl, n.nbc);
  
  if (m.nbc != n.nbl)
    err("MATRIX : Multiplication with incompatible dimensions.");
  
  for (int l = 0; l < m.nbl; l++)
    for (int c = 0; c < n.nbc; c++)
      p.w(l,c, lc_product(m, l, n, c) );
  return p;
}
template<class T> matrix_gsp<T> operator *(const matrix_gsp<T> &m, T x)
{ matrix_gsp<T> p = m; p *= x; return p; }
template<class T> matrix_gsp<T> operator *(T x, const matrix_gsp<T> &m)
{ matrix_gsp<T> p = m; p *= x; return p; }
template<class T> matrix_gsp<T> operator /(const matrix_gsp<T>& m,
						      const matrix_gsp<T>& n)
{ matrix_gsp<T> p = n; p.inv(); return m * p; }
template<class T> matrix_gsp<T> operator /(const matrix_gsp<T>& m, T x)
{ matrix_gsp<T> p = m; p /= x; return p; }
template<class T> matrix_gsp<T> operator /(T x, const matrix_gsp<T>& m)
{ matrix_gsp<T> p = m; p.inv(); return p * x; }


/*********  intermediary structure for r/w operation.	*******************/

template<class T> class ref_elt_matrix_gsp /* ref. on an matrix element.  */
{ 
    matrix_gsp<T> *pm;
    int l,c;

  public :

    operator T() const { return (*pm).r(l,c); }
    ref_elt_matrix_gsp(matrix_gsp<T> *p, int ll, int cc)
      { pm = p; l = ll; c = cc; }
    inline ref_elt_matrix_gsp operator =(T v)
      { (*pm).w(l,c,v); return *this; }
    inline ref_elt_matrix_gsp operator +=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) + v); return *this; }
    inline ref_elt_matrix_gsp operator -=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) - v); return *this; }
    inline ref_elt_matrix_gsp operator /=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) / v); return *this; }
    inline ref_elt_matrix_gsp operator *=(T v)
      { (*pm).w(l,c,(*pm).r(l,c) * v); return *this; }
    inline ref_elt_matrix_gsp operator =(const ref_elt_matrix_gsp &re)
      { *this = T(re); }
    friend T operator +(const ref_elt_matrix_gsp<T> &re);
    friend T operator -(const ref_elt_matrix_gsp<T> &re);
    friend T operator +(const ref_elt_matrix_gsp<T> &re, T v);
    friend T operator +(T v, const ref_elt_matrix_gsp<T> &re);
    friend T operator -(const ref_elt_matrix_gsp<T> &re, T v);
    friend T operator -(T v, const ref_elt_matrix_gsp<T> &re);
    friend T operator *(const ref_elt_matrix_gsp<T> &re, T v);
    friend T operator *(T v, const ref_elt_matrix_gsp<T> &re);
    friend T operator /(const ref_elt_matrix_gsp<T> &re, T v);
    friend T operator /(T v, const ref_elt_matrix_gsp<T> &re);
};  

template<class T> T operator +(const ref_elt_matrix_gsp<T> &re)
{ return T(re); }
template<class T> T operator -(const ref_elt_matrix_gsp<T> &re)
{ return -T(re); }
template<class T> T operator +(const ref_elt_matrix_gsp<T> &re, T v)
{ return T(re)+ v; }
template<class T> T operator +(T v, const ref_elt_matrix_gsp<T> &re)
{ return v+ T(re); }
template<class T> T operator -(const ref_elt_matrix_gsp<T> &re, T v)
{ return T(re)- v; }
template<class T> T operator -(T v, const ref_elt_matrix_gsp<T> &re)
{ return v- T(re); }
template<class T> T operator *(const ref_elt_matrix_gsp<T> &re, T v)
{ return T(re)* v; }
template<class T> T operator *(T v, const ref_elt_matrix_gsp<T> &re)
{ return v* T(re); }
template<class T> T operator /(const ref_elt_matrix_gsp<T> &re, T v)
{ return T(re)/ v; }
template<class T> T operator /(T v, const ref_elt_matrix_gsp<T> &re)
{ return v/ T(re); }

#endif /* __MATRIX_GSP_H */
