/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_matrix.h : small matrices.                             */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
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


#ifndef __BGEOT_MATRIX_H
#define __BGEOT_MATRIX_H

#include <bgeot_vector.h>

namespace bgeot
{

  /* ******************************************************************** */
  /*		Class fsmatrix<T,n>: fixed size matrix.	        	  */
  /* ******************************************************************** */

  template<class T, int N> class fsmatrix : public fsvector<T, N*N>
  {
    public :

      typedef T value_type;
      typedef typename fsvector<T, N*N>::size_type size_type;
      typedef typename fsvector<T, N*N>::iterator iterator;
      typedef typename fsvector<T, N*N>::const_iterator const_iterator;

      inline const T& operator ()(size_type l, size_type c) const
      {
        #ifdef __GETFEM_VERIFY
	if (l >= N || c >= N) out_of_range_error();
        #endif
	return *(begin() + c*N+l);
      }
      inline T& operator ()(size_type l, size_type c)
      {
        #ifdef __GETFEM_VERIFY
	if (l >= N || c >= N) out_of_range_error();
        #endif
        return *(begin() + c*N+l);
      }

      void clear_line(size_type i);
      void line_exchange(size_type, size_type);
      void transpose(void);
      void fill(T a, T b = T(0));
      inline size_type nrows(void) const { return N; }
      inline size_type ncols(void) const { return N; }
      void l_mul(size_type, const T&);
      fsmatrix &operator *=(const fsmatrix<T,N>&);
      fsmatrix &operator *=(const T &a)
      { fsvector<T, N*N>::operator *=(a); return *this; }
     
      fsmatrix(size_type l, size_type c) {
	if ((N!=l) || (N!=c)) DAL_THROW(std::invalid_argument,"bad arguments");
      }
      fsmatrix(void) {}
  };

  template<class T, int N> void fsmatrix<T,N>::fill(T a, T b)
  { 
    fsvector<T,N*N>::fill(b);
    iterator p = begin(), e = end();
    while (p < e) { *p = a; p += N+1; }
  }

  template<class T, int N> void fsmatrix<T,N>::clear_line(size_type i)
  { iterator p = begin()+i, e = end(); while (p < e) { *p = T(0); p += N; } }

  template<class T, int N>
    void fsmatrix<T,N>::line_exchange(size_type i, size_type j)
  {
    iterator p1 = begin(), e = end(), p2 = p1; p1 += i; p2 += j;
    while (p1 < e) { std::swap(*p1, *p2); p1 += N; p2 += N; }
  }

  template<class T, int N> void fsmatrix<T,N>::transpose(void)
  {
    for (size_type l = 0; l < N; l++)
    {
      iterator p1 = begin()+l+l*N, e = end(), p2 = p1; p1 += N; p2++;
      while (p1 < e) { std::swap(*p1, *p2); p1 += N; p2++; }
    }
  }

  template<class T, int N> void fsmatrix<T,N>::l_mul(size_type i, const T &x)
  { iterator p = begin()+i, e = end(); while (p1 < e) { *p *= x; p += N; } }

  template<class T, int N>
    fsmatrix<T,N>& fsmatrix<T,N>::operator *=(const fsmatrix<T,N>& m)
  {
    if (begin() != m.begin())
    {
      fsvector<T,N> tmp;
      for (size_type l = 0; l < N; l++)
      {
	iterator p1 = begin()+l, e = end(), p2 = tmp.begin(), te = tmp.end();
	const_iterator p3 = p1;
	while (p3 < e) { *p2++ = *p3; p3 += N; }
	for (size_type c = 0; c < N; c++)
	{
	  *p1 = T(0);
	  p2 = tmp.begin(); p3 = m.begin() + N * c;
	  while (p2 < te) *p1 += (*p2++) * (*p3++);
	  p1 += N;
	}
      }
    }
    else
    {
      fsmatrix<T,N> q = *this;
      for (size_type l = 0; l < N; l++)
      {
	iterator p1 = begin() + l, e = q.end(), p2, p3;
	for (c = 0; c < N; c++)
	{
	  *p1 = T(0);
	  p2 = q.begin()+l; p3 = q.begin() + N * c;
	  while (p2 < e) { *p1 += (*p2) * (*p3++); p2 += N; }
	  p1 += N;
	}
      }

    }
    return *this;
  }

  template<class T, int N>
    fsvector<T,N>& operator *=(fsvector<T,N>& v, const fsmatrix<T,N>& m)
  {
    fsvector<T,N> tmp = v;
    typename fsvector<T,N>::iterator p1 = v.begin(), te = tmp.end(), p2;
    typename fsmatrix<T,N>::const_iterator p3;
    for (typename vsvector<T>::size_type l = 0; l < N; l++)
    {
      *p1 = T(0); p2 = tmp.begin(); p3 = m.begin() + l;
      while (p2 < te) { *p1 += (*p2++) * (*p3); p3 += N; }
      p1++;
    }
    return v;
  }

  template<class T, int N>
    fsvector<T,N> operator *(const fsmatrix<T,N>& m, const fsvector<T,N>& v)
  {
    fsvector<T,N> res;
    typename fsvector<T,N>::iterator p1 = res.begin();
    typename fsvector<T,N>::const_iterator te = v.end(), p2;
    typename fsmatrix<T,N>::const_iterator p3;

    for (typename fsvector<T,N>::size_type l = 0; l < N; l++, p1++)
    {
      *p1 = T(0); p2 = v.begin(); p3 = m.begin() + l;
      while (p2 < te) { *p1 += (*p2++) * (*p3); p3 += m.nrows(); }
    }
  
    return res;
  }

  template<class T, int N>
    fsmatrix<T,N> operator *(const fsmatrix<T,N>& m, const fsmatrix<T,N>& n)
  {
    fsmatrix<T,N> res;
    typename fsmatrix<T,N>::iterator p1 = res.begin();
    typename fsmatrix<T,N>::const_iterator p2, p3, e;
    
    for (typename fsmatrix<T,N>::size_type c = 0; c < N; c++)
      for (typename fsmatrix<T,N>::size_type l = 0; l < N; l++, p1++)
      {
	*p1 = T(0);
	p2 = m.begin() + l; e = m.end(); p3 = n.begin() + c * n.nrows();
	while (p2 < e) { *p1 += (*p2) * (*p3++); p2 += m.nrows(); }
      }
    
    return res;
  }

  /* ******************************************************************** */
  /*		Class vsmatrix<T>: variable size matrix.     		  */
  /* ******************************************************************** */

  template<class T> class vsmatrix : public vsvector<T>
  {
    public:

      typedef typename vsvector<T>::size_type size_type;
      typedef typename vsvector<T>::iterator iterator;
      typedef typename vsvector<T>::const_iterator const_iterator;

    protected:
      size_type nbc, nbl;

    public:

      typedef vsvector<T> vector_type;

      inline const T& operator ()(size_type l, size_type c) const
      {
        #ifdef __GETFEM_VERIFY
	if (l >= nbl || c >= nbc) out_of_range_error();
        #endif
	return *(begin() + c*nbl+l);
      }
      inline T& operator ()(size_type l, size_type c)
      {
        #ifdef __GETFEM_VERIFY
	if (l >= nbl || c >= nbc) out_of_range_error();
        #endif
	return *(begin() + c*nbl+l);
      }

      void resize(size_type l, size_type c)
      { if (c * l != size()) vsvector<T>::resize(c*l); nbl = l; nbc = c; }
      void clear_line(size_type i); 
      void line_exchange(size_type, size_type);
      void transpose(void);
      void fill(T a, T b = T(0));
      inline size_type nrows(void) const { return nbl; }
      inline size_type ncols(void) const { return nbc; }
      void l_mul(size_type, const T&);
      vsmatrix<T>& operator *=(const vsmatrix<T>&);

      /// Add matrix w to current matrix.
      vsmatrix<T>& operator +=(const vsmatrix<T>& w)
      { this->vsvector<T>::operator +=(w); return *this; }
      /// Substract matrix w to current matrix v.
      vsmatrix<T>& operator -=(const vsmatrix<T>& w)
      { this->vsvector<T>::operator -=(w); return *this; }
      /// Multiply the current matrix with the scalar x.
      vsmatrix<T>& operator *=(const T &x)
      { this->vsvector<T>::operator *=(x); return *this; }
      /// Divide the current matrix with the scalar x.
      vsmatrix<T>& operator /=(const T &x)
      { this->vsvector<T>::operator /=(x); return *this; }
      

      vsmatrix(size_type l, size_type c) : vsvector<T>(c*l) { nbl=l; nbc=c; }
      vsmatrix(void) : vsvector<T>(1) { nbl = nbc = 1; }

  };

  template<class T> void vsmatrix<T>::fill(T a, T b)
  { 
    vsvector<T>::fill(b);
    iterator p = begin(), e = end();
    while (p < e) { *p = a; p += nbl+1; }
  }

  template<class T> void vsmatrix<T>::clear_line(size_type i)
  { iterator p = begin()+i, e = end(); while (p < e) { *p = T(0); p += nbl; } }

  template<class T> void vsmatrix<T>::line_exchange(size_type i, size_type j)
  {
    iterator p1 = begin(), e = end(), p2 = p1; p1 += i; p2 += j;
    while (p1 < e) { std::swap(*p1, *p2); p1 += nbl; p2 += nbl; }
  }

  template<class T> void vsmatrix<T>::transpose(void)
  { /* à optimiser. */
    if (nbl == nbc)
    {
      for (size_type l = 0; l < this->nbl; l++)
      {
	iterator p1 = begin()+l+l*nbl, e = end(), p2 = p1; p1 += nbl; p2++;
	while (p1 < e) { std::swap(*p1, *p2); p1 += nbl; p2++; }
      }
    }
    else
    {
      vsmatrix<T> m(nbc, nbl);
      for (size_type l = 0; l < nbl; l++)
	for (size_type c = 0; c < nbc; c++)
	  m(c,l) = (*this)(l,c);
      *this = m;
    }
  }

  template<class T> void vsmatrix<T>::l_mul(size_type i, const T &x)
  { iterator p = begin()+i, e = end(); while (p1 < e) { *p *= x; p += nbl; } }

  template<class T>  vsmatrix<T>& vsmatrix<T>::operator *=(const vsmatrix<T>& m)
  {
    if (nbc != m.nbl || nbc != m.nbc)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    if (begin() != m.begin())
    {
      vsvector<T> tmp(nbc);
      for (size_type l = 0; l < nbl; l++)
      {
	iterator p1 = begin()+l, e = end(), p2 = tmp.begin(), te = tmp.end();
	const_iterator p3 = p1;
	while (p3 < e) { *p2++ = *p3; p3 += nbl; }
	for (size_type c = 0; c < nbc; c++)
	{
	  *p1 = T(0);
	  p2 = tmp.begin(); p3 = m.begin() + nbc * c;
	  while (p2 < te) *p1 += (*p2++) * (*p3++);
	  p1 += nbl;
	}
      }
    }
    else
    {
      vsmatrix<T> q = *this;
      for (size_type l = 0; l < nbl; l++)
      {
	iterator p1 = begin() + l, e = q.end(), p2, p3;
	for (size_type c = 0; c < nbc; c++)
	{
	  *p1 = T(0);
	  p2 = q.begin()+l; p3 = q.begin() + c * nbl;
	  while (p2 < e) { *p1 += (*p2) * (*p3++); p2 += nbl; }
	  p1 += nbl;
	}
      }

    }
    return *this;
  }


  template<class T>  vsvector<T>& operator *=(vsvector<T>& v, const vsmatrix<T>& m)
  {
    if (v.size() != m.nrows() || v.size() != m.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");
  
    vsvector<T> tmp = v;
    typename vsvector<T>::iterator p1 = v.begin(), te = tmp.end(), p2;
    typename vsmatrix<T>::const_iterator p3;
    for (typename vsvector<T>::size_type l = 0; l < m.nrows(); l++, p1++)
    {
      *p1 = T(0); p2 = tmp.begin(); p3 = m.begin() + l;
      while (p2 < te) { *p1 += (*p2++) * (*p3); p3 += m.nrows(); }
    }
  
    return v;
  }

  template<class T>
    vsvector<T> operator *(const vsmatrix<T>& m, const vsvector<T>& v)
  {

    if (v.size() != m.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");

    vsvector<T> res(m.nrows());
    typename vsvector<T>::iterator p1 = res.begin();
    typename vsvector<T>::const_iterator te = v.end(), p2;
    typename vsmatrix<T>::const_iterator p3;

    for (typename vsvector<T>::size_type l = 0; l < m.nrows(); l++, p1++)
    {
      *p1 = T(0); p2 = v.begin(); p3 = m.begin() + l;
      while (p2 < te) { *p1 += (*p2++) * (*p3); p3 += m.nrows(); }
    }
  
    return res;
  }

  template<class T>
    vsmatrix<T> operator *(const vsmatrix<T>& m, const vsmatrix<T>& n)
  {

    if (m.ncols() != n.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");

    vsmatrix<T> res(m.nrows(), n.ncols());
    typename vsmatrix<T>::iterator p1 = res.begin();
    typename vsmatrix<T>::const_iterator p2, p3, e;
    
    for (typename vsmatrix<T>::size_type c = 0; c < p.ncols(); c++)
      for (typename vsmatrix<T>::size_type l = 0; l < p.nrows(); l++, p1++)
      {
	*p1 = T(0);
	p2 = m.begin() + l; e = m.end(); p3 = n.begin() + c * n.nrows();
	while (p2 < e) { *p1 += (*p2) * (*p3++); p2 += m.nrows(); }
      }
    
    return res;
  }


  /* ******************************************************************** */
  /*		Generic functions on matrix.            		  */
  /* ******************************************************************** */

  template<class MAT>
    typename MAT::value_type lc_product(const MAT& m1,
				     typename MAT::size_type i,
				     const MAT& m2, typename MAT::size_type j)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (m1.ncols() != m2.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");
    typename MAT::const_iterator p1 = m1.begin() + i, e = m1.end();
    typename MAT::const_iterator p2 = m2.begin() + j * m2.nrows();
    while (p1 < e) { res += (*p1) * (*p2++); p1 += m1.nrows(); }
    return res;
  }

  template<class MAT> 
    typename MAT::value_type partial_lc_product(const MAT& m1,
			 typename MAT::size_type i, const MAT& m2,
			 typename MAT::size_type j, typename MAT::size_type k)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
        
    if (m1.ncols() < k || m2.nrows() < k)
      DAL_THROW(dimension_error, "dimensions mismatch");

    typename MAT::const_iterator p1 = m1.begin() + i;
    typename MAT::const_iterator e  = m1.begin() + k * m1.nrows();
    typename MAT::const_iterator p2 = m2.begin() + j * m2.nrows();
    while (p1 < e) { res += (*p1) * (*p2++); p1 += m1.nrows(); }
    return res;
  }

  template<class MAT>
    typename MAT::value_type cc_product(const MAT& m1,
				     typename MAT::size_type i,
				     const MAT& m2, typename MAT::size_type j)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (m1.nrows() != m2.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");

    typename MAT::const_iterator p1 = m1.begin() + i * m1.nrows(), 
                                 e = p1 + m1.nrows();
    typename MAT::const_iterator p2 = m2.begin() + j * m2.nrows();
    while (p1 < e) res += (*p1++) * (*p2++);
    return res;
  }

  template<class MAT, class VEC>
    typename MAT::value_type lv_product(const MAT& m1,
				      typename MAT::size_type i, const VEC& v)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (m1.ncols() != v.size())
      DAL_THROW(dimension_error, "dimensions mismatch");

    typename MAT::const_iterator p1 = m1.begin() + i, e = m1.end();
    typename VEC::const_iterator p2 = v.begin();
    while (p1 < e) { res += (*p1) * (*p2++); p1 += m1.nrows(); }
    return res;
  }

  template<class MAT, class VEC> 
    typename MAT::value_type partial_lv_product(const MAT& m1,
			typename MAT::size_type i, const VEC& v,
                        typename MAT::size_type i1, typename MAT::size_type i2)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (!(m1.ncols() >= i2 && v.size() >= i2))
      DAL_THROW(dimension_error, "dimensions mismatch");
  
    typename MAT::const_iterator p1 = m1.begin() + i + i1 * m1.nrows();
    typename VEC::const_iterator p2 = v.begin() + i1, e = v.begin() + i2;
    while (p2 < e) { res += (*p1) * (*p2++); p1 += m1.nrows(); }
    return res;
  }

  template<class MAT, class VEC>
    typename MAT::value_type cv_product(const MAT& m1,
				       typename MAT::size_type i, const VEC& v)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (m1.nrows() != v.size())
      DAL_THROW(dimension_error, "dimensions mismatch");

    typename MAT::const_iterator p1 = m1.begin() + i * m1.nrows();
    typename VEC::const_iterator p2 = v.begin(), e = v.end();
    while (p2 < e) res += (*p1++) * (*p2++);
    return res;
  }
  
  template<class MAT, class VEC> 
    typename MAT::value_type partial_cv_product(const MAT& m1,
		       typename MAT::size_type i, const VEC& v,
                       typename MAT::size_type i1, typename MAT::size_type i2)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (!(m1.nrows() >= i2 && v.size() >= i2))
      DAL_THROW(dimension_error, "dimensions mismatch");
      
    typename MAT::const_iterator p1 = m1.begin() + i1 + i * m1.nrows();
    typename VEC::const_iterator p2 = v.begin() + i1, e = v.begin() + i2;
    while (p2 < e) res += (*p1++) * (*p2++);
    return res;
  }

  template<class MAT>
    typename MAT::value_type ll_product(const MAT& m1,
				     typename MAT::size_type i,
			       	     const MAT& m2, typename MAT::size_type j)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (m1.nrows() != m2.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");
   
    typename MAT::const_iterator p1 = m1.begin() + i, e = m1.end();
    typename MAT::const_iterator p2 = m2.begin() + j;
    while (p1 < e) { res += (*p1)*(*p2); p1 += m1.nrows(); p2 += m2.nrows(); }
    return res;
  }

  template<class MAT> 
    typename MAT::value_type partial_ll_product(const MAT& m1,
			 typename MAT::size_type i, const MAT& m2,
			 typename MAT::size_type j, typename MAT::size_type k)
  {
    typename MAT::value_type res = typename MAT::value_type(0);
    if (!(m1.nrows() >= k && m2.nrows() >= k))
      DAL_THROW(dimension_error, "dimensions mismatch");

    typename MAT::const_iterator p1 = m1.begin() + i;
    typename MAT::const_iterator e  = m1.begin() + k * m1.nrows();
    typename MAT::const_iterator p2 = m2.begin() + j;
    while (p1 < e) { res += (*p1)*(*p2); p1 += m1.nrows(); p2 += m2.nrows(); }
    return res;
  }
  
  template<class MAT1, class MAT2, class MATR>
    void mat_product(const MAT1 &m1, const MAT2 &m2, MATR &mr)
  { // mr = m1 * m2; optimisable.
   
    if (m1.ncols() != m2.nrows() || mr.nrows() != m1.nrows()
	  || mr.ncols() != m2.ncols() || ((void *)(&mr) == (void *)(&m1))
	  || ((void *)(&mr) == (void *)(&m2)))
      DAL_THROW(dimension_error, "dimensions mismatch");

    for (size_t i = 0; i < mr.nrows(); ++i)
      for (size_t j = 0; j < mr.ncols(); ++j)
	mr(i,j) = lc_product(m1, i, m2, j);
  }

  template<class MAT, class VECT, class VECTR>
    void mat_vect_product(const MAT &m, const VECT &v, VECTR &vr)
  { // vr = m * v; optimisable.
    if (m.ncols() != v.size() || vr.size() != m.nrows()
	|| ((void *)(&vr) == (void *)(&v)))
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    for (size_t i = 0; i < vr.size(); ++i) vr[i] = lv_product(m, i, v);
  }

  template<class MAT1, class MAT2, class MATR>
    void mat_product_tt(const MAT1 &m1, const MAT2 &m2, MATR &mr) 
  { // mr = transp(m1) * trans(m2); optimisable.
    if (m1.nrows() != m2.ncols() || mr.nrows() != m1.ncols()
	  || mr.ncols() != m2.nrows() || ((void *)(&mr) == (void *)(&m1))
	  || ((void *)(&mr) == (void *)(&m2)))
      DAL_THROW(dimension_error, "dimensions mismatch");

    for (size_t i = 0; i < mr.nrows(); ++i)
      for (size_t j = 0; j < mr.ncols(); ++j)
	mr(i,j) = lc_product(m2, j, m1, i);
  }

  template<class MAT, class VECT, class VECTR>
    void mat_vect_product_t(const MAT &m, const VECT &v, VECTR &vr)
  { // vr = transp(m) * v; optimisable.
    if (m.ncols() != vr.size() || v.size() != m.nrows()
	|| ((void *)(&vr) == (void *)(&v)))
      DAL_THROW(dimension_error, "dimensions mismatch");
  
    for (size_t i = 0; i < vr.size(); ++i) vr[i] = cv_product(m, i, v);
  }

  template<class MAT1, class MAT2, class MATR>
    void mat_product_tn(const MAT1 &m1, const MAT2 &m2, MATR &mr) 
  { // mr = trans(m1) * m2; optimisable.
    if (m1.nrows() != m2.nrows() || mr.nrows() != m1.ncols()
	  || mr.ncols() != m2.ncols() || ((void *)(&mr) == (void *)(&m1))
	  || ((void *)(&mr) == (void *)(&m2)))
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    for (size_t i = 0; i < mr.nrows(); ++i)
      for (size_t j = 0; j < mr.ncols(); ++j)
	mr(i,j) = cc_product(m1, i, m2, j);
  }

  template<class MAT1, class MAT2, class MATR>
    void mat_add_product_tn(const MAT1 &m1, const MAT2 &m2, MATR &mr) 
  { // mr += trans(m1) * m2; optimisable.
    if (m1.nrows() != m2.nrows() || mr.nrows() != m1.nrows()
	  || mr.ncols() != m2.ncols() || ((void *)(&mr) == (void *)(&m1))
	  || ((void *)(&mr) == (void *)(&m2)))
      DAL_THROW(dimension_error, "dimensions mismatch");

    for (size_t i = 0; i < mr.nrows(); ++i)
      for (size_t j = 0; j < mr.ncols(); ++j)
	mr(i,j) += cc_product(m1, i, m2, j);
  }

  template<class MAT> inline void mat_out(std::ostream &o, const MAT &m)
  {
    for (typename MAT::size_type l = 0; l < m.nrows(); l++ )
    {
      for (typename MAT::size_type c = 0; c < m.ncols(); c++ )
	o << m(l,c) << ' ';
      o << '\n';
    }
  }

  template<class T>  std::ostream &operator <<(std::ostream &o, const vsmatrix<T>& m)
  { mat_out(o,m); return o; }

  template<class T, int N>
    std::ostream &operator <<(std::ostream &o, const fsmatrix<T,N>& m)
  { mat_out(o,m); return o; }


  /* ******************************************************************** */
  /*	resolution d'un systeme lineaire par pivot de gauss.		  */
  /* ******************************************************************** */

  template<class MAT, class VECT>
    void mat_gauss_solve(const MAT &m, const VECT &b, VECT &x, MAT &tmp,
			 double EPS = 1E-12)
  {
    x = b; tmp = m;
    if ((m.nrows() != m.ncols()) || (m.ncols() != b.size())
	|| (m.ncols() != x.size()))
      DAL_THROW(dimension_error, "dimensions mismatch");
 
    typename MAT::value_type vlp;
    typename MAT::size_type l, c, pc, l2, nbl = m.nrows();
    
    /* Triangulation superieure du système.                               */
    for (l = 0; l < nbl; l++)
    {
      /* Choix du pivot.                                                  */  
      for (c = l+1, pc = l, vlp = dal::abs(tmp(pc,l)); c < nbl; c++)
	if ( dal::abs(tmp(c,l)) > vlp ) { pc = c; vlp = dal::abs(tmp(pc,l)); }
      
      /* Echange pour avoir le pivot en bonne place.                      */
      if (l != pc) { tmp.line_exchange(l, pc); std::swap(x[l], x[pc]); }
      
      /* triangulation.                                                   */
      vlp = tmp(l,l);

      if (dal::abs(vlp) < EPS)
	DAL_THROW(failure_error, "Non invertible matrix");

      for (c = l; c < nbl; c++) tmp(l, c) /= vlp;
      x[l] /= vlp;
      for (l2 = l + 1; l2 < nbl; l2++ )
      {
	vlp = tmp(l2,l);
	for (c = l + 1; c < nbl; c++) tmp(l2, c) -= tmp(l,c) * vlp;
	x[l2] -= x[l] * vlp;
      }
    }  
    
    /* Resolution.                                                        */
    for (c = nbl - 1; c != typename MAT::size_type(-1); c--)
      for (l = c - 1; l != typename MAT::size_type(-1); l--)
	x[l] -=  x[c] * tmp(l,c);
  }

  template<class MAT, class VECT> inline
    void mat_gauss_solve(const MAT &m, const VECT &b, VECT &x,
			 double EPS = 1E-12)
  { MAT tmp; mat_gauss_solve(m, b, x, tmp, EPS); }


  /* ******************************************************************** */
  /*	Inversion d'une matrice par pivot de gauss.			  */
  /* ******************************************************************** */

  template<class MAT>
    void mat_gauss_inverse(MAT &m, MAT &tmp, double EPS = 1E-12)
  {
    if (m.nrows() != m.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");
      
    typename MAT::value_type vlp;
    typename MAT::size_type l, c, pc, l2, nbl = m.nrows();
    tmp = m;
    m.fill(typename MAT::value_type(1), typename MAT::value_type(0));
   
    /* Triangulation Inférieure du système.               */
    for (l = 0; l < nbl; l++)
    {
      /* Choix du pivot.                                  */  
      for (c = l+1, pc = l; c < nbl; c++)
	if ( dal::abs(tmp(c,l)) > dal::abs(tmp(pc,l)) ) pc = c;
      
      /* Echange pour avoir le pivot en bonne place.      */
      if (l != pc) { tmp.line_exchange(l, pc); m.line_exchange(l, pc); }
      
      /* triangulation.                                   */
      vlp = tmp(l,l); tmp(l,l) = 1;
      if (dal::abs(vlp) < EPS)
	DAL_THROW(failure_error, "Non invertible matrix");
	
      for (c = l + 1; c < nbl; c++) tmp(l, c) /= vlp;
      for (c = 0; c < nbl; c++) m(l,c) /= vlp;
      for (l2 = l+1; l2 < nbl; l2++ )
      {
	vlp = tmp(l2,l);
	for (c = l + 1; c < nbl; c++) tmp(l2, c) -= tmp(l,c) * vlp;
	for (c = 0; c < nbl; c++) m(l2,c) -= m(l,c) * vlp;
      }
    }  
    
    /* Triangulation superieure du systeme.              */

    for (l = nbl - 1; l  != typename MAT::size_type(-1); l--)
      for (l2 = l - 1; l2  != typename MAT::size_type(-1); l2--)
      { vlp = tmp(l2,l); for (c = 0; c < nbl; c++) m(l2, c) -=  m(l,c)*vlp; }

  }

  template<class MAT> inline void mat_gauss_inverse(MAT &m, double EPS = 1E-12)
  { MAT tmp; mat_gauss_inverse(m, tmp, EPS); }

 
  /* ******************************************************************** */
  /*	Calcul du determinant par pivot de gauss			  */
  /* ******************************************************************** */

  template<class MAT> 
    typename MAT::value_type mat_gauss_det(const MAT &m, MAT &tmp,
					                double EPS = 1E-12)
  {
    typename MAT::size_type nbl = m.nrows(), nbc = m.ncols();
    if (nbl != nbc)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    tmp = m;
    typename MAT::size_type l, c, pc, l2;
    typename MAT::value_type vlp, vlpr, res = typename MAT::value_type(1);

    /* Triangulation Inferieure du systeme.                               */
    for (l = 0; l < nbl; l++)
    {
      /* Choix du pivot.                                                  */  
      for (c = l+1, pc = l; c < nbl; c++)
	if ( dal::abs(tmp(c,l)) > dal::abs(tmp(pc,l)) ) pc = c;
      
      /* Echange pour avoir le pivot en bonne place.                      */
      if (l != pc) tmp.line_exchange(l, pc);
      
      /* triangulation.                                                   */
      vlpr = tmp(l,l); res *= vlpr;
      if (dal::abs(vlpr) < EPS) return typename MAT::value_type(0);
      for (l2 = l+1; l2 < nbl; l2++ )
      {
	vlp = tmp(l2,l);
	for (c = l; c < nbl; c++) tmp(l2, c) = tmp(l2, c) - tmp(l,c)*vlp/vlpr;
      }
    }
    return res;
  }

  template<class MAT> inline typename MAT::value_type mat_gauss_det(MAT &m,
							   double EPS = 1E-12)
  { MAT tmp; return mat_gauss_det(m, tmp, EPS); }

  
  template<class MAT>
    typename MAT::value_type mat_gauss_sub_det(const MAT &m,
			 typename MAT::size_type i, typename MAT::size_type j,
	              			         MAT &tmp, double EPS = 1E-12)
  {
    typename MAT::size_type nbl = m.nrows(), nbc = m.ncols();
    if (!(nbl == nbc && i < nbl && j < nbl))
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    if (nbl == 1) return typename MAT::value_type(1);
    tmp = m; nbl--;
    typename MAT::size_type l, c, pc, l2;
    
    for (l = 0; l < i; l++)
      for (c = j; c < nbl; c++) tmp(l,c) = tmp(l, c+1);
    for (l = i; l < nbl; l++)
      for (c = 0; c < j; c++) tmp(l,c) = tmp(l+1, c);
    for (l = i; l < nbl; l++)
      for (c = 0; c < j; c++) tmp(l,c) = tmp(l+1, c+1);
    
    
    typename MAT::value_type vlp, vlpr, res = typename MAT::value_type(1);
    
    /* Triangulation Inferieure du systeme.                               */
    for (l = 0; l < nbl; l++)
    {
      /* Choix du pivot.                                                  */  
      for (c = l+1, pc = l; c < nbl; c++)
	if ( dal::abs(tmp(c,l)) > dal::abs(tmp(pc,l)) ) pc = c;
      
      /* Echange pour avoir le pivot en bonne place.                      */
      if (l != pc) tmp.line_exchange(l, pc);
      
      /* triangulation.                                                   */
      vlpr = tmp(l,l); res *= vlpr;
      if (dal::abs(vlpr) < EPS) return typename MAT::value_type(0);
      for (l2 = l+1; l2 < nbl; l2++ )
      {
	vlp = tmp(l2,l);
	for (c = l; c < nbl; c++)
	  tmp(l2, c) = tmp(l2, c) - tmp(l,c) * vlp / vlpr;
      }
    }

    return res;
  }

  template<class MAT> inline typename MAT::value_type mat_gauss_sub_det(MAT &m,
                         typename MAT::size_type i, typename MAT::size_type j,
							    double EPS = 1E-12)
  { MAT tmp; return mat_gauss_sub_det(m, i, j, tmp, EPS); }


  /* ******************************************************************** */
  /*	Decomposition L U.						  */
  /* ******************************************************************** */


  template <class MAT> void mat_decomposition_lu(const MAT &a, MAT &l, MAT &u)
  {
    typename MAT::size_type n = a.nrows(), i, j;
    if (l.nrows() != n || l.ncols() != n) l = a;
    if (u.nrows() != n || u.ncols() != n) u = a;

    if (n !=  a.ncols || n != l.ncols || n != u.ncols)
      DAL_THROW(dimension_error, "dimensions mismatch");

    l.fill(typename MAT::value_type(0)); u.fill(typename MAT::value_type(0));

    for (i = 0; i < n; i++)
    {
      typename MAT::value_type e = a(i,i) - partial_lc_product(l,i,u,i,i);
      typename MAT::value_type lii = sqrt(dal::abs(e)), uii;
      if (e == typename MAT::value_type(0))
	DAL_THROW(failure_error, "Non invertible matrix ?");
      if (e > typename MAT::value_type(0)) uii = lii; else uii = -lii;

      l(i,i) = lii; u(i,i) = uii;
      
      for (j = i+1; j < n; j++)
      { 
	typename MAT::value_type aux;
	if ((aux = (a(j,i) - partial_lc_product(l,j,u,i,i)) / uii ) != 0) 
	  l(j,i) = aux;
	if ((aux = (a(i,j) - partial_lc_product(l,i,u,j,i)) / lii ) != 0)
	  u(i,j) = aux;
      }
    }
  }

  /* ******************************************************************** */
  /*	Decomposition L Lt.						  */
  /* ******************************************************************** */

  template <class MAT> void mat_decomposition_llt(const MAT &a, MAT &l)
  {
    typename MAT::size_type n = a.nrows(), i, j;
    if (l.nrows() != n || l.ncols() != n) l = a;
    if (n != a.ncols() || n != l.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    l.fill((typename MAT::value_type)(0));
    typename MAT::value_type lii;
    for (i = 0; i < n; i++)
    {
      lii = a(i,i) - partial_ll_product(l,i,l,i,i);
      if (lii == typename MAT::value_type(0))
	DAL_THROW(failure_error, "Non invertible matrix ?");
      
      lii = sqrt(lii); l(i,i) = lii;
      for (j = i+1; j < n; j++)
      { l(j,i) = (a(i,j) - partial_ll_product(l,j,l,i,i)) / lii;}
    }
  }
  
  template <class MAT> void mat_decomposition_llt(MAT &a)
  {
    typename MAT::size_type n = a.nrows(), i, j;
    if (n != a.ncols())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    typename MAT::value_type lii;
    for (i = 0; i < n; i++)
    {
      lii = a(i,i) - partial_ll_product(a,i,a,i,i);
      if (lii == typename MAT::value_type(0))
	DAL_THROW(failure_error, "Non invertible matrix ?");

      lii = sqrt(lii); a(i,i) = lii;
      for (j = i+1; j < n; j++)
      { 
	a(j,i) = (a(i,j) - partial_ll_product(a,j,a,i,i)) / lii;
      }
    }
  }
  
  template <class MAT, class VECT>
    void mat_solve_llt(const MAT &l, const VECT &b, VECT &x)
  {
    typename MAT::size_type n = l.nrows(), j; 
    if (b.size() != x.size()) x = b;
    if ((n != l.ncols() || n != b.size() || n != x.size()))
      DAL_THROW(dimension_error, "dimensions mismatch");
     
    for (j = 0; j < n; j++)
    { x[j] = (b[j] - partial_lv_product(l, j, x, 0, j)) / l(j,j); }
    for (j = n-1; j != typename MAT::size_type(-1); j--)
    { x[j] = (x[j] - partial_cv_product(l, j, x, j+1, n)) / l(j,j); }
  }

  /* ******************************************************************** */
  /*    Inversion par cholesky d'une matrice S.D.P.	        	  */
  /* ******************************************************************** */

  template<class MAT> void mat_inv_triang_inf(MAT &m)
  {
    typename MAT::size_type n = m.nrows(), j, i;
    for(j = 0; j < n; j++)
    {
      m(j,j) = 1.0 / m(j,j);
      for(i = j+1; i < n; i++)
	m(i,j) = - partial_lc_product(m, i, m, j, i) / m(i,i);
    }
  }

  template<class MAT>
    typename MAT::value_type mat_inv_cholesky(MAT &m, MAT &tmp)
  {
    typename MAT::size_type n = m.nrows(), j;
    mat_decomposition_llt(m,tmp);
    typename MAT::value_type d = 1.0;
    for (j = 0; j < n; j++) d *= tmp(j,j);
    if (d != 0.0)
    { mat_inv_triang_inf(tmp); m = tmp; m.transpose(); m *= tmp; }
    return d*d;
  }
  
  template<class MAT> inline
    typename MAT::value_type mat_inv_cholesky(MAT &m)
  { MAT tmp; return mat_inv_cholesky(m, tmp); }
  
  typedef vsmatrix<scalar_type> base_matrix;


}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_MATRIX_H */
