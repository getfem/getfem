/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_tensor.h : plain tensors                               */
/*     									   */
/*                                                                         */
/* Date : October 09, 2000.                                                */
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


#ifndef __BGEOT_TENSOR_H
#define __BGEOT_TENSOR_H

#include <bgeot_matrix.h>

namespace bgeot
{

  /* ********************************************************************* */
  /*		Class tensor<T>.     		                           */
  /* ********************************************************************* */

  typedef size_t size_type;
  typedef dal::uint16_type short_type;

  class multi_index : public std::vector<short_type>
  {
    public :

      void incrementation(const multi_index &m)
      { /* a compiler ... */
	iterator it = begin(), ite = end();
	const_iterator itm = m.begin();
	
	++(*it);
	while (*it >= *itm && it != (ite-1)) { *it = 0; ++it; ++itm; ++(*it);}
      }

      void reset(void) { std::fill(begin(), end(), 0); }

      inline bool finished(const multi_index &m)
      { return ((*this)[size()-1] >= m[size()-1]); }

      multi_index(size_t n) : std::vector<short_type>(n)
      { std::fill(begin(), end(), 0); }

      multi_index(void) {}

  };

  inline STD_NEEDED ostream &operator <<(STD_NEEDED ostream &o,
					 const multi_index& mi)
  {  /* a compiler ... */
    multi_index::const_iterator it = mi.begin(), ite = mi.end();
    bool f = true;
    o << "(";
    for ( ; it != ite; ++it) 
      { if (!f) o << ", "; o << *it; f = false; }
    o << ")";
    return o;
  }

  template<class T> class tensor : public vsvector<T>
  {
    protected:

      multi_index _sizes, coeff;
      
    public:

      typedef typename vsvector<T>::size_type size_type;
      typedef typename vsvector<T>::iterator iterator;
      typedef typename vsvector<T>::const_iterator const_iterator;

      template<class CONT> inline const T& operator ()(const CONT &c) const
      {
	typename CONT::const_iterator it = c.begin();
	multi_index::const_iterator q = coeff.begin(), e = coeff.end();
	#ifdef __GETFEM_VERIFY
	  multi_index::const_iterator qv = _sizes.begin();
	#endif
	size_type d = 0;
	for ( ; q != e; ++q, ++it)
	{ 
	  d += (*q) * (*it);
	  #ifdef __GETFEM_VERIFY
            assert(*it < *qv); ++qv;
          #endif
	}
	return *(begin() + d);
      }

      template<class CONT> inline T& operator ()(const CONT &c)
      {
	typename CONT::const_iterator it = c.begin();
	multi_index::iterator q = coeff.begin(), e = coeff.end();
	size_type d = 0;
	for ( ; q != e; ++q, ++it) d += (*q) * (*it);
	#ifdef __GETFEM_VERIFY
          assert(d < size());
        #endif
	return *(begin() + d);
      }

      inline size_type size(void) const { return vsvector<T>::size(); }
      inline size_type size(int i) const { return _sizes[i]; }
      inline const multi_index &sizes(void) const { return _sizes; }
      inline size_type order(void) const { return _sizes.size(); }

      void init(const multi_index &c)
      {
	multi_index::const_iterator it = c.begin();
	size_type d = 1;
	_sizes = c; coeff.resize(c.size());
	multi_index::iterator p = coeff.begin(), pe = coeff.end();
	for ( ; p != pe; ++p, ++it) { *p = d; d *= *it; }
	resize(d);
      }

      void adjust_sizes(const multi_index &mi)
      {
	if ((mi.size() != sizes().size())
	    || !(std::equal(mi.begin(), mi.end(), sizes().begin())))
	  init(mi);
      }

      tensor(const multi_index &c) { init(c); }
      tensor(void) {}

      void mat_transp_reduction(const tensor &t, const vsmatrix<T> &m, int ni)
      { 
	/* reduction du tenseur t par son indice ni et la matrice          */
	/* transposee de m.                                                */

	static std::vector<T> *tmp;
	static multi_index *mi;
	static bool isinit = false;
	if (!isinit) {
	  tmp = new std::vector<T>(3); mi = new multi_index(); isinit = true;
	}

	*mi = t.sizes();
	size_type dimt = (*mi)[ni], dim = m.nrows();
	#ifdef __GETFEM_VERIFY
          assert(dimt == m.ncols());
        #endif
	(*mi)[ni] = dim;
	if (tmp->size() < dimt) tmp->resize(dimt);
	adjust_sizes(*mi);
	const_iterator pft = t.begin();
	iterator pf = begin();
	size_type dd  =   coeff[ni]*(  sizes()[ni]-1)-1, co  =   coeff[ni];
	size_type ddt = t.coeff[ni]*(t.sizes()[ni]-1)-1, cot = t.coeff[ni];
	std::fill(mi->begin(), mi->end(), 0);
	for (;!mi->finished(sizes()); mi->incrementation(sizes()), ++pf, ++pft)
	  if ((*mi)[ni] != 0)
	  { 
	    for (short_type k = 0; k <= ni; ++k) (*mi)[k] = sizes()[k] - 1;
	    pf += dd; pft += ddt;
	  }
	  else
	  {
	    const_iterator pl = pft; iterator pt = tmp->begin();
	    for(size_type k = 0; k < dimt; ++k, pl += cot, ++pt) *pt = *pl;
	    
	    iterator pff = pf; pl = m.begin();
	    for (size_type k = 0; k < dim; ++k, pff += co)
	    {
	      *pff = T(0); pt = tmp->begin();
	      for (size_type l = 0; l < dimt; ++l, ++pt, ++pl)
		*pff += (*pl) * (*pt);
	    }
	  }
      }

      void mat_reduction(const tensor &t, const vsmatrix<T> &m, int ni)
      {
	/* reduction du tenseur t par son indice ni et la matrice m.       */
	static std::vector<T> *tmp;
	static multi_index *mi;
	static bool isinit = false;
	if (!isinit) {
	  tmp = new std::vector<T>(3); mi = new multi_index(); isinit = true;
	}
	*mi = t.sizes();
	size_type dimt = (*mi)[ni], dim = m.ncols();
	#ifdef __GETFEM_VERIFY
          assert(dimt == m.nrows());
        #endif
	(*mi)[ni] = dim;
	if (tmp->size() < dimt) tmp->resize(dimt);
	adjust_sizes(*mi);
	const_iterator pft = t.begin();
	iterator pf = begin();
	size_type dd  =   coeff[ni]*(  sizes()[ni]-1)-1, co  =   coeff[ni];
	size_type ddt = t.coeff[ni]*(t.sizes()[ni]-1)-1, cot = t.coeff[ni];
	std::fill(mi->begin(), mi->end(), 0);
	for (;!mi->finished(sizes()); mi->incrementation(sizes()), ++pf, ++pft)
	  if ((*mi)[ni] != 0)
	  { 
	    for (short_type k = 0; k <= ni; ++k) (*mi)[k] = sizes()[k] - 1;
	    pf += dd; pft += ddt;
	  }
	  else
	  {
	    const_iterator pl = pft; iterator pt = tmp->begin();
	    for(size_type k = 0; k < dimt; ++k, pl += cot, ++pt) *pt = *pl;
	    
	    iterator pff = pf;
	    for (size_type k = 0; k < dim; ++k, pff += co)
	    {
	      *pff = T(0); pt = tmp->begin(); pl = m.begin() + k;
	      for (size_type l = 0; l < dimt; ++l, ++pt, pl += dim)
		*pff += (*pl) * (*pt);
	    }
	  }
      }
  };

  template<class T> STD_NEEDED ostream &operator <<(STD_NEEDED ostream &o,
						    const tensor<T>& t)
  { // a ameliorer ...
    o << "sizes " << t.sizes() << endl;
    o << *((const vsvector<T> *)(&t));
    return o;
  }

}  /* end of namespace bgeot.                                              */


#endif  /* __BGEOT_TENSOR_H */
