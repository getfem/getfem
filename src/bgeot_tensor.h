// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_tensor.h : plain tensors
//           
// Date    : October 09, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#ifndef BGEOT_TENSOR_H__
#define BGEOT_TENSOR_H__

#include <bgeot_vector.h>

namespace bgeot
{

  /* ********************************************************************* */
  /*		Class tensor<T>.     		                           */
  /* ********************************************************************* */

  typedef size_t size_type;
  typedef dal::uint16_type short_type;

  class multi_index : public std::vector<short_type> {
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
    
    multi_index(size_type i, size_type j)
      : std::vector<short_type>(2) {
      (*this)[0] = i; (*this)[1] = j; 
    } 
    multi_index(size_type i, size_type j, size_type k, size_type l)
      : std::vector<short_type>(4) {
      (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; (*this)[3] = l; 
    } 

    multi_index(void) {}
    
    size_type memsize() const {
      return std::vector<short_type>::capacity()*sizeof(short_type) + 
	sizeof(multi_index);
    }
  };

  inline std::ostream &operator <<(std::ostream &o,
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

      multi_index sizes_, coeff;
      
    public:

      typedef typename vsvector<T>::size_type size_type;
      typedef typename vsvector<T>::iterator iterator;
      typedef typename vsvector<T>::const_iterator const_iterator;

      template<class CONT> inline const T& operator ()(const CONT &c) const
      {
	typename CONT::const_iterator it = c.begin();
	multi_index::const_iterator q = coeff.begin(), e = coeff.end();
	#ifdef GETFEM_VERIFY
	  multi_index::const_iterator qv = sizes_.begin();
	#endif
	size_type d = 0;
	for ( ; q != e; ++q, ++it)
	{ 
	  d += (*q) * (*it);
	  #ifdef GETFEM_VERIFY
	    if (*it >= *qv) DAL_THROW(std::out_of_range, "index out of range");
	    ++qv;
          #endif
	}
	return *(this->begin() + d);
      }

      inline T& operator ()(size_type i, size_type j, size_type k,
			    size_type l) {
	if (order() != 4)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k + coeff[3]*l;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }
    
      inline T& operator ()(size_type i, size_type j, size_type k) {
	if (order() != 3)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }
    
      inline T& operator ()(size_type i, size_type j) {
	if (order() != 2)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }

      inline const T& operator ()(size_type i, size_type j, size_type k,
			    size_type l) const {
	if (order() != 4)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k + coeff[3]*l;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }
    
      inline const T& operator ()(size_type i, size_type j,
				  size_type k) const {
	if (order() != 3)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }
    
      inline const T& operator ()(size_type i, size_type j) const {
	if (order() != 2)
	  DAL_THROW(std::out_of_range, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j;
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }


    

      template<class CONT> inline T& operator ()(const CONT &c) {
	typename CONT::const_iterator it = c.begin();
	multi_index::iterator q = coeff.begin(), e = coeff.end();
	size_type d = 0;
	for ( ; q != e; ++q, ++it) d += (*q) * (*it);
	
	if (d >= size()) DAL_THROW(std::out_of_range, "index out of range");
	return *(this->begin() + d);
      }

      inline size_type size(void) const { return vsvector<T>::size(); }
      inline size_type size(int i) const { return sizes_[i]; }
      inline const multi_index &sizes(void) const { return sizes_; }
      inline size_type order(void) const { return sizes_.size(); }

      void init(const multi_index &c)
      {
	multi_index::const_iterator it = c.begin();
	size_type d = 1;
	sizes_ = c; coeff.resize(c.size());
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
    tensor(size_type i, size_type j, size_type k, size_type l)
    { init(multi_index(i, j, k, l)); }
    tensor(void) {}

      void mat_transp_reduction(const tensor &t, const gmm::dense_matrix<T> &m, int ni);

      void mat_reduction(const tensor &t, const gmm::dense_matrix<T> &m, int ni);
    size_type memsize() const { return vsvector<T>::memsize() +
				  sizes_.memsize() + coeff.memsize(); }
  };

  template<class T> void tensor<T>::mat_transp_reduction (const tensor &t,
					      const gmm::dense_matrix<T> &m, int ni) { 
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
    
    if (dimt != m.ncols()) DAL_THROW(dimension_error, "dimensions mismatch");
    if (&t == this)
      DAL_THROW(std::invalid_argument,
		"does not work when t and *this are the same");

    (*mi)[ni] = dim;
    if (tmp->size() < dimt) tmp->resize(dimt);
    adjust_sizes(*mi);
    const_iterator pft = t.begin();
    iterator pf = this->begin();
    size_type dd  =   coeff[ni]*(  sizes()[ni]-1)-1, co  =   coeff[ni];
    size_type ddt = t.coeff[ni]*(t.sizes()[ni]-1)-1, cot = t.coeff[ni];
    std::fill(mi->begin(), mi->end(), 0);
    for (;!mi->finished(sizes()); mi->incrementation(sizes()), ++pf, ++pft)
      if ((*mi)[ni] != 0) { 
	for (short_type k = 0; k <= ni; ++k) (*mi)[k] = sizes()[k] - 1;
	pf += dd; pft += ddt;
      }
      else {
	const_iterator pl = pft; iterator pt = tmp->begin();
	for(size_type k = 0; k < dimt; ++k, pl += cot, ++pt) *pt = *pl;
	
	iterator pff = pf;
	for (size_type k = 0; k < dim; ++k, pff += co) {
	  *pff = T(0); pt = tmp->begin(); pl = m.begin() + k;
	  for (size_type l = 0; l < dimt; ++l, ++pt, pl += dim)
	    *pff += (*pl) * (*pt);
	}
      }
  }
  
  template<class T> void tensor<T>::mat_reduction(const tensor &t,
					 const gmm::dense_matrix<T> &m, int ni) {
    /* reduction du tenseur t par son indice ni et la matrice m.       */
    static std::vector<T> *tmp;
    static multi_index *mi;
    static bool isinit = false;
    if (!isinit) {
      tmp = new std::vector<T>(3); mi = new multi_index(); isinit = true;
    }
    *mi = t.sizes();
    size_type dimt = (*mi)[ni], dim = m.ncols();
    if (dimt != m.nrows())
      DAL_THROW(dimension_error, "dimensions mismatch");
    if (&t == this)
      DAL_THROW(std::invalid_argument,
		"does not work when t and *this are the same");
    
    (*mi)[ni] = dim;
    if (tmp->size() < dimt) tmp->resize(dimt);
    adjust_sizes(*mi);
    const_iterator pft = t.begin();
    iterator pf = this->begin();
    size_type dd  =   coeff[ni]*(  sizes()[ni]-1)-1, co  =   coeff[ni];
    size_type ddt = t.coeff[ni]*(t.sizes()[ni]-1)-1, cot = t.coeff[ni];
    std::fill(mi->begin(), mi->end(), 0);
    for (;!mi->finished(sizes()); mi->incrementation(sizes()), ++pf, ++pft)
      if ((*mi)[ni] != 0) { 
	for (short_type k = 0; k <= ni; ++k) (*mi)[k] = sizes()[k] - 1;
	pf += dd; pft += ddt;
      }
      else {
	const_iterator pl = pft; iterator pt = tmp->begin();
	for(size_type k = 0; k < dimt; ++k, pl += cot, ++pt) *pt = *pl;
	
	iterator pff = pf; pl = m.begin();
	for (size_type k = 0; k < dim; ++k, pff += co) {
	  *pff = T(0); pt = tmp->begin();
	  for (size_type l = 0; l < dimt; ++l, ++pt, ++pl)
	    *pff += (*pl) * (*pt);
	}
      }
  }
  

  template<class T> std::ostream &operator <<(std::ostream &o,
					      const tensor<T>& t)
  { // a ameliorer ...
    o << "sizes " << t.sizes() << endl;
    o << *((const vsvector<T> *)(&t));
    return o;
  }

  typedef tensor<scalar_type> base_tensor;


}  /* end of namespace bgeot.                                              */


#endif  /* BGEOT_TENSOR_H__ */
