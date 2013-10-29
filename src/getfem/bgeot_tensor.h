/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2000-2013 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file bgeot_tensor.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date October 09, 2000.
   @brief tensor class, used in mat_elem computations.
*/
#ifndef BGEOT_TENSOR_H__
#define BGEOT_TENSOR_H__

#include "bgeot_vector.h"
#include "getfem/getfem_omp.h"


namespace bgeot {

  /* ********************************************************************* */
  /*		Class tensor<T>.     		                           */
  /* ********************************************************************* */

  typedef size_t size_type;
  typedef gmm::uint16_type short_type;

  class multi_index : public std::vector<size_type> {
  public :
    
    void incrementation(const multi_index &m) { /* a compiler ... */
      iterator it = begin(), ite = end();
      const_iterator itm = m.begin();
      
      ++(*it);
      while (*it >= *itm && it != (ite-1)) { *it = 0; ++it; ++itm; ++(*it);}
    }
    
    void reset(void) { std::fill(begin(), end(), 0); }
    
    inline bool finished(const multi_index &m)
    { return ((*this)[size()-1] >= m[size()-1]); }
    
    multi_index(size_t n) : std::vector<size_type>(n)
    { std::fill(begin(), end(), size_type(0)); }
    multi_index(size_type i, size_type j)
      : std::vector<size_type>(2)
    { (*this)[0] = i; (*this)[1] = j; }
    multi_index(size_type i, size_type j, size_type k)
      : std::vector<size_type>(3)
    { (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; }
    multi_index(size_type i, size_type j, size_type k, size_type l)
      : std::vector<size_type>(4)
    { (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; (*this)[3] = l; } 

    multi_index(void) {}
    
    size_type memsize() const {
      return std::vector<size_type>::capacity()*sizeof(size_type) + 
	sizeof(multi_index);
    }
  };

  inline std::ostream &operator <<(std::ostream &o,
				   const multi_index& mi) { /* a compiler ...*/
    multi_index::const_iterator it = mi.begin(), ite = mi.end();
    bool f = true;
    o << "(";
    for ( ; it != ite; ++it) 
      { if (!f) o << ", "; o << *it; f = false; }
    o << ")";
    return o;
  }

  template<class T> class tensor : public std::vector<T> {
    protected:

    multi_index sizes_, coeff;
      
    public:

    typedef typename std::vector<T>::size_type size_type;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;

    template<class CONT> inline const T& operator ()(const CONT &c) const
    {
      typename CONT::const_iterator it = c.begin();
      multi_index::const_iterator q = coeff.begin(), e = coeff.end();
#ifndef NDEBUG
      multi_index::const_iterator qv = sizes_.begin();
#endif
      size_type d = 0;
      for ( ; q != e; ++q, ++it) { 
        d += (*q) * (*it);
        GMM_ASSERT2(*it < *qv++, "Index out of range.");
      }
      return *(this->begin() + d);
    }
    
    inline T& operator ()(size_type i, size_type j, size_type k,
                          size_type l) {
      GMM_ASSERT2(order() == 4, "Bad tensor order.");
      size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k + coeff[3]*l;
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }
    
    inline T& operator ()(size_type i, size_type j, size_type k) {
      GMM_ASSERT2(order() == 3, "Bad tensor order.");
      size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k;
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }
    
    inline T& operator ()(size_type i, size_type j) {
      GMM_ASSERT2(order() == 2, "Bad tensor order");
	size_type d = coeff[0]*i + coeff[1]*j;
	GMM_ASSERT2(d < size(), "Index out of range.");
	return *(this->begin() + d);
    }
    
    inline const T& operator ()(size_type i, size_type j, size_type k,
                                size_type l) const {
      GMM_ASSERT2(order() == 4, "Bad tensor order.");
      size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k + coeff[3]*l;
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }
    
    inline const T& operator ()(size_type i, size_type j,
                                size_type k) const {
      GMM_ASSERT2(order() == 3, "Bad tensor order.");
      size_type d = coeff[0]*i + coeff[1]*j + coeff[2]*k;
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }
    
    inline const T& operator ()(size_type i, size_type j) const {
      GMM_ASSERT2(order() == 2, "Bad tensor order.");
      size_type d = coeff[0]*i + coeff[1]*j;
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }

    template<class CONT> inline T& operator ()(const CONT &c) {
      typename CONT::const_iterator it = c.begin();
      multi_index::iterator q = coeff.begin(), e = coeff.end();
      size_type d = 0;
      for ( ; q != e; ++q, ++it) d += (*q) * (*it);
      
      GMM_ASSERT2(d < size(), "Index out of range.");
      return *(this->begin() + d);
    }
    
    inline size_type size(void) const { return std::vector<T>::size(); }
    inline size_type size(size_type i) const { return sizes_[i]; }
    inline const multi_index &sizes(void) const { return sizes_; }
    inline size_type order(void) const { return sizes_.size(); }
    
    void init(const multi_index &c) {
      multi_index::const_iterator it = c.begin();
      size_type d = 1;
      sizes_ = c; coeff.resize(c.size());
      multi_index::iterator p = coeff.begin(), pe = coeff.end();
      for ( ; p != pe; ++p, ++it) { *p = d; d *= *it; }
      this->resize(d);
    }
    
    void adjust_sizes(const multi_index &mi) {
      if (!mi.size() || (mi.size() != sizes().size())
          || !(std::equal(mi.begin(), mi.end(), sizes().begin())))
        init(mi);
    }
    
    tensor(const multi_index &c) { init(c); }
    tensor(size_type i, size_type j, size_type k, size_type l)
    { init(multi_index(i, j, k, l)); }
    tensor(void) {}

    /** reduction of tensor t with respect to index ni with matrix m:
     *  t(...,j,...) <-- t(...,i,..) m(i, j)
     */
    void mat_reduction(const tensor &t, const gmm::dense_matrix<T> &m, int ni);
    void mat_transp_reduction(const tensor &t, const gmm::dense_matrix<T> &m,
			      int ni);
    /** mm(i,j) = t(i,j,k,l) * m(k,l); For order four tensor. */
    void mat_mult(const gmm::dense_matrix<T> &m, gmm::dense_matrix<T> &mm);
    

    size_type memsize() const {
      return sizeof(T) * this->size()
	+ sizeof(*this) + sizes_.memsize() + coeff.memsize();
    }

    std::vector<T> &as_vector(void) { return *this; }
    const std::vector<T> &as_vector(void) const { return *this; }
    

    tensor<T>& operator +=(const tensor<T>& w)
    { gmm::add(w.as_vector(), this->as_vector()); return *this; }
    
    tensor<T>& operator -=(const tensor<T>& w) {
      gmm::add(gmm::scaled(w.as_vector(), T(-1)), this->as_vector());
      return *this;
    }
    
    tensor<T>& operator *=(const scalar_type w)
    { gmm::scale(this->as_vector(), w); return *this; }

    tensor<T>& operator /=(const scalar_type w)
    { gmm::scale(this->as_vector(), scalar_type(1)/w); return *this; }
  };

  template<class T> void tensor<T>::mat_transp_reduction
  (const tensor &t, const gmm::dense_matrix<T> &m, int ni) { 
    /* reduction du tenseur t par son indice ni et la matrice          */
    /* transposee de m.                                                */
    
	DEFINE_STATIC_THREAD_LOCAL(std::vector<T>*,tmp);
	DEFINE_STATIC_THREAD_LOCAL(multi_index*,mi);
	DEFINE_STATIC_THREAD_LOCAL_INITIALIZED(bool,isinit,false);

    if (!isinit) {
      tmp = new std::vector<T>(3); mi = new multi_index(); isinit = true;
    }
    
    *mi = t.sizes();
    size_type dimt = (*mi)[ni], dim = m.nrows();
    
    GMM_ASSERT2(dimt == m.ncols(), "Dimensions mismatch.");
    GMM_ASSERT2(&t != this, "Does not work when t and *this are the same.");

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
	for (size_type k = 0; k <= size_type(ni); ++k)
	  (*mi)[k] = size_type(sizes()[k] - 1);
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

  template<class T> void tensor<T>::mat_mult(const gmm::dense_matrix<T> &m,
                                             gmm::dense_matrix<T> &mm) {
    GMM_ASSERT2(order() == 4,
                "This operation is for order four tensors only.");
    GMM_ASSERT2(sizes_[2] == gmm::mat_nrows(m) &&
                sizes_[3] == gmm::mat_ncols(m), "Dimensions mismatch.");
    gmm::resize(mm, sizes_[0], sizes_[1]);
    gmm::clear(mm);

    const_iterator pt = this->begin();
    const_iterator pm = m.begin();
    for (size_type l = 0; l < sizes_[3]; ++l)
      for (size_type k = 0; k < sizes_[2]; ++k) {
        iterator pmm = mm.begin();
        for (size_type j = 0; j < sizes_[1]; ++j)
          for (size_type i = 0; i < sizes_[0]; ++i)
            *pmm++ += *pt++ * (*pm);
        ++pm;
      }
  }
  
  template<class T> void tensor<T>::mat_reduction
  (const tensor &t, const gmm::dense_matrix<T> &m, int ni) {
    /* reduction du tenseur t par son indice ni et la matrice m.       */
	DEFINE_STATIC_THREAD_LOCAL(std::vector<T>*,tmp);
	DEFINE_STATIC_THREAD_LOCAL(multi_index*,mi);
	DEFINE_STATIC_THREAD_LOCAL_INITIALIZED(bool,isinit,false);
    if (!isinit) {
      tmp = new std::vector<T>(3); mi = new multi_index(); isinit = true;
    }
    *mi = t.sizes();
    size_type dimt = (*mi)[ni], dim = m.ncols();
    GMM_ASSERT2(dimt == m.nrows(), "Dimensions mismatch.");
    GMM_ASSERT2(&t != this, "Does not work when t and *this are the same.");
    
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
	for (size_type k = 0; k <= size_type(ni); ++k)
	  (*mi)[k] = size_type(sizes()[k] - 1);
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
  

  template<class T> std::ostream &operator <<
    (std::ostream &o, const tensor<T>& t) { // a ameliorer ...
    o << "sizes " << t.sizes() << endl;
    o << t.as_vector();
    return o;
  }

  typedef tensor<scalar_type> base_tensor;


}  /* end of namespace bgeot.                                              */


#endif  /* BGEOT_TENSOR_H */
