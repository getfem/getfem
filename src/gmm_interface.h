/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_interface.h : generic algorithms on linear algebra       */
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

#ifndef __GMM_INTERFACE_H
#define __GMM_INTERFACE_H

namespace bgeot {
  template <class T, int N> class fsvector;
  template <class T> class vsvector;
  template <class VECT> class PT;
  template <class T, int N> class fsmatrix;
  template <class T> class vsmatrix;
}

namespace gmm {

  /* ********************************************************************* */
  /*                                                                       */
  /* What is needed for a Vector type :                                    */
  /*   Vector v(n) defines a vector with n components.                     */
  /*   v[i] allows to access to the ith component of v.                    */
  /*   linalg_traits<Vector> should be filled with appropriate definitions */
  /*                                                                       */
  /* What is needed for a Matrix type :                                    */
  /*   Matrix m(n, m) defines a matrix with n rows and m columns.          */
  /*   m(i, j) allows to access to the element at row i and column j.      */
  /*   linalg_traits<Matrix> should be filled with appropriate definitions */
  /*                                                                       */
  /* What is needed for an iterator on plain vector                        */
  /*    to be standard ramdom access iterator                              */
  /*                                                                       */
  /* What is needed for an iterator on a sparse vector                     */
  /*    to be a standard forward iterator                                  */
  /*    it.index() gives the index of the non-zero element.                */
  /*    if it.index() is equal to size_type(-1) this means that the        */
  /*       element is not present (eliminated : usefull for sub vectors).  */
  /*       then (*it) should return the value 0.                           */
  /*                                                                       */
  /* Remark : If original iterators are not convenient, they could be      */
  /*   redefined and interfaced in linalg_traits<Vector> without changin   */
  /*   the original Vector type.                                           */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*		Simple references on vectors            		   */
  /* ********************************************************************* */

  template <class PT> class simple_vector_ref {
    protected :
      PT l;

    public :
      typedef typename std::iterator_traits<PT>::value_type V;
      typedef typename std::iterator_traits<PT>::reference ref_V;
      typedef typename linalg_traits<V>::value_type value_type;
      typedef typename vect_ref_type<PT,  V>::access_type access_type;
      simple_vector_ref(ref_V v) : l(&v) {}
      ref_V deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      access_type operator[](size_type i) { return (*l)[i]; }
      value_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class PT> struct linalg_traits<simple_vector_ref<PT> > {
    typedef simple_vector_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename vect_ref_type<PT,  V>::iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    size_type size(const this_type &v) { return vect_size(v.deref()); }
    iterator begin(this_type &v) { return vect_begin(v.deref()); }
    const_iterator const_begin(const this_type &v) 
      { return vect_begin(v.deref()); }
    iterator end(this_type &v) { return vect_end(v.deref()); }
    const_iterator const_end(const this_type &v) 
      { return vect_end(v.deref()); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { clear(v.deref()); }
  };

  /* ********************************************************************* */
  /*		Scaled references on vectors            		   */
  /* ********************************************************************* */

  template <class V> struct scaled_vector_const_ref {

      const V *l;
      typename linalg_traits<V>::value_type r;
      scaled_vector_const_ref(const V &v,
			      typename linalg_traits<V>::value_type x)
	: l(&v), r(x) {}
      typedef typename linalg_traits<V>::value_type value_type;
      const V &deref(void) const { return *l; }
      value_type operator[](size_type i) const { return (*l)[i] * r; }
  };

  template <class IT> struct scaled_const_iterator {
    typedef typename std::iterator_traits<IT>::value_type      value_type;
    typedef typename std::iterator_traits<IT>::pointer         pointer;
    typedef typename std::iterator_traits<IT>::reference       reference;
    typedef typename std::iterator_traits<IT>::difference_type difference_type;
    typedef typename std::iterator_traits<IT>::iterator_category
    iterator_category;

    IT it;
    value_type r;
    
    scaled_const_iterator(void) {}
    scaled_const_iterator(const IT &i, value_type x) : it(i), r(x) {}
    
    inline size_type index(void) const { return it.index(); }
    scaled_const_iterator operator ++(int)
    { scaled_const_iterator tmp = *this; ++it; return tmp; }
    scaled_const_iterator operator --(int) 
    { scaled_const_iterator tmp = *this; --it; return tmp; }
    scaled_const_iterator &operator ++() { ++it; return *this; }
    scaled_const_iterator &operator --() { --it; return *this; }
    scaled_const_iterator &operator +=(difference_type i)
      { it += i; return *this; }
    scaled_const_iterator &operator -=(difference_type i)
      { it -= i; return *this; }
    scaled_const_iterator operator +(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb += i); }
    scaled_const_iterator operator -(difference_type i) const
      { scaled_const_iterator itb = *this; return (itb -= i); }
    difference_type operator -(const scaled_const_iterator &i) const
      { return difference_type(it - i.it); }
    
    value_type operator  *() const { return (*it) * r; }
    value_type operator [](size_type ii) const { return it[in+ii] * r; }
    
    bool operator ==(const scaled_const_iterator &i) const
      { return (i.it == it); }
    bool operator !=(const scaled_const_iterator &i) const
      { return (i.it != it); }
    bool operator < (const scaled_const_iterator &i) const
      { return (it < i.it); }
  };

  template <class V> struct linalg_traits<scaled_vector_const_ref<V> > {
    typedef scaled_vector_const_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef scaled_const_iterator<typename linalg_traits<V>::const_iterator>
    iterator;
    typedef scaled_const_iterator<typename linalg_traits<V>::const_iterator>
    const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    size_type size(const this_type &v)  { return vect_size(v.deref()); }
    iterator begin(this_type &v)
    { return iterator(vect_begin(v.deref()), v.r); }
    const_iterator const_begin(const this_type &v)
    { return iterator(vect_begin(v.deref()), v.r); }
    iterator end(this_type &v)
    { return iterator(vect_end(v.deref()), v.r); }
    const_iterator const_end(const this_type &v)
    { return iterator(vect_end(v.deref()), v.r); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };

  /* ********************************************************************* */
  /*		transposed reference                    		   */
  /* ********************************************************************* */
  
  template <class PT> class  transposed_ref {
    protected :
      PT l;

    public :
      typedef typename std::iterator_traits<PT>::value_type V;
      typedef typename std::iterator_traits<PT>::reference ref_V;
      typedef typename linalg_traits<V>::value_type value_type;
      typedef typename mat_ref_type<PT,  V>::access_type access_type;
      transposed_ref(ref_V v) : l(&v) {}
      ref_V deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      access_type operator()(size_type i, size_type j) { return (*l)(j,i); }
      value_type operator()(size_type i, size_type j) const
      { return (*l)(j,i); }
  };

  template <class PT> struct linalg_traits<transposed_ref<PT> > {
    typedef transposed_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::sub_col_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_col_type const_sub_row_type;
    typedef typename linalg_traits<V>::sub_row_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_row_type
    const_sub_col_type;
    typedef typename transposed_type<typename 
      linalg_traits<V>::sub_orientation>::t_type sub_orientation;
    size_type nrows(const this_type &v) { return mat_ncols(v.deref()); }
    size_type ncols(const this_type &v) { return mat_nrows(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return mat_const_col(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return mat_const_row(v.deref(), i); }
    sub_row_type row(this_type &v, size_type i)
    { return mat_col(v.deref(), i); }
    sub_col_type col(this_type &v, size_type i)
    { return mat_row(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { clear(v.deref()); }
  };

  template <class PT> class  vect_transposed_ref {
    protected :
      PT l;

    public :
      typedef typename std::iterator_traits<PT>::value_type V;
      typedef typename std::iterator_traits<PT>::reference ref_V;
      typedef typename linalg_traits<V>::value_type value_type;
      typedef typename vect_ref_type<PT,  V>::access_type access_type;
      vect_transposed_ref(ref_V v) : l(&v) {}
      ref_V deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      access_type operator()(size_type i, size_type j) { return (*l)[j]; }
      value_type operator()(size_type i, size_type j) const
      { return (*l)[j]; }
  };

  template <class PT> struct linalg_traits<vect_transposed_ref<PT> > {
    typedef vect_transposed_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef simple_vector_ref<V *> sub_row_type;
    typedef simple_vector_ref<const V *> const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &v) { return 1; }
    size_type ncols(const this_type &v) { return vect_size(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type)
    { return const_sub_row_type(v.deref()); }
    sub_row_type row(this_type &v, size_type)
    { return sub_row_type(v.deref()); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { clear(v.deref()); }
  };

  /* ********************************************************************* */
  /*		                                         		   */
  /*		Traits for S.T.L. object                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class T> struct linalg_traits<std::vector<T> > {
    typedef std::vector<T> this_type;
    typedef linalg_false is_reference;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { std::fill(v.begin(), v.end(), T(0)); }
  };

  // to be done :  std::valarray<T> ...

  /* ********************************************************************* */
  /*		                                         		   */
  /*		Traits for dal objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class IT>
  class tab_ref_with_origin : public dal::tab_ref<IT> {
  protected :
    const void *pt;
    
  public :
    typedef dal::tab_ref<IT>                          base_tab;
    typedef typename base_tab::value_type             value_type;
    typedef typename base_tab::pointer                pointer;
    typedef typename base_tab::const_pointer          const_pointer;
    typedef typename base_tab::reference              reference;
    typedef typename base_tab::const_reference        const_reference;
    typedef typename base_tab::difference_type        difference_type;
    typedef typename base_tab::iterator               iterator;
    typedef typename base_tab::const_iterator         const_iterator;
    typedef typename base_tab::const_reverse_iterator const_reverse_iterator;
    typedef typename base_tab::reverse_iterator       reverse_iterator;
    typedef typename base_tab::size_type              size_type;

    tab_ref_with_origin(void) {}
    tab_ref_with_origin(const IT &b, const IT &e, const void *p)
      : dal::tab_ref<IT>(b,e), pt(p) {}

    const void *origin(void) const { return pt; }
  };

  template <class IT> struct linalg_traits<tab_ref_with_origin<IT> > {
    typedef tab_ref_with_origin<IT> this_type;
    typedef linalg_true is_reference;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::value_type& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v)
    { std::fill(v.begin(), v.end(), value_type(0)); }
  };

  template <class IT>
  class tab_ref_reg_spaced_with_origin : public dal::tab_ref_reg_spaced<IT> {
  protected :
    const void *pt;
    
  public :
    typedef dal::tab_ref_reg_spaced<IT>               base_tab;
    typedef typename base_tab::value_type             value_type;
    typedef typename base_tab::pointer                pointer;
    typedef typename base_tab::const_pointer          const_pointer;
    typedef typename base_tab::reference              reference;
    typedef typename base_tab::const_reference        const_reference;
    typedef typename base_tab::difference_type        difference_type;
    typedef typename base_tab::iterator               iterator;
    typedef typename base_tab::const_iterator         const_iterator;
    typedef typename base_tab::const_reverse_iterator const_reverse_iterator;
    typedef typename base_tab::reverse_iterator       reverse_iterator;
    typedef typename base_tab::size_type              size_type;

    tab_ref_reg_spaced_with_origin(void) {}
    tab_ref_reg_spaced_with_origin(const IT &b, const IT &e, size_type n,
				   const void *p) 
      : dal::tab_ref_reg_spaced<IT>(b,e,n), pt(p) {}

    const void *origin(void) const { return pt; }
  };

  template <class IT>
  struct linalg_traits<tab_ref_reg_spaced_with_origin<IT> > {
    typedef tab_ref_reg_spaced_with_origin<IT> this_type;
    typedef linalg_true is_reference;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::value_type& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v)
    { std::fill(v.begin(), v.end(), value_type(0)); }
  };

  template <class IT, class ITINDEX>
  class tab_ref_index_ref_with_origin 
    : public dal::tab_ref_index_ref<IT, ITINDEX> {
  protected :
    const void *pt;
    
  public :
    typedef dal::tab_ref_index_ref<IT, ITINDEX>       base_tab;
    typedef typename base_tab::value_type             value_type;
    typedef typename base_tab::pointer                pointer;
    typedef typename base_tab::const_pointer          const_pointer;
    typedef typename base_tab::reference              reference;
    typedef typename base_tab::const_reference        const_reference;
    typedef typename base_tab::difference_type        difference_type;
    typedef typename base_tab::iterator               iterator;
    typedef typename base_tab::const_iterator         const_iterator;
    typedef typename base_tab::const_reverse_iterator const_reverse_iterator;
    typedef typename base_tab::reverse_iterator       reverse_iterator;
    typedef typename base_tab::size_type              size_type;

    tab_ref_index_ref_with_origin(void) {}
    tab_ref_index_ref_with_origin(const IT &b, const ITINDEX &bi,
				  const ITINDEX &ei, const void *p)
      : dal::tab_ref_index_ref<IT, ITINDEX>(b, bi, ei), pt(p) {}

    const void *origin(void) const { return pt; }
  };

  template <class IT, class ITINDEX>
  struct linalg_traits<tab_ref_index_ref_with_origin<IT, ITINDEX> > {
    typedef tab_ref_index_ref_with_origin<IT, ITINDEX> this_type;
    typedef linalg_true is_reference;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::value_type& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v)
    { std::fill(v.begin(), v.end(), value_type(0)); }
  };

  /* ********************************************************************* */
  /*		                                         	 	   */
  /*		Traits for bgeot objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class T, int N> struct linalg_traits<bgeot::fsvector<T, N> > {
    typedef bgeot::fsvector<T, N> this_type;
    typedef linalg_false is_reference;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &) { return N; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };

  

  template <class T> struct linalg_traits<bgeot::vsvector<T> > {
    typedef bgeot::vsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };

  

  template <class VECT> struct linalg_traits<bgeot::PT<VECT> > {
    typedef bgeot::PT<VECT> this_type;
    typedef linalg_false is_reference;
    typedef typename linalg_traits<VECT>::value_type value_type;
    typedef typename linalg_traits<VECT>::reference_type reference_type;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<VECT>::iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return linalg_traits<VECT>().begin(v); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<VECT>().const_begin(v); }
    iterator end(this_type &v) { return linalg_traits<VECT>().end(v); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<VECT>().const_end(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };

  

  template <class T, int N> struct linalg_traits<bgeot::fsmatrix<T, N> > {
    typedef bgeot::fsmatrix<T, N> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type::iterator>
    sub_row_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type
    ::const_iterator> const_sub_row_type;
    typedef tab_ref_with_origin<typename this_type::iterator> sub_col_type;
    typedef tab_ref_with_origin<typename this_type::const_iterator>
    const_sub_col_type;
    typedef col_and_row sub_orientation;
    size_type nrows(const this_type &m) { return N; }
    size_type ncols(const this_type &m) { return N; } 
    const_sub_row_type row(const this_type &m, size_type i)
    {
      typename this_type::const_iterator b = m.begin();
      return const_sub_row_type(b + i, b + i + mat_ncols(m)* N, N, &m);
    }
    const_sub_col_type col(const this_type &m, size_type i) {
      typename this_type::const_iterator b = m.begin();
      return const_sub_col_type(b + i * N, b + (i+1) * N, &m);
    }
    sub_row_type row(this_type &m, size_type i)
    {
      typename this_type::iterator b = m.begin();
      return sub_row_type(b + i, b + i + mat_ncols(m)* N, N, &m);
    }
    sub_col_type col(this_type &m, size_type i) { 
      typename this_type::iterator b = m.begin();
      return sub_col_type(b + i * N, b + (i+1) * N, &m);
    }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };

 

  template <class T> struct linalg_traits<bgeot::vsmatrix<T> > {
    typedef bgeot::vsmatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type::iterator>
    sub_row_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type
    ::const_iterator> const_sub_row_type;
    typedef tab_ref_with_origin<typename this_type::iterator> sub_col_type;
    typedef tab_ref_with_origin<typename this_type::const_iterator>
    const_sub_col_type;
    typedef col_and_row sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); } 
    const_sub_row_type row(const this_type &m, size_type i)
    {
      typename this_type::const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return const_sub_row_type(b + i, b + i + mat_ncols(m)* nr, nr, &m);
    }
    const_sub_col_type col(const this_type &m, size_type i) {
      typename this_type::const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return const_sub_col_type(b + i * nr, b + (i+1) * nr, &m);
    }
    sub_row_type row(this_type &m, size_type i)
    {
      typename this_type::iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return sub_row_type(b + i, b + i + mat_ncols(m)* nr, nr, &m);
    }
    sub_col_type col(this_type &m, size_type i) { 
      typename this_type::iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return sub_col_type(b + i * nr, b + (i+1) * nr, &m);
    }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };


  /* ******************************************************************** */
  /*		Reference on a compressed sparse column matrix            */
  /* ******************************************************************** */

  template <class PT1, class PT2, int shift = 0>
  struct cs_vector_ref_iterator {
    PT1 pr;
    PT2 ir;

    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef PT1 pointer;
    typedef typename std::iterator_traits<PT1>::reference  reference;
    typedef size_t        size_type;
    typedef ptrdiff_t     difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef cs_vector_ref_iterator<PT1, PT2, shift> iterator;
    
    cs_vector_ref_iterator(void) {}
    cs_vector_ref_iterator(PT1 p1, PT2 p2) : pr(p1), ir(p2) {}

    inline size_type index(void) const { return (*ir) - shift; }
    iterator &operator ++() { ++pr; ++ir; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    
    reference operator  *() const { return *pr; }
    pointer   operator ->() const { return pr; }
    
    bool operator ==(const iterator &i) const { return (i.pr==pr);}
    bool operator !=(const iterator &i) const { return (i.pr!=pr);}
  };
    
  // read only reference on a compressed sparse vector
  template <class PT1, class PT2, int shift = 0> struct cs_vector_ref {
    PT1 pr;
    PT2 ir;
    size_type n, _size;
    
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef cs_vector_ref_iterator<PT1, PT2, shift>  iterator;
    typedef cs_vector_ref_iterator<typename const_pointer<PT1>::pointer,
      typename const_pointer<PT2>::pointer, shift> const_iterator;

    cs_vector_ref(PT1 pt1, PT2 pt2, size_type nnz, size_type ns)
      : pr(pt1), ir(pt2), n(nnz), _size(ns) {}
    cs_vector_ref(void) {}

    size_type size(void) const { return _size; }
    
    iterator begin(void) { return iterator(pr, ir); }
    const_iterator begin(void) const { return const_iterator(pr, ir); }
    iterator end(void) { return iterator(pr+n, ir+n); }
    const_iterator end(void) const { return const_iterator(pr+n, ir+n); }
    
    value_type operator[](size_type i) const {
      static value_type zero(0);
      if (n == 0) return zero;
      PT2 p = std::lower_bound(ir, ir+n, i+shift);
      return (*p == i+shift) ? *p : zero;
    }
  };

  template <class PT1, class PT2, int shift>
  struct linalg_traits<cs_vector_ref<PT1, PT2, shift> > {
    typedef cs_vector_ref<PT1, PT2, shift> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::value_type reference_type;
    typedef cs_vector_ref_iterator<PT1, PT2, shift>  iterator;
    typedef cs_vector_ref_iterator<typename const_pointer<PT1>::pointer,
      typename const_pointer<PT2>::pointer, shift> const_iterator;
    typedef abstract_sparse storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.pr; }
  };


  template <class PT1, class PT2, class PT3, int shift = 0>
  struct csc_matrix_ref {
    PT1 pr;
    PT2 ir;
    PT3 jc;
    size_type nc, nr;
    
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::reference access_type;
    csc_matrix_ref(PT1 pt1, PT2 pt2, PT3 pt3, size_type nrr, size_type ncc)
      : pr(pt1), ir(pt2), jc(pt3), nc(ncc), nr(nrr) {}
    csc_matrix_ref(void) {}

    // void do_clear(void) ... to be done
    
    size_type nrows(void) const { return nr; }
    size_type ncols(void) const { return nc; }
   
    // access_type operator()(size_type i, size_type j) // to be done
    //  { return mat_col(*this, j)[i]; }
    value_type operator()(size_type i, size_type j) const
      { return mat_col(*this, j)[i]; }
  };
  
  template <class PT1, class PT2, class PT3, int shift>
    struct linalg_traits<csc_matrix_ref<PT1, PT2, PT3, shift> > {
    typedef csc_matrix_ref<PT1, PT2, PT3, shift> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::value_type reference_type;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef cs_vector_ref<typename const_pointer<PT1>::pointer,
      typename const_pointer<PT2>::pointer, shift> sub_col_type;
    typedef cs_vector_ref<typename const_pointer<PT1>::pointer,
      typename const_pointer<PT2>::pointer, shift> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) {
      return const_sub_col_type(m.pr + m.jc[i]-shift, m.ir + m.jc[i] - shift,
				m.jc[i+1] - m.jc[i], m.nr);
    }
    const void* origin(const this_type &m) { return m.pr; }
    void do_clear(this_type &m) { m.do_clear(); }
  };

}


#endif //  __GMM_INTERFACE_H
