/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_linalg_interface.h : generic algorithms on linear      */
/*                                      algebra                            */
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

#ifndef __BGEOT_LINALG_INTERFACE_H
#define __BGEOT_LINALG_INTERFACE_H

namespace bgeot {

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

  template <class V> struct linalg_traits;

  template <> struct linalg_traits<abstract_null_type> {
    typedef abstract_null_type is_reference;
    typedef abstract_null_type this_type;
    typedef abstract_null_type linalg_type;
    typedef abstract_null_type value_type;
    typedef abstract_null_type reference_type;
    typedef void * iterator;
    typedef const void * const_iterator;
    typedef abstract_null_type storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type sub_orientation;
    size_type size(const this_type &) { return 0; }
    size_type nrows(const this_type &) { return 0; }
    size_type ncols(const this_type &) { return 0; }
    iterator begin(this_type &) { return 0; }
    const_iterator const_begin(const this_type &) { return 0; }
    iterator end(this_type &) { return 0; }
    const_iterator const_end(const this_type &) { return 0; }
    const_sub_row_type row(const this_type &, size_type)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &, size_type)
      { DAL_THROW(failure_error,"Columns inaccessible for this object"); }
    sub_row_type row(this_type &, size_type)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &, size_type) 
      { DAL_THROW(failure_error,"Columns inaccessible for this object"); }
    const void* origin(const this_type &) { return 0; }
    void do_clear(this_type &) { }
  };

  /* ********************************************************************* */
  /*		Simple references on vectors            		   */
  /* ********************************************************************* */
  
  template <class V> class simple_vector_ref {
    protected :
      V *l;

    public :
      typedef typename linalg_traits<V>::value_type value_type;
      typedef typename linalg_traits<V>::reference_type reference_type;
      simple_vector_ref(V &v) : l(&v) {}
      V &deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      reference_type operator[](size_type i) { return (*l)[i]; }
      value_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_ref<V> > {
    typedef simple_vector_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::sub_row_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type const_sub_row_type;
    typedef typename linalg_traits<V>::sub_col_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_col_type
    const_sub_col_type;
    typedef typename linalg_traits<V>::sub_orientation sub_orientation;
    size_type size(const this_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const this_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    size_type ncols(const this_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    iterator begin(this_type &v)
    { return linalg_traits<V>().begin(v.deref()); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(this_type &v)
    { return linalg_traits<V>().end(v.deref()); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    sub_row_type row(this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    sub_col_type col(this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { clear(v.deref()); }
  };

  template <class V> class simple_vector_const_ref {
    protected :
      const V *l;

    public :
      simple_vector_const_ref(const simple_vector_ref<V> &v)
	: l(&(v.deref())) {}
      simple_vector_const_ref(const V &v) : l(&v) {}
      typedef typename linalg_traits<V>::value_type value_type;
      const V &deref(void) const { return *l; }
      value_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_const_ref<V> > {
    typedef simple_vector_const_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::const_iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::const_sub_row_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type const_sub_row_type;
    typedef typename linalg_traits<V>::const_sub_col_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_col_type const_sub_col_type;
    typedef typename linalg_traits<V>::sub_orientation sub_orientation;
    size_type size(const this_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const this_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    size_type ncols(const this_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    iterator begin(this_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(this_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
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
    typedef typename linalg_traits<V>::const_sub_row_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type const_sub_row_type;
    typedef this_type sub_col_type;
    typedef this_type const_sub_col_type;
    typedef typename linalg_traits<V>::sub_orientation sub_orientation;
    size_type size(const this_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const this_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    size_type ncols(const this_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    iterator begin(this_type &v) {
      return iterator(vect_begin(v.deref()), v.r);
    }
    const_iterator const_begin(const this_type &v) {
      return iterator(vect_begin(v.deref()), v.r);
    }
    iterator end(this_type &v) {
      return iterator(vect_end(v.deref()), v.r);
    }
    const_iterator const_end(const this_type &v) {
      return iterator(vect_end(v.deref()), v.r);
    }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return v; }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };

  /* ********************************************************************* */
  /*		transposed reference                    		   */
  /* ********************************************************************* */
  
  template <class V> class transposed_ref {
    protected :
      V *l;

    public :
      typedef typename linalg_traits<V>::value_type value_type;
      typedef typename linalg_traits<V>::reference_type reference_type;
      transposed_ref(V &v) : l(&v) {}
      V &deref(void) { return *l; }
      const V &deref(void) const { return *l; }
    // value_type &operator[](size_type i) { return (*l)[i]; }
    // value_type operator[](size_type i) const { return (*l)[i]; }
      reference_type operator()(size_type i, size_type j) { return (*l)(j,i); }
      value_type operator()(size_type i, size_type j) const 
      { return (*l)(j,i); }      
  };

  template <class V> struct linalg_traits<transposed_ref<V> > {
    typedef transposed_ref<V> this_type;
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
    size_type size(const this_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const this_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    size_type ncols(const this_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    sub_row_type row(this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    sub_col_type col(this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { linalg_traits<V>().clear(v.deref()); }
  };

  template <class V> class transposed_const_ref {
    protected :
      const V *l;

    public :
      typedef typename linalg_traits<V>::value_type value_type;
      transposed_const_ref(const V &v) : l(&v) {}
      transposed_const_ref(const transposed_ref<V> &v) : l(&(v.deref())) {}
      const V &deref(void) const { return *l; }
    // value_type &operator[](size_type i) { return (*l)[i]; }
    // value_type operator[](size_type i) const { return (*l)[i]; }
      value_type operator()(size_type i, size_type j) const 
      { return (*l)(j,i); }      
  };

  template <class V> struct linalg_traits<transposed_const_ref<V> > {
    typedef transposed_const_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::const_sub_col_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_col_type const_sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_row_type
    const_sub_col_type;
    typedef typename transposed_type<typename 
      linalg_traits<V>::sub_orientation>::t_type sub_orientation;
    size_type size(const this_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const this_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    size_type ncols(const this_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
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
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
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
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &v) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v) { std::fill(v.begin(), v.end(), value_type(0)); }
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
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &v) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v) { std::fill(v.begin(), v.end(), value_type(0)); }
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
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &v) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return v.origin(); }
    void do_clear(this_type &v) { std::fill(v.begin(), v.end(), value_type(0)); }
  };

  /* ********************************************************************* */
  /*		                                         	 	   */
  /*		Traits for bgeot objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  // to be done fsmatrix<T>

  template <class T, int N> struct linalg_traits<fsvector<T, N> > {
    typedef fsvector<T, N> this_type;
    typedef linalg_false is_reference;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<fsvector<T, N> > sub_col_type;
    typedef simple_vector_const_ref<fsvector<T, N> >
    const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &) { return N; }
    size_type nrows(const this_type &) { return N; }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };


  template <class T> struct linalg_traits<vsvector<T> > {
    typedef vsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };

  template <class VECT> struct linalg_traits<PT<VECT> > {
    typedef PT<VECT> this_type;
    typedef linalg_false is_reference;
    typedef typename linalg_traits<VECT>::value_type value_type;
    typedef typename linalg_traits<VECT>::reference_type reference_type;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<VECT>::iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef typename linalg_traits<VECT>::sub_col_type sub_col_type;
    typedef typename linalg_traits<VECT>::sub_row_type const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v) { return linalg_traits<VECT>().begin(v); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<VECT>().const_begin(v); }
    iterator end(this_type &v) { return linalg_traits<VECT>().end(v); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<VECT>().const_end(v); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(this_type &v, size_type i) { return sub_col_type(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.fill(T(0)); }
  };

  template <class T> struct linalg_traits<svector<T> > {
    typedef svector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef ref_elt_svector<T> reference_type;
    typedef svector_iterator<T>  iterator;
    typedef svector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<svector<T> > sub_col_type;
    typedef simple_vector_const_ref<svector<T> > const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return v.size(); }
    size_type nrows(const this_type &v) { return v.size(); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator const_begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator const_end(const this_type &v) { return v.end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_col_type col(this_type &v, size_type i) 
    { return sub_col_type(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<vsmatrix<T> > {
    typedef vsmatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_reg_spaced_with_origin<typename vsmatrix<T>::iterator>
    sub_row_type;
    typedef tab_ref_reg_spaced_with_origin<typename vsmatrix<T>
    ::const_iterator> const_sub_row_type;
    typedef tab_ref_with_origin<typename vsmatrix<T>::iterator> sub_col_type;
    typedef tab_ref_with_origin<typename vsmatrix<T>::const_iterator>
    const_sub_col_type;
    typedef col_and_row sub_orientation;
    size_type size(const this_type &m) { return m.size(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); } 
    const_sub_row_type row(const this_type &m, size_type i)
    {
      typename vsmatrix<T>::const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return const_sub_row_type(b + i, b + i + mat_ncols(m)* nr, nr, &m);
    }
    const_sub_col_type col(const this_type &m, size_type i) {
      typename vsmatrix<T>::const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return const_sub_col_type(b + i * nr, b + (i+1) * nr, &m);
    }
    sub_row_type row(this_type &m, size_type i)
    {
      typename vsmatrix<T>::iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return sub_row_type(b + i, b + i + mat_ncols(m)* nr, nr, &m);
    }
    sub_col_type col(this_type &m, size_type i) { 
      typename vsmatrix<T>::iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return sub_col_type(b + i * nr, b + (i+1) * nr, &m);
    }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<smatrix<T> > {
    typedef smatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef ref_elt_smatrix<T> reference_type;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef simple_vector_ref<svector<T> > sub_row_type;
    typedef simple_vector_const_ref<svector<T> > const_sub_row_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.size(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_row_type row(const this_type &m, size_type i)
    { return simple_vector_const_ref<svector<T> >(m.row(i)); }
    sub_col_type col(this_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_row_type row(this_type &m, size_type i) 
    { return simple_vector_ref<svector<T> >(m.row(i)); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { v.clear(); }
  };


  /* ********************************************************************* */
  /*		sparse sub-vectors                                         */
  /* ********************************************************************* */

  template <class IT> struct sparse_sub_vector_iterator {

    const reverse_index *_r_i;
    IT itb;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef const typename traits_type::value_type *pointer;
    typedef typename traits_type::reference         reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef typename traits_type::size_type         size_type;
    typedef sparse_sub_vector_iterator<IT>          iterator;

    size_type index(void) const { return (*_r_i)[it.index()]; }
    iterator &operator ++() { ++itb; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    reference operator *() const
    {
      static value_type zero(0);
      return (index() == size_type(-1)) ? zero : *it;
    }

    bool operator ==(const iterator &i) const { return itb == i.itb; }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return itb < i.itb; }

    sparse_sub_vector_iterator(void) {}
    sparse_sub_vector_iterator(const IT &it, const reverse_index &ri)
      : itb(it), _r_i(&ri) {}
    
  };

  template <class V, class IT> class sparse_sub_vector {
  protected :
    const reverse_index *_r_i;
    IT begin, end;
    V *v;

  public :
    const reverse_index &rindex(void) const { return *_r_i; }
    V& deref(void) const { return *v; }
    typedef typename linalg_traits<V>::reference_type reference_type;
    reference_type operator[](size_type i) { return (*v)[begin[i]]; }

    sparse_sub_vector(V& w, const IT &b,
		      const IT &e, const reverse_index &ri)
      : _r_i(&ri), begin(b), end(e), v(&w) {}
  };

  template <class IT> struct const_sparse_sub_vector_iterator {

    const reverse_index *_r_i;
    IT itb;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef const typename traits_type::value_type *pointer;
    typedef typename traits_type::value_type        reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef typename traits_type::size_type         size_type;
    typedef sparse_sub_vector_iterator<IT>          iterator;

    size_type index(void) const { return (*_r_i)[it.index()]; }
    iterator &operator ++() { ++itb; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    reference operator *() const
    { return (index() == size_type(-1)) ? value_type(0) : *it; }

    bool operator ==(const iterator &i) const { return itb == i.itb; }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return itb < i.itb; }

    const_sparse_sub_vector_iterator(void) {}
    const_sparse_sub_vector_iterator(const IT &it, const reverse_index &ri)
      : itb(it), _r_i(&ri) {}
    
  };

  template <class V, class IT> class const_sparse_sub_vector {
  protected :
    const reverse_index *_r_i;
    IT begin, end;
    const V *v;

  public :
    const reverse_index &rindex(void) const { return *_r_i; }
    const V& deref(void) const { return *v; }
    typedef typename linalg_traits<V>::value_type value_type;
    value_type operator[](size_type i) const { return (*v)[begin[i]]; }

    const_sparse_sub_vector(const V& w, const IT &b,
			    const IT &e, const reverse_index &ri)
      : _r_i(&ri), begin(b), end(e), v(&w) {}

  };

  template <class V, class IT>
  struct linalg_traits<sparse_sub_vector<V, IT> > {
    typedef sparse_sub_vector<V, IT> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef sparse_sub_vector_iterator<typename linalg_traits<V>::iterator>
    iterator;
    typedef const_sparse_sub_vector_iterator<typename linalg_traits<V>
    ::const_iterator> const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return vect_size(v.deref()); }
    size_type nrows(const this_type &v) { return vect_size(v.deref()); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v)
      { return iterator(vect_begin(v.deref()), v.rbegin()); }
    const_iterator const_begin(const this_type &v)
      { return const_iterator(vect_begin(v.deref()), v.rbegin()); }
    iterator end(this_type &v)
      { return iterator(vect_end(v.deref()), v.rbegin()); }
    const_iterator const_end(const this_type &v)
      { return const_iterator(vect_end(v.deref()), v.rbegin()); }
    const_sub_row_type row(const this_type &, size_type )
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &, size_type )
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_col_type col(this_type &v, size_type i) 
    { return sub_col_type(v); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { std::fill(begin(), end(), value_type(0)); }
  };


  template <class V, class IT>
  struct linalg_traits<const_sparse_sub_vector<V, IT> > {
    typedef const_sparse_sub_vector<V, IT> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::value_type reference_type;
    typedef const_sparse_sub_vector_iterator<typename linalg_traits<V>
    ::const_iterator> iterator;
    typedef const_sparse_sub_vector_iterator<typename linalg_traits<V>
    ::const_iterator> const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_const_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &v) { return vect_size(v.deref()); }
    size_type nrows(const this_type &v) { return vect_size(v.deref()); }
    size_type ncols(const this_type &) { return 1; }
    iterator begin(this_type &v)
      { return iterator(vect_begin(v.deref()), v.rbegin()); }
    const_iterator const_begin(const this_type &v)
      { return const_iterator(vect_begin(v.deref()), v.rbegin()); }
    iterator end(this_type &v)
      { return iterator(vect_end(v.deref()), v.rbegin()); }
    const_iterator const_end(const this_type &v)
      { return const_iterator(vect_end(v.deref()), v.rbegin()); }
    const_sub_row_type row(const this_type &, size_type )
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &, size_type )
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_col_type col(this_type &v, size_type i) 
    { return sub_col_type(v); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { DAL_THROW(failure_error,"impossible"); }
  };

  /* ********************************************************************* */
  /*		sparse sub-matrices                                        */
  /* ********************************************************************* */

  template <class M, class IT1, class IT2> struct sparse_row_sub_matrix {
    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    M *m;

    M &deref(void) const { return *m; }
    const reverse_index &rindex(void) const { return *_r_i; }
    typedef typename linalg_traits<M>::reference_type reference_type;
    reference_type operator()(size_type i, size_type j) const
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    sparse_row_sub_matrix(M &mm, const IT1 &it1, const IT1 &e1,
	         const IT1 &it2, const IT1 &e2, const reverse_index &rindex1,
			  const reverse_index &)
      : _r_i(&rindex1), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}


  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<sparse_row_sub_matrix<M, IT1, IT2> > {
    typedef sparse_row_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_sparse storage_type;
    typedef sparse_sub_vector<typename linalg_traits<M>
    ::sub_row_type, IT1> sub_row_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_row_type, IT1> const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return const_sub_row_type(mat_row(m.deref(), (m.begin2)[i]),
				m.begin1, m.end1, m.rindex());
    }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(mat_row(m.deref(), (m.begin2)[i]),
			  m.begin1, m.end1, m.rindex());
    }
    sub_col_type col(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < nrows(); ++i) clear(row(m,i)); }
  };

  template <class M, class IT1, class IT2> struct const_sparse_row_sub_matrix {

    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    const M *m;

    const reverse_index &rindex(void) const { return *_r_i; }
    const M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::value_type value_type;
    value_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    const_sparse_row_sub_matrix(const M &mm, const IT1 &it1, const IT1 &e1,
		  const IT1 &it2, const IT1 &e2, const reverse_index &rindex1,
				const reverse_index &)
      : _r_i(&rindex1), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}


  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<const_sparse_row_sub_matrix<M, IT1, IT2> > {
    typedef const_sparse_row_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_sparse storage_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_row_type, IT1> sub_row_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_row_type, IT1> const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return const_sub_row_type(mat_row(m.deref(), (m.begin2)[i]),
				m.begin1, m.end1, m.rindex());
    }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(mat_row(m.deref(), (m.begin2)[i]),
			  m.begin1, m.end1, m.rindex());
    }
    sub_col_type col(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { DAL_THROW(failure_error,"impossible"); }
  };

  template <class M, class IT1, class IT2> struct sparse_col_sub_matrix {

    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    M *m;

    const reverse_index &rindex(void) const { return *_r_i; }
    M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::reference_type reference_type;
    reference_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    sparse_col_sub_matrix(M &mm, const IT1 &it1, const IT1 &e1,
		 const IT1 &it2, const IT1 &e2, const reverse_index &,
			  const reverse_index &rindex2)
      : _r_i(&rindex2), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}


  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<sparse_col_sub_matrix<M, IT1, IT2> > {
    typedef sparse_col_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_sparse storage_type;
    typedef sparse_sub_vector<typename linalg_traits<M>
    ::sub_col_type, IT2> sub_col_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_col_type, IT2> const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(mat_col(m.deref(), (m.begin1)[i]),
				m.begin2, m.end2, m.rindex());
    }
    const_sub_row_type row(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_col_type col(this_type &m, size_type i)
    { return sub_col_type(mat_row(m.deref(), (m.begin1)[i]),
			  m.begin2, m.end2, m.rindex());
    }
    sub_row_type row(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < ncols(); ++i) clear(col(m,i)); }
  };

  template <class M, class IT1, class IT2> struct const_sparse_col_sub_matrix {
  
    const reverse_index *_r_i;
    IT1 begin1, end1;
    IT2 begin2, end2;
    const M *m;

    const reverse_index &rindex(void) const { return *_r_i; }
    const M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::value_type value_type;
    value_type operator()(size_type i, size_type j) const
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    const_sparse_col_sub_matrix(const M &mm, const IT1 &it1, const IT1 &e1,
		 const IT1 &it2, const IT1 &e2, const reverse_index &,
				const reverse_index &rindex2)
      : _r_i(&rindex2), begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}

  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<const_sparse_col_sub_matrix<M, IT1, IT2> > {
    typedef const_sparse_col_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_sparse storage_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_col_type, IT2> sub_col_type;
    typedef const_sparse_sub_vector<typename linalg_traits<M>
    ::const_sub_col_type, IT2> const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(mat_col(m.deref(), (m.begin1)[i]),
				m.begin2, m.end2, m.rindex());
    }
    const_sub_row_type row(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_col_type col(this_type &m, size_type i)
    { return sub_col_type(mat_row(m.deref(), (m.begin1)[i]),
			  m.begin2, m.end2, m.rindex());
    }
    sub_row_type row(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m) { DAL_THROW(failure_error,"impossible"); }
  };

  /* ********************************************************************* */
  /*		plain sub-matrices                                         */
  /* ********************************************************************* */

  template <class M, class IT1, class IT2> struct plain_row_sub_matrix {
    IT1 begin1, end1;
    IT2 begin2, end2;
    M *m;

    M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::reference_type reference_type;
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }
    reference_type operator()(size_type i, size_type j) const
    { 
      cout << "access (" << i << "," << j << ") "; 
      cout << "which is (" << begin1[i] << "," << begin2[j]<< ")\n";
      cout << "nrows() = " << nrows() << "ncols() = " << ncols() << endl;
      return (*m)(begin1[i], begin2[j]);
    }
    plain_row_sub_matrix(M &mm, const IT1 &it1, const IT1 &e1,
			 const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}

  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<plain_row_sub_matrix<M, IT1, IT2> > {
    typedef plain_row_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::sub_row_type>::iterator, IT1> sub_row_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_row_type>::const_iterator, IT1>
    const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return
	const_sub_row_type(vect_const_begin(mat_const_row(m.deref(), (m.begin2)[i])),
			   m.begin1, m.end1, linalg_origin(m.deref()));
    }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(vect_begin(mat_row(m.deref(), (m.begin2)[i])),
			  m.begin1, m.end1, linalg_origin(m.deref()));
    }
    sub_col_type col(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < nrows(); ++i) clear(row(m,i)); }
  };

  template <class M, class IT1, class IT2> struct const_plain_row_sub_matrix {

    IT1 begin1, end1;
    IT2 begin2, end2;
    const M *m;

    const M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::value_type value_type;
    value_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    const_plain_row_sub_matrix(const M &mm, const IT1 &it1, const IT1 &e1,
			       const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}

  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<const_plain_row_sub_matrix<M, IT1, IT2> > {
    typedef const_plain_row_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::reference_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_row_type>::const_iterator, IT1> sub_row_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_row_type>::const_iterator, IT1>
    const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i) { 
      return const_sub_row_type(vect_begin(mat_row(m.deref(), (m.begin2)[i])),
				m.begin1, m.end1, linalg_origin(m.deref()));
    }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type row(this_type &m, size_type i)
    { return sub_row_type(vect_begin(mat_row(m.deref(), (m.begin2)[i])),
			  m.begin1, m.end1, linalg_origin(m.deref()));
    }
    sub_col_type col(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { DAL_THROW(failure_error,"impossible"); }
  };

  template <class M, class IT1, class IT2> struct plain_col_sub_matrix {

    IT1 begin1, end1;
    IT2 begin2, end2;
    M *m;

    M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::reference_type reference_type;
    reference_type operator()(size_type i, size_type j) const 
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    plain_col_sub_matrix(M &mm, const IT1 &it1, const IT1 &e1,
			 const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}


  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<plain_col_sub_matrix<M, IT1, IT2> > {
    typedef plain_col_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::sub_col_type>::iterator, IT1> sub_col_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_col_type>::const_iterator, IT1>
    const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
				m.begin2, m.end2, linalg_origin(m.deref()));
    }
    const_sub_row_type row(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type col(this_type &m, size_type i)
    { return sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
			  m.begin2, m.end2, linalg_origin(m.deref()));
    }
    sub_row_type row(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m)
      { for (size_type i = 0; i < ncols(); ++i) clear(col(m,i)); }
  };

  template <class M, class IT1, class IT2> struct const_plain_col_sub_matrix {
  
    IT1 begin1, end1;
    IT2 begin2, end2;
    const M *m;

    const M &deref(void) const { return *m; }
    typedef typename linalg_traits<M>::value_type value_type;
    value_type operator()(size_type i, size_type j) const
    { return (*m)(begin1[i], begin2[j]); }
    size_type nrows(void) const { return end1 - begin1; }
    size_type ncols(void) const { return end2 - begin2; }

    const_plain_col_sub_matrix(const M &mm, const IT1 &it1, const IT1 &e1,
			       const IT1 &it2, const IT1 &e2)
      : begin1(it1), end1(e1), begin2(it2), end2(e2), m(&mm) {}

  };

  template <class M, class IT1, class IT2>
  struct linalg_traits<const_plain_col_sub_matrix<M, IT1, IT2> > {
    typedef const_plain_col_sub_matrix<M, IT1, IT2> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename linalg_traits<M>::value_type reference_type;
    typedef abstract_plain storage_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_col_type>::const_iterator, IT1> sub_col_type;
    typedef tab_ref_index_ref_with_origin<typename linalg_traits<typename
    linalg_traits<M>::const_sub_col_type>::const_iterator, IT1>
    const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type size(const this_type &m) { return m.nrows() * m.ncols(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i) { 
      return const_sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
				m.begin2, m.end2, linalg_origin(m.deref()));
    }
    const_sub_row_type row(const this_type &, size_type)
    { DAL_THROW(failure_error,"inaccessible"); }
    sub_row_type col(this_type &m, size_type i)
    { return sub_col_type(vect_begin(mat_col(m.deref(), (m.begin1)[i])),
			  m.begin2, m.end2, linalg_origin(m.deref()));
    }
    sub_row_type row(this_type &, size_type ) 
    { DAL_THROW(failure_error,"inaccessible");  }
    const void* origin(const this_type &m) { return linalg_origin(m.deref()); }
    void do_clear(this_type &m) { DAL_THROW(failure_error,"impossible"); }
  };

}


#endif //  __BGEOT_LINALG_INTERFACE_H
