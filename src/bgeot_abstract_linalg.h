/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_abstract_linalg.h : generic algorithms on linear       */
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

//
// A faire
//
//   . donner l'accés aux lignes et colonnes dans le cas de la matrice
//     vsmatrix<T>
//   . donner l'accés aux colonnes pour vsmatrix<T> et smatrix<T>
//   . mult : optimisable dans certains cas.
//   . faire scale et scaled sur les matrices aussi
//   . donner une origine correcte pour les dal::tab_ref ...
//   . extraction de sous matrices.
//

// Inspired from M.T.L. (http://www.osl.iu.edu/research/mtl)

#ifndef __BGEOT_ABSTRACT_LINALG_H
#define __BGEOT_ABSTRACT_LINALG_H

#include <dal_ref.h>
#include <bgeot_matrix.h>
#include <bgeot_smatrix.h>

namespace bgeot {

  /* ******************************************************************** */
  /*		Specifier types                             		  */
  /* ******************************************************************** */

  struct abstract_null_type {}; // specify an information lake.

  struct linalg_true {};
  struct linalg_false {};

  struct abstract_vector {};
  struct abstract_matrix {};
  
  struct abstract_sparse {};    // sparse matrix or vector
  struct abstract_plain {};     // plain matrix or vector
  struct abstract_indirect {};  // matrix given by the product with a vector

  struct row_major {};          // matrix with a row access.
  struct col_major {};       // matrix with a column access
  struct row_and_col {};     // both accesses but row preference
  struct col_and_row {};     // both accesses but column preference

  template <class T> struct transposed_type;
  template<> struct transposed_type<row_major>   {typedef col_major   t_type;};
  template<> struct transposed_type<col_major>   {typedef row_major   t_type;};
  template<> struct transposed_type<row_and_col> {typedef col_and_row t_type;};
  template<> struct transposed_type<col_and_row> {typedef row_and_col t_type;};
  

  template <class V> struct linalg_traits;

  template <> struct linalg_traits<abstract_null_type> {
    typedef abstract_null_type is_reference;
    typedef abstract_null_type this_type;
    typedef abstract_null_type linalg_type;
    typedef abstract_null_type base_type;
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
    void clear(this_type &) { }
  };

  /* ******************************************************************** */
  /*		Operations on scalars                         		  */
  /* ******************************************************************** */

  template <class T> T conj_product(T a, T b) { return a * b; }
  // à définir sur les complexes
  template <class T> T modulus(T a) { return dal::abs(a); }
  // à définir sur les complexes

  /* ******************************************************************** */
  /*		Simple references on vectors            		  */
  /* ******************************************************************** */
  
  template <class V> class simple_vector_ref {
    protected :
      V *l;

    public :
      typedef typename linalg_traits<V>::base_type base_type;
      typedef typename linalg_traits<V>::reference_type reference_type;
      simple_vector_ref(V &v) : l(&v) {}
      V &deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      reference_type operator[](size_type i) { return (*l)[i]; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_ref<V> > {
    typedef simple_vector_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
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
    void clear(this_type &v) { linalg_traits<V>().clear(v.deref()); }
  };

  template <class V> class simple_vector_const_ref {
    protected :
      const V *l;

    public :
      simple_vector_const_ref(const simple_vector_ref<V> &v)
	: l(&(v.deref())) {}
      simple_vector_const_ref(const V &v) : l(&v) {}
      typedef typename linalg_traits<V>::base_type base_type;
      const V &deref(void) const { return *l; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_const_ref<V> > {
    typedef simple_vector_const_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
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
    void clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };

  /* ******************************************************************** */
  /*		Scaled references on vectors            		  */
  /* ******************************************************************** */
  
  template <class V> struct scaled_vector_const_ref {

      const V *l;
      typename linalg_traits<V>::base_type r;
      scaled_vector_const_ref(const V &v,
			      typename linalg_traits<V>::base_type x)
	: l(&v), r(x) {}
      typedef typename linalg_traits<V>::base_type base_type;
      const V &deref(void) const { return *l; }
      base_type operator[](size_type i) const { return (*l)[i] * r; }
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
    typedef typename linalg_traits<V>::base_type base_type;
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
    void clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };



  /* ******************************************************************** */
  /*		transposed reference                    		  */
  /* ******************************************************************** */
  
  template <class V> class transposed_ref {
    protected :
      V *l;

    public :
      typedef typename linalg_traits<V>::base_type base_type;
      typedef typename linalg_traits<V>::reference_type reference_type;
      transposed_ref(V &v) : l(&v) {}
      V &deref(void) { return *l; }
      const V &deref(void) const { return *l; }
    // base_type &operator[](size_type i) { return (*l)[i]; }
    // base_type operator[](size_type i) const { return (*l)[i]; }
      reference_type operator()(size_type i, size_type j) { return (*l)(j,i); }
      base_type operator()(size_type i, size_type j) const 
      { return (*l)(j,i); }      
  };

  template <class V> struct linalg_traits<transposed_ref<V> > {
    typedef transposed_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
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
    iterator begin(this_type &v)
    { return linalg_traits<V>().begin(v.deref()); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(this_type &v)
    { return linalg_traits<V>().end(v.deref()); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    sub_row_type row(this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    sub_col_type col(this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void clear(this_type &v) { linalg_traits<V>().clear(v.deref()); }
  };

  template <class V> class transposed_const_ref {
    protected :
      const V *l;

    public :
      typedef typename linalg_traits<V>::base_type base_type;
      transposed_const_ref(const V &v) : l(&v) {}
      transposed_const_ref(const transposed_ref<V> &v) : l(&(v.deref())) {}
      const V &deref(void) const { return *l; }
    // base_type &operator[](size_type i) { return (*l)[i]; }
    // base_type operator[](size_type i) const { return (*l)[i]; }
      base_type operator()(size_type i, size_type j) const 
      { return (*l)(j,i); }      
  };

  template <class V> struct linalg_traits<transposed_const_ref<V> > {
    typedef transposed_const_ref<V> this_type;
    typedef linalg_true is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::const_iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
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
    iterator begin(this_type &v)
    { return linalg_traits<V>().begin(v.deref()); }
    const_iterator const_begin(const this_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(this_type &v)
    { return linalg_traits<V>().end(v.deref()); }
    const_iterator const_end(const this_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const this_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void clear(this_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };

  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for S.T.L. object                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

    template <class T> struct linalg_traits<std::vector<T> > {
    typedef std::vector<T> this_type;
    typedef linalg_false is_reference;
    typedef T base_type;
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
    void clear(this_type &v) { v.fill(T(0)); }
  };

  // to be done :  std::valarray<T> ...

  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for dal objects                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

  template <class IT> struct linalg_traits<dal::tab_ref<IT> > {
    typedef dal::tab_ref<IT> this_type;
    typedef linalg_false is_reference;
    typedef typename std::iterator_traits<IT>::value_type base_type;
    typedef typename std::iterator_traits<IT>::value_type& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    size_type size(const this_type &v) { return v.end() - v.begin(); }
    size_type nrows(const this_type &v) { return v.end() - v.begin(); }
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
    const void* origin(const this_type &v) { return &v; } /* faux ... */
    void clear(this_type &v) { std::fill(v.begin(), v.end(), base_type(0)); }
  };

    template <class IT> struct linalg_traits<dal::tab_ref_reg_spaced<IT> > {
    typedef dal::tab_ref_reg_spaced<IT> this_type;
    typedef linalg_false is_reference;
    typedef typename std::iterator_traits<IT>::value_type base_type;
    typedef typename std::iterator_traits<IT>::value_type& reference_type;
    typedef abstract_vector linalg_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<this_type> sub_col_type;
    typedef simple_vector_const_ref<this_type> const_sub_col_type;
    size_type size(const this_type &v) { return v.end() - v.begin(); }
    size_type nrows(const this_type &v) { return v.end() - v.begin(); }
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
    const void* origin(const this_type &v) { return &v; } /* faux ... */
    void clear(this_type &v) { std::fill(v.begin(), v.end(), base_type(0)); }
  };

  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for bgeot objects                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

  // to be done fsmatrix<T>

  template <class T, int N> struct linalg_traits<fsvector<T, N> > {
    typedef fsvector<T, N> this_type;
    typedef linalg_false is_reference;
    typedef T base_type;
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
    void clear(this_type &v) { v.fill(T(0)); }
  };


  template <class T> struct linalg_traits<vsvector<T> > {
    typedef vsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef T base_type;
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
    void clear(this_type &v) { v.fill(T(0)); }
  };

  template <class VECT> struct linalg_traits<PT<VECT> > {
    typedef PT<VECT> this_type;
    typedef linalg_false is_reference;
    typedef typename linalg_traits<VECT>::base_type base_type;
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
    void clear(this_type &v) { v.fill(T(0)); }
  };

  template <class T> struct linalg_traits<svector<T> > {
    typedef svector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T base_type;
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
    iterator begin(this_type &v) { return v.tas_begin(); }
    const_iterator const_begin(const this_type &v) { return v.tas_begin(); }
    iterator end(this_type &v) { return v.tas_end(); }
    const_iterator const_end(const this_type &v) { return v.tas_end(); }
    const_sub_row_type row(const this_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_col_type col(const this_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(this_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_col_type col(this_type &v, size_type i) 
    { return sub_col_type(v); }
    const void* origin(const this_type &v) { return &v; }
    void clear(this_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<vsmatrix<T> > {
    typedef vsmatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T base_type;
    typedef T& reference_type;
    typedef typename this_type::iterator  iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef dal::tab_ref_reg_spaced<iterator> sub_row_type;
    typedef dal::tab_ref_reg_spaced<const_iterator> const_sub_row_type;
    typedef dal::tab_ref<iterator> sub_col_type;
    typedef dal::tab_ref<const_iterator> const_sub_col_type;
    typedef col_and_row sub_orientation;
    size_type size(const this_type &m) { return m.size(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    iterator begin(this_type &m) { return m.begin(); }
    const_iterator const_begin(const this_type &m) { return m.begin(); }
    iterator end(this_type &m) { return m.end(); }
    const_iterator const_end(const this_type &m) { return m.end(); }
    const_sub_row_type row(const this_type &m, size_type i)
    {
      const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return dal::tab_ref_reg_spaced<const_iterator>
	(b + i, b + i + mat_ncols(m)* nr, nr);
    }
    const_sub_col_type col(const this_type &m, size_type i) {
      const_iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return dal::tab_ref<const_iterator>(b + i * nr, b + (i+1) * nr);
    }
    sub_row_type row(this_type &m, size_type i)
    {
      iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return dal::tab_ref_reg_spaced<iterator>
	(b + i, b + i + mat_ncols(m)* nr, nr);
    }
    sub_col_type col(this_type &m, size_type i) { 
      iterator b = m.begin();
      size_type nr = mat_nrows(m);
      return dal::tab_ref<iterator>(b + i * nr, b + (i+1) * nr);
    }
    const void* origin(const this_type &v) { return &v; }
    void clear(this_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<smatrix<T> > {
    typedef smatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T base_type;
    typedef ref_elt_smatrix<T> reference_type;
    typedef svector_iterator<T>  iterator;
    typedef svector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef simple_vector_ref<svector<T> > sub_row_type;
    typedef simple_vector_const_ref<svector<T> > const_sub_row_type;
    typedef row_major sub_orientation;
    size_type size(const this_type &m) { return m.size(); }
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    iterator begin(this_type &m) { return m.row(0).tas_begin(); }
    const_iterator const_begin(const this_type &m)
    { return m.row(0).tas_begin(); }
    iterator end(this_type &m) { return m.row(0).tas_end(); }
    const_iterator const_end(const this_type &m)
    { return m.row(0).tas_end(); }
    const_sub_col_type col(const this_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_row_type row(const this_type &m, size_type i)
    { return simple_vector_const_ref<svector<T> >(m.row(i)); }
    sub_col_type col(this_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_row_type row(this_type &m, size_type i) 
    { return simple_vector_ref<svector<T> >(m.row(i)); }
    const void* origin(const this_type &v) { return &v; }
    void clear(this_type &v) { v.clear(); }
  };


  /* ******************************************************************** */
  /*		                                         		  */
  /*		Generic algorithm                           		  */
  /*		                                         		  */
  /* ******************************************************************** */

  template <class V> inline size_type vect_size(const V &v)
  { return linalg_traits<V>().size(v); }

  template <class MAT> inline size_type mat_nrows(const MAT &m)
  { return linalg_traits<MAT>().nrows(m); }

  template <class MAT> inline size_type mat_ncols(const MAT &m)
  { return linalg_traits<MAT>().ncols(m); }

  template <class L>
  inline const void *linalg_origin(const L &l)
  { return linalg_traits<L>().origin(l); }  

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_begin(const V &v)
  { return linalg_traits<V>().const_begin(v); }

  template <class V>
  inline typename linalg_traits<V>::const_iterator vect_end(const V &v)
  { return linalg_traits<V>().const_end(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_begin(V &v)
  { return linalg_traits<V>().begin(v); }

  template <class V>
  inline typename linalg_traits<V>::iterator vect_end(V &v)
  { return linalg_traits<V>().end(v); }

  template <class MAT> inline 
    typename linalg_traits<MAT>::const_sub_row_type
    mat_row(const MAT &m, size_type i)
  { return linalg_traits<MAT>().row(m, i); }

  template <class MAT> inline  
    typename linalg_traits<MAT>::const_sub_col_type
    mat_col(const MAT &m, size_type i)
  { return linalg_traits<MAT>().col(m, i); }

  template <class MAT> inline 
    typename linalg_traits<MAT>::sub_row_type
    mat_row(MAT &m, size_type i)
  { return linalg_traits<MAT>().row(m, i); }

  template <class MAT> inline  
    typename linalg_traits<MAT>::sub_col_type
    mat_col(MAT &m, size_type i)
  { return linalg_traits<MAT>().col(m, i); }

  template <class L> inline void clear(L &l)
  { return linalg_traits<L>().clear(l); }

  template <class L> inline
    scaled_vector_const_ref<L> scaled(const L &l,
				      typename linalg_traits<L>::base_type x)
  { return scaled_vector_const_ref<L>(l, x); }
  

  template <class L> inline transposed_ref<L> transposed(L &l)
  { return transposed_ref<L>(l); }

  template <class L> inline transposed_const_ref<L> transposed(const L &l)
  { return transposed_const_ref<L>(l); }

  /* ******************************************************************** */
  /*		Scalar product                             		  */
  /* ******************************************************************** */

  template <class V1, class V2>
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    return vect_sp(v1, v2,
		   typename linalg_traits<V1>::storage_type(), 
		   typename linalg_traits<V2>::storage_type());
  }


  template <class IT1, class IT2>
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_plain(IT1 it, IT1 ite, IT2 it2) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it, ++it2) res += conj_product(*it, *it2);
    return res;
  }

  template <class IT1, class V>
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_sparse(IT1 it, IT1 ite, const V &v) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it) res += conj_product(*it, v[it.index()]);
    return res;
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2, abstract_plain, abstract_plain) {
    return _vect_sp_plain(vect_begin(v1), vect_end(v1), vect_begin(v2));
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_plain) {
    return _vect_sp_sparse(vect_begin(v1), vect_end(v1), v2);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2, abstract_plain,abstract_sparse) {
    return _vect_sp_sparse(vect_begin(v2), vect_end(v2), v1);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_sparse) {
    return _vect_sp_sparse(vect_begin(v1), vect_end(v1), v2);
  }

  /* ******************************************************************** */
  /*		Euclidian norm                             		  */
  /* ******************************************************************** */

   template <class V>
    typename linalg_traits<V>::base_type norm2(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res += dal::sqr(modulus(*it));
    return sqrt(res);
  }

  /* ******************************************************************** */
  /*		Inifity norm                              		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::base_type norminf(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res = std::max(res, modulus(*it));
    return res;
  }
  
  /* ******************************************************************** */
  /*		norm1                                    		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::base_type norm1(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = vect_begin(v), ite = vect_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res += modulus(*it);
    return res;
  }

  /* ******************************************************************** */
  /*		Copy                                    		  */
  /* ******************************************************************** */

  template <class L1, class L2> inline
  void copy(const L1& l1, L2& l2) { 
    #ifdef __GETFEM_VERIFY
      if (linalg_origin(l1) == linalg_origin(l2))
	cerr << "Warning : a conflict is possible in vector copy\n";
    #endif
    if ((const void *)(&l1) != (const void *)(&l2))
      copy(l1, l2, typename linalg_traits<L1>::linalg_type(),
	   typename linalg_traits<L2>::linalg_type());
  }

  template <class L1, class L2> inline
  void copy(const L1& l1, const L2& l2) {
    copy_ref(l1, l2, typename linalg_traits<L2>::is_reference());
  }

  template <class L1, class L2> inline
  void copy_ref(const L1& l1, const L2& l2, linalg_true) {
    L2 temp(l2); copy(l1, temp);
  }

  template <class L1, class L2>
  void copy(const L1& l1, L2& l2, abstract_vector, abstract_vector) {
    if (vect_size(l1) != vect_size(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    copy_vect(l1, l2, typename linalg_traits<L1>::storage_type(),
	      typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2>
  void copy(const L1& l1, L2& l2, abstract_matrix, abstract_matrix) {
    if (mat_ncols(l1) != mat_ncols(l2) || mat_nrows(l1) != mat_nrows(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    copy_mat(l1, l2, typename linalg_traits<L1>::sub_orientation(),
	     typename linalg_traits<L2>::sub_orientation());
  }

  template <class L1, class L2>
  void copy_mat_by_row(const L1& l1, L2& l2)
  {
    size_type nbr = mat_nrows(l1);
    for (size_type i = 0; i < nbr; ++i) {
      typename linalg_traits<L2>::sub_row_type r = mat_row(l2, i);
      copy_vect(mat_row(l1, i), r,
		typename linalg_traits<typename linalg_traits<L1>
		 ::const_sub_row_type>::storage_type(),
		typename linalg_traits<typename linalg_traits<L2>
		 ::sub_row_type>::storage_type());
    }
  }

  template <class L1, class L2>
  void copy_mat_by_col(const L1 &l1, L2 &l2) {
    size_type nbc = mat_ncols(l1);
    for (size_type i = 0; i < nbc; ++i) {
      typename linalg_traits<L2>::sub_col_type c = mat_col(l2, i);
      copy_vect(mat_col(l1, i), c,
		typename linalg_traits<typename linalg_traits<L1>
		::const_sub_col_type>::storage_type(),
		typename linalg_traits<typename linalg_traits<L2>
		::sub_col_type>::storage_type());
    }
  }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, row_major)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_major, col_and_row)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, row_and_col)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, col_and_row)
  { copy_mat_by_row(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, col_and_row)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, row_and_col)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, row_and_col, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, col_major)
  { copy_mat_by_col(l1, l2); }

  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_and_row, col_and_row)
  { copy_mat_by_col(l1, l2); }
  
  template <class L1, class L2> inline
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i) {
    copy_mat_mixed_rc(l1, l2, i, typename linalg_traits<L1>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it) l2(i, it.index()) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(i, j) = *it;
  }

  template <class L1, class L2> inline
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i) {
    copy_mat_mixed_cr(l1, l2, i, typename linalg_traits<L1>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it) l2(it.index(), i) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(j, i) = *it;
  }

  template <class L1, class L2>
  void copy_mat(const L1& l1, L2& l2, row_major, col_major) {
    clear(l2);
    size_type nbr = mat_nrows(l1);
    for (size_type i = 0; i < nbr; ++i)
      copy_mat_mixed_rc(mat_row(l1, i), l2, i);
  }
  
  template <class L1, class L2>
  void copy_mat(const L1& l1, L2& l2, col_major, row_major) {
    clear(l2);
    size_type nbc = mat_ncols(l1);
    for (size_type i = 0; i < nbc; ++i)
      copy_mat_mixed_cr(mat_col(l1, i), l2, i);
  }
  
  template <class L1, class L2> inline
  void copy_vect(const L1 &l1, L2 &l2, abstract_plain, abstract_plain) {
    std::copy(vect_begin(l1), vect_end(l1), vect_begin(l2));
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_plain) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it) l2[it.index()] = *it;
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (; it != ite; ++it)
      if (*it != (typename linalg_traits<L1>::base_type)(0))
	l2[it.index()] = *it;
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_plain, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = vect_begin(l1), ite = vect_end(l1);
    for (size_type i = 0; it != ite; ++it, ++i)
      if (*it != (typename linalg_traits<L1>::base_type)(0))
	l2[i] = *it;
  }

  /* ******************************************************************** */
  /*		Vector Addition                                    	  */
  /*   algorithms are build in order to avoid some conflicts whith        */
  /*   repeated arguments or with overlapping part of a same object.      */
  /*   In the latter case, conflicts are still possible.                  */
  /* ******************************************************************** */
  
  template <class L1, class L2> inline
    void add(const L1& l1, L2& l2) {
      add_spec(l1, l2, typename linalg_traits<L2>::linalg_type());
  }

  template <class L1, class L2> inline
  void add(const L1& l1, const L2& l2) {
    add_ref(l1, l2, typename linalg_traits<L2>::is_reference());
  }

  template <class L1, class L2> inline
  void add_ref(const L1& l1, const L2& l2, linalg_true) {
    L2 temp(l2); add(l1, temp);
  }

  template <class L1, class L2>
    void add_spec(const L1& l1, L2& l2, abstract_vector) {
    if (vect_size(l1) != vect_size(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (mat_ncols(l1) != 1 || mat_ncols(l2) != 1)
      DAL_THROW(to_be_done_error,"to be done.");
    add(l1, l2, typename linalg_traits<L1>::storage_type(),
	typename linalg_traits<L2>::storage_type());
  }

  template <class L1, class L2, class L3>
    void add(const L1& l1, const L2& l2, L3& l3) {
    if (vect_size(l1) != vect_size(l2) || vect_size(l1) != vect_size(l3))
      DAL_THROW(dimension_error,"dimensions mismatch"); 
    if (mat_ncols(l1) != 1 || mat_ncols(l2) != 1 || mat_ncols(l3) != 1)
      DAL_THROW(to_be_done_error,"to be done.");
    if ((const void *)(&l1) == (const void *)(&l3))
      add(l2, l3);
    else if ((const void *)(&l2) == (const void *)(&l3))
      add(l1, l3);
    else
      add(l1, l2, l3, typename linalg_traits<L1>::storage_type(),
	  typename linalg_traits<L2>::storage_type(),
	  typename linalg_traits<L3>::storage_type());
  }

  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, const L3& l3) {
    add_ref(l1, l2, l3, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3> inline
  void add_ref(const L1& l1, const L2& l2, const L3& l3, linalg_true) {
    L3 temp(l3); add(l1, l2, temp);
  }

  template <class IT1, class IT2, class IT3>
    void _add_full(IT1 it1, IT2 it2, IT3 it3, IT3 ite) {
    for (; it3 != ite; ++it3, ++it2, ++it1) *it3 = *it1 + *it2;
  }

  template <class IT1, class IT2, class IT3>
    void _add_almost_full(IT1 it1, IT1 ite1, IT2 it2, IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it, ++it2) *it = *it2;
    for (; it1 != ite1; ++it1) *(it3 + it1.index()) += *it1;
  }

  template <class IT1, class IT2, class IT3>
  void _add_to_full(IT1 it1, IT1 ite1, IT2 it2, IT2 ite2,
		    IT3 it3, IT3 ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it) *it = 0;
    for (; it1 != ite1; ++it1) *(it3 + it1.index()) = *it1;
    for (; it2 != ite2; ++it2) *(it3 + it2.index()) += *it2;    
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_plain, abstract_plain, abstract_plain) {
    _add_full(vect_begin(l1), vect_begin(l2),
	      vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_plain, abstract_plain) {
    _add_almost_full(vect_begin(l1), vect_end(l1), vect_begin(l2),
		     vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_plain, abstract_sparse, abstract_plain) {
    _add_almost_full(vect_begin(l2), vect_end(l2), vect_begin(l1),
		     vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_plain) {
    _add_to_full(vect_begin(l1), vect_end(l1),
		 vect_begin(l2), vect_end(l2),
		 vect_begin(l3), vect_end(l3));
  }
  
  template <class L1, class L2, class L3>
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    typename linalg_traits<L2>::const_iterator
      it2 = vect_begin(l2), ite2 = vect_end(l2);
    clear(l3);
    while (it1 != ite1 && it2 != ite2) {
      ptrdiff_t d = it1.index() - it2.index();
      if (d < 0)
	{ l3[it1.index()] += *it1; ++it1; }
      else if (d > 0)
	{ l3[it2.index()] += *it2; ++it2; }
      else
	{ l3[it1.index()] = *it1 + *it2; ++it1; ++it2; }
    }
    for (; it1 != ite1; ++it1) l3[it1.index()] += *it1;
    for (; it2 != ite2; ++it2) l3[it2.index()] += *it2;   
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_plain, abstract_sparse, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_plain, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  template <class L1, class L2, class L3> inline
  void add(const L1&, const L2&, L3&,
	   abstract_sparse, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); }
  
  template <class L1, class L2>
  void add(const L1& l1, L2& l2,
	   abstract_plain, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1); 
    typename linalg_traits<L2>::iterator
      it2 = vect_begin(l2), ite = vect_end(l2);
    for (; it2 != ite; ++it2, ++it1) *it2 += *it1;
  }
  
  template <class L1, class L2>
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    for (; it1 != ite1; ++it1) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2>
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = vect_begin(l1), ite1 = vect_end(l1);
    for (; it1 != ite1; ++it1) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1&, L2&, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); } 

  /* ******************************************************************** */
  /*		scale                                    	          */
  /* ******************************************************************** */

  template <class L> void scale(L& l, typename linalg_traits<L>::base_type a)
  {
    typename linalg_traits<L>::iterator it = vect_begin(l), ite = vect_end(l);
    for ( ; it != ite; ++it) *it *= a;
  }

  template <class L> inline
  void scale(const L& l, typename linalg_traits<L>::base_type a)
  { scale_const(l, a, typename linalg_traits<L>::is_reference()); }

  template <class L> inline
  void scale_const(const L& l, linalg_true)
  { L temp(l); scale_const(temp, a); }


  /* ******************************************************************** */
  /*		Matrix-vector mult                                    	  */
  /* ******************************************************************** */

  template <class L1, class L2, class L3>
  void mult(const L1& l1, const L2& l2, L3& l3) {
    if (mat_ncols(l1) != vect_size(l2) || mat_nrows(l1) != vect_size(l3))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (linalg_origin(l2) != linalg_origin(l3))
      mult_spec(l1, l2, l3, typename linalg_traits<L1>::sub_orientation());
    else {
      #ifdef __GETFEM_VERIFY
        cerr << "Warning, A temporary is used for mult\n";
      #endif
      L3 temp(vect_size(l3));
      mult_spec(l1, l2, temp, typename linalg_traits<L1>::sub_orientation());
      copy(temp, l3);
    }
  }

  template <class L1, class L2, class L3> inline
  void mult(const L1& l1, const L2& l2, const L3& l3) {
    mult_const(l1, l2, l3, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3> inline
  void mult_const(const L1& l1, const L2& l2, const L3& l3, linalg_true)
  { L3 temp(l3); mult(l1, l2, temp); }

  template <class L1, class L2, class L3>
  void mult_by_row(const L1& l1, const L2& l2, L3& l3, abstract_sparse) {
    clear(l3);
    size_type nr = mat_nrows(l1);
    for (size_type i = 0; i < nr; ++i) {
      typename linalg_traits<L1>::base_type aux = vect_sp(mat_row(l1, i), l2);
      if (aux != 0) l3[i] = aux;
    }
  }

  template <class L1, class L2, class L3>
  void mult_by_row(const L1& l1, const L2& l2, L3& l3, abstract_plain) {
    typename linalg_traits<L3>::iterator
      it = vect_begin(l3), ite = vect_end(l3);
    for (size_type i = 0; it != ite; ++it, ++i)
      *it = vect_sp(mat_row(l1, i), l2);
  }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, row_major)
  { mult_by_row(l1, l2, l3, typename linalg_traits<L3>::storage_type()); }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, row_and_col)
  { mult_by_row(l1, l2, l3, typename linalg_traits<L3>::storage_type()); }
  
  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, col_and_row) 
  { mult_by_row(l1, l2, l3, typename linalg_traits<L3>::storage_type()); }
  

  template <class L1, class L2, class L3>
  void mult_spec(const L1& l1, const L2& l2, L3& l3, col_major) {
    clear(l3);
    size_type nc = mat_ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l3);
  }

  template <class L1, class L2, class L3> inline
  void mult_spec(const L1& l1, const L2& l2, L3& l3, abstract_null_type)
  { mult_ind(l1, l2, l3, typename linalg_traits<L1>::storage_type()); }

  template <class L1, class L2, class L3>
  void mult_ind(const L1& l1, const L2& l2, L3& l3, abstract_indirect) {
    DAL_THROW(failure_error,
	  "You have to define the mult(m, v1, v2) for this kind of matrix");
  }

  template <class L1, class L2, class L3, class L4>
  void mult(const L1& l1, const L2& l2, const L3& l3, L4& l4) {
    if (mat_ncols(l1) != vect_size(l2) || mat_nrows(l1) != vect_size(l3)
	|| mat_nrows(l1) != vect_size(l4))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (linalg_origin(l2) != linalg_origin(l4))
      mult_spec(l1, l2, l3, l4, typename linalg_traits<L1>::sub_orientation());
    else {
      #ifdef __GETFEM_VERIFY
        cerr << "Warning, A temporary is used for mult\n";
      #endif
      L3 temp(vect_size(l3));
      mult_spec(l1,l2,l3, temp, typename linalg_traits<L1>::sub_orientation());
      copy(temp, l4);
    }
  }
  
  template <class L1, class L2, class L3, class L4> inline
  void mult(const L1& l1, const L2& l2, const L3& l3, const L4& l4) {
    mult_const(l1, l2, l3, l4, typename linalg_traits<L3>::is_reference());
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_const(const L1& l1, const L2& l2, const L3& l3,
		  const L4& l4, linalg_true)
  { L4 temp(l4); mult(l1, l2, l3, temp); }

  template <class L1, class L2, class L3, class L4>
  void mult_by_row(const L1& l1, const L2& l2, const L3& l3,
		   L4& l4, abstract_sparse) {
    if ((const void *)(&l3) != (const void *)(&l4)) copy(l3, l4);
    size_type nr = nrows(l1);
    for (size_type i = 0; i < nr; ++i) {
      typename linalg_traits<L1>::base_type
	aux = vect_sp(mat_row(l1, i), l2);
      if (aux != 0) l4[i] += aux;
    }
  }

  template <class L1, class L2, class L3, class L4>
  void mult_by_row(const L1& l1, const L2& l2, const L3& l3, L4& l4,
		   abstract_plain) {
    if ((const void *)(&l3) != (const void *)(&l4)) copy(l3, l4);
    typename linalg_traits<L4>::iterator
      it = vect_begin(l4), ite = vect_end(l4);
    for (size_type i = 0; it != ite; ++it, ++i)
      *it += vect_sp(mat_row(l1, i), l2);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, row_major)
  { mult_by_row(l1, l2, l3, l4, typename linalg_traits<L4>::storage_type()); }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, row_and_col)
  { mult_by_row(l1, l2, l3, l4, typename linalg_traits<L4>::storage_type()); }
  
  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, col_and_row)
  { mult_by_row(l1, l2, l3, l4, typename linalg_traits<L4>::storage_type()); }
  

  template <class L1, class L2, class L3, class L4>
  void mult_spec(const L1& l1, const L2& l2, const L3& l3, L4& l4, col_major) {
    copy(l3, l4);
    clear(l4);
    size_type nc = ncols(l1);
    for (size_type i = 0; i < nc; ++i)
      add(scaled(mat_col(l1, i), l2[i]), l4);
  }

  template <class L1, class L2, class L3, class L4> inline
  void mult_spec(const L1& l1, const L2& l2, const L3& l3,
		 L4& l4, abstract_null_type)
  { mult_ind(l1, l2, l3, l4, typename linalg_traits<L1>::storage_type()); }

  template <class L1, class L2, class L3, class L4>
  void mult_ind(const L1& l1, const L2& l2, const L3& l3,
		L4& l4, abstract_indirect) {
    DAL_THROW(failure_error,
	  "You have to define the mult(m, v1, v2) for this kind of matrix");
  }

}


#endif //  __BGEOT_ABSTRACT_LINALG_H
