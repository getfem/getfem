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
//   . un programme de test qui mélange un peu tout les cas
//


#ifndef __BGEOT_ABSTRACT_LINALG_H
#define __BGEOT_ABSTRACT_LINALG_H

#include <bgeot_matrix.h>
#include <bgeot_smatrix.h>

namespace bgeot {

  /* ******************************************************************** */
  /*		Specifier types                             		  */
  /* ******************************************************************** */

  struct abstract_null_type {}; // specify an information lake.
  
  struct abstract_sparse {};    // sparse matrix or vector
  struct abstract_plain {};     // plain matrix or vector
  struct abstract_indirect {};  // matrix given by the product with a vector

  struct row_major {};          // matrix with a row access.
  struct col_major {};       // matrix with a column access
  struct row_and_col {};     // both accesses but row preference
  struct col_and_row {};     // both accesses but column preference

  template <class V> struct linalg_traits;

  template <> struct linalg_traits<abstract_null_type> {
    typedef abstract_null_type linalg_type;
    typedef abstract_null_type base_type;
    typedef void * iterator;
    typedef const void * const_iterator;
    typedef abstract_null_type storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type sub_orientation;
    size_type size(const linalg_type &) { return 0; }
    size_type nrows(const linalg_type &) { return 0; }
    size_type ncols(const linalg_type &) { return 0; }
    iterator begin(linalg_type &) { return 0; }
    const_iterator const_begin(const linalg_type &) { return 0; }
    iterator end(linalg_type &) { return 0; }
    const_iterator const_end(const linalg_type &) { return 0; }
    const_sub_row_type row(const linalg_type &, size_type)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const linalg_type &, size_type)
      { DAL_THROW(failure_error,"Columns inaccessible for this object"); }
    sub_row_type row(linalg_type &, size_type)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(linalg_type &, size_type) 
      { DAL_THROW(failure_error,"Columns inaccessible for this object"); }
    void clear(linalg_type &) { }
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
      simple_vector_ref(V &v) : l(&v) {}
      V &deref(void) { return *l; }
      const V &deref(void) const { return *l; }
      base_type &operator[](size_type i) { return (*l)[i]; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_ref<V> > {
    typedef simple_vector_ref<V> linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
    typedef typename linalg_traits<V>::iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::sub_row_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type const_sub_row_type;
    typedef typename linalg_traits<V>::sub_col_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_col_type
    const_sub_col_type;
    typedef typename linalg_traits<V>::sub_orientation sub_orientation;
    size_type size(const linalg_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const linalg_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    size_type ncols(const linalg_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    iterator begin(linalg_type &v)
    { return linalg_traits<V>().begin(v.deref()); }
    const_iterator const_begin(const linalg_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(linalg_type &v)
    { return linalg_traits<V>().end(v.deref()); }
    const_iterator const_end(const linalg_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const linalg_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const_sub_col_type col(const linalg_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    sub_row_type row(linalg_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    sub_col_type col(linalg_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    void clear(linalg_type &v) { linalg_traits<V>().clear(v.deref()); }
  };

  template <class V> class simple_vector_const_ref {
    protected :
      const V *l;

    public :
      simple_vector_const_ref(simple_vector_ref<V> &v) : l(&(v.deref())) {}
      simple_vector_const_ref(const V &v) : l(&v) {}
      typedef typename linalg_traits<V>::base_type base_type;
      const V &deref(void) const { return *l; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class V> struct linalg_traits<simple_vector_const_ref<V> > {
    typedef simple_vector_const_ref<V> linalg_type;
    typedef typename linalg_traits<V>::base_type base_type;
    typedef typename linalg_traits<V>::const_iterator  iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::const_sub_row_type sub_row_type;
    typedef typename linalg_traits<V>::const_sub_row_type const_sub_row_type;
    typedef typename linalg_traits<V>::const_sub_col_type sub_col_type;
    typedef typename linalg_traits<V>::const_sub_col_type const_sub_col_type;
    typedef typename linalg_traits<V>::sub_orientation sub_orientation;
    size_type size(const linalg_type &v)
    { return linalg_traits<V>().size(v.deref()); }
    size_type nrows(const linalg_type &v)
    { return linalg_traits<V>().nrows(v.deref()); }
    size_type ncols(const linalg_type &v)
    { return linalg_traits<V>().ncols(v.deref()); }
    iterator begin(linalg_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    const_iterator const_begin(const linalg_type &v)
    { return linalg_traits<V>().const_begin(v.deref()); }
    iterator end(linalg_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_iterator const_end(const linalg_type &v)
    { return linalg_traits<V>().const_end(v.deref()); }
    const_sub_row_type row(const linalg_type &v, size_type i)
    { return linalg_traits<V>().row(v.deref(), i); }
    const_sub_col_type col(const linalg_type &v, size_type i)
    { return linalg_traits<V>().col(v.deref(), i); }
    void clear(linalg_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
  };

  /* ******************************************************************** */
  /*		standard extern references for plain vectors              */
  /* ******************************************************************** */
  
  template <class V, class T> class extern_plain_vector_ref {
    protected :
    V *l;
    
    public :
    typedef T base_type;
    V &deref(void) { return *l; }
    const V &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_plain_vector_ref(V &v) : l(&v) {}
  };
  
  template <class V, class T> class extern_plain_vector_const_ref {
    protected :
    const V *l;
    
    public :
    typedef T base_type;
    const V &deref(void) const { return *l; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_plain_vector_const_ref(const V &v) : l(&v) {}
    extern_plain_vector_const_ref(const extern_plain_vector_ref<V, T> &v)
    { l = &(v.deref()); }
  };

  template <class V, class T>
    struct linalg_traits<extern_plain_vector_ref<V, T> > {
      typedef extern_plain_vector_ref<V, T> linalg_type;
      typedef T base_type;
      typedef typename V::iterator  iterator;
      typedef typename V::const_iterator const_iterator;
      typedef abstract_plain storage_type;
      typedef T sub_row_type;
      typedef T const_sub_row_type;
      typedef linalg_type sub_col_type;
      typedef extern_plain_vector_const_ref<V, T> const_sub_col_type;
      typedef col_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      const_iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      const_iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_row_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_col_type col(const linalg_type &v, size_type i)
      { return v; }
      sub_row_type row(linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      sub_col_type col(linalg_type &v, size_type i) { return v; }
      void clear(linalg_type &v)
      { std::fill(v.deref().begin(), v.deref().end(), T(0)); }
    };

    template <class V, class T>
    struct linalg_traits<extern_plain_vector_const_ref<V, T> > {
      typedef extern_plain_vector_const_ref<V, T> linalg_type;
      typedef T base_type;
      typedef typename V::iterator  iterator;
      typedef typename V::const_iterator const_iterator;
      typedef abstract_plain storage_type;
      typedef T sub_row_type;
      typedef T const_sub_row_type;
      typedef linalg_type sub_col_type;
      typedef linalg_type const_sub_col_type;
      typedef col_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      const_iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      const_iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_col_type col(const linalg_type &v, size_type i)
	{ return v; }
      const_sub_row_type row(const linalg_type &, size_type)
	{ DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      void clear(linalg_type &)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
    };

  /* ******************************************************************** */
  /*		standard extern references for sparse vectors             */
  /* ******************************************************************** */
  
  template <class V, class T> class extern_sparse_vector_ref {
    protected :
    V *l;
    
    public :
    typedef T base_type;
    V &deref(void) { return *l; }
    const V &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_sparse_vector_ref(V &v) : l(&v) {}
  };

  template <class V, class T> class extern_sparse_vector_const_ref {
    protected :
    const V *l;
    
    public :
    typedef T base_type;
    V &deref(void) { return *l; }
    const V &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_sparse_vector_const_ref(const V &v) : l(&v) {}
    extern_sparse_vector_const_ref(const extern_sparse_vector_ref<V, T> &v)
    { l = &(v.deref()); }
  };

  template <class V, class T>
    struct linalg_traits<extern_sparse_vector_ref<V, T> > {
      typedef extern_sparse_vector_ref<V, T> linalg_type;
      typedef T base_type;
      typedef typename V::iterator  iterator;
      typedef typename V::const_iterator const_iterator;
      typedef abstract_sparse storage_type;
      typedef T sub_row_type;
      typedef T const_sub_row_type;
      typedef linalg_type sub_col_type;
      typedef extern_sparse_vector_const_ref<V, T> const_sub_col_type;
      typedef col_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      const_iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      const_iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_row_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_col_type col(const linalg_type &v, size_type i)
	{ return v; }
      sub_row_type row(linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      sub_col_type col(linalg_type &v, size_type i) { return v; }
      void clear(linalg_type &v) { v.deref().clear(); }
    };

  template <class V, class T>
    struct linalg_traits<extern_sparse_vector_const_ref<V, T> > {
      typedef extern_sparse_vector_const_ref<V, T> linalg_type;
      typedef T base_type;
      typedef typename V::iterator  iterator;
      typedef typename V::const_iterator const_iterator;
      typedef abstract_sparse storage_type;
      typedef T sub_row_type;
      typedef T const_sub_row_type;
      typedef linalg_type sub_col_type;
      typedef linalg_type const_sub_col_type;
      typedef col_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      const_iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      const_iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_row_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_col_type col(const linalg_type &v, size_type i) { return v; }
      void clear(linalg_type &v)
      { DAL_THROW(failure_error,"Impossible to clear a constant object");}
    };


  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for S.T.L. object                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

  // to be done : std::vector<T> and std::valarray<T> ...

  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for bgeot objects                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

  // to be done fsmatrix<T>

  template <class T, int N> struct linalg_traits<fsvector<T, N> > {
    typedef fsvector<T, N> linalg_type;
    typedef T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<fsvector<T, N> > sub_col_type;
    typedef simple_vector_const_ref<fsvector<T, N> >
    const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const linalg_type &v) { return N; }
    size_type nrows(const linalg_type &v) { return N; }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.begin(); }
    const_iterator const_begin(const linalg_type &v) { return v.begin(); }
    iterator end(linalg_type &v) { return v.end(); }
    const_iterator const_end(const linalg_type &v) { return v.end(); }
    const_sub_row_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const linalg_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(linalg_type &v, size_type i) { return sub_col_type(v); }
    void clear(linalg_type &v) { v.fill(T(0)); }
  };


  template <class T> struct linalg_traits<vsvector<T> > {
    typedef vsvector<T> linalg_type;
    typedef T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<linalg_type> sub_col_type;
    typedef simple_vector_const_ref<linalg_type> const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const linalg_type &v) { return v.size(); }
    size_type nrows(const linalg_type &v) { return v.size(); }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.begin(); }
    const_iterator const_begin(const linalg_type &v) { return v.begin(); }
    iterator end(linalg_type &v) { return v.end(); }
    const_iterator const_end(const linalg_type &v) { return v.end(); }
    const_sub_row_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const linalg_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(linalg_type &v, size_type i) { return sub_col_type(v); }
    void clear(linalg_type &v) { v.fill(T(0)); }
  };

  template <class T> struct linalg_traits<svector<T> > {
    typedef svector<T> linalg_type;
    typedef T base_type;
    typedef svector_iterator<T>  iterator;
    typedef svector_const_iterator<T> const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<svector<T> > sub_col_type;
    typedef simple_vector_const_ref<svector<T> > const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const linalg_type &v) { return v.size(); }
    size_type nrows(const linalg_type &v) { return v.size(); }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.tas_begin(); }
    const_iterator const_begin(const linalg_type &v) { return v.tas_begin(); }
    iterator end(linalg_type &v) { return v.tas_end(); }
    const_iterator const_end(const linalg_type &v) { return v.tas_end(); }
    const_sub_row_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_col_type col(const linalg_type &v, size_type i)
    { return const_sub_col_type(v); }
    sub_row_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_col_type col(linalg_type &v, size_type i) 
    { return sub_col_type(v); }
    void clear(linalg_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<vsmatrix<T> > {
    typedef vsmatrix<T> linalg_type;
    typedef T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<vsvector<T> > sub_col_type;
    typedef simple_vector_const_ref<vsvector<T> > const_sub_col_type;
    typedef col_major sub_orientation;
    size_type size(const linalg_type &m) { return m.size(); }
    size_type nrows(const linalg_type &m) { return m.nrows(); }
    size_type ncols(const linalg_type &m) { return m.ncols(); }
    iterator begin(linalg_type &m) { return m.begin(); }
    const_iterator const_begin(const linalg_type &m) { return m.begin(); }
    iterator end(linalg_type &m) { return m.end(); }
    const_iterator const_end(const linalg_type &m) { return m.end(); }
    const_sub_row_type row(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_row_type row(linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(linalg_type &m, size_type i) 
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    void clear(linalg_type &v) { v.clear(); }
  };

  template <class T> struct linalg_traits<smatrix<T> > {
    typedef smatrix<T> linalg_type;
    typedef T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef simple_vector_ref<svector<T> > sub_col_type;
    typedef simple_vector_const_ref<svector<T> > const_sub_col_type;
    typedef row_major sub_orientation;
    size_type size(const linalg_type &m) { return m.size(); }
    size_type nrows(const linalg_type &m) { return m.nrows(); }
    size_type ncols(const linalg_type &m) { return m.ncols(); }
    iterator begin(linalg_type &m) { return m.begin(); }
    const_iterator const_begin(const linalg_type &m) { return m.begin(); }
    iterator end(linalg_type &m) { return m.end(); }
    const_iterator const_end(const linalg_type &m) { return m.end(); }
    const_sub_row_type row(const linalg_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_col_type col(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_row_type row(linalg_type &, size_type)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_col_type col(linalg_type &m, size_type i) 
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    void clear(linalg_type &v) { v.clear(); }
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

  /* ******************************************************************** */
  /*		Scalar product                             		  */
  /* ******************************************************************** */

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
    return _vect_sp_plain(linalg_traits<V1>().const_begin(v1),
			  linalg_traits<V1>().const_end(v1),
			  linalg_traits<V2>().const_begin(v2));
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_plain) {
    return _vect_sp_sparse(linalg_traits<V1>().const_begin(v1),
			   linalg_traits<V1>().const_end(v1), v2);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2, abstract_plain,abstract_sparse) {
    return _vect_sp_sparse(linalg_traits<V2>().const_begin(v2),
			   linalg_traits<V2>().const_end(v2), v1);
  }

  template <class V1, class V2> inline
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2,abstract_sparse,abstract_sparse) {
    return _vect_sp_sparse(linalg_traits<V1>().const_begin(v1),
			   linalg_traits<V1>().const_end(v1), v2);
  }

  template <class V1, class V2>
    typename linalg_traits<V1>::base_type
    vect_sp(const V1 &v1, const V2 &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    return vect_sp(v1, v2,
		   typename linalg_traits<V1>::storage_type(), 
		   typename linalg_traits<V2>::storage_type());
  }

  /* ******************************************************************** */
  /*		Euclidian norm                             		  */
  /* ******************************************************************** */

   template <class V>
    typename linalg_traits<V>::base_type vect_norm_new(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = linalg_traits<V>().const_begin(v),
      ite = linalg_traits<V>().const_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res += dal::sqr(modulus(*it));
    return sqrt(res);
  }

  /* ******************************************************************** */
  /*		Inifity norm                              		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::base_type vect_norminf_new(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = linalg_traits<V>().const_begin(v),
      ite = linalg_traits<V>().const_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res = std::max(res, modulus(*it));
    return res;
  }
  
  /* ******************************************************************** */
  /*		norm1                                    		  */
  /* ******************************************************************** */

  template <class V>
    typename linalg_traits<V>::base_type vect_norm1_new(const V &v) {
    typename linalg_traits<V>::const_iterator
      it = linalg_traits<V>().const_begin(v),
      ite = linalg_traits<V>().const_end(v);
    typename linalg_traits<V>::base_type res(0);
    for (; it != ite; ++it) res += modulus(*it);
    return res;
  }

  /* ******************************************************************** */
  /*		Copy                                    		  */
  /* ******************************************************************** */

  template <class L1, class L2>
  void copy(const L1& l1, L2& l2) {
    if (vect_size(l1) != vect_size(l2) || mat_ncols(l1) != mat_ncols(l2)
	|| mat_nrows(l1) != mat_nrows(l2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if (mat_ncols(l1) != 1)
      copy_mat(l1, l2, typename linalg_traits<L1>::sub_orientation(),
	 typename linalg_traits<L2>::sub_orientation());
    else
      copy_vect(l1, l2, typename linalg_traits<L1>::storage_type(),
		typename linalg_traits<L2>::storage_type());
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
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (; it != ite; ++it) l2(i, it.index()) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_rc(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(i, j) = *it;
  }

  template <class L1, class L2> inline
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i) {
    copy_mat_mixed_cr(l1, l2, i, typename linalg_traits<L1>::storage_type());
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (; it != ite; ++it) l2(it.index(), i) = *it;
  }

  template <class L1, class L2>
  void copy_mat_mixed_cr(const L1& l1, L2& l2, size_type i, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (size_type j = 0; it != ite; ++it, ++j) l2(j, i) = *it;
  }

  template <class L1, class L2>
  void copy_mat(const L1& l1, L2& l2, row_major, col_major) {
    clear(l2);
    size_type nbr = mat_nrows(l1);
    for (size_type i = 0; i < nbr; ++i)
      copy_mat_mixed_rc(mat_row(l1, i), l2, i);
  }
  
  template <class L1, class L2> inline
  void copy_mat(const L1& l1, L2& l2, col_major, row_major) {
    clear(l2);
    size_type nbc = mat_ncols(l1);
    for (size_type i = 0; i < nbr; ++i)
      copy_mat_mixed_cr(mat_col(l1, i), l2, i);
  }
  
  template <class L1, class L2> inline
  void copy_vect(const L1 &l1, L2 &l2,
		 abstract_plain, abstract_plain) {
    std::copy(linalg_traits<L1>().const_begin(l1),
	      linalg_traits<L1>().const_end(l1),
	      linalg_traits<L2>().begin(l2));
  }
  
  template <class L1, class L2>
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_plain) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (; it != ite; ++it) l2[it.index()] = *it;
  }
  
  template <class L1, class L2> inline
  void copy_vect(const L1& l1, L2& l2,
		 abstract_sparse, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
    for (; it != ite; ++it) l2[it.index()] = *it;
  }
  
  template <class L1, class L2> inline
  void copy_vect(const L1& l1, L2& l2,
		 abstract_plain, abstract_sparse) {
    clear(l2);
    typename linalg_traits<L1>::const_iterator
      it  = linalg_traits<L1>().const_begin(l1),
      ite = linalg_traits<L1>().const_end(l1);
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

  template <class L1, class L2>
    void add(const L1& l1, L2& l2) {
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
    _add_full(linalg_traits<L1>().const_begin(l1),
	      linalg_traits<L2>().const_begin(l2),
	      linalg_traits<L3>().begin(l3),
	      linalg_traits<L3>().end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_plain, abstract_plain) {
    _add_almost_full(linalg_traits<L1>().const_begin(l1),
		     linalg_traits<L1>().const_end(l1),
		     linalg_traits<L2>().const_begin(l2),
		     linalg_traits<L3>().begin(l3),
		     linalg_traits<L3>().end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_plain, abstract_sparse, abstract_plain) {
    _add_almost_full(linalg_traits<L2>().const_begin(l2),
		     linalg_traits<L2>().const_end(l2),
		     linalg_traits<L1>().const_begin(l1),
		     linalg_traits<L3>().begin(l3),
		     linalg_traits<L3>().end(l3));
  }
  
  template <class L1, class L2, class L3> inline
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_plain) {
    _add_almost_full(linalg_traits<L1>().const_begin(l1),
		     linalg_traits<L1>().const_end(l1),
		     linalg_traits<L2>().const_begin(l2),
		     linalg_traits<L2>().const_end(l2),
		     linalg_traits<L3>().begin(l3),
		     linalg_traits<L3>().end(l3));
  }
  
  template <class L1, class L2, class L3>
  void add(const L1& l1, const L2& l2, L3& l3,
	   abstract_sparse, abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = linalg_traits<L1>().const_begin(l1), 
      ite1 = linalg_traits<L1>().const_end(l1);
    typename linalg_traits<L2>::const_iterator
      it2 = linalg_traits<L2>().const_begin(l2), 
      ite2 = linalg_traits<L2>().const_end(l2);
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
      it1 = linalg_traits<L1>().const_begin(l1); 
    typename linalg_traits<L2>::iterator
      it2 = linalg_traits<L2>().begin(l2), 
      ite = linalg_traits<L2>().end(l2);
    for (; it2 != ite; ++it2, ++it1) *it2 += *it1;
  }
  
  template <class L1, class L2>
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_plain) {
    typename linalg_traits<L1>::const_iterator
      it1 = linalg_traits<L1>().const_begin(l1), 
      ite1 = linalg_traits<L1>().const_end(l1);
    for (; it1 != ite1; ++it1) l2[it1.index()] = *it1;
  }
  
  template <class L1, class L2>
  void add(const L1& l1, L2& l2,
	   abstract_sparse, abstract_sparse) {
    typename linalg_traits<L1>::const_iterator
      it1 = linalg_traits<L1>().const_begin(l1), 
      ite1 = linalg_traits<L1>().const_end(l1);
    for (; it1 != ite1; ++it1) l2[it1.index()] += *it1;
  }
  
  template <class L1, class L2> inline
  void add(const L1&, L2&, abstract_plain, abstract_sparse)
  { DAL_THROW(failure_error,"Unauthorized addition"); } 



  /* ******************************************************************** */
  /*		conjugate gradient  (unpreconditionned        		  */
  /* ******************************************************************** */

  template < class Matrix, class Vector>
  int cg_new(const Matrix& A, Vector& x, const Vector& b, int itemax, 
	 double residu, bool noisy = true) {
    typedef typename linalg_traits<Vector>::base_type base_type;
    base_type rho(0), rho_1(0), alpha(0), beta(0);
    Vector p(x.size()), q(x.size()), r(x.size());
    int iter = 0;
    mult(A, scaled(x, -1.0), b, r);
    rho = vect_sp(r,r);
    
    while (::sqrt(rho) > residu) {
      if (iter == 0) copy(r, p);		  
      else { beta = rho / rho_1; add(r, scaled(p, beta), p); }
      
      mult(A, p, q);
      alpha = rho / vect_sp(p, q);
      
      add(scaled(p, alpha), x);
      add(scaled(q, -alpha), r);
      
      rho_1 = rho;
      
      ++iter;
      rho = vect_sp(r, r);
      if (noisy) cout << "iter " << iter << " residu " << ::sqrt(rho) << endl;
      if (iter >= itemax) return 1;
    }
    return 0;
  }
  
}


#endif //  __BGEOT_ABSTRACT_LINALG_H
