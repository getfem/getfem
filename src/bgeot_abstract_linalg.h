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
  struct column_major {};       // matrix with a column access
                                // (default for vector)
  struct row_and_column {};     // both accesses.

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
  
  template <class VECT> class simple_vector_ref {
    protected :
      VECT *l;

    public :
      typedef typename linalg_traits<VECT>::base_type base_type;
      VECT &deref(void) { return *l; }
      const VECT &deref(void) const { return *l; }
      base_type &operator[](size_type i) { return (*l)[i]; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class VECT> struct linalg_traits<simple_vector_ref<VECT> > {
    typedef typename simple_vector_ref<VECT> linalg_type;
    typedef typename linalg_traits<VECT>::base_type base_type;
    typedef typename linalg_traits<VECT>::iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    typedef typename linalg_traits<VECT>::sub_type sub_type;
    typedef typename linalg_traits<VECT>::const_sub_type const_sub_type;
    typedef typename linalg_traits<VECT>::sub_orientation sub_orientation;
    size_type size(const linalg_type &v)
    { return linalg_traits<VECT>::size(v.deref()); }
    size_type nrows(const linalg_type &v)
    { return linalg_traits<VECT>::nrows(v.deref()); }
    size_type ncols(const linalg_type &v)
    { return linalg_traits<VECT>::ncols(v.deref()); }
    iterator begin(linalg_type &v)
    { return linalg_traits<VECT>::begin(v.deref()); }
    iterator const_begin(const linalg_type &v)
    { return linalg_traits<VECT>::const_begin(v.deref()); }
    iterator end(linalg_type &v)
    { return linalg_traits<VECT>::end(v.deref()); }
    iterator const_end(const linalg_type &v)
    { return linalg_traits<VECT>::const_end(v.deref()); }
    const_sub_type row(const linalg_type &v, size_type i)
    { return linalg_traits<VECT>::row(v.deref(), i); }
    const_sub_type column(const linalg_type &v, size_type i)
    { return linalg_traits<VECT>::column(v.deref(), i); }
    sub_type row(linalg_type &v, size_type i)
    { return linalg_traits<VECT>::row(v.deref(), i); }
    sub_type column(linalg_type &v, size_type i)
    { return linalg_traits<VECT>::column(v.deref(), i); }
  };

  template <class VECT> class simple_vector_const_ref {
    protected :
      const VECT *l;

    public :
      typedef typename linalg_traits<VECT>::base_type base_type;
      const VECT &deref(void) const { return *l; }
      base_type operator[](size_type i) const { return (*l)[i]; }
  };

  template <class VECT> struct linalg_traits<simple_vector_const_ref<VECT> > {
    typedef typename simple_vector_const_ref<VECT> linalg_type;
    typedef typename linalg_traits<VECT>::base_type base_type;
    typedef typename linalg_traits<VECT>::const_iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    typedef typename linalg_traits<VECT>::const_sub_type sub_type;
    typedef typename linalg_traits<VECT>::const_sub_type sub_type;
    typedef typename linalg_traits<VECT>::sub_orientation sub_orientation;
    size_type size(const linalg_type &v)
    { return linalg_traits<VECT>::size(v.deref()); }
    size_type nrows(const linalg_type &v)
    { return linalg_traits<VECT>::nrows(v.deref()); }
    size_type ncols(const linalg_type &v)
    { return linalg_traits<VECT>::ncols(v.deref()); }
    iterator begin(linalg_type &v)
    { return linalg_traits<VECT>::const_begin(v.deref()); }
    iterator const_begin(const linalg_type &v)
    { return linalg_traits<VECT>::const_begin(v.deref()); }
    iterator end(linalg_type &v)
    { return linalg_traits<VECT>::const_end(v.deref()); }
    iterator const_end(const linalg_type &v)
    { return linalg_traits<VECT>::const_end(v.deref()); }
    const_sub_type row(const linalg_type &v, size_type i)
    { return linalg_traits<VECT>::row(v.deref(), i); }
    const_sub_type column(const linalg_type &v, size_type i)
    { return linalg_traits<VECT>::column(v.deref(), i); }
  };

  /* ******************************************************************** */
  /*		standard extern references for plain vectors              */
  /* ******************************************************************** */
  
  template <class VECT, class T> class extern_plain_vector_ref {
    protected :
    VECT *l;
    
    public :
    typedef T base_type;
    VECT &deref(void) { return *l; }
    const VECT &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_plain_vector_ref(VECT &v) : l(&v) {}
  };
  
  template <class VECT, class T> class extern_plain_vector_const_ref {
    protected :
    const VECT *l;
    
    public :
    typedef T base_type;
    const VECT &deref(void) const { return *l; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_plain_vector_const_ref(const VECT &v) : l(&v) {}
    extern_plain_vector_const_ref(const extern_plain_vector_ref<VECT, T> &v)
    { l = &(v.deref()); }
  };

  template <class VECT, class T>
    struct linalg_traits<extern_plain_vector_ref<VECT, T> > {
      typedef typename extern_plain_vector_ref<VECT, T> linalg_type;
      typedef typename T base_type;
      typedef typename VECT::iterator  iterator;
      typedef typename VECT::const_iterator const_iterator;
      typedef typename abstract_plain storage_type;
      typedef typename linalg_type sub_type;
      typedef typename extern_plain_vector_const_ref<VECT, T> const_sub_type;
      typedef typename column_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_type column(const linalg_type &v, size_type i) { return v; }
      sub_type row(linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      sub_type column(linalg_type &v, size_type i) { return v; }
    };

    template <class VECT, class T>
    struct linalg_traits<extern_plain_vector_const_ref<VECT, T> > {
      typedef typename extern_plain_vector_const_ref<VECT, T> linalg_type;
      typedef typename T base_type;
      typedef typename VECT::iterator  iterator;
      typedef typename VECT::const_iterator const_iterator;
      typedef typename abstract_plain storage_type;
      typedef typename linalg_type sub_type;
      typedef typename linalg_type const_sub_type;
      typedef typename column_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_type column(const linalg_type &v, size_type i) { return v; }
    };

  /* ******************************************************************** */
  /*		standard extern references for sparse vectors             */
  /* ******************************************************************** */
  
  template <class VECT, class T> class extern_sparse_vector_ref {
    protected :
    VECT *l;
    
    public :
    typedef T base_type;
    VECT &deref(void) { return *l; }
    const VECT &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_sparse_vector_ref(VECT &v) : l(&v) {}
  };

  template <class VECT, class T> class extern_sparse_vector_const_ref {
    protected :
    const VECT *l;
    
    public :
    typedef T base_type;
    VECT &deref(void) { return *l; }
    const VECT &deref(void) const { return *l; }
    base_type &operator[](size_type i) { return (*l)[i]; }
    base_type operator[](size_type i) const { return (*l)[i]; }
    extern_sparse_vector_const_ref(const VECT &v) : l(&v) {}
    extern_sparse_vector_const_ref(const extern_sparse_vector_ref<VECT, T> &v)
    { l = &(v.deref()); }
  };

  template <class VECT, class T>
    struct linalg_traits<extern_sparse_vector_ref<VECT, T> > {
      typedef typename extern_sparse_vector_ref<VECT, T> linalg_type;
      typedef typename T base_type;
      typedef typename VECT::iterator  iterator;
      typedef typename VECT::const_iterator const_iterator;
      typedef typename abstract_sparse storage_type;
      typedef typename linalg_type sub_type;
      typedef typename extern_sparse_vector_const_ref<VECT, T> const_sub_type;
      typedef typename column_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_type column(const linalg_type &v, size_type i) { return v; }
      sub_type row(linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      sub_type column(linalg_type &v, size_type i) { return v; }
    };

  template <class VECT, class T>
    struct linalg_traits<extern_sparse_vector_const_ref<VECT, T> > {
      typedef typename extern_sparse_vector_const_ref<VECT, T> linalg_type;
      typedef typename T base_type;
      typedef typename VECT::iterator  iterator;
      typedef typename VECT::const_iterator const_iterator;
      typedef typename abstract_sparse storage_type;
      typedef typename linalg_type sub_type;
      typedef typename linalg_type const_sub_type;
      typedef typename column_major sub_orientation;
      size_type size(const linalg_type &v) { return v.size(); }
      size_type nrows(const linalg_type &v) { return v.size(); }
      size_type ncols(const linalg_type &v) { return 1; }
      iterator begin(linalg_type &v) { return v.begin(); }
      iterator const_begin(const linalg_type &v) { return v.begin(); }
      iterator end(linalg_type &v) { return v.end(); }
      iterator const_end(const linalg_type &v) { return v.end(); }
      const_sub_type row(const linalg_type &v, size_type i)
      { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
      const_sub_type column(const linalg_type &v, size_type i) { return v; }
    };


  /* ******************************************************************** */
  /*		                                         		  */
  /*		Traits for bgeot objects                     		  */
  /*		                                         		  */
  /* ******************************************************************** */

  template <class T, int N> struct linalg_traits<fsvector<T, N> > {
    typedef typename fsvector<T> linalg_type;
    typedef typename T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef typename abstract_plain storage_type;
    typedef typename simple_vector_ref<fsvector<T, N> > sub_type;
    typedef typename simple_vector_const_ref<fsvector<T, N> > const_sub_type;
    typedef typename column_major sub_orientation;
    size_type size(const linalg_type &v) { return N; }
    size_type nrows(const linalg_type &v) { return N; }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.begin(); }
    iterator const_begin(const linalg_type &v) { return v.begin(); }
    iterator end(linalg_type &v) { return v.end(); }
    iterator const_end(const linalg_type &v) { return v.end(); }
    const_sub_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_type column(const linalg_type &v, size_type i)
    { return const_sub_type(v); }
    sub_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_type column(linalg_type &v, size_type i) { return sub_type(v); }
  };


  template <class T> struct linalg_traits<vsvector<T> > {
    typedef typename vsvector<T> linalg_type;
    typedef typename T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef typename abstract_plain storage_type;
    typedef typename simple_vector_ref<vsvector<T> > sub_type;
    typedef typename simple_vector_const_ref<vsvector<T> > const_sub_type;
    typedef typename column_major sub_orientation;
    size_type size(const linalg_type &v) { return v.size(); }
    size_type nrows(const linalg_type &v) { return v.size(); }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.begin(); }
    iterator const_begin(const linalg_type &v) { return v.begin(); }
    iterator end(linalg_type &v) { return v.end(); }
    iterator const_end(const linalg_type &v) { return v.end(); }
    const_sub_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    const_sub_type column(const linalg_type &v, size_type i)
    { return const_sub_type(v); }
    sub_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Rows inaccessible for this object"); }
    sub_type column(linalg_type &v, size_type i) { return sub_type(v); }
  };

  template <class T> struct linalg_traits<svector<T> > {
    typedef typename svector<T> linalg_type;
    typedef typename T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef typename abstract_sparse storage_type;
    typedef typename simple_vector_ref<svector<T> > sub_type;
    typedef typename simple_vector_const_ref<svector<T> > const_sub_type;
    typedef typename column_major sub_orientation;
    size_type size(const linalg_type &v) { return v.size(); }
    size_type nrows(const linalg_type &v) { return v.size(); }
    size_type ncols(const linalg_type &v) { return 1; }
    iterator begin(linalg_type &v) { return v.begin(); }
    iterator const_begin(const linalg_type &v) { return v.begin(); }
    iterator end(linalg_type &v) { return v.end(); }
    iterator const_end(const linalg_type &v) { return v.end(); }
    const_sub_type row(const linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_type column(const linalg_type &v, size_type i)
    { return const_sub_type(v); }
    sub_type row(linalg_type &v, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_type column(linalg_type &v, size_type i) 
    { return sub_type(v); }
  };

  template <class T> struct linalg_traits<vsmatrix<T> > {
    typedef typename vsmatrix<T> linalg_type;
    typedef typename T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef typename abstract_plain storage_type;
    typedef typename simple_vector_ref<vsvector<T> > sub_type;
    typedef typename simple_vector_const_ref<vsvector<T> > const_sub_type;
    typedef typename column_major sub_orientation;
    size_type size(const linalg_type &m) { return m.size(); }
    size_type nrows(const linalg_type &m) { return m.nrows(); }
    size_type ncols(const linalg_type &m) { return m.ncols(); }
    iterator begin(linalg_type &m) { return m.begin(); }
    iterator const_begin(const linalg_type &m) { return m.begin(); }
    iterator end(linalg_type &m) { return m.end(); }
    iterator const_end(const linalg_type &m) { return m.end(); }
    const_sub_type row(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_type column(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_type row(linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_type column(linalg_type &m, size_type i) 
    { DAL_THROW(failure_error,"Sorry, to be done"); }
  };

  template <class T> struct linalg_traits<smatrix<T> > {
    typedef typename smatrix<T> linalg_type;
    typedef typename T base_type;
    typedef typename linalg_type::iterator  iterator;
    typedef typename linalg_type::const_iterator const_iterator;
    typedef typename abstract_sparse storage_type;
    typedef typename simple_vector_ref<svector<T> > sub_type;
    typedef typename simple_vector_const_ref<svector<T> > const_sub_type;
    typedef typename row_major sub_orientation;
    size_type size(const linalg_type &m) { return m.size(); }
    size_type nrows(const linalg_type &m) { return m.nrows(); }
    size_type ncols(const linalg_type &m) { return m.ncols(); }
    iterator begin(linalg_type &m) { return m.begin(); }
    iterator const_begin(const linalg_type &m) { return m.begin(); }
    iterator end(linalg_type &m) { return m.end(); }
    iterator const_end(const linalg_type &m) { return m.end(); }
    const_sub_type row(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    const_sub_type column(const linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_type row(linalg_type &m, size_type i)
    { DAL_THROW(failure_error,"Sorry, to be done"); }
    sub_type column(linalg_type &m, size_type i) 
    { DAL_THROW(failure_error,"Sorry, to be done"); }
  };


  /* ******************************************************************** */
  /*		                                         		  */
  /*		Generic algorithm                           		  */
  /*		                                         		  */
  /* ******************************************************************** */

  template <class VECT> inline size_type vect_size(const VECT &v)
    { return typename linalg_traits<VECT>::size(v); }

  template <class MAT> inline size_type mat_nrows(const MAT &m)
    { return typename linalg_traits<MAT>::nrows(v); }

  template <class MAT> inline size_type mat_ncols(const MAT &m)
    { return typename linalg_traits<MAT>::ncols(v); }

  /* ******************************************************************** */
  /*		Scalar product                             		  */
  /* ******************************************************************** */

  template <class VECT1, class VECT2>
    typename linalg_traits<VECT1>::base_type
    vect_sp(const VECT1 &v1, const VECT2 &v2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    return vect_sp(v1, v2,
		   (typename linalg_traits<VECT1>::storage_type)(), 
		   (typename linalg_traits<VECT2>::storage_type)());
  }

  template <class IT1, class IT2>
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_plain(IT1 &it, IT1 &ite, IT2 &it2) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it, ++it2) res += conj_product(*it, *it2);
    return res;
  }

  template <class IT1, class VECT>
    typename std::iterator_traits<IT1>::value_type
    _vect_sp_sparse(IT1 &it, IT1 &ite, const VECT &v) {
    typename std::iterator_traits<IT1>::value_type res(0);
    for (; it != ite; ++it) res += conj_product(*it, v[it->index()]);
    return res;
  }

  template <class VECT1, class VECT2> inline
    typename linalg_traits<VECT1>::base_type
    vect_sp(const VECT1 &v1, const VECT2 &v2, abstract_plain, abstract_plain) {
    return _vect_sp_plain(linalg_traits<VECT1>::const_begin(v1),
			  linalg_traits<VECT1>::const_end(v1),
			  linalg_traits<VECT2>::const_begin(v2));
  }

  template <class VECT1, class VECT2> inline
    typename linalg_traits<VECT1>::base_type
    vect_sp(const VECT1 &v1, const VECT2 &v2,abstract_sparse,abstract_plain) {
    return _vect_sp_sparse(linalg_traits<VECT1>::const_begin(v1),
			   linalg_traits<VECT1>::const_end(v1), v2);
  }

  template <class VECT1, class VECT2> inline
    typename linalg_traits<VECT1>::base_type
    vect_sp(const VECT1 &v1, const VECT2 &v2, abstract_plain,abstract_sparse) {
    return _vect_sp_sparse(linalg_traits<VECT2>::const_begin(v2),
			   linalg_traits<VECT2>::const_end(v2), v1);
  }

  template <class VECT1, class VECT2> inline
    typename linalg_traits<VECT1>::base_type
    vect_sp(const VECT1 &v1, const VECT2 &v2,abstract_sparse,abstract_sparse) {
    return _vect_sp_sparse(linalg_traits<VECT1>::const_begin(v1),
			   linalg_traits<VECT1>::const_end(v1), v2);
  }

  /* ******************************************************************** */
  /*		Euclidian norm                             		  */
  /* ******************************************************************** */

   template <class VECT>
    typename linalg_traits<VECT>::base_type vect_norminf(const VECT &v) {
    linalg_traits<VECT>::const_iterator
      it = linalg_traits<VECT>::const_begin(v),
      ite = linalg_traits<VECT>::const_end(v);
    typename linalg_traits<VECT>::base_type res(0);
    for (; it != ite; ++it) res += dal::sqr(modulus(*it));
    return sqrt(res);
  }

  /* ******************************************************************** */
  /*		Inifity norm                              		  */
  /* ******************************************************************** */

  template <class VECT>
    typename linalg_traits<VECT>::base_type vect_norminf(const VECT &v) {
    linalg_traits<VECT>::const_iterator
      it = linalg_traits<VECT>::const_begin(v),
      ite = linalg_traits<VECT>::const_end(v);
    typename linalg_traits<VECT>::base_type res(0);
    for (; it != ite; ++it) res = std::max(res, modulus(*it));
    return res;
  }
  
  /* ******************************************************************** */
  /*		norm1                                    		  */
  /* ******************************************************************** */

  template <class VECT>
    typename linalg_traits<VECT>::base_type vect_norm1(const VECT &v) {
    linalg_traits<VECT>::const_iterator
      it = linalg_traits<VECT>::const_begin(v),
      ite = linalg_traits<VECT>::const_end(v);
    typename linalg_traits<VECT>::base_type res(0);
    for (; it != ite; ++it) res += modulus(*it);
    return res;
  }

  /* ******************************************************************** */
  /*		addition                                    		  */
  /* ******************************************************************** */

  template <class LINALG1, class LINALG2>
    void add(const LINALG1& l1, LINALG2& l2) {
    if (vect_size(v1) != vect_size(v2))
      DAL_THROW(dimension_error,"dimensions mismatch");
    add(l1, l2, (typename linalg_traits<LINALG1>::storage_type)(),
	(typename linalg_traits<LINALG2>::storage_type)());
  }

  template <class LINALG1, class LINALG2, class LINALG3>
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3) {
    if (vect_size(v1) != vect_size(v2) || vect_size(v1) != vect_size(v3))
      DAL_THROW(dimension_error,"dimensions mismatch");
    if ((const void *)(&l1) == (const void *)(&l3))
      add(l2, l3);
    else if ((const void *)(&l2) == (const void *)(&l3))
      add(l1, l3);
    else
      add(l1, l2, l3, (typename linalg_traits<LINALG1>::storage_type)(),
	  (typename linalg_traits<LINALG2>::storage_type)(),
	  (typename linalg_traits<LINALG3>::storage_type)());
  }

  template <class IT1, class IT2, class IT3>
    _add_full(IT1 &it1, IT2 &it2, IT3 &it3, IT3 &ite) {
    for (; it3 != ite; ++it3, ++it2, ++it1) *it3 = *it1 + *it2;
  }

  template <class IT1, class IT2, class IT3>
    _add_almost_full(IT1 &it1, IT1 &ite1, IT2 &it2, IT3 &it3, IT3 &ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it, ++it2) *it = *it2;
    for (; it1 != ite1; ++it1) *(it3 + it1->index()) += *it1;
  }

  template <class IT1, class IT2, class IT3>
    _add_to_full(IT1 &it1, IT1 &ite1, IT2 &it2, IT2 &ite2,IT3 &it3,IT3 &ite3) {
    IT3 it = it3;
    for (; it != ite3; ++it) *it = 0;
    for (; it1 != ite1; ++it1) *(it3 + it1->index()) += *it1;
    for (; it2 != ite2; ++it2) *(it3 + it2->index()) += *it2;    
  }
  
  template <class LINALG1, class LINALG2, class LINALG3> inline
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3,
	     abstract_plain, abstract_plain, abstract_plain) {
    _add_full(linalg_traits<LINALG1>::const_begin(l1),
	      linalg_traits<LINALG2>::const_begin(l2),
	      linalg_traits<LINALG3>::begin(l3),
	      linalg_traits<LINALG3>::end(l3));
  }
  
  template <class LINALG1, class LINALG2, class LINALG3> inline
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3,
	     abstract_sparse, abstract_plain, abstract_plain) {
    _add_almost_full(linalg_traits<LINALG1>::const_begin(l1),
		     linalg_traits<LINALG1>::const_end(l1),
		     linalg_traits<LINALG2>::const_begin(l2),
		     linalg_traits<LINALG3>::begin(l3),
		     linalg_traits<LINALG3>::end(l3));
  }
  
  template <class LINALG1, class LINALG2, class LINALG3> inline
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3,
	     abstract_plain, abstract_sparse, abstract_plain) {
    _add_almost_full(linalg_traits<LINALG2>::const_begin(l2),
		     linalg_traits<LINALG2>::const_end(l2),
		     linalg_traits<LINALG1>::const_begin(l1),
		     linalg_traits<LINALG3>::begin(l3),
		     linalg_traits<LINALG3>::end(l3));
  }

  template <class LINALG1, class LINALG2, class LINALG3> inline
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3,
	     abstract_sparse, abstract_sparse, abstract_plain) {
    _add_almost_full(linalg_traits<LINALG1>::const_begin(l1),
		     linalg_traits<LINALG1>::const_end(l1),
		     linalg_traits<LINALG2>::const_begin(l2),
		     linalg_traits<LINALG2>::const_end(l2),
		     linalg_traits<LINALG3>::begin(l3),
		     linalg_traits<LINALG3>::end(l3));
  }

  template <class LINALG1, class LINALG2, class LINALG3> inline
    void add(const LINALG1& l1, const LINALG2& l2, LINALG3& l3,
	     abstract_sparse, abstract_sparse, abstract_sparse) {
    typename linalg_traits<LINALG1>::const_iterator
      it1 = linalg_traits<LINALG1>::const_begin(l1), 
      ite1 = linalg_traits<LINALG1>::const_end(l1);
    typename linalg_traits<LINALG2>::const_iterator
      it2 = linalg_traits<LINALG2>::const_begin(l2), 
      ite2 = linalg_traits<LINALG2>::const_end(l2);
    linalg_traits<LINALG3>::clear(l3);
    
  }
  


}


#endif //  __BGEOT_ABSTRACT_LINALG_H
