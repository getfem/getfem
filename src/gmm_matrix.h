/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_matrix.h : matrices.                                     */
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

#ifndef __GMM_MATRIX_H
#define __GMM_MATRIX_H

#include <gmm.h>

namespace gmm
{

  /* ******************************************************************** */
  /*		Identity matrix                         		  */
  /* ******************************************************************** */

  struct identity_matrix {};

  template <class V1, class V2> inline
  void mult(identity_matrix, const V1 &v1, V2 &v2) { copy(v1, v2); }
  template <class V1, class V2> inline
  void mult(identity_matrix, const V1 &v1, const V2 &v2) { copy(v1, v2); }
  template <class V1, class V2, class V3> inline
  void mult(identity_matrix, const V1 &v1, const V2 &v2, V3 &v3)
  { add(v1, v2, v3); }
  template <class V1, class V2, class V3> inline
  void mult(const identity_matrix&, const V1 &v1, const V2 &v2, const V3 &v3)
  { add(v1, v2, v3); }
  template <class M> void copy(const identity_matrix&, M &m) {
    size_type i = 0, n = std::max(mat_nrows(m), mat_ncols(m)); clear(m);
    for (; i < n; ++i) m(i,i) = typename linalg_traits<M>::value_type(1);
  }
  template <class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp(const identity_matrix &, const V1 &v1, const V2 &v2)
  { return vect_sp(v1, v2); }
  template <class M> inline void copy(const identity_matrix&, const M &m)
  { copy_ident(m, linalg_traits<M>::is_reference()); }
  template <class M> inline void copy_ident(const M &m, linalg_true)
  { copy(identity_matrix(), const_cast<M &>(m)); }
  template<class M> inline bool is_identity(const M&) { return false; }
  inline bool is_identity(const identity_matrix&) { return true; }

  /* ******************************************************************** */
  /*		Row matrix                                   		  */
  /* ******************************************************************** */


  template<class V> class row_matrix {
  protected :
    std::vector<V> li; /* array of rows.                                   */
    
  public :
    
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::value_type value_type;
    
    row_matrix(size_type r, size_type c) : li(r)
    { for (size_type i = 0; i < r; ++i) li[i] = V(c); }
    row_matrix(void) {}
    reference_type operator ()(size_type l, size_type c)
    { return li[l][c]; }
    value_type operator ()(size_type l, size_type c) const
    { return li[l][c]; }

    void clear_row(size_type i) { clear(li[i]); }
    void clear() { for (size_type i=0; i < nrows(); ++i) clear_row(i); }
    
    V& row(size_type i) { return li[i]; }
    const V& row(size_type i) const { return li[i]; }
    
    inline size_type nrows(void) const { return li.size(); }
    inline size_type ncols(void) const
    { return (nrows() == 0) ? 0 : li[0].size(); }
  };

  template <class V> struct linalg_traits<row_matrix<V> > {
    typedef row_matrix<V> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef simple_vector_ref<V *> sub_row_type;
    typedef simple_vector_ref<const V *> const_sub_row_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const this_type &m, size_type i)
    { return const_sub_row_type(m.row(i)); }
    sub_row_type row(this_type &m, size_type i) 
    { return sub_row_type(m.row(i)); }
    const void* origin(const this_type &m) { return &m; }
    void do_clear(this_type &m) { m.clear(); }
  };

  /* ******************************************************************** */
  /*		Column matrix                                		  */
  /* ******************************************************************** */

  template<class V> class col_matrix {
  protected :
    std::vector<V> li; /* array of columns.                               */
    
  public :
    
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::value_type value_type;
    
    col_matrix(size_type r, size_type c) : li(c)
    { for (size_type i = 0; i < c; ++i) li[i] = V(r); }
    col_matrix(void) {}
    reference_type operator ()(size_type l, size_type c)
    { return li[c][l]; }
    value_type operator ()(size_type l, size_type c) const
    { return li[c][l]; }

    void clear_col(size_type i) { clear(li[i]); }
    void clear() { for (size_type i=0; i < ncols(); ++i) clear_col(i); }
    
    V& col(size_type i) { return li[i]; }
    const V& col(size_type i) const { return li[i]; }
    
    inline size_type ncols(void) const { return li.size(); }
    inline size_type nrows(void) const
    { return (nrows() == 0) ? 0 : li[0].size(); }
  };

  template <class V> struct linalg_traits<col_matrix<V> > {
    typedef col_matrix<V> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef simple_vector_ref<V *> sub_col_type;
    typedef simple_vector_ref<const V *> const_sub_col_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_col_type col(const this_type &m, size_type i)
    { return const_sub_col_type(m.col(i)); }
    sub_col_type col(this_type &m, size_type i) 
    { return sub_col_type(m.col(i)); }
    const void* origin(const this_type &m) { return &m; }
    void do_clear(this_type &m) { m.clear(); }
  };


}

#endif /* __GMM_MATRIX_H */
