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
