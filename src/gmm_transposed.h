/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_transposed.h : generic transposed matrices.              */
/*     									   */
/* Date : November 10, 2002.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef __GMM_TRANSPOSED_H
#define __GMM_TRANSPOSED_H

namespace gmm {

  /* ********************************************************************* */
  /*		transposed reference                    		   */
  /* ********************************************************************* */
  
  template <class PT> struct  transposed_row_ref {
    
    typedef transposed_row_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef M * CPT;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<this_type>::col_iterator iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<M>::access_type access_type;

    iterator _begin, _end;
    const void *origin; 

    transposed_row_ref(ref_M m) : _begin(mat_row_begin(m)), 
      _end(mat_row_end(m)),
      origin(linalg_origin(m)) {}

    transposed_row_ref(const transposed_row_ref<CPT> &cr) :
      _begin(cr._begin),_end(cr._end),origin(cr.origin) {}

    reference operator()(size_type i, size_type j) const
    { return access_type()(_begin+j, i); }
  };

  template <class PT> struct transposed_row_matrix_access {
    typedef transposed_row_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::col_iterator iterator;
    typedef typename linalg_traits<M>::access_type access_type;
    
    reference operator()(const iterator &itcol, size_type i)
    { return access_type()(itcol, i); }
  };

  template <class PT> struct linalg_traits<transposed_row_ref<PT> > {
    typedef transposed_row_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename select_return<value_type,
            typename linalg_traits<M>::reference, PT>::return_type reference;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type row_iterator;
    typedef abstract_null_type const_row_iterator;
    typedef typename linalg_traits<M>::const_sub_row_type const_sub_col_type;
    typedef typename select_return<const_sub_col_type, typename
            linalg_traits<M>::sub_row_type, PT>::return_type sub_col_type;
    typedef typename linalg_traits<M>::const_row_iterator const_col_iterator;
    typedef typename select_return<const_col_iterator, typename
            linalg_traits<M>::row_iterator, PT>::return_type col_iterator;
    typedef col_major sub_orientation;
    typedef transposed_row_matrix_access<PT> access_type;
    static size_type ncols(const this_type &v) { return v._end - v._begin; }
    static size_type nrows(const this_type &v)
    { return (ncols(v) == 0) ? 0 : vect_size(const_mat_col(v, 0)); }
    static const_sub_col_type col(const const_col_iterator &it)
    { return linalg_traits<M>::row(it); }
    static sub_col_type col(col_iterator &it)
    { return linalg_traits<M>::row(it); }
    static col_iterator col_begin(this_type &m) { return m._begin; }
    static col_iterator col_end(this_type &m) { return m._end; }
    static const_col_iterator col_begin(const this_type &m)
    { return m._begin; }
    static const_col_iterator col_end(const this_type &m) { return m._end; }
    static const void* origin(const this_type &v) { return v.origin; }
    static void do_clear(this_type &v);
  };
  
  template <class PT> 
  void linalg_traits<transposed_row_ref<PT> >::do_clear(this_type &v) { 
    col_iterator it = mat_col_begin(v), ite = mat_col_end(v);
    for (; it != ite; ++it) clear(col(it));
  }
  

  // for GCC 2.95
  template <class PT> struct linalg_traits<const transposed_row_ref<PT> > 
  : public linalg_traits<transposed_row_ref<PT> > {}; 

  template<class PT> std::ostream &operator <<
  (std::ostream &o, const transposed_row_ref<PT>& m)
  { gmm::write(o,m); return o; }

  template <class PT> struct  transposed_col_ref {
    
    typedef transposed_col_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef M * CPT;
    typedef typename std::iterator_traits<PT>::reference ref_M;
    typedef typename linalg_traits<this_type>::row_iterator iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<M>::access_type access_type;
    
    iterator _begin, _end;
    const void *origin; 

    transposed_col_ref(ref_M m) : _begin(mat_col_begin(m)),
				  _end(mat_col_end(m)),
				  origin(linalg_origin(m)) {}

    transposed_col_ref(const transposed_col_ref<CPT> &cr) :
      _begin(cr._begin),_end(cr._end),origin(cr.origin) {}

    reference operator()(size_type i, size_type j) const
    { return access_type()(_begin+i, j); }
  };

  template <class PT> struct transposed_col_matrix_access {
    typedef transposed_col_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::row_iterator iterator;
    typedef typename linalg_traits<M>::access_type access_type;
    
    reference operator()(const iterator &itcol, size_type j)
    { return access_type()(itcol, j); }
  };

  template <class PT> struct linalg_traits<transposed_col_ref<PT> > {
    typedef transposed_col_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type M;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<M>::value_type value_type;
    typedef typename select_return<value_type,
            typename linalg_traits<M>::reference, PT>::return_type reference;
    typedef typename linalg_traits<M>::storage_type storage_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type col_iterator;
    typedef abstract_null_type const_col_iterator;
    typedef typename linalg_traits<M>::const_sub_col_type const_sub_row_type;
    typedef typename select_return<const_sub_row_type, typename
            linalg_traits<M>::sub_col_type, PT>::return_type sub_row_type;
    typedef typename linalg_traits<M>::const_col_iterator const_row_iterator;
    typedef typename select_return<const_row_iterator, typename
            linalg_traits<M>::col_iterator, PT>::return_type row_iterator;
    typedef row_major sub_orientation;
    typedef transposed_col_matrix_access<PT> access_type;
    static size_type nrows(const this_type &v) { return v._end - v._begin; }
    static size_type ncols(const this_type &v)
    { return (nrows(v) == 0) ? 0 : vect_size(const_mat_row(v, 0)); }
    static const_sub_row_type row(const const_row_iterator &it)
    { return linalg_traits<M>::col(it); }
    static sub_row_type row(row_iterator &it)
    { return linalg_traits<M>::col(it); }
    static row_iterator row_begin(this_type &m) { return m._begin; }
    static row_iterator row_end(this_type &m) { return m._end; }
    static const_row_iterator row_begin(const this_type &m) { return m._begin; }
    static const_row_iterator row_end(const this_type &m) { return m._end; }
    static const void* origin(const this_type &m) { return m.origin; }
    static void do_clear(this_type &m);
  };

  template <class PT> 
  void linalg_traits<transposed_col_ref<PT> >::do_clear(this_type &v) { 
    row_iterator it = mat_row_begin(v), ite = mat_row_end(v);
    for (; it != ite; ++it) clear(row(it));
  }


  // for GCC 2.95
  template <class PT> struct linalg_traits<const transposed_col_ref<PT> > 
  : public linalg_traits<transposed_col_ref<PT> > {}; 

  template<class PT> std::ostream &operator <<
  (std::ostream &o, const transposed_col_ref<PT>& m)
  { gmm::write(o,m); return o; }

  template <class TYPE, class PT> struct _transposed_return;
  template <class PT> struct _transposed_return<row_major, PT> {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename select_return<transposed_row_ref<const L *>,
            transposed_row_ref< L *>, PT>::return_type return_type;
  };
  template <class PT> struct _transposed_return<col_major, PT> {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename select_return<transposed_col_ref<const L *>,
            transposed_col_ref< L *>, PT>::return_type return_type;
  };
  template <class PT> struct transposed_return {
    typedef typename std::iterator_traits<PT>::value_type L;
    typedef typename _transposed_return<typename principal_orientation_type<
            typename linalg_traits<L>::sub_orientation>::potype,
	    PT>::return_type return_type;
  };

  template <class L> inline 
  typename transposed_return<const L *>::return_type transposed(const L &l)
  { return typename transposed_return<const L *>::return_type(linalg_cast(l));}

  template <class L> inline 
  typename transposed_return<L *>::return_type transposed(L &l)
  { return typename transposed_return<L *>::return_type(linalg_cast(l)); }

}

#endif //  __GMM_TRANSPOSED_H
