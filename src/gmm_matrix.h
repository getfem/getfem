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

#ifndef __GMM_MATRIX_H
#define __GMM_MATRIX_H

namespace gmm
{

  /* ******************************************************************** */
  /*		Identity matrix                         		  */
  /* ******************************************************************** */

  struct identity_matrix {};

  template <class V1, class V2> inline
  void mult(const identity_matrix&, const V1 &v1, V2 &v2)
  { copy(v1, v2); }
  template <class V1, class V2> inline
  void mult(const identity_matrix&, const V1 &v1, const V2 &v2) 
  { copy(v1, v2); }
  template <class V1, class V2, class V3> inline
  void mult(const identity_matrix&, const V1 &v1, const V2 &v2, V3 &v3)
  { add(v1, v2, v3); }
  template <class V1, class V2, class V3> inline
  void mult(const identity_matrix&, const V1 &v1, const V2 &v2, const V3 &v3)
  { add(v1, v2, v3); }
  template <class M> void copy(const identity_matrix&, M &m) {
    size_type i = 0, n = std::max(mat_nrows(m), mat_ncols(m)); clear(m);
    for (; i < n; ++i) m(i,i) = typename linalg_traits<M>::value_type(1);
  } 
  template <class M> inline void copy(const identity_matrix &, const M &m)
  { copy_ident(identity_matrix(), linalg_const_cast(m)); }
  template <class V1, class V2> inline
  typename linalg_traits<V1>::value_type
  vect_sp(const identity_matrix &, const V1 &v1, const V2 &v2)
  { return vect_sp(v1, v2); }
  template<class M> inline bool is_identity(const M&) { return false; }
  inline bool is_identity(const identity_matrix&) { return true; }

  /* ******************************************************************** */
  /*		Row matrix                                   		  */
  /* ******************************************************************** */

  template<class V> class row_matrix {
  protected :
    std::vector<V> li; /* array of rows.                                   */
    
  public :
    
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::value_type value_type;
    
    row_matrix(size_type r, size_type c) : li(r)
    { for (size_type i = 0; i < r; ++i) li[i] = V(c); }
    row_matrix(void) {}
    reference operator ()(size_type l, size_type c)
    { return li[l][c]; }
    value_type operator ()(size_type l, size_type c) const
    { return li[l][c]; }

    void clear_row(size_type i) { clear(li[i]); }
    void clear_mat();
    void resize(size_type i) { li.resize(i); }

    typename std::vector<V>::iterator begin(void)
    { return li.begin(); }
    typename std::vector<V>::iterator end(void)  
    { return li.end(); }
    typename std::vector<V>::const_iterator begin(void) const
    { return li.begin(); }
    typename std::vector<V>::const_iterator end(void) const
    { return li.end(); }
    
    
    V& row(size_type i) { return li[i]; }
    const V& row(size_type i) const { return li[i]; }
    
    inline size_type nrows(void) const { return li.size(); }
    inline size_type ncols(void) const
    { return (nrows() == 0) ? 0 : li[0].size(); }
  };

  template<class V> void row_matrix<V>::clear_mat()
  { for (size_type i=0; i < nrows(); ++i) clear_row(i); }

  template <class V> struct row_matrix_access {
    typedef typename linalg_traits<row_matrix<V> >::reference reference;
    typedef typename linalg_traits<row_matrix<V> >::row_iterator iterator;
    typedef typename linalg_traits<row_matrix<V> >::value_type value_type;
    typedef typename linalg_traits<row_matrix<V> >::const_row_iterator
          const_iterator;
    
    reference operator()(const iterator &itrow, size_type j)
    { return (*itrow)[j]; }
    value_type operator()(const const_iterator &itrow, size_type j)
    { return (*itrow)[j]; }
  };

  template <class V> struct linalg_traits<row_matrix<V> > {
    typedef row_matrix<V> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef simple_vector_ref<V *> sub_row_type;
    typedef simple_vector_ref<const V *> const_sub_row_type;
    typedef typename std::vector<V>::iterator row_iterator;
    typedef typename std::vector<V>::const_iterator const_row_iterator;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type col_iterator;
    typedef abstract_null_type const_col_iterator;
    typedef row_major sub_orientation;
    typedef row_matrix_access<V> access_type;
    static size_type nrows(const this_type &m) { return m.nrows(); }
    static size_type ncols(const this_type &m) { return m.ncols(); }
    static row_iterator row_begin(this_type &m) { return m.begin(); }
    static row_iterator row_end(this_type &m) { return m.end(); }
    static const_row_iterator row_begin(const this_type &m) { return m.begin(); }
    static const_row_iterator row_end(const this_type &m) { return m.end(); }
    static const_sub_row_type row(const const_row_iterator &it)
    { return const_sub_row_type(*it); }
    static sub_row_type row(const row_iterator &it) 
    { return sub_row_type(*it); }
    static const void* origin(const this_type &m) { return &m; }
    static void do_clear(this_type &m) { m.clear_mat(); }
  };

#ifdef USING_BROKEN_GCC295
  template <class V> struct linalg_traits<const row_matrix<V> >
    : public linalg_traits<row_matrix<V> > {};
#endif

   template<class V> std::ostream &operator <<
  (std::ostream &o, const row_matrix<V>& m) { gmm::write(o,m); return o; }

  /* ******************************************************************** */
  /*		Column matrix                                		  */
  /* ******************************************************************** */

  template<class V> class col_matrix {
  protected :
    std::vector<V> li; /* array of columns.                               */
    
  public :
    
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::value_type value_type;
    
    col_matrix(size_type r, size_type c) : li(c)
    { for (size_type i = 0; i < c; ++i) li[i] = V(r); }
    col_matrix(void) {}
    reference operator ()(size_type l, size_type c)
    { return li[c][l]; }
    value_type operator ()(size_type l, size_type c) const
    { return li[c][l]; }

    void clear_col(size_type i) { clear(li[i]); }
    void clear_mat();
    void resize(size_type i) { li.resize(i); }

    V& col(size_type i) { return li[i]; }
    const V& col(size_type i) const { return li[i]; }

    typename std::vector<V>::iterator begin(void)
    { return li.begin(); }
    typename std::vector<V>::iterator end(void)  
    { return li.end(); }
    typename std::vector<V>::const_iterator begin(void) const
    { return li.begin(); }
    typename std::vector<V>::const_iterator end(void) const
    { return li.end(); }
    
    inline size_type ncols(void) const { return li.size(); }
    inline size_type nrows(void) const
    { return (ncols() == 0) ? 0 : li[0].size(); }
  };

  template<class V> void col_matrix<V>::clear_mat()
  { for (size_type i=0; i < ncols(); ++i) clear_col(i); }

  template <class V> struct col_matrix_access {
    typedef typename linalg_traits<col_matrix<V> >::reference reference;
    typedef typename linalg_traits<col_matrix<V> >::col_iterator iterator;
    typedef typename linalg_traits<col_matrix<V> >::value_type value_type;
    typedef typename linalg_traits<col_matrix<V> >::const_col_iterator
          const_iterator;
    
    reference operator()(const iterator &itcol, size_type j)
    { return (*itcol)[j]; }
    value_type operator()(const const_iterator &itcol, size_type j)
    { return (*itcol)[j]; }
  };

  template <class V> struct linalg_traits<col_matrix<V> > {
    typedef col_matrix<V> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference reference;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef simple_vector_ref<V *> sub_col_type;
    typedef simple_vector_ref<const V *> const_sub_col_type;
    typedef typename std::vector<V>::iterator col_iterator;
    typedef typename std::vector<V>::const_iterator const_col_iterator;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type row_iterator;
    typedef abstract_null_type const_row_iterator;
    typedef col_major sub_orientation;
    typedef col_matrix_access<V> access_type;
    static size_type nrows(const this_type &m) { return m.nrows(); }
    static size_type ncols(const this_type &m) { return m.ncols(); }
    static col_iterator col_begin(this_type &m) { return m.begin(); }
    static col_iterator col_end(this_type &m) { return m.end(); }
    static const_col_iterator col_begin(const this_type &m) { return m.begin(); }
    static const_col_iterator col_end(const this_type &m) { return m.end(); }
    static const_sub_col_type col(const const_col_iterator &it)
    { return const_sub_col_type(*it); }
    static sub_col_type col(const col_iterator &it) 
    { return sub_col_type(*it); }
    static const void* origin(const this_type &m) { return &m; }
    static void do_clear(this_type &m) { m.clear_mat(); }
  };

  template<class V> std::ostream &operator <<
  (std::ostream &o, const col_matrix<V>& m) { gmm::write(o,m); return o; }

#ifdef USING_BROKEN_GCC295
  template <class V> struct linalg_traits<const col_matrix<V> >
    : public linalg_traits<col_matrix<V> > {};
#endif

  /* ******************************************************************** */
  /*		Block matrix                                		  */
  /* ******************************************************************** */

  template <class MAT> class block_matrix {
  protected :
    std::vector<MAT> blocks;
    size_type _nrowblocks;
    size_type _ncolblocks;
    std::vector<sub_interval> introw, intcol;

  public :
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename linalg_traits<MAT>::reference reference;

    size_type nrows(void) const { return introw[_nrowblocks-1].max + 1; }
    size_type ncols(void) const { return intcol[_ncolblocks-1].max + 1; }
    size_type nrowblocks(void) const { return _nrowblocks; }
    size_type ncolblocks(void) const { return _ncolblocks; }
    const sub_interval &subrowinterval(size_type i) const { return introw[i]; }
    const sub_interval &subcolinterval(size_type i) const { return intcol[i]; }
    const MAT &block(size_type i, size_type j) const 
    { return blocks[j*_ncolblocks+i]; }
    MAT &block(size_type i, size_type j)
    { return blocks[j*_ncolblocks+i]; }
    void do_clear(void);
    // to be done : read and write access to a component
    value_type operator() (size_type i, size_type j) const {
      size_type k, l;
      for (k = 0; k < _nrowblocks; ++k)
	if (i >= introw[k].min && i <=  introw[k].max) break;
      for (l = 0; l < _nrowblocks; ++l)
	if (j >= introw[l].min && j <=  introw[l].max) break;
      return (block(k, l))(i - introw[k].min, j - introw[l].min);
    }
    reference operator() (size_type i, size_type j) {
      size_type k, l;
      for (k = 0; k < _nrowblocks; ++k)
	if (i >= introw[k].min && i <=  introw[k].max) break;
      for (l = 0; l < _nrowblocks; ++l)
	if (j >= introw[l].min && j <=  introw[l].max) break;
      return (block(k, l))(i - introw[k].min, j - introw[l].min);
    }
    
    template <class CONT> void resize(const CONT &c1, const CONT &c2);
    template <class CONT> block_matrix(const CONT &c1, const CONT &c2)
    { resize(c1, c2); }
    block_matrix(void) {}

  };

  template <class MAT> struct linalg_traits<block_matrix<MAT> > {
    typedef block_matrix<MAT> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename linalg_traits<MAT>::value_type value_type;
    typedef typename linalg_traits<MAT>::reference reference;
    typedef typename linalg_traits<MAT>::storage_type storage_type;
    typedef abstract_null_type sub_row_type; // to be done ...
    typedef abstract_null_type const_sub_row_type; // to be done ...
    typedef abstract_null_type row_iterator; // to be done ...
    typedef abstract_null_type const_row_iterator; // to be done ...
    typedef abstract_null_type sub_col_type; // to be done ...
    typedef abstract_null_type const_sub_col_type; // to be done ...
    typedef abstract_null_type col_iterator; // to be done ...
    typedef abstract_null_type const_col_iterator; // to be done ...
    typedef abstract_null_type sub_orientation; // to be done ...
    typedef abstract_null_type access_type; // to be done ...
    static size_type nrows(const this_type &m) { return m.nrows(); }
    static size_type ncols(const this_type &m) { return m.ncols(); }
    static const void* origin(const this_type &m) { return &m; }
    static void do_clear(this_type &m) { m.do_clear(); }
  };

  template <class MAT> void block_matrix<MAT>::do_clear(void) { 
    for (size_type j = 0, l = 0; j < _ncolblocks; ++j)
      for (size_type i = 0, k = 0; i < _nrowblocks; ++i)
	clear(block(i,j));
  }

  template <class MAT> template <class CONT>
  void block_matrix<MAT>::resize(const CONT &c1, const CONT &c2) {
    _nrowblocks = c1.size(); _ncolblocks = c2.size();
    blocks.resize(_nrowblocks * _ncolblocks);
    intcol.resize(_ncolblocks);
    introw.resize(_nrowblocks);
    for (size_type j = 0, l = 0; j < _ncolblocks; ++j) {
      intcol[j] = sub_interval(l, c2[j]); l += c2[j];
      for (size_type i = 0, k = 0; i < _nrowblocks; ++i) {
	if (j == 0) { introw[i] = sub_interval(k, c1[i]); k += c1[i]; }
	block(i, j) = MAT(c1[i], c2[j]);
      }
    }
  }

  template <class M1, class M2>
  void copy(const block_matrix<M1> &m1, M2 &m2) {
    for (size_type j = 0; j < m1.ncolblocks(); ++j)
      for (size_type i = 0; i < m1.nrowblocks(); ++i)
	copy(m1.block(i,j), sub_matrix(m2, m1.subrowinterval(i), 
				       m1.subcolinterval(j)));
  }

  template <class M1, class M2>
  void copy(const block_matrix<M1> &m1, const M2 &m2)
  { copy(m1, linalg_const_cast(m2)); }
  

  template <class MAT, class V1, class V2>
  void mult(const block_matrix<MAT> &m, const V1 &v1, V2 &v2) {
    clear(v2);
    typename sub_vector_type<V2 *, sub_interval>::vector_type sv;
    for (size_type i = 0; i < m.nrowblocks() ; ++i)
      for (size_type j = 0; j < m.ncolblocks() ; ++j) {
	sv = sub_vector(v2, m.subrowinterval(i));
	mult(m.block(i,j),
	     sub_vector(v1, m.subcolinterval(j)), sv, sv);
      }
  }

  template <class MAT, class V1, class V2, class V3>
  void mult(const block_matrix<MAT> &m, const V1 &v1, const V2 &v2, V3 &v3) {
    typename sub_vector_type<V3 *, sub_interval>::vector_type sv;
    for (size_type i = 0; i < m.nrowblocks() ; ++i)
      for (size_type j = 0; j < m.ncolblocks() ; ++j) {
	sv = sub_vector(v3, m.subrowinterval(i));
	if (j == 0)
	  mult(m.block(i,j),
	       sub_vector(v1, m.subcolinterval(j)),
	       sub_vector(v2, m.subrowinterval(i)), sv);
	else
	  mult(m.block(i,j),
	       sub_vector(v1, m.subcolinterval(j)), sv, sv);
      }
    
  }

  template <class MAT, class V1, class V2>
  void mult(const block_matrix<MAT> &m, const V1 &v1, const V2 &v2)
  { mult(m, v1, linalg_const_cast(v2)); }

  template <class MAT, class V1, class V2, class V3>
  void mult(const block_matrix<MAT> &m, const V1 &v1, const V2 &v2, 
	    const V3 &v3)
  { mult_const(m, v1, v2, linalg_const_cast(v3)); }

}

#endif /* __GMM_MATRIX_H */
