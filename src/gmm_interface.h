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
  /*   for a plain vector : the minimum is two random iterators (begin and */
  /*                        end) and a const void * pointer to a valid     */
  /*                        origin.                                        */
  /*   for a sparse vector : the minimum is two forward iterators, with    */
  /*                         a method it.index() which gives the index of  */
  /*                         a non zero element, an interface object       */
  /*                         should descibe the method to add new non      */
  /*                         zero element, and  a const void * pointer to  */
  /*                         a valid origin.                               */
  /*                                                                       */
  /* What is needed for a Matrix type :                                    */
  /*   Matrix m(n, m) defines a matrix with n rows and m columns.          */
  /*   m(i, j) allows to access to the element at row i and column j.      */
  /*   linalg_traits<Matrix> should be filled with appropriate definitions */
  /*                                                                       */
  /* What is needed for an iterator on plain vector                        */
  /*    to be standard random access iterator                              */
  /*                                                                       */
  /* What is needed for an iterator on a sparse vector                     */
  /*    to be a standard forward iterator                                  */
  /*    it.index() gives the index of the non-zero element.                */
  /*                                                                       */
  /* Remark : If original iterators are not convenient, they could be      */
  /*   redefined and interfaced in linalg_traits<Vector> without changin   */
  /*   the original Vector type.                                           */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*		Simple references on vectors            		   */
  /* ********************************************************************* */

  template <class PT> struct simple_vector_ref {
    typedef simple_vector_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename std::iterator_traits<PT>::reference ref_V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<this_type>::iterator iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<V>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    size_type _size;

    simple_vector_ref(ref_V v)
      : _begin(vect_begin(v)), _end(vect_end(v)), origin(linalg_origin(v)),
	_size(vect_size(v)) {}

    reference operator[](size_type i) const
    { return access_type()(origin, _begin, _end, i); }
  };

  template <class PT> struct linalg_traits<simple_vector_ref<PT> > {
    typedef simple_vector_ref<PT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename select_return<value_type, typename
            linalg_traits<V>::reference, PT>::return_type reference;
    typedef typename select_return<typename linalg_traits<V>::const_iterator,
	    typename linalg_traits<V>::iterator, PT>::return_type iterator;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::storage_type storage_type;
    typedef typename linalg_traits<V>::access_type access_type;
    typedef typename linalg_traits<V>::clear_type  clear_type;
    size_type size(const this_type &v) { return v._size; }
    iterator begin(this_type &v) { return v._begin; }
    const_iterator begin(const this_type &v) { return v._begin; }
    iterator end(this_type &v) { return v._end; }
    const_iterator end(const this_type &v) { return v._end; }
    const void* origin(const this_type &v) { return v.origin; }
    void do_clear(this_type &v) { clear_type()(v.origin, v._begin, v._end); }
  };

  // for GCC 2.95
  template <class PT> struct linalg_traits<const simple_vector_ref<PT> >
  : public linalg_traits<simple_vector_ref<PT> > {};

  /* ********************************************************************* */
  /*		                                         		   */
  /*		Traits for S.T.L. object                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class T> struct linalg_traits<std::vector<T> > {
    typedef std::vector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };
}
namespace std {
  template <class T> ostream &operator <<
  (std::ostream &o, const vector<T>& m)
  { gmm::write(o,m); return o; }
}
namespace gmm {
  // for GCC 2.95
  template <class T> struct linalg_traits<const std::vector<T> > 
    : public linalg_traits<std::vector<T> > {};

  // to be done :  std::valarray<T> ...

  /* ********************************************************************* */
  /*		                                         		   */
  /*		Traits for dal objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class IT> struct tab_ref_with_origin : public dal::tab_ref<IT> {

    const void *origin;
   
    tab_ref_with_origin(void) {}
    tab_ref_with_origin(const IT &b, const IT &e, const void *p)
      : dal::tab_ref<IT>(b,e), origin(p) {}
    template <class V> tab_ref_with_origin(const V &v, const sub_interval &si)
      : dal::tab_ref<IT>(vect_begin(v)+si.min, vect_begin(v)+si.max+1),
        origin(linalg_origin(v)) {}
    template <class V> tab_ref_with_origin(V &v, const sub_interval &si)
      : dal::tab_ref<IT>(vect_begin(v)+si.min, vect_begin(v)+si.max+1),
        origin(linalg_origin(v)) {}
  };

  template <class IT> struct linalg_traits<tab_ref_with_origin<IT> > {
    typedef typename std::iterator_traits<IT>::pointer PT;
    typedef tab_ref_with_origin<IT> this_type;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::reference reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin; }
    void do_clear(this_type &v) { clear_type()(v.origin, v.begin(), v.end()); }
  };

  template <class IT> std::ostream &operator <<
  (std::ostream &o, const tab_ref_with_origin<IT>& m)
  { gmm::write(o,m); return o; }

  // for GCC 2.95
  template <class IT> struct linalg_traits<const tab_ref_with_origin<IT> >
    : public linalg_traits<tab_ref_with_origin<IT> >{};

  template <class IT>
  struct tab_ref_reg_spaced_with_origin : public dal::tab_ref_reg_spaced<IT> {

    const void *origin;
    
    tab_ref_reg_spaced_with_origin(void) {}
    tab_ref_reg_spaced_with_origin(const IT &b, const IT &e, size_type n,
				   const void *p) 
      : dal::tab_ref_reg_spaced<IT>(b,e,n), origin(p) {}
    template <class V> tab_ref_reg_spaced_with_origin(const V &v,
						      const sub_slice &si) :
       dal::tab_ref_reg_spaced<IT>(vect_begin(v) + si.min,
				   vect_begin(v) + si.max+si.N,
				   si.N), origin(linalg_origin(v)) {}
    template <class V> tab_ref_reg_spaced_with_origin(V &v,
						      const sub_slice &si) :
       dal::tab_ref_reg_spaced<IT>(vect_begin(v) + si.min,
				   vect_begin(v) + si.max+si.N,
				   si.N), origin(linalg_origin(v)) {}

  };

  template <class IT> 
    struct linalg_traits<tab_ref_reg_spaced_with_origin<IT> > {
    typedef typename std::iterator_traits<IT>::pointer PT;
    typedef tab_ref_reg_spaced_with_origin<IT> this_type;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::reference reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin; }
    void do_clear(this_type &v) { clear_type()(v.origin, v.begin(), v.end()); }
  };
  
  template <class IT> std::ostream &operator <<
  (std::ostream &o, const tab_ref_reg_spaced_with_origin<IT>& m)
  { gmm::write(o,m); return o; }

  template <class IT> // for gcc 2.95
  struct linalg_traits<const tab_ref_reg_spaced_with_origin<IT> >
    : public linalg_traits<tab_ref_reg_spaced_with_origin<IT> > {};

  template <class IT, class ITINDEX>
  struct tab_ref_index_ref_with_origin 
    : public dal::tab_ref_index_ref<IT, ITINDEX> {

    const void *origin;

    tab_ref_index_ref_with_origin(void) {}
    tab_ref_index_ref_with_origin(const IT &b, const ITINDEX &bi,
				  const ITINDEX &ei, const void *p)
      : dal::tab_ref_index_ref<IT, ITINDEX>(b, bi, ei), origin(p) {}

    template <class V> tab_ref_index_ref_with_origin(const V &v,
						     const sub_index &si) :
      dal::tab_ref_index_ref<IT, ITINDEX>(vect_begin(v),
					  si.ind.begin(), si.ind.end()),
                                          origin(linalg_origin(v)) {}
    template <class V> tab_ref_index_ref_with_origin(V &v,
						     const sub_index &si) :
      dal::tab_ref_index_ref<IT, ITINDEX>(vect_begin(v),
					  si.ind.begin(), si.ind.end()),
                                          origin(linalg_origin(v)) {}
  };

  template <class IT, class ITINDEX>
  struct linalg_traits<tab_ref_index_ref_with_origin<IT, ITINDEX> > {
    typedef typename std::iterator_traits<IT>::pointer PT;
    typedef tab_ref_index_ref_with_origin<IT, ITINDEX> this_type;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename std::iterator_traits<IT>::value_type value_type;
    typedef typename std::iterator_traits<IT>::reference reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.origin; }
    void do_clear(this_type &v) { clear_type()(v.origin, v.begin(), v.end()); }
  };

  template <class IT, class ITINDEX> std::ostream &operator <<
  (std::ostream &o, const tab_ref_index_ref_with_origin<IT, ITINDEX>& m)
  { gmm::write(o,m); return o; }

  // for GCC 2.95
  template <class IT, class ITINDEX>
  struct linalg_traits<const tab_ref_index_ref_with_origin<IT, ITINDEX> >
    : public  linalg_traits<tab_ref_index_ref_with_origin<IT, ITINDEX> > {};

  /* ********************************************************************* */
  /*		                                         	 	   */
  /*		Traits for bgeot objects                     		   */
  /*		                                         		   */
  /* ********************************************************************* */

  template <class T, int N> struct linalg_traits<bgeot::fsvector<T, N> > {
    typedef bgeot::fsvector<T, N> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

  // for GCC 2.95
  template <class T, int N> struct linalg_traits<const bgeot::fsvector<T, N> >
    : public linalg_traits<bgeot::fsvector<T, N> > {};

  template <class T> struct linalg_traits<bgeot::vsvector<T> > {
    typedef bgeot::vsvector<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_plain storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

  // for GCC 2.95
  template <class T> struct linalg_traits<const bgeot::vsvector<T> > 
    : public linalg_traits<bgeot::vsvector<T> > {};

  template <class VECT> struct linalg_traits<bgeot::PT<VECT> > {
    typedef bgeot::PT<VECT> this_type;
    typedef linalg_false is_reference;
    typedef typename linalg_traits<VECT>::value_type value_type;
    typedef typename linalg_traits<VECT>::reference reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<VECT>::iterator  iterator;
    typedef typename linalg_traits<VECT>::const_iterator const_iterator;
    typedef typename linalg_traits<VECT>::storage_type storage_type;
    typedef plain_access<this_type> access_type;
    typedef plain_clear<this_type> clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return vect_begin(v); }
    const_iterator begin(const this_type &v) { return vect_begin(v); }
    iterator end(this_type &v) { return vect_end(v); }
    const_iterator end(const this_type &v) { return vect_end(v); }
    const void* origin(const this_type &v) { return &v; }
    void do_clear(this_type &v) { clear_type()(origin(v), begin(v), end(v)); }
  };

  // for GCC 2.95
  template <class VECT> struct linalg_traits<const bgeot::PT<VECT> >
  : public linalg_traits<bgeot::PT<VECT> > {};

  template<class ITER, class MIT> struct plain_compressed_iterator
  {
    typedef ITER value_type;
    typedef ITER *pointer;
    typedef ITER &reference;
    typedef ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t size_type;
    typedef plain_compressed_iterator<ITER, MIT> iterator;

    ITER it;
    size_type N, nrows, ncols;
    const void *origin;
    
    iterator operator ++(int) { iterator tmp = *this; it += N; return tmp; }
    iterator operator --(int) { iterator tmp = *this; it -= N; return tmp; }
    iterator &operator ++()   { it += N; return *this; }
    iterator &operator --()   { it -= N; return *this; }
    iterator &operator +=(difference_type i) { it += i * N; return *this; }
    iterator &operator -=(difference_type i) { it -= i * N; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const
    { return (it - i.it) / N; }

    ITER operator *() const { return it; }
    ITER operator [](int ii) const { return it + ii * N; }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    plain_compressed_iterator(void) {}
    plain_compressed_iterator(const plain_compressed_iterator<MIT, MIT> &ii)
    : it(ii.it), N(ii.N), nrows(ii.nrows),ncols(ii.ncols),origin(ii.origin) {}
    plain_compressed_iterator(const ITER &iter, size_type n, size_type r,
			      size_type c, const void *o)
      : it(iter), N(n), nrows(r), ncols(c), origin(o) { }
    
  };

  template <class T, int N> struct bgeot_fsmatrix_access {
    typedef typename linalg_traits<bgeot::fsmatrix<T, N> >::reference
            reference;
    typedef typename linalg_traits<bgeot::fsmatrix<T, N> >::col_iterator
            iterator;
    
    reference operator()(const iterator &itcol, size_type j)
    { return (*itcol)[j]; }
  };

  // row and col iterators could be more specialized to take into account 
  // the fixed size.
  template <class T, int N> struct linalg_traits<bgeot::fsmatrix<T, N> > {
    typedef bgeot::fsmatrix<T, N> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef abstract_plain storage_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type::iterator>
            sub_row_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type
            ::const_iterator> const_sub_row_type;
    typedef plain_compressed_iterator<typename this_type::iterator,
            typename this_type::iterator> row_iterator;
    typedef plain_compressed_iterator<typename this_type::const_iterator, 
            typename this_type::iterator> const_row_iterator;
    typedef tab_ref_with_origin<typename this_type::iterator> sub_col_type;
    typedef tab_ref_with_origin<typename this_type::const_iterator>
            const_sub_col_type;
    typedef plain_compressed_iterator<typename this_type::iterator,
            typename this_type::iterator> col_iterator;
    typedef plain_compressed_iterator<typename this_type::const_iterator, 
            typename this_type::iterator> const_col_iterator;
    typedef col_and_row sub_orientation;
    typedef bgeot_fsmatrix_access<T, N> access_type;
    size_type nrows(const this_type &m) { return N; }
    size_type ncols(const this_type &m) { return N; } 
    const_sub_row_type row(const const_row_iterator &it)
    { return const_sub_row_type(it.it, it.it + N * N, N, it.origin);  }
    const_sub_col_type col(const const_col_iterator &it)
    { return const_sub_col_type(it.it, it.it + N, it.origin); }
    sub_row_type row(const row_iterator &it) 
    { return sub_row_type(it.it, it.it + N * N, it.nrows, it.origin); }
    sub_col_type col(const col_iterator &it)
    { return sub_col_type(it.it, it.it + N, it.origin); }
    row_iterator row_begin(this_type &m)
    { return row_iterator(m.begin(), 1, N, N, &m); }
    row_iterator row_end(this_type &m)
    { return row_iterator(m.begin() + N, 1, N, N, &m); }
    const_row_iterator row_begin(const this_type &m)
    { return const_row_iterator(m.begin(), 1, N, N, &m); }
    const_row_iterator row_end(const this_type &m)
    { return const_row_iterator(m.begin() + N, 1, N, N, &m); }
    col_iterator col_begin(this_type &m)
    { return col_iterator(m.begin(), N, N, N, &m); }
    col_iterator col_end(this_type &m)
    { return col_iterator(m.end(), N, N, N, &m); }
    const_col_iterator col_begin(const this_type &m)
    { return const_col_iterator(m.begin(), N, N, N, &m); }
    const_col_iterator col_end(const this_type &m)
    { return const_col_iterator(m.end(), N, N, N, &m); }
    const void* origin(const this_type &m) { return &m; }
    void do_clear(this_type &m);
  };
  template <class T, int N>
  void linalg_traits<bgeot::fsmatrix<T, N> >::do_clear(this_type &m)
  { std::fill(m.begin(), m.end(), value_type(0)); }


  // for GCC 2.95
  template <class T, int N> struct linalg_traits<const bgeot::fsmatrix<T, N> >
    : public linalg_traits<bgeot::fsmatrix<T, N> > {};

  template <class T> struct bgeot_vsmatrix_access {
    typedef typename linalg_traits<bgeot::vsmatrix<T> >::reference reference;
    typedef typename linalg_traits<bgeot::vsmatrix<T> >::col_iterator iterator;
    
    reference operator()(const iterator &itcol, size_type j)
    { return (*itcol)[j]; }
  };

  template <class T> struct linalg_traits<bgeot::vsmatrix<T> > {
    typedef bgeot::vsmatrix<T> this_type;
    typedef linalg_false is_reference;
    typedef abstract_matrix linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef abstract_plain storage_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type::iterator>
            sub_row_type;
    typedef tab_ref_reg_spaced_with_origin<typename this_type
            ::const_iterator> const_sub_row_type;
    typedef plain_compressed_iterator<typename this_type::iterator,
	    typename this_type::iterator> row_iterator;
    typedef plain_compressed_iterator<typename this_type::const_iterator,
	    typename this_type::iterator> const_row_iterator;
    typedef tab_ref_with_origin<typename this_type::iterator> sub_col_type;
    typedef tab_ref_with_origin<typename this_type::const_iterator>
            const_sub_col_type;
    typedef plain_compressed_iterator<typename this_type::iterator,
	    typename this_type::iterator> col_iterator;
    typedef plain_compressed_iterator<typename this_type::const_iterator,
	    typename this_type::iterator> const_col_iterator;
    typedef col_and_row sub_orientation;
    typedef bgeot_vsmatrix_access<T> access_type;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_sub_row_type row(const const_row_iterator &it) {
      return const_sub_row_type(it.it, it.it + it.ncols * it.nrows,
				it.nrows, it.origin); 
    }
    const_sub_col_type col(const const_col_iterator &it)
    { return const_sub_col_type(it.it, it.it + it.nrows, it.origin); }
    sub_row_type row(const row_iterator &it) {
      return sub_row_type(it.it, it.it + it.ncols * it.nrows,
			  it.nrows, it.origin);
    }
    sub_col_type col(const col_iterator &it)
    { return sub_col_type(it.it, it.it + it.nrows, it.origin); }
    row_iterator row_begin(this_type &m)
    { return row_iterator(m.begin(), 1, m.nrows(), m.ncols(), &m); }
    row_iterator row_end(this_type &m)
    { return row_iterator(m.begin()+m.nrows(), 1, m.nrows(), m.ncols(), &m); }
    const_row_iterator row_begin(const this_type &m)
    { return const_row_iterator(m.begin(), 1, m.nrows(), m.ncols(), &m); }
    const_row_iterator row_end(const this_type &m) {
      return const_row_iterator(m.begin()+m.nrows(), 1, m.nrows(),
				m.ncols(), &m);
    }
    col_iterator col_begin(this_type &m)
    { return col_iterator(m.begin(), m.nrows(), m.nrows(), m.ncols(), &m); }
    col_iterator col_end(this_type &m)
    { return col_iterator(m.end(), m.nrows(), m.nrows(), m.ncols(), &m); }
    const_col_iterator col_begin(const this_type &m)
    { return const_col_iterator(m.begin(),m.nrows(),m.nrows(),m.ncols(),&m); }
    const_col_iterator col_end(const this_type &m)
    { return const_col_iterator(m.end(), m.nrows(),m.nrows(),m.ncols(), &m); }
    const void* origin(const this_type &m) { return &m; }
    void do_clear(this_type &m) { m.fill(value_type(0)); }
  };

  // for GCC 2.95
  template <class T> struct linalg_traits<const bgeot::vsmatrix<T> >
    : public linalg_traits<bgeot::vsmatrix<T> > {};


  /* ******************************************************************** */
  /*	    Read only reference on a compressed sparse vector             */
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
    
  template <class PT1, class PT2, int shift = 0> struct cs_vector_ref {
    PT1 pr;
    PT2 ir;
    size_type n, _size;

    typedef cs_vector_ref<PT1, PT2, shift> this_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename linalg_traits<this_type>::const_iterator const_iterator;

    cs_vector_ref(PT1 pt1, PT2 pt2, size_type nnz, size_type ns)
      : pr(pt1), ir(pt2), n(nnz), _size(ns) {}
    cs_vector_ref(void) {}

    size_type size(void) const { return _size; }
    
    const_iterator begin(void) const { return const_iterator(pr, ir); }
    const_iterator end(void) const { return const_iterator(pr+n, ir+n); }
    
    value_type operator[](size_type i) const
    { typename linalg_traits<this_type>::access_type()(pr, begin(), end(),i); }
  };

  template <class PT1, class PT2, int shift> struct cs_vector_access {
    typedef cs_vector_ref<PT1, PT2, shift> V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::const_iterator const_iterator;
    
    value_type operator()(const void *, const const_iterator &b,
			  const const_iterator &e, size_type i);
  };

  template <class PT1, class PT2, int shift>
  typename linalg_traits<cs_vector_ref<PT1, PT2, shift> >::value_type
  cs_vector_access<PT1, PT2, shift>::operator()(const void *,
                          const const_iterator &b,
			  const const_iterator &e, size_type i) {
    static value_type zero(0);
    if (n == 0) return zero;
    PT2 p = std::lower_bound(b.ir, e.ir, i+shift);
    return (p != b.ir && *p == i+shift) ? b.pr[p-b.ir] : zero;
  }

  template <class PT1, class PT2, int shift>
  struct linalg_traits<cs_vector_ref<PT1, PT2, shift> > {
    typedef cs_vector_ref<PT1, PT2, shift> this_type;
    typedef linalg_const is_reference;
    typedef abstract_vector linalg_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::value_type reference;
    typedef cs_vector_ref_iterator<typename const_pointer<PT1>::pointer,
	    typename const_pointer<PT2>::pointer, shift>  iterator;
    typedef iterator const_iterator;
    typedef abstract_sparse storage_type;
    typedef cs_vector_access<PT1, PT2, shift> access_type;
    typedef abstract_null_type clear_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v) { return v.begin(); }
    const_iterator begin(const this_type &v) { return v.begin(); }
    iterator end(this_type &v) { return v.end(); }
    const_iterator end(const this_type &v) { return v.end(); }
    const void* origin(const this_type &v) { return v.pr; }
  };

  template <class PT1, class PT2, int shift>
  std::ostream &operator <<
  (std::ostream &o, const cs_vector_ref<PT1, PT2, shift>& m)
  { gmm::write(o,m); return o; }

  /* ******************************************************************** */
  /*	    Read only reference on a compressed sparse column matrix      */
  /* ******************************************************************** */

    template <class PT1, class PT2, class PT3, int shift = 0>
  struct sparse_compressed_iterator {
    typedef PT1 value_type;
    typedef PT1 *pointer;
    typedef PT1 &reference;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef sparse_compressed_iterator<PT1, PT2, PT3, shift> iterator;

    PT1 pr;
    PT2 ir;
    PT3 jc;
    size_type n;
    const void *origin;
    
    iterator operator ++(int) { iterator tmp = *this; jc++; return tmp; }
    iterator operator --(int) { iterator tmp = *this; jc--; return tmp; }
    iterator &operator ++()   { jc++; return *this; }
    iterator &operator --()   { jc--; return *this; }
    iterator &operator +=(difference_type i) { jc += i; return *this; }
    iterator &operator -=(difference_type i) { jc -= i; return *this; }
    iterator operator +(difference_type i) const 
    { iterator itt = *this; return (itt += i); }
    iterator operator -(difference_type i) const
    { iterator itt = *this; return (itt -= i); }
    difference_type operator -(const iterator &i) const { return jc - i.jc; }

    reference operator *() const { return pr + *jc - shift; }
    reference operator [](int ii) { return pr + *(jc+ii) - shift; }

    bool operator ==(const iterator &i) const { return (jc == i.jc); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (jc < i.jc); }

    sparse_compressed_iterator(void) {}
    sparse_compressed_iterator(PT1 p1,PT2 p2,PT3 p3,size_type nn,const void *o)
      : pr(p1), ir(p2), jc(p3), n(nn), origin(o) { }
    
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
    
    size_type nrows(void) const { return nr; }
    size_type ncols(void) const { return nc; }
   
    value_type operator()(size_type i, size_type j) const
      { return mat_col(*this, j)[i]; }
  };

  template <class PT1, class PT2, class PT3, int shift>
  struct csc_matrix_access {
    typedef csc_matrix_ref<PT1, PT2, PT3, shift> this_type;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::col_iterator iterator;
    
    reference operator()(const iterator &itcol, size_type j)
    { return linalg_traits<this_type>().col(itcol)[j]; }
  };

  template <class PT1, class PT2, class PT3, int shift>
  struct linalg_traits<csc_matrix_ref<PT1, PT2, PT3, shift> > {
    typedef csc_matrix_ref<PT1, PT2, PT3, shift> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::value_type reference;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_row_type;
    typedef abstract_null_type const_sub_row_type;
    typedef abstract_null_type row_iterator;
    typedef abstract_null_type const_row_iterator;
    typedef cs_vector_ref<typename const_pointer<PT1>::pointer,
	    typename const_pointer<PT2>::pointer, shift> sub_col_type;
    typedef cs_vector_ref<typename const_pointer<PT1>::pointer,
            typename const_pointer<PT2>::pointer, shift> const_sub_col_type;
    typedef sparse_compressed_iterator<typename const_pointer<PT1>::pointer,
				       typename const_pointer<PT2>::pointer,
				       typename const_pointer<PT3>::pointer,
				       shift>  col_iterator;
    typedef col_iterator const_col_iterator;
    typedef csc_matrix_access<PT1, PT2, PT3, shift> access_type;
    typedef col_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_col_iterator col_begin(const this_type &m)
    { return const_col_iterator(m.pr, m.ir, m.jc, m.nr, m.pr); }
    const_col_iterator col_end(const this_type &m)
    { return const_col_iterator(m.pr, m.ir, m.jc + m.nc, m.nr, m.pr); }
    const_sub_col_type col(const const_col_iterator &it) {
      return const_sub_col_type(it.pr + *(it.jc) - shift,
	     it.ir + *(it.jc) - shift, *(it.jc + 1) - *(it.jc), it.n);
    }
    const void* origin(const this_type &m) { return m.pr; }
    void do_clear(this_type &m) { m.do_clear(); }
  };

  template <class PT1, class PT2, class PT3, int shift>
  std::ostream &operator <<
  (std::ostream &o, const csc_matrix_ref<PT1, PT2, PT3, shift>& m)
  { gmm::write(o,m); return o; }

  /* ******************************************************************** */
  /*	   Read only reference on a compressed sparse row matrix          */
  /* ******************************************************************** */

  template <class PT1, class PT2, class PT3, int shift = 0>
  struct csr_matrix_ref {
    PT1 pr;
    PT2 ir;
    PT3 jc;
    size_type nc, nr;
    
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::reference access_type;
    csr_matrix_ref(PT1 pt1, PT2 pt2, PT3 pt3, size_type nrr, size_type ncc)
      : pr(pt1), ir(pt2), jc(pt3), nc(ncc), nr(nrr) {}
    csr_matrix_ref(void) {}
    
    size_type nrows(void) const { return nr; }
    size_type ncols(void) const { return nc; }
   
    value_type operator()(size_type i, size_type j) const
      { return mat_col(*this, i)[j]; }
  };

  template <class PT1, class PT2, class PT3, int shift>
  struct csr_matrix_access {
    typedef csr_matrix_ref<PT1, PT2, PT3, shift> this_type;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::row_iterator iterator;
    
    reference operator()(const iterator &itrow, size_type i)
    { return linalg_traits<this_type>().row(itrow)[i]; }
  };
  
  template <class PT1, class PT2, class PT3, int shift>
  struct linalg_traits<csr_matrix_ref<PT1, PT2, PT3, shift> > {
    typedef csr_matrix_ref<PT1, PT2, PT3, shift> this_type;
    typedef linalg_const is_reference;
    typedef abstract_matrix linalg_type;
    typedef typename std::iterator_traits<PT1>::value_type value_type;
    typedef typename std::iterator_traits<PT1>::value_type reference;
    typedef abstract_sparse storage_type;
    typedef abstract_null_type sub_col_type;
    typedef abstract_null_type const_sub_col_type;
    typedef abstract_null_type col_iterator;
    typedef abstract_null_type const_col_iterator;
    typedef cs_vector_ref<typename const_pointer<PT1>::pointer,
      typename const_pointer<PT2>::pointer, shift> sub_row_type;
    typedef sub_row_type const_sub_row_type;
    typedef sparse_compressed_iterator<typename const_pointer<PT1>::pointer,
				       typename const_pointer<PT2>::pointer,
				       typename const_pointer<PT3>::pointer,
				       shift>  row_iterator;
    typedef row_iterator const_row_iterator;
    typedef row_major sub_orientation;
    size_type nrows(const this_type &m) { return m.nrows(); }
    size_type ncols(const this_type &m) { return m.ncols(); }
    const_row_iterator row_begin(const this_type &m)
    { return const_row_iterator(m.pr, m.ir, m.jc, m.nc, m.pr); }
    const_row_iterator row_end(const this_type &m)
    { return const_row_iterator(m.pr, m.ir, m.jc + m.nr, m.nc, m.pr); }
    const_sub_row_type row(const const_row_iterator &it) {
      return const_sub_row_type(it.pr + *(it.jc) - shift,
	     it.ir + *(it.jc) - shift, *(it.jc + 1) - *(it.jc), it.n);
    }
    const void* origin(const this_type &m) { return m.pr; }
    void do_clear(this_type &m) { m.do_clear(); }
  };

  template <class PT1, class PT2, class PT3, int shift>
  std::ostream &operator <<
  (std::ostream &o, const csr_matrix_ref<PT1, PT2, PT3, shift>& m)
  { gmm::write(o,m); return o; }


}


#endif //  __GMM_INTERFACE_H
