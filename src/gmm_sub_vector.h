/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_sub_vector.h : generic sub vectors.                      */
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

#ifndef __GMM_SUB_VECTOR_H
#define __GMM_SUB_VECTOR_H

namespace gmm {

  /* ********************************************************************* */
  /*		sparse sub-vectors                                         */
  /* ********************************************************************* */

  template <class IT> struct sparse_sub_vector_iterator {

    const reverse_index *_r_i;
    IT itb;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef typename traits_type::pointer           pointer;
    typedef typename traits_type::reference         reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef size_t                                  size_type;
    typedef sparse_sub_vector_iterator<IT>          iterator;

    size_type index(void) const { return (*_r_i)[itb.index()]; }
    iterator &operator ++() { ++itb; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    reference operator *() const {
      static value_type zero(0);
      if (index() == size_type(-1)) return zero; else return *itb;
    }

    bool operator ==(const iterator &i) const { return itb == i.itb; }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return itb < i.itb; }

    sparse_sub_vector_iterator(void) {}
    sparse_sub_vector_iterator(const IT &it, const reverse_index &ri)
      : itb(it), _r_i(&ri) {}
    
  };

  template <class PT, class IT> class sparse_sub_vector {
  protected :
    const reverse_index *_r_i;
    IT begin, end;
    PT v;

  public :
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename std::iterator_traits<PT>::reference ref_V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename vect_ref_type<PT,  V>::access_type access_type;
    size_type size(void) const { return end - begin; }
    const reverse_index &rindex(void) const { return *_r_i; }
    ref_V deref(void) const { return *v; }
    access_type operator[](size_type i) { return (*v)[begin[i]]; }
    value_type operator[](size_type i) const { return (*v)[begin[i]]; }

    sparse_sub_vector(ref_V w, const IT &b,
		      const IT &e, const reverse_index &ri)
      : _r_i(&ri), begin(b), end(e), v(&w) {}
    sparse_sub_vector() {}
  };

  template <class PT, class IT>
  struct linalg_traits<sparse_sub_vector<PT, IT> > {
    typedef sparse_sub_vector<PT, IT> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef sparse_sub_vector_iterator<typename vect_ref_type<PT,  V>::iterator>
    iterator;
    typedef sparse_sub_vector_iterator<typename linalg_traits<V>
    ::const_iterator> const_iterator;
    typedef abstract_sparse storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v)
      { return iterator(vect_begin(v.deref()), v.rindex()); }
    const_iterator const_begin(const this_type &v)
      { return const_iterator(vect_begin(v.deref()), v.rindex()); }
    iterator end(this_type &v)
      { return iterator(vect_end(v.deref()), v.rindex()); }
    const_iterator const_end(const this_type &v)
      { return const_iterator(vect_end(v.deref()), v.rindex()); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { std::fill(begin(), end(), value_type(0)); }
  };

  /* ******************************************************************** */
  /*		sub vector with an array of indexes.                      */
  /* ******************************************************************** */

  template <class PT, class IT, class st_type> struct svrt_ir;

  template <class PT, class IT>
  struct svrt_ir<PT, IT, abstract_plain> {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename vect_ref_type<PT,  V>::iterator iterator;
    typedef tab_ref_index_ref_with_origin<iterator, IT> vector_type;
  }; 

  template <class PT, class IT>
  struct svrt_ir<PT, IT, abstract_sparse> {
    typedef sparse_sub_vector<PT, IT> vector_type;
  };

  template <class PT, class IT>
  struct sub_vector_type {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename svrt_ir<PT, IT,
      typename linalg_traits<V>::storage_type>::vector_type vector_type;
  };
  
  template <class V, class IT> inline
  typename sub_vector_type<const V *, IT>::vector_type
  sub_vector(const V &v, const IT &it, const IT &e) {
    return sub_vector_stc(v, it, e,
			  typename linalg_traits<V>::storage_type());
  }

  template <class V, class IT>  inline
  typename sub_vector_type<V *, IT>::vector_type
  sub_vector(V &v, const IT &it, const IT &e) {
    return sub_vector_st(v, it, e,
			 typename linalg_traits<V>::storage_type());
  }

  template <class V, class IT> inline
  typename sub_vector_type<const V *, IT>::vector_type
  sub_vector(const V &v, const IT &it, const IT &e, 
	     const reverse_index &rindex) {
    return sub_vector_stc(v, it, e, 
			  typename linalg_traits<V>::storage_type(), &rindex);
  }

  template <class V, class IT>  inline
  typename sub_vector_type<V *, IT>::vector_type
  sub_vector(V &v, const IT &it, const IT &e, const reverse_index &rindex) {
    return sub_vector_st(v, it, e, 
			 typename linalg_traits<V>::storage_type(), &rindex);
  }

  template <class V, class IT> inline
  typename sub_vector_type<const V *, IT>::vector_type
  sub_vector_stc(const V &v, const IT &it, const IT &e,
		 abstract_plain, const reverse_index * = 0) {
    return  typename sub_vector_type<const V *, IT>
      ::vector_type(vect_begin(v), it, e, linalg_origin(v));
  }

  template <class V, class IT> inline
  typename sub_vector_type<const V *, IT>::vector_type
  sub_vector_stc(const V &v, const IT &it, const IT &e,
		 abstract_sparse, const reverse_index *rindex = 0) {
    return typename sub_vector_type<const V *, IT>
      ::vector_type(v, it, e, *rindex);
  }

  template <class V, class IT> inline
  typename sub_vector_type<V *, IT>::vector_type
  sub_vector_st(V &v, const IT &it, const IT &e,
		abstract_plain, const reverse_index * = 0) {
    return typename sub_vector_type<V *, IT>
      ::vector_type(vect_begin(v), it, e, linalg_origin(v));
  }

  template <class V, class IT> inline
  typename sub_vector_type<V *, IT>::vector_type
  sub_vector_st(V &v, const IT &it, const IT &e,
		abstract_sparse, const reverse_index *rindex = 0) {
    return typename sub_vector_type<V *, IT>
      ::vector_type(v, it, e, *rindex);
  }

}

#endif //  __GMM_SUB_VECTOR_H
