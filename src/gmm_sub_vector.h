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

  /* ******************************************************************** */
  /*		reverse indexes                             		  */
  /* ******************************************************************** */

  struct sub_index {

    typedef std::vector<size_t> base_type;
    typedef base_type::const_iterator const_iterator;

    std::vector<size_t> ind;
    std::vector<size_t> rind;

    size_type size(void) const { return ind.size(); }
    size_type index(size_type i) const { return ind[i]; }
    size_type rindex(size_type i) const { return rind[i]; }
   
    const_iterator  begin(void) const { return  ind.begin(); }
    const_iterator    end(void) const { return  ind.end();   }
    const_iterator rbegin(void) const { return rind.begin(); }
    const_iterator   rend(void) const { return rind.end();   }

    sub_index(void) {}
    template <class IT> void init(IT, IT, size_type, abstract_sparse);
    template <class IT> void init(IT, IT, size_type, abstract_plain);
    template <class IT, class L> sub_index(IT it, IT ite,
					       size_type n, const L &)
    { init(it, ite, n, typename linalg_traits<L>::storage_type()); }

  };

  template <class IT>
  void sub_index::init(IT it, IT ite, size_type n, abstract_sparse) {
      rind.resize(n); std::fill(rind.begin(), rind.end(), size_t(-1));
      for (size_type i = 0; it != ite; ++it, ++i) rind[*it] = i;
      ind.resize(ite - it); std::fill(it, ite, ind.begin());
  }

  template <class IT>
  void sub_index::init(IT it, IT ite, size_type n, abstract_plain)
  { ind.resize(ite - it); std::copy(it, ite, ind.begin()); }


  struct sub_interval {
    size_type min, max; 

    size_type size(void) const { return max - min; }
    size_type index(size_type i) const { return min + i; }
    size_type rindex(size_type i) const
    { if (i >= min && i <= max) return i - min; return size_type(-1); }
    sub_interval(size_type mi, size_type ma) : min(mi), max(ma) {}
    sub_interval() {}
  };

  struct sub_slice {
    size_type min, max, N; 

    size_type size(void) const { return (max - min) / N; }
    size_type index(size_type i) const { return min + N * i; }
    size_type rindex(size_type i) const { 
      if (i >= min && i <= max)
	{ size_type j = (i - min); if (j % N == 0) return j / N; }
      return size_type(-1);
    }
    sub_slice(size_type mi, size_type ma, size_type n)
      : min(mi), max(ma), N(n) {}
    sub_slice(void) {}
  };

  /* ********************************************************************* */
  /*		sparse sub-vectors                                         */
  /* ********************************************************************* */

  template <class IT, class SUBI> struct sparse_sub_vector_iterator {

    const SUBI *psi;
    IT itb;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef typename traits_type::pointer           pointer;
    typedef typename traits_type::reference         reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef size_t                                  size_type;
    typedef sparse_sub_vector_iterator<IT, SUBI>    iterator;

    size_type index(void) const { return (*_r_i)[itb.index()]; }
    iterator &operator ++() { ++itb; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    reference operator *() const {
      static value_type zero(0);
      return (index() == size_type(-1)) ? zero : *it;
    }

    bool operator ==(const iterator &i) const { return itb == i.itb; }
    bool operator !=(const iterator &i) const { return !(i == *this); }

    sparse_sub_vector_iterator(void) {}
    sparse_sub_vector_iterator(const IT &it, const SUBI &si)
      : itb(it), psi(&si) {}
  };

  template <class PT, class SUBI> class sparse_sub_vector {
  protected :
    const SUBI *psi;
    PT pv; // garder uniquement un itérateur !!

  public :
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename std::iterator_traits<PT>::reference ref_V;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename vect_ref_type<PT,  V>::access_type access_type;
    size_type size(void) const { si->size(); }
    const SUBI &sindex(void) const { return *psi; }
    ref_V deref(void) const { return *pv; }
    access_type operator[](size_type i) { return (*pv)[psi->index(i)]; }
    value_type operator[](size_type i) const { return (*pv)[psi->index(i)]; }

    sparse_sub_vector(ref_V w, const SUBI &si) : psi(&si), pv(&w) {}
    sparse_sub_vector() {}
  };

  template <class PT, class SUBI>
  struct linalg_traits<sparse_sub_vector<PT, SUBI> > {
    typedef sparse_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef linalg_true is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename linalg_traits<V>::reference_type reference_type;
    typedef sparse_sub_vector_iterator<typename vect_ref_type<PT,  V>::iterator,
	    SUBI> iterator;
    typedef sparse_sub_vector_iterator<typename linalg_traits<V>
    ::const_iterator, SUBI> const_iterator;
    typedef abstract_sparse storage_type;
    size_type size(const this_type &v) { return v.size(); }
    iterator begin(this_type &v)
      { return iterator(vect_begin(v.deref()), v.sindex()); }
    const_iterator begin(const this_type &v)
      { return const_iterator(vect_begin(v.deref()), v.sindex()); }
    iterator end(this_type &v)
      { return iterator(vect_end(v.deref()), v.sindex()); }
    const_iterator end(const this_type &v)
      { return const_iterator(vect_end(v.deref()), v.sindex()); }
    const void* origin(const this_type &v) { return linalg_origin(v.deref()); }
    void do_clear(this_type &v) { std::fill(begin(v), end(v), value_type(0)); }
  };

  /* ******************************************************************** */
  /*		sub vector.                                               */
  /* ******************************************************************** */
  /* sub_vector_type<PT, SUBI>::vector_type is the sub vector type        */
  /* returned by sub_vector(v, sub_index)                                 */
  /************************************************************************/

  template <class PT, class SUBI, class st_type> struct svrt_ir;

  template <class PT>
  struct svrt_ir<PT, sub_index, abstract_plain> {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename vect_ref_type<PT,  V>::iterator iterator;
    typedef tab_ref_index_ref_with_origin<iterator,
      sub_index::const_iterator> vector_type;
  }; 

  template <class PT>
  struct svrt_ir<PT, sub_interval, abstract_plain> {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename vect_ref_type<PT,  V>::iterator iterator;
    typedef tab_ref_with_origin<iterator> vector_type;
  }; 

  template <class PT>
  struct svrt_ir<PT, sub_slice, abstract_plain> {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename vect_ref_type<PT,  V>::iterator iterator;
    typedef tab_ref_reg_spaced_with_origin<iterator> vector_type;
  }; 

  template <class PT, class SUBI>
  struct svrt_ir<PT, SUBI, abstract_sparse> {
    typedef sparse_sub_vector<PT, SUBI> vector_type;
  };

  // -------------- Principal interface ------------------------
  template <class PT, class SUBI>
  struct sub_vector_type {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename svrt_ir<PT, SUBI,
      typename linalg_traits<V>::storage_type>::vector_type vector_type;
  };
  
  // à spécifier suivant le type de la référence

  template <class V, class SUBI> inline
  typename sub_vector_type<const V *, SUBI>::vector_type
  sub_vector(const V &v, const SUBI &si) {
    return sub_vector_stc(v, si, typename linalg_traits<V>::storage_type());
  }

  template <class V, class SUBI>  inline
  typename sub_vector_type<V *, SUBI>::vector_type
  sub_vector(V &v, const SUBI &si) {
    return sub_vector_st(v, si, typename linalg_traits<V>::storage_type());
  }
  // -----------------------------------------------------------

  template <class V> inline
  typename sub_vector_type<const V *, sub_index>::vector_type
  sub_vector_stc(const V &v, const sub_index &si, abstract_plain) {
    return  typename sub_vector_type<const V *, sub_index >
      ::vector_type(vect_begin(v), si.begin(), si.end(), linalg_origin(v));
  }

  template <class V> inline
  typename sub_vector_type<V *, sub_index>::vector_type
  sub_vector_st(V &v, const sub_index &si, abstract_plain) {
    return  typename sub_vector_type<V *, sub_index>
      ::vector_type(vect_begin(v), si.begin(), si.end(), linalg_origin(v));
  }

  template <class V> inline
  typename sub_vector_type<const V *, sub_interval>::vector_type
  sub_vector_stc(const V &v, const sub_interval &si, abstract_plain) {
    return  typename sub_vector_type<const V *, sub_interval>
      ::vector_type(vect_begin(v) + si.min, vect_begin(v) + si.max+1,
		    linalg_origin(v));
  }

  template <class V> inline
  typename sub_vector_type<V *, sub_interval>::vector_type
  sub_vector_st(V &v, const sub_interval &si, abstract_plain) {
    return  typename sub_vector_type<V *, sub_interval>
      ::vector_type(vect_begin(v) + si.min, vect_begin(v) + si.max+1,
		    linalg_origin(v));
  }

  template <class V> inline
  typename sub_vector_type<const V *, sub_slice>::vector_type
  sub_vector_stc(const V &v, const sub_slice &si, abstract_plain) {
    return  typename sub_vector_type<const V *, sub_slice>
      ::vector_type(vect_begin(v), v.begin() + si.min, v.begin() + si.max+1,
		    si.N, linalg_origin(v));
  }

  template <class V> inline
  typename sub_vector_type<V *, sub_slice>::vector_type
  sub_vector_st(V &v, const sub_slice &si, abstract_plain) {
    return  typename sub_vector_type<V *, sub_slice>
      ::vector_type(vect_begin(v), v.begin() + si.min, v.begin() + si.max+1,
		    si.N, linalg_origin(v));
  }

  template <class V, class SUBI> inline
  typename sub_vector_type<const V *, SUBI>::vector_type
  sub_vector_stc(const V &v, const SUBI &si, abstract_sparse) {
    return typename sub_vector_type<const V *, SUBI>::vector_type(v, si);
  }

  template <class V, class SUBI> inline
  typename sub_vector_type<V *, SUBI>::vector_type
  sub_vector_st(V &v, const SUBI &si, abstract_sparse) {
    return typename sub_vector_type<V *, SUBI>::vector_type(v, si);
  }

}

#endif //  __GMM_SUB_VECTOR_H
