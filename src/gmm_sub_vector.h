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

#ifndef __GMM_SUB_VECTOR_H
#define __GMM_SUB_VECTOR_H

namespace gmm {

  /* ********************************************************************* */
  /*		sparse sub-vectors                                         */
  /* ********************************************************************* */

  template <class IT, class MIT, class SUBI>
  struct sparse_sub_vector_iterator {

    IT itb, itbe;
    SUBI si;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef typename traits_type::pointer           pointer;
    typedef typename traits_type::reference         reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef size_t                                  size_type;
    typedef sparse_sub_vector_iterator<IT, MIT, SUBI>    iterator;

    size_type index(void) const { return si.rindex(itb.index()); }
    void forward(void);
    iterator &operator ++()
    { ++itb; forward(); return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    reference operator *() const { return *itb; }

    bool operator ==(const iterator &i) const { return itb == i.itb; }
    bool operator !=(const iterator &i) const { return !(i == *this); }

    sparse_sub_vector_iterator(void) {}
    sparse_sub_vector_iterator(const IT &it, const IT &ite, const SUBI &s)
      : itb(it), itbe(ite), si(s) { forward(); }
    sparse_sub_vector_iterator(const sparse_sub_vector_iterator<MIT, MIT,
	 SUBI> &it) : itb(it.itb), itbe(it.itbe), si(it.si) {}
  };
  
  template <class IT, class MIT, class SUBI>
  void  sparse_sub_vector_iterator<IT, MIT, SUBI>::forward(void)
  { while(itb!=itbe && index()==size_type(-1)) ++itb; }

  template <class PT, class SUBI> struct sparse_sub_vector {
    typedef sparse_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef V * CPT;
    typedef typename select_return<typename linalg_traits<V>::const_iterator,
            typename linalg_traits<V>::iterator, PT>::return_type iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<V>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    SUBI si;

    size_type size(void) const { return si.size(); }
   
    reference operator[](size_type i) const
    { return access_type()(origin, _begin, _end, si.index(i)); }

    sparse_sub_vector(V &v, const SUBI &s) : _begin(vect_begin(v)),
       _end(vect_end(v)), origin(linalg_origin(v)), si(s) {}
    sparse_sub_vector(const V &v, const SUBI &s) 
      : _begin(vect_begin(const_cast<V &>(v))),
       _end(vect_end(const_cast<V &>(v))),
	origin(linalg_origin(const_cast<V &>(v))), si(s) {}
    sparse_sub_vector() {}
    sparse_sub_vector(const sparse_sub_vector<CPT, SUBI> &cr)
      : _begin(cr._begin),_end(cr._end),origin(cr.origin), si(cr.si) {}
  };

  template <class PT, class SUBI> struct sparse_sub_vector_access {
    typedef sparse_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename linalg_traits<this_type>::value_type value_type;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::iterator iterator;
    typedef typename linalg_traits<this_type>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::access_type access_type;
    
    reference operator()(const void *o, const iterator &it,
			 const iterator &ite, size_type i)
    { return access_type()(o, it.itb, ite.itb, it.si.index(i)); }
    
    value_type operator()(const void *o, const const_iterator &it,
			 const const_iterator &ite, size_type i)
    { return access_type()(o, it.itb, ite.itb, it.si.index(i)); }
  };

  template <class PT, class SUBI> struct sparse_sub_vector_clear {
    typedef sparse_sub_vector<PT, SUBI> this_type;
    typedef typename linalg_traits<this_type>::iterator iterator;
    typedef typename linalg_traits<this_type>::value_type value_type;
    typedef typename linalg_traits<this_type>::access_type access_type;
    
    void operator()(const void *o,const iterator &_begin,const iterator &_end);
  };

  template <class PT, class SUBI> void
  sparse_sub_vector_clear<PT, SUBI>::operator()(const void *o,
		      const iterator &_begin,const iterator &_end) {
    std::vector<size_type> tab(100);
    iterator it = _begin; size_type i = 0;
    for (; it != _end; ++it)
      { tab[i++] = it.index(); if (i >= tab.size()) tab.resize(i + 100); }
    for (size_type j = 0; j < i; ++j)
      access_type()(o, _begin, _end, tab[j]) = value_type(0);
  }

  template <class PT, class SUBI>
  struct linalg_traits<sparse_sub_vector<PT, SUBI> > {
    typedef sparse_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename select_return<value_type, typename
            linalg_traits<V>::reference, PT>::return_type reference;
    typedef typename select_return<typename linalg_traits<V>::const_iterator,
	    typename linalg_traits<V>::iterator, PT>::return_type pre_iterator;
    typedef typename select_return<abstract_null_type, 
	    sparse_sub_vector_iterator<pre_iterator, pre_iterator, SUBI>,
	    PT>::return_type iterator;
    typedef sparse_sub_vector_iterator<typename linalg_traits<V>
            ::const_iterator, pre_iterator, SUBI> const_iterator;
    typedef abstract_sparse storage_type;
    typedef sparse_sub_vector_access<PT, SUBI> access_type;
    typedef sparse_sub_vector_clear<PT, SUBI> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v)
    { return iterator(v._begin, v._end, v.si); }
    static const_iterator begin(const this_type &v)
    { return const_iterator(v._begin, v._end, v.si); }
    static iterator end(this_type &v)
    { return iterator(v._end, v._end, v.si); }
    static const_iterator end(const this_type &v)
    { return const_iterator(v._end, v._end, v.si); }
    static const void* origin(const this_type &v) { return v.origin; }
    static void do_clear(this_type &v)
      { clear_type()(v.origin, begin(v), end(v)); }
  };

  template <class PT, class SUBI> std::ostream &operator <<
  (std::ostream &o, const sparse_sub_vector<PT, SUBI>& m)
  { gmm::write(o,m); return o; }

#ifdef USING_BROKEN_GCC295
  template <class PT, class SUBI>
  struct linalg_traits<const sparse_sub_vector<PT, SUBI> >
    : public linalg_traits<sparse_sub_vector<PT, SUBI> > {};
#endif

  /* ********************************************************************* */
  /*		skyline sub-vectors                                        */
  /* ********************************************************************* */

  template <class IT, class MIT>
  struct skyline_sub_vector_iterator {

    IT it;
    int N;
    size_type id;

    typedef std::iterator_traits<IT>                traits_type;
    typedef typename traits_type::value_type        value_type;
    typedef typename traits_type::pointer           pointer;
    typedef typename traits_type::reference         reference;
    typedef typename traits_type::difference_type   difference_type;
    typedef std::random_access_iterator_tag         iterator_category;
    //    typedef size_t                                  size_type;
    typedef skyline_sub_vector_iterator<IT, MIT>    iterator;

    size_type index(void) const { return id; }
    
    iterator &operator ++() { it += N; id++; return *this; }
    iterator operator ++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator &operator --() { it -= N; id--; return *this; }
    iterator operator --(int) { iterator tmp = *this; --(*this); return tmp; }
    iterator &operator +=(difference_type i)
    { it += N * i; id += i; return *this; }
    iterator &operator -=(difference_type i)
    { it -= N * i; id -= i; return *this; }
    iterator operator +(difference_type i) const
    { iterator ii = *this; return (ii += i); }
    iterator operator -(difference_type i) const
    { iterator ii = *this; return (ii -= i); }
    difference_type operator -(const iterator &i) const
    { return (it - i.it) / N; }

    reference operator *() const { return *it; }
    reference operator [](int ii) { return *(it + ii * N); }

    bool operator ==(const iterator &i) const { return it == i.it;    }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return it  < i.it;    }
    skyline_sub_vector_iterator(void) {}
    skyline_sub_vector_iterator(const IT &i, int NN, size_type j)
      : it(i), N(NN), id(j) { }
    skyline_sub_vector_iterator(const skyline_sub_vector_iterator<MIT, MIT>
				&i) : it(i.it), N(i.N), id(i.id) {}
  };

  template <class PT, class SUBI> struct skyline_sub_vector {
    typedef skyline_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef V * CPT;
    typedef typename select_return<typename linalg_traits<V>::const_iterator,
            typename linalg_traits<V>::iterator, PT>::return_type iterator;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<V>::access_type access_type;

    iterator _begin, _end;
    const void *origin;
    SUBI si;

    size_type size(void) const { return si.size(); }
   
    reference operator[](size_type i) const
    { return access_type()(origin, _begin, _end, si.index(i)); }

    skyline_sub_vector(V &v, const SUBI &s) : _begin(vect_begin(v)),
       _end(vect_end(v)), origin(linalg_origin(v)), si(s) {}
    skyline_sub_vector(const V &v, const SUBI &s)
      : _begin(vect_begin(const_cast<V &>(v))),
	_end(vect_end(const_cast<V &>(v))),
	origin(linalg_origin(const_cast<V &>(v))), si(s) {}
    skyline_sub_vector() {}
    skyline_sub_vector(const skyline_sub_vector<CPT, SUBI> &cr)
      : _begin(cr._begin),_end(cr._end),origin(cr.origin), si(cr.si) {}
  };

  template <class PT, class SUBI> struct skyline_sub_vector_access {
    typedef skyline_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename linalg_traits<this_type>::value_type value_type;
    typedef typename linalg_traits<this_type>::reference reference;
    typedef typename linalg_traits<this_type>::iterator iterator;
    typedef typename linalg_traits<this_type>::const_iterator const_iterator;
    typedef typename linalg_traits<V>::access_type access_type;
    
    reference operator()(const void *o, const iterator &it,
			 const iterator &ite, size_type i)
    { return access_type()(o, it.it, ite.it, it.si.index(i)); }
    
    value_type operator()(const void *o, const const_iterator &it,
			 const const_iterator &ite, size_type i)
    { return access_type()(o, it.it, ite.it, it.si.index(i)); }
  };

  template <class PT, class SUBI> struct skyline_sub_vector_clear {
    typedef skyline_sub_vector<PT, SUBI> this_type;
    typedef typename linalg_traits<this_type>::iterator iterator;
    typedef typename linalg_traits<this_type>::value_type value_type;
    typedef typename linalg_traits<this_type>::access_type access_type;
    
    void operator()(const void */*o*/,const iterator &_begin,const iterator &_end)
      { std::fill(_begin, _end, value_type(0)); }
  };

  template <class PT, class SUBI>
  struct linalg_traits<skyline_sub_vector<PT, SUBI> > {
    typedef skyline_sub_vector<PT, SUBI> this_type;
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename which_reference<PT>::is_reference is_reference;
    typedef abstract_vector linalg_type;
    typedef typename linalg_traits<V>::value_type value_type;
    typedef typename select_return<value_type, typename
            linalg_traits<V>::reference, PT>::return_type reference;
    typedef typename select_return<typename linalg_traits<V>::const_iterator,
	    typename linalg_traits<V>::iterator, PT>::return_type pre_iterator;
    typedef typename select_return<abstract_null_type, 
	    skyline_sub_vector_iterator<pre_iterator, pre_iterator>,
	    PT>::return_type iterator;
    typedef skyline_sub_vector_iterator<typename linalg_traits<V>
            ::const_iterator, pre_iterator> const_iterator;
    typedef abstract_skyline storage_type;
    typedef skyline_sub_vector_access<PT, SUBI> access_type;
    typedef skyline_sub_vector_clear<PT, SUBI> clear_type;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator end(this_type &v) {
      size_type i;
      if (v.si.max + v.si.step() < v._end.index()) {
	i = v.si.max + v.si.step() - v._begin.index();
	size_type j = (v._begin.index() + i - v.si.min) / v.si.step();
	return iterator(v._begin + i, v.si.step(), j);
      }
      else {
	i = (v.si.max + v.si.step() + v.si.step() * v._end.index()
	     - v._end.index()) % v.si.step();
	size_type j = (v._end.index() + i - v.si.min) / v.si.step();
	return iterator(v._end + i, v.si.step(), j);
      }
    }
    static const_iterator end(const this_type &v) {
      size_type i;
      if (v.si.max + v.si.step() < v._end.index()) {
	i = v.si.max + v.si.step() - v._begin.index();
	size_type j = (v._begin.index() + i - v.si.min) / v.si.step();
	return const_iterator(v._begin + i, v.si.step(), j);
      }
      else {
	i = (v.si.max + v.si.step() + v.si.step() * v._end.index()
	     - v._end.index()) % v.si.step();
	size_type j = (v._end.index() + i - v.si.min) / v.si.step();
	return const_iterator(v._end + i, v.si.step(), j);
      }
    }
    static iterator begin(this_type &v) {
      size_type i;
      if (v.si.min >= v._end.index()) return end(v);
      if (v.si.min > v._begin.index()) i = v.si.min - v._begin.index();
      else i = (v.si.min + v.si.step() * v._begin.index()
		- v._begin.index()) % v.si.step();
      if (v._begin.index() + i > v.si.max) return end(v);
      size_type j = (v._begin.index() + i - v.si.min) / v.si.step();
      return iterator(v._begin + i, v.si.step(), j);
    }
    static const_iterator begin(const this_type &v) {
      size_type i;
      if (v.si.min >= v._end.index()) return end(v);
      if (v.si.min > v._begin.index()) i = v.si.min - v._begin.index();
      else i = (v.si.min + v.si.step() * v._begin.index()
		- v._begin.index()) % v.si.step();
      if (v._begin.index() + i > v.si.max) return end(v);
      size_type j = (v._begin.index() + i - v.si.min) / v.si.step();
      return const_iterator(v._begin+i, v.si.step(), j);
    }
    static const void* origin(const this_type &v) { return v.origin; }
    static void do_clear(this_type &v)
      { clear_type()(v.origin, begin(v), end(v)); }
  };

  template <class PT, class SUBI> std::ostream &operator <<
  (std::ostream &o, const skyline_sub_vector<PT, SUBI>& m)
  { gmm::write(o,m); return o; }

#ifdef USING_BROKEN_GCC295
  template <class PT, class SUBI>
  struct linalg_traits<const skyline_sub_vector<PT, SUBI> >
    : public linalg_traits<skyline_sub_vector<PT, SUBI> > {};
#endif

  /* ******************************************************************** */
  /*		sub vector.                                               */
  /* ******************************************************************** */
  /* sub_vector_type<PT, SUBI>::vector_type is the sub vector type        */
  /* returned by sub_vector(v, sub_index)                                 */
  /************************************************************************/

  template <class PT, class SUBI, class st_type> struct svrt_ir {
    typedef abstract_null_type vector_type;
  };

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
  struct svrt_ir<PT, SUBI, abstract_skyline> {
    typedef skyline_sub_vector<PT, SUBI> vector_type;
  };

  template <class PT>
  struct svrt_ir<PT, sub_index, abstract_skyline> {
    typedef sparse_sub_vector<PT, sub_index> vector_type;
  };


  template <class PT, class SUBI>
  struct svrt_ir<PT, SUBI, abstract_sparse> {
    typedef sparse_sub_vector<PT, SUBI> vector_type;
  };

  template <class PT, class SUBI>
  struct sub_vector_type {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename svrt_ir<PT, SUBI,
      typename linalg_traits<V>::storage_type>::vector_type vector_type;
  };

  template <class V, class SUBI>
  typename select_return<
    typename sub_vector_type<const V *, SUBI>::vector_type,
    typename sub_vector_type<V *, SUBI>::vector_type, const V *>::return_type
  sub_vector(const V &v, const SUBI &si) {
    return typename select_return<
      typename sub_vector_type<const V *, SUBI>::vector_type,
      typename sub_vector_type<V *, SUBI>::vector_type, const V *>::return_type
      (linalg_cast(v), si);
  }

  template <class V, class SUBI>
  typename select_return<
    typename sub_vector_type<const V *, SUBI>::vector_type,
    typename sub_vector_type<V *, SUBI>::vector_type, V *>::return_type
  sub_vector(V &v, const SUBI &si) {
    return  typename select_return<
      typename sub_vector_type<const V *, SUBI>::vector_type,
      typename sub_vector_type<V *, SUBI>::vector_type, V *>::return_type
      (linalg_cast(v), si);
  }

}

#endif //  __GMM_SUB_VECTOR_H
