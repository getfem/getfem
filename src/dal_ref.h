/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_ref.h : Structures which refere to containers.           */
/*     									   */
/*                                                                         */
/* Date : August 26, 2000.                                                 */
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


#ifndef __DAL_REF_H
#define __DAL_REF_H

/* *********************************************************************** */
/* WARNING : modifiying the container infirm the validity of references.   */
/* *********************************************************************** */


#include <iterator>
#include <dal_basic.h>

namespace dal
{

  /* ********************************************************************* */
  /* Simple reference.                                                     */
  /* ********************************************************************* */

  template<class ITER> class tab_ref
  {
    protected :

      ITER _begin, _end;

    public :

      typedef typename std::iterator_traits<ITER>::value_type  value_type;
      typedef typename std::iterator_traits<ITER>::pointer     pointer;
      typedef typename std::iterator_traits<ITER>::pointer     const_pointer;
      typedef typename std::iterator_traits<ITER>::reference   reference;
      typedef typename std::iterator_traits<ITER>::reference   const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type
	                                                       difference_type;
      typedef ITER                            iterator;
      typedef ITER                            const_iterator;
      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;
      typedef size_t size_type;
    
      bool empty(void) const { return _begin == _end; }
      size_type size(void) const { return _end - _begin; }

      const iterator &begin(void) { return _begin; }
      const const_iterator &begin(void) const { return _begin; }
      const iterator &end(void) { return _end; }
      const const_iterator &end(void) const { return _end; }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

      reference front(void) { return *begin(); }
      const_reference front(void) const { return *begin(); }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
      void pop_front(void) { ++_begin; }

      const_reference operator [](size_type ii) const { return *(_begin + ii);}
      reference operator [](size_type ii) { return *(_begin + ii); }

      tab_ref(void) {}
      tab_ref(const ITER &b, const ITER &e) : _begin(b), _end(e) {}

  };


  /* ********************************************************************* */
  /* Reference with index.                                                 */
  /* ********************************************************************* */

  template<class ITER> struct _tab_ref_index_iterator
    : public dynamic_array<size_t>::const_iterator
  {
    typedef typename std::iterator_traits<ITER>::value_type  value_type;
    typedef typename std::iterator_traits<ITER>::pointer     pointer;
    typedef typename std::iterator_traits<ITER>::reference   reference;
    typedef typename std::iterator_traits<ITER>::difference_type  
    difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t size_type;
    typedef dynamic_array<size_type>::const_iterator _dnas_iterator;
    typedef _tab_ref_index_iterator<ITER> iterator;
    

    ITER piter;
    
    iterator operator ++(int)
    { iterator tmp = *this; ++(*((_dnas_iterator *)(this))); return tmp; }
    iterator operator --(int)
    { iterator tmp = *this; --(*((_dnas_iterator *)(this))); return tmp; }
    iterator &operator ++()
    { ++(*((_dnas_iterator *)(this))); return *this; }
    iterator &operator --()
    { --(*((_dnas_iterator *)(this))); return *this; }
    iterator &operator +=(difference_type i)
    { (*((_dnas_iterator *)(this))) += i; return *this; }
    iterator &operator -=(difference_type i)
    { (*((_dnas_iterator *)(this))) -= i; return *this; }
    iterator operator +(difference_type i) const
    { iterator it = *this; return (it += i); }
    iterator operator -(difference_type i) const
    { iterator it = *this; return (it -= i); }
    difference_type operator -(const iterator &i) const
    { return *((_dnas_iterator *)(this)) - *((_dnas_iterator *)(&i)); }
	
    reference operator *() const
    { return *(piter + *((*((_dnas_iterator *)(this))))); }
    reference operator [](int ii)
    { return *(piter + *((*((_dnas_iterator *)(this+ii))))); }
    
    bool operator ==(const iterator &i) const
    { 
      return ((piter) == ((i.piter))
       && *((_dnas_iterator *)(this)) == *((*((_dnas_iterator *)(this)))));
    }
    bool operator !=(const iterator &i) const
    { return !(i == *this); }
    bool operator < (const iterator &i) const
    { 
      return ((piter) == ((i.piter))
	 && *((_dnas_iterator *)(this)) < *((*((_dnas_iterator *)(this)))));
    }

    _tab_ref_index_iterator(void) {}
    _tab_ref_index_iterator(const ITER &iter, const _dnas_iterator &dnas_iter)
      : _dnas_iterator(dnas_iter), piter(iter) {}
  };


  template<class ITER> class tab_ref_index
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type
	                                                       difference_type;
      typedef size_t size_type; 
      typedef _tab_ref_index_iterator<ITER> iterator;
      typedef iterator                          const_iterator;
      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;
    
    protected :

      ITER _begin;
      dynamic_array<size_type> _index;

    public :

      bool empty(void) const { return _index.empty(); }
      size_type size(void) const { return _index.size(); }


      iterator begin(void) { return iterator(_begin, _index.begin()); }
      const_iterator begin(void) const
      { return iterator(_begin, _index.begin()); }
      iterator end(void) { return iterator(_begin, _index.end()); }
      const_iterator end(void) const { return iterator(_begin, _index.end()); }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }


      reference front(void) { return *(_begin +_index[0]); }
      const_reference front(void) const { return *(_begin +_index[0]); }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
   
      tab_ref_index(void) {}
      tab_ref_index(const ITER &b, const dynamic_array<size_type> &ind)
      { _begin = b; _index = ind; }


      const_reference operator [](size_type ii) const
      { return *(_begin + _index[ii]);}
      reference operator [](size_type ii) { return *(_begin + _index[ii]); }

  };


  /* ********************************************************************* */
  /* Reference with reference on index.                                    */
  /* ********************************************************************* */

  template<class ITER, class ITER_INDEX>
    struct _tab_ref_index_ref_iterator
    {
      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::difference_type
                                                              difference_type;
      typedef std::random_access_iterator_tag iterator_category;
      typedef _tab_ref_index_ref_iterator<ITER, ITER_INDEX> iterator;
      typedef size_t size_type;

      ITER piter;
      ITER_INDEX iter_index;
      
      iterator operator ++(int)
      { iterator tmp = *this; ++iter_index; return tmp; }
      iterator operator --(int)
      { iterator tmp = *this; --iter_index; return tmp; }
      iterator &operator ++() { ++iter_index; return *this; }
      iterator &operator --() { --iter_index; return *this; }
      iterator &operator +=(difference_type i)
      { iter_index += i; return *this; }
      iterator &operator -=(difference_type i)
      { iter_index -= i; return *this; }
      iterator operator +(difference_type i) const
      { iterator it = *this; return (it += i); }
      iterator operator -(difference_type i) const
      { iterator it = *this; return (it -= i); }
      difference_type operator -(const iterator &i) const
      { return iter_index - i.iter_index; }
	
      reference operator *() const
      { return *(piter + *iter_index); }
      reference operator [](int ii) const
      { return *(piter + *(iter_index+ii)); }
      
      bool operator ==(const iterator &i) const
      { return ((piter) == ((i.piter)) && iter_index == i.iter_index); }
      bool operator !=(const iterator &i) const { return !(i == *this); }
      bool operator < (const iterator &i) const
      { return ((piter) == ((i.piter)) && iter_index < i.iter_index); }

      _tab_ref_index_ref_iterator(void) {}
      _tab_ref_index_ref_iterator(const ITER &iter, 
				  const ITER_INDEX &dnas_iter)
	: piter(iter), iter_index(dnas_iter) {}
      
    };


  template<class ITER, class ITER_INDEX> class tab_ref_index_ref
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type  
	                                                       difference_type;
      typedef size_t size_type;
      typedef _tab_ref_index_ref_iterator<ITER, ITER_INDEX> iterator;
      typedef iterator                          const_iterator;
      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;
    
    protected :

      ITER _begin;
      ITER_INDEX _index_begin, _index_end;

    public :

      bool empty(void) const { return _index_begin == _index_end; }
      size_type size(void) const { return _index_end - _index_begin; }

      iterator begin(void) { return iterator(_begin, _index_begin); }
      const_iterator begin(void) const
      { return iterator(_begin, _index_begin); }
      iterator end(void) { return iterator(_begin, _index_end); }
      const_iterator end(void) const { return iterator(_begin, _index_end); }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

      reference front(void) { return *(_begin + *_index_begin); }
      const_reference front(void) const { return *(_begin + *_index_begin); }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
      void pop_front(void) { ++_index_begin; }

      tab_ref_index_ref(void) {}
      tab_ref_index_ref(const ITER &b, const ITER_INDEX &bi,
			               const ITER_INDEX &ei)
	: _begin(b), _index_begin(bi), _index_end(ei) {}

      const_reference operator [](size_type ii) const
      { return *(_begin + _index_begin[ii]);}
      reference operator [](size_type ii)
      { return *(_begin + _index_begin[ii]); }

  };


  /* ********************************************************************* */
  /* Reference on regularly spaced elements.                               */
  /* ********************************************************************* */

  template<class ITER> struct _tab_ref_reg_spaced_iterator
  {
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::difference_type
                                                            difference_type;
    typedef typename std::iterator_traits<ITER>::iterator_category
                                                            iterator_category;
    typedef size_t size_type;
    typedef _tab_ref_reg_spaced_iterator<ITER> iterator;

    ITER it;
    size_type N;
    
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

    reference operator *() const { return *it; }
    reference operator [](int ii) { return *(it + ii * N); }

    bool operator ==(const iterator &i) const { return (it == i.it); }
    bool operator !=(const iterator &i) const { return !(i == *this); }
    bool operator < (const iterator &i) const { return (it < i.it); }

    _tab_ref_reg_spaced_iterator(void) {}
    _tab_ref_reg_spaced_iterator(const ITER &iter, size_type n)
      : it(iter), N(n) { }
    
  };

  template<class ITER> class tab_ref_reg_spaced
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef typename std::iterator_traits<ITER>::difference_type
	                                                       difference_type;
      typedef size_t size_type;
      typedef _tab_ref_reg_spaced_iterator<ITER> iterator;
      typedef iterator                          const_iterator;
      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;
    
    protected :

      ITER _begin, _end;
      size_type N;

    public :

      bool empty(void) const { return _begin == _end; }
      size_type size(void) const { return (_end - _begin) / N; }

      iterator begin(void) { return iterator(_begin, N); }
      const_iterator begin(void) const { return iterator(_begin, N); }
      iterator end(void) { return iterator(_end, N); }
      const_iterator end(void) const { return iterator(_end, N); }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }

      reference front(void) { return *_begin; }
      const_reference front(void) const { return *_begin; }
      reference back(void) { return *(--(end())); }
      const_reference back(void) const { return *(--(end())); }
      void pop_front(void) { _begin += N; }

      tab_ref_reg_spaced(void) {}
      tab_ref_reg_spaced(const ITER &b, const ITER &e, size_type n)
	: _begin(b), _end(e), N(n) {}


      const_reference operator [](size_type ii) const
      { return *(_begin + ii * N);}
      reference operator [](size_type ii) { return *(_begin + ii * N); }

  };

  /* ********************************************************************* */
  /* Reference elements selected with a condition.                         */
  /* ********************************************************************* */

  template<class ITER, class COND> 
    struct _tab_ref_with_selection_iterator : public ITER
  {
    typedef typename std::iterator_traits<ITER>::value_type value_type;
    typedef typename std::iterator_traits<ITER>::pointer    pointer;
    typedef typename std::iterator_traits<ITER>::reference  reference;
    typedef typename std::iterator_traits<ITER>::difference_type
                                                              difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef _tab_ref_with_selection_iterator<ITER, COND> iterator;
    const COND *cond;
    
    void forward(void) { while (!(*cond)(*this)) (*((ITER *)(this)))++; }
    iterator &operator ++()
    { (*((ITER *)(this)))++; forward(); return *this; }
    iterator operator ++(int)
    { iterator tmp = *this; ++(*this); return tmp; }
    void begin(const ITER &iter, const COND *c)
    { *((ITER *)(this)) = iter; cond = c; forward(); }
    
    _tab_ref_with_selection_iterator(void) {}
    _tab_ref_with_selection_iterator(const ITER &iter, const COND *c)
      : ITER(iter) { cond = c; }
    
  };

  template<class ITER, class COND> class tab_ref_with_selection
  {
    public :

      typedef typename std::iterator_traits<ITER>::value_type value_type;
      typedef typename std::iterator_traits<ITER>::pointer    pointer;
      typedef typename std::iterator_traits<ITER>::pointer    const_pointer;
      typedef typename std::iterator_traits<ITER>::reference  reference;
      typedef typename std::iterator_traits<ITER>::reference  const_reference;
      typedef size_t  size_type;
      typedef _tab_ref_with_selection_iterator<ITER, COND> iterator;
      typedef iterator   const_iterator;
    
    protected :

      ITER _begin, _end;
      COND cond;

    public :

      iterator begin(void)
      { iterator it; it.begin(_begin, &cond); return it; }
      const_iterator begin(void) const { return iterator(_begin, &cond); }
      iterator end(void) { return iterator(_end, &cond); }
      const_iterator end(void) const { return iterator(_end, &cond); }
      bool empty(void) const { return _begin == _end; }

      reference front(void) { return *begin(); }
      const_reference front(void) const { return *begin(); }
      void pop_front(void) { ++_begin; _begin = begin(); }

      COND &condition(void) { return cond; }
      const COND &condition(void) const { return cond; }
   
      tab_ref_with_selection(void) {}
      tab_ref_with_selection(const ITER &b, const ITER &e, const COND &c)
	: _begin(b), _end(e), cond(c) { _begin = begin(); }

  };

}

#endif /* __DAL_REF_H  */
