/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_bit_vector.h : bit vector based on a dynamic array.      */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
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


#ifndef __DAL_BIT_VECTOR_H
#define __DAL_BIT_VECTOR_H

/* *********************************************************************** */
/* Remarks                                                                 */
/* 									   */
/* - As a convention, the default value of a bit is false.                 */
/* 									   */
/* - Some operations can be optimized (card, first_true ... for instance)  */
/*                                                                         */
/* *********************************************************************** */

#include <dal_basic.h>
#include <limits.h>


namespace dal
{
  typedef unsigned int bit_support;
  static const bit_support WD_BIT = bit_support(CHAR_BIT*sizeof(bit_support));
  static const bit_support WD_MASK = WD_BIT - 1;
  typedef dynamic_array<bit_support, 4> bit_container;

  class bit_vector;

  struct bit_reference
  {
    typedef size_t  size_type;

    bit_support* p;
    bit_support mask;
    size_type ind;
    bit_vector* bv;
    
    bit_reference(bit_support* x, bit_support m, size_type y, bit_vector* z)
    { p = x; ind = y; bv = z; mask = m; }
    bit_reference(void) {}
    operator bool(void) const { return (*p & mask) != 0; }
    bit_reference& operator = (bool x);
    bit_reference& operator=(const bit_reference& x)
    { return *this = bool(x); }
    bool operator==(const bit_reference& x) const
    { return bool(*this) == bool(x); }
    bool operator<(const bit_reference& x) const
      { return bool(*this) < bool(x); }
    void flip(void) { if (bool(*this)) *this = false; else *this = true; }
  };

  struct bit_iterator
  {
    typedef bool             value_type;
    typedef bit_reference    reference;
    typedef bit_reference*   pointer;
    typedef size_t           size_type;
    typedef ptrdiff_t        difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    size_type ind;
    bit_support mask;
    bit_container::iterator p;
    bit_vector* bv;
    
    inline void bump_up()
    { ind++; if (!(mask <<= 1)) { ++p; mask = 1;} }
    inline void bump_down()
    { ind--; if (!(mask >>= 1)) { --p; mask = 1; mask <<= WD_MASK; } }
    
    bit_iterator(void) {}
    bit_iterator(bit_vector &b, size_type i);
    reference operator*() const
    { return reference(&(*p), mask, ind, bv); }
    bit_iterator& operator++() { bump_up(); return *this; }
    bit_iterator operator++(int)
    { bit_iterator tmp=*this; bump_up(); return tmp; }
    bit_iterator& operator--() { bump_down(); return *this; }
    bit_iterator operator--(int)
    { bit_iterator tmp = *this; bump_down(); return tmp; }
    bit_iterator& operator+=(difference_type i);
    bit_iterator& operator-=(difference_type i)
    { *this += -i; return *this; }
    bit_iterator operator+(difference_type i) const
    { bit_iterator tmp = *this; return tmp += i; }
    bit_iterator operator-(difference_type i) const
    { bit_iterator tmp = *this; return tmp -= i; }
    difference_type operator-(bit_iterator x) const { return ind - x.ind; }
    reference operator[](difference_type i) { return *(*this + i); }
    size_type index(void) const { return ind; }
    bool operator==(const bit_iterator& x) const { return ind == x.ind; }
    bool operator!=(const bit_iterator& x) const { return ind != x.ind; }
    bool operator<(bit_iterator x) const { return ind < x.ind; }
  };

  struct bit_const_iterator
  {
    typedef bool             value_type;
    typedef bool             reference;
    typedef const bool*      pointer;
    typedef size_t           size_type;
    typedef ptrdiff_t        difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    
    size_type ind;
    bit_support mask;
    bit_container::const_iterator p;
    const bit_vector* bv;
    
    inline void bump_up()
    { ind++; if (!(mask <<= 1)) { ++p; mask = 1;} }
    inline void bump_down()
    { ind--; if (!(mask >>= 1)) { --p; mask = 1; mask <<= WD_MASK; } }
    
    bit_const_iterator() {}
    bit_const_iterator(const bit_vector &b, size_type i);
    bit_const_iterator(const bit_iterator& x)
      : p(x.p), ind(x.ind), bv(x.bv), mask(x.mask) {}
    reference operator*() const { return (*p & mask) != 0; }
    bit_const_iterator& operator++() { bump_up();  return *this; }
    bit_const_iterator operator++(int)
    { bit_const_iterator tmp = *this; bump_up(); return tmp; }
    bit_const_iterator& operator--() { bump_down(); return *this; }
    bit_const_iterator operator--(int)
    { bit_const_iterator tmp = *this; bump_down(); return tmp; }
    bit_const_iterator& operator+=(difference_type i);
    bit_const_iterator& operator-=(difference_type i)
    { *this += -i; return *this; }
    bit_const_iterator operator+(difference_type i) const
    { bit_const_iterator tmp = *this; return tmp += i; }
    bit_const_iterator operator-(difference_type i) const
    { bit_const_iterator tmp = *this; return tmp -= i; }
    difference_type operator-(bit_const_iterator x) const { return ind-x.ind; }
    reference operator[](difference_type i) { return *(*this + i); }
    size_type index(void) const { return ind; }
    bool operator==(const bit_const_iterator& x) const { return ind == x.ind; }
    bool operator!=(const bit_const_iterator& x) const { return ind != x.ind; }
    bool operator<(bit_const_iterator x) const { return ind < x.ind; }
  };

  class bit_vector : public bit_container
  {
    public :
      
      typedef bool         value_type;
      typedef size_t       size_type;
      typedef ptrdiff_t    difference_type;
      typedef bool         const_reference;
      typedef const bool*  const_pointer;
      typedef bit_reference reference;
      typedef bit_reference*   pointer;
      typedef bit_iterator iterator;
      typedef bit_const_iterator const_iterator;

    protected :

      mutable size_type ifirst_true, ilast_true;
      mutable size_type ifirst_false, ilast_false;
      mutable size_type icard;
      mutable bool icard_valid;

      void fill_false(size_type i1, size_type i2);
 
   public : 
      
      void change_for_true(size_type i)
      {
	ifirst_true = std::min(ifirst_true, i);
	ilast_true = std::max(ilast_true, i);
	++icard;
      }
      void change_for_false(size_type i)
      {
	ifirst_false = std::min(ifirst_false, i);
	ilast_false = std::max(ilast_false, i);
	--icard;
      }

      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;
      size_type size(void) const { return std::max(ilast_true, ilast_false)+1;}
      
      iterator begin(void) { return iterator(*this, 0); }
      const_iterator begin(void) const { return const_iterator(*this, 0); }
      iterator end(void) { return iterator(*this, size()); }
      const_iterator end(void) const { return const_iterator(*this, size()); }
      reverse_iterator rbegin(void) { return reverse_iterator(end()); }
      const_reverse_iterator rbegin(void) const
      { return const_reverse_iterator(end()); }
      reverse_iterator rend(void) { return reverse_iterator(begin()); }
      const_reverse_iterator rend(void) const
      { return const_reverse_iterator(begin()); }
      
      size_type capacity(void) const
      { return bit_container::capacity() * WD_BIT; }
      size_type max_size(void) const { return (size_type(-1)); }
      // bool empty(void) const { return card() == 0; } /* ?? */
      reference front(void) { return *begin(); }
      const_reference front(void) const { return *begin(); }
      reference back(void) { return *(end() - 1); }
      const_reference back(void) const { return *(end() - 1); }
      

      const_reference operator [](size_type ii) const
      { return (ii >= size()) ? false : *const_iterator(*this, ii); }
      reference operator [](size_type ii)
      { if (ii >= size()) fill_false(size(),ii); return *iterator(*this, ii); }

      void swap(bit_vector &da);
      void clear(void);
      void swap(size_type i1, size_type i2);
      size_type memsize(void) const
      {	return bit_container::memsize() + sizeof(bit_vector); }
      size_type card(void) const;
      size_type first_true(void) const;
      size_type first_false(void) const;
      size_type last_true(void) const;
      size_type last_false(void) const;
      bit_vector &operator |=(const bit_vector &bv);
      bit_vector &operator &=(const bit_vector &bv);

      bit_vector operator |(const bit_vector &bv) const
      { bit_vector r(*this); r |= bv; return r; }
      bit_vector operator &(const bit_vector &bv) const
      { bit_vector r(*this); r &= bv; return r; }
      bool operator ==(const bit_vector &bv) const;
      bool operator !=(const bit_vector &bv) const
      { return !((*this) == bv); }
     
      bit_vector(void) { clear(); }

  /* ********************************************************************* */
  /*									   */
  /*	     Adaptation for old structure int_set.                         */
  /*									   */
  /* ********************************************************************* */
  
    public : 

      bool is_in(size_type i) const { return (*this)[i]; }
      void add(size_type i) { (*this)[i] = true; }
      void add(size_type i, size_type j);
      void sup(size_type i) { (*this)[i] = false; }
      void sup(size_type i, size_type j);
      int first(void) const { return (card() == 0) ? -1 : int(first_true()); }
      int last(void) const { return (card() == 0) ? -1 : int(last_true()); }
      inline int take_first(void)
      { int res = first(); if (res >= 0) sup(res); return res; }
      inline int take_last(void)
      { int res = last(); if (res >= 0) sup(res); return res; }
  };

  inline int &operator << (int &i, bit_vector &s)
  { i = s.take_first(); return i; }
  inline const int &operator >> (const int &i, bit_vector &s)
  { s.add(i); return i; }

  inline size_t &operator << (size_t &i, bit_vector &s)
  { i = s.take_first(); return i; }
  inline const size_t &operator >> (const size_t &i, bit_vector &s)
  { s.add(i); return i; }

  STD_NEEDED ostream &operator <<(STD_NEEDED ostream &o, const bit_vector &s);

}

#endif /* __DAL_BIT_VECTOR_H */
