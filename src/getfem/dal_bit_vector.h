/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 1995-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/


#ifndef DAL_BIT_VECTOR_H__
#define DAL_BIT_VECTOR_H__


/** @file dal_bit_vector.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 01, 1995.
    @brief Provide a dynamic bit container.
    
    Provide a dynamic bit container, which can also be considered as a
    set of integers.
    
    As a convention, the default value of a bit is false.  The main
    member functions are dal::bit_vector::is_in,
    dal::bit_vector::add, dal::bit_vector::sup. Iterate over the
    bit_vector with dal::bv_visitor
*/

#include "dal_basic.h"
#include <limits.h>
#include <bitset>

namespace dal {

  typedef unsigned int bit_support;
  static const bit_support WD_BIT = bit_support(CHAR_BIT*sizeof(bit_support));
  static const bit_support WD_MASK = WD_BIT - 1;
  typedef dynamic_array<bit_support, 4> bit_container;

  class bit_vector;

  struct APIDECL bit_reference {
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

  struct APIDECL bit_iterator {
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

  struct APIDECL bit_const_iterator {
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
      : ind(x.ind),  mask(x.mask), p(x.p), bv(x.bv) {}
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

  ///Dynamic bit container. 
  class APIDECL bit_vector : public bit_container {
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
    
    void change_for_true(size_type i) {
      ifirst_true = std::min(ifirst_true, i);
      ilast_true = std::max(ilast_true, i);
      ++icard;
    }
    void change_for_false(size_type i) {
      ifirst_false = std::min(ifirst_false, i);
      ilast_false = std::max(ilast_false, i);
      --icard;
    }
    
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
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
    { if (ii >= size()) fill_false(size(),ii); return *iterator(*this, ii);}
    
    void swap(bit_vector &da);
    
    void clear(void) {
      icard = 0; icard_valid = true;
      ifirst_false = ilast_false = ifirst_true = ilast_true = 0;
      fill_false(0,0); 
    }
    void swap(size_type i1, size_type i2) {
      if (i1 != i2) {
	reference r1 = (*this)[i1], r2 = (*this)[i2];
	bool tmp = r1; r1 = r2; r2 = tmp;
      }
    }
    size_type memsize(void) const {
      return bit_container::memsize() + sizeof(bit_vector) 
	- sizeof(bit_container);
    }
    size_type card(void) const;
    /// index of first non-zero entry (size_type(-1) for an empty bit_vector)
    size_type first_true(void) const;
    /// index of first zero entry (size_type(0) for an empty bit_vector)
    size_type first_false(void) const;
      /// index of last non-zero entry (size_type(-1) for an empty bit_vector)
    size_type last_true(void) const;
    /// index of last zero entry (size_type(0) for an empty bit_vector)
    size_type last_false(void) const;
    /// remove all elements found in bv
    bit_vector &setminus(const bit_vector &bv);
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
    template <size_t N> bit_vector(const std::bitset<N> &bs) {
      clear();
      for (size_type i=0; i < bs.size(); ++i) { if (bs[i]) add(i); }
    }
    
    /// merges the integer values of the supplied container into the bit_vector
    template <typename ICONT> dal::bit_vector& merge_from(const ICONT& c) {
      for (typename ICONT::const_iterator it = c.begin(); it != c.end(); ++it)
	add(*it);
      return *this;
    }
    /** merges the integer values of the supplied iterator range into
     * the bit_vector */
    template <typename IT> dal::bit_vector& merge_from(IT b, IT e) {
      while (b != e) { add(*b++); }
      return *this;
    }
    /** return true if the supplied bit_vector is a subset of the current
     * bit_vector */
    bool contains(const dal::bit_vector &other) const;
    
  public : 
    /// return true if (*this)[i] == true
    bool is_in(size_type i) const { 
      if (i < ifirst_true || i > ilast_true) return false;
      else return (((*(const bit_container*)(this))[i / WD_BIT]) & 
		   (bit_support(1) << (i & WD_MASK))) ? true : false; }
    void add(size_type i) { (*this)[i] = true; }
    /** set the interval [i .. i+nb-1] to true */
    void add(size_type i, size_type nb);
    void sup(size_type i) { (*this)[i] = false; } /* deprecated ...*/
    void del(size_type i) { (*this)[i] = false; }
    /** set the interval [i .. i+nb-1] to false */
    void sup(size_type i, size_type nb); /* deprecated ...*/
    void del(size_type i, size_type nb);
    int first(void) const { return (card() == 0) ? -1 : int(first_true()); }
    int last(void) const { return (card() == 0) ? -1 : int(last_true()); }
    inline int take_first(void)
    { int res = first(); if (res >= 0) sup(res); return res; }
    inline int take_last(void)
    { int res = last(); if (res >= 0) sup(res); return res; }
  };

  /**
  Iterator class for bv_visitor and bv_visitor_c.
  This iterator class enables the use of c++11 range-based loop
  feature.
  example:
  @code
  for (auto i : bv_visitor(v)) {
  .... (use i as an unsigned int)
  }
  */
  template<typename VISITOR>
  class const_visitor_iterator
  {
    typedef dal::bit_vector::size_type size_type;

  public:
    const_visitor_iterator(const VISITOR* p_visitor, size_type pos)
      : p_visitor_(const_cast<VISITOR*>(p_visitor)), pos_(pos)
    {}

    bool operator!= (const const_visitor_iterator &other) const{
      return pos_ < other.pos_;
    }

    size_type operator *() const{
      return dal::bit_vector::size_type(*p_visitor_);
    }

    const const_visitor_iterator &operator++(){
      ++*p_visitor_;
      pos_ = *p_visitor_;
      return *this;
    }

  private:
    VISITOR *p_visitor_;
    size_type pos_;
  };
  
  /**
     if you are only interested in indexes of true values of a bit_vector
     (i.e. if you use it as an int set), use bv_visitor instead of
     bit_vector::const_iterator (much faster)

     example:
     @code
     for (bv_visitor i(v); !i.finished(); ++i) {
       .... (use i as an unsigned int)
     }
     @endcode
     CAUTION: use bv_visitor_c instead of bv_visitor if the class bv_visitor
     need to store a copy of the bit_vector 
     (if the original is destroyed just after the creation...)
  */
  class APIDECL bv_visitor {
    typedef dal::bit_vector::size_type size_type;
    bit_container::const_iterator it;
    size_type ilast,ind;
    bit_support v;
  public:
    bv_visitor(const dal::bit_vector& b) : 
      it(((const bit_container&)b).begin()+b.first()/WD_BIT),
      ilast(b.last()+1), ind(b.first()), v(0) {
      if (ind < ilast) { v = *it; v >>= (ind&WD_MASK); }
    }
    bool finished() const { return ind >= ilast; }
    bool operator++();
    operator size_type() const { return ind; }

    size_type get_last_index() const { return ilast;}

    const_visitor_iterator<bv_visitor> begin() const{
      return const_visitor_iterator<bv_visitor>(this, *this);
    }

    const_visitor_iterator<bv_visitor> end() const{
      return const_visitor_iterator<bv_visitor>(this, ilast);
    }
  };

  /**
    bv_visitor with local copy of the bit_vector
  */
  class APIDECL bv_visitor_c {
    bit_vector bv;
    bv_visitor v; // no inheritance since v must be init after bv
  public:
    bv_visitor_c(const dal::bit_vector& b) : bv(b), v(bv) {}
    bool finished() const { return v.finished(); }
    bool operator++() { return ++v; }
    operator dal::bit_vector::size_type() const
    { return dal::bit_vector::size_type(v); }

    const_visitor_iterator<bv_visitor_c> begin() const{
      return const_visitor_iterator<bv_visitor_c>(this, *this);
    }
    const_visitor_iterator<bv_visitor_c> end() const{
      return const_visitor_iterator<bv_visitor_c>(this, v.get_last_index());
    }
  };

  /// extract index of first entry in the bit_vector
  inline int APIDECL &operator << (int &i, bit_vector &s)
  { i = s.take_first(); return i; }
  inline const int APIDECL &operator >> (const int &i, bit_vector &s)
  { s.add(i); return i; }

  inline size_t APIDECL &operator << (size_t &i, bit_vector &s)
  { i = s.take_first(); return i; }
  inline const size_t &operator >> (const size_t &i, bit_vector &s)
  { s.add(i); return i; }

  std::ostream APIDECL &operator <<(std::ostream &o, const bit_vector &s);

}

#endif /* DAL_BIT_VECTOR_H__ */
