// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_bit_vector.cc : bit vector based on a dynamic array.
//           
// Date    : June 01, 1995.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1995-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#include <dal_bit_vector.h>

namespace dal
{

  bit_reference& bit_reference::operator = (bool x) {
      if (x) { if (!(*p & mask)) { *p |= mask; bv->change_for_true(ind); } } 
      else  { if (*p & mask) { *p &= ~mask; bv->change_for_false(ind); } }
      return *this;
  }

  bit_iterator::bit_iterator(bit_vector &b, size_type i) : p(b, i / WD_BIT)
  { ind = i; bv = &b, mask = bit_support(1) << (i & WD_MASK); }

  bit_iterator& bit_iterator::operator+=(difference_type i){
      ind+=i; mask = bit_support(1) << (ind & WD_MASK); 
      p=bit_container::iterator(*bv, ind/WD_BIT);
      return *this;
  }

  bit_const_iterator::bit_const_iterator(const bit_vector &b, size_type i)
    : p(b, i / WD_BIT)
  { ind = i; bv = &b, mask = bit_support(1) << (i & WD_MASK); }

  bit_const_iterator& bit_const_iterator::operator+=(difference_type i) {
    ind += i; mask = bit_support(1) << (ind & WD_MASK); 
    p = bit_container::const_iterator(*bv, ind/WD_BIT);
    return *this;
  }


  bit_vector::size_type bit_vector::card(void) const {
    if (!icard_valid) {
      bit_container::const_iterator itb = bit_container::begin() + ifirst_true/WD_BIT, ite = bit_container::end();
      bit_support x;
      icard = 0;
      for (; itb != ite; ++itb) { 
        if ((x = *itb)) /* fast count of the nb of bits of a integer (faster than shifts) */  
          do icard++; while ((x &= x-1));
      }
      /*      const_iterator itb = begin(), ite = end();
      icard = 0;
      while (itb != ite) { if (*itb) ++icard; ++itb; }
      icard_valid = true;*/
    }
    return icard;
  }

  bit_vector::size_type bit_vector::first_true(void) const {
    const_iterator itx = begin(), ite = end(); itx += ifirst_true;
    while (itx != ite && !*itx ) { ++itx; ++(ifirst_true); }
    return ifirst_true;
  }
  
  bit_vector::size_type bit_vector::first_false(void) const {
    const_iterator itx = begin(), ite = end(); itx += ifirst_false;
    while (itx != ite && *itx) { ++itx; ++(ifirst_false); }
    return ifirst_false;
  }

  bit_vector::size_type bit_vector::last_true(void) const {
    const_iterator itb = begin(), itx = itb; itx += ilast_true;
    while (itx != itb && !*itx) { --itx; --(ilast_true); }
    return ilast_true;
  }
  
  bit_vector::size_type bit_vector::last_false(void) const {
    const_iterator itb = begin(), itx = itb; itx += ilast_false;
    while (itx != itb && *itx) { --itx; --(ilast_false); }
    return ilast_false;
  }
  
  bit_vector &bit_vector::setminus(const bit_vector& b) {
    for (bv_visitor i(b); !i.finished(); ++i) sup(i); 
    return *this;
  }

  bit_vector &bit_vector::operator |=(const bit_vector &bv) {
    for (bv_visitor i(bv); !i.finished(); ++i) add(i);
    return *this;
  }
  
  bit_vector &bit_vector::operator &=(const bit_vector &bv) {
    bit_container::iterator it1b = bit_container::begin();
    bit_container::iterator it1e = bit_container::end();
    bit_container::const_iterator it2b = bv.bit_container::begin();
    bit_container::const_iterator it2e = bv.bit_container::end();
    
    while (it1b != it1e && it2b != it2e) { *it1b++ &= *it2b++; }
    while (it1b != it1e) { *it1b++ = 0; }
    icard_valid = false;
    ilast_true = std::min(ilast_true, bv.ilast_true);
    ifirst_true = std::max(ifirst_true, bv.ifirst_true);
    ilast_false = std::min(size()-1, std::max(ilast_false,bv.ilast_false));
    ifirst_false = std::min(ifirst_false, bv.ifirst_false);
    return *this;
  }

  bool bit_vector::operator ==(const bit_vector &bv) const {
    bit_container::const_iterator it1b = bit_container::begin();
    bit_container::const_iterator it1e = bit_container::end();
    bit_container::const_iterator it2b = bv.bit_container::begin();
    bit_container::const_iterator it2e = bv.bit_container::end();
    
    while (it1b!=it1e && it2b!=it2e) if (*it1b++ != *it2b++) return false;
    while (it1b != it1e) if (*it1b++ != 0) return false;
    while (it2b != it2e) if (*it2b++ != 0) return false;
    return true;
  }

  void bit_vector::add(size_type i, size_type nb) { 
    if (nb) {
      add(i+nb-1);
      std::fill(this->begin()+i, this->begin()+(i+nb), true);
    }
  }

  void bit_vector::sup(size_type i, size_type nb) {
    if (nb) {
      sup(i+nb-1);
      std::fill(this->begin()+i, this->begin()+(i+nb), false);
    }
  }

  bool bit_vector::contains(const dal::bit_vector& other) const {
    for (dal::bv_visitor i(other); !i.finished(); ++i) {
      if (!this->is_in(i)) return false;
    }
    return true;
  }

  std::ostream &operator <<(std::ostream &o, const bit_vector &s) {
    bool first = true;
    o << "[";
    for (bv_visitor i(s); !i.finished(); ++i) {
      if (!first) o << " ";
      o << bit_vector::size_type(i);
      first = false;
    }
    o << "]";
    return o;
  }
}

