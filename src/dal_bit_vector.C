/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_bit_vector.C : bit vector based on a dynamic array.      */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1995-2002  Yves Renard.                                   */
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
      const_iterator itb = begin(), ite = end();
      icard = 0;
      while (itb != ite) { if (*itb) ++icard; ++itb; }
      icard_valid = true;
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
  
  bit_vector &bit_vector::operator |=(const bit_vector &bv) {
    bit_container::iterator it1b = bit_container::begin();
    bit_container::iterator it1e = bit_container::end();
    bit_container::const_iterator it2b = bv.bit_container::begin();
    bit_container::const_iterator it2x = it2b;
    bit_container::const_iterator it2e = bv.bit_container::end();
    
    while (it1b != it1e && it2x != it2e) { *it1b++ |= *it2x++; }
    while (it2x != it2e) (*(bit_container *)(this))[it2x - it2b] = *it2x++;
    icard_valid = false;
    ilast_false = std::min(ilast_false, bv.ilast_false);
    ifirst_false = std::max(ifirst_false, bv.ifirst_false);
    ilast_true = std::max(ilast_true, bv.ilast_true);
    ifirst_true = std::min(ifirst_true, bv.ifirst_true);
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

  void bit_vector::add(size_type i, size_type j) { 
    if (j < i) std::swap(i,j); add(j);
    std::fill(this->begin()+i, this->begin()+j+1, true);
  }

  void bit_vector::sup(size_type i, size_type j) {
    if (j < i) std::swap(i,j); sup(j);
    std::fill(this->begin()+i, this->begin()+j+1, false);
  }

  std::ostream &operator <<(std::ostream &o, const bit_vector &s) {
    bit_vector u = s; int i;
    o << "int set ";
    for (i << u; i != int(-1); i << u) o << " : " << i;
    o << " ";
    return o;
  }
}

