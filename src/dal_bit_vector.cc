/*===========================================================================
 
 Copyright (C) 1995-2015 Yves Renard
 
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
 
===========================================================================*/


#include "getfem/dal_bit_vector.h"

namespace dal {

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

  void bit_vector::fill_false(size_type i1, size_type i2) {
    size_type f = i1 / WD_BIT, r = i1 & (WD_BIT-1), l = i2 / WD_BIT;
    (*((bit_container *)(this)))[l];
    if (r != 0) f++; l++;
    if (f < l) std::fill(dal::bit_container::begin()+f,
			 dal::bit_container::begin()+l, 0);
    ilast_false = i2;
  }

  void bit_vector::swap(bit_vector &da) {
    ((bit_container *)(this))->swap(da);
    std::swap(ifirst_true, da.ifirst_true);
    std::swap(ifirst_false, da.ifirst_false);
    std::swap(ilast_true, da.ilast_true);
    std::swap(ilast_false, da.ilast_false);
    std::swap(icard, da.icard);
    std::swap(icard_valid, da.icard_valid);
  }

  bit_vector::size_type bit_vector::card(void) const {
    if (!icard_valid) {
//       bit_container::const_iterator itb = bit_container::begin()
// 	+ ifirst_true/WD_BIT, ite = bit_container::end();
//       bit_support x;
      icard = 0;
//       for (; itb != ite; ++itb) { 
//         if ((x = *itb)) // fast count of the nb of bits
// 	                // of an integer (faster than shifts)  
//           do icard++; while ((x &= x-1));
//       }
//       icard_valid = true;

      const_iterator itb = begin(), ite = end();
      icard = 0;
      while (itb != ite) { if (*itb) ++icard; ++itb; }
      icard_valid = true;
    }
    return icard;
  }

  bit_vector::size_type bit_vector::first_true(void) const {
    assert(ifirst_true <= ilast_true);
    const_iterator itx = begin(), ite = end(); itx += ifirst_true;
    while (itx != ite && !*itx ) { ++itx; ++(ifirst_true); }
    if (is_in(ifirst_true)) return ifirst_true;
    else { ifirst_true = ilast_true = 0; return size_type(-1); }
  }
  
  bit_vector::size_type bit_vector::first_false(void) const {
    const_iterator itx = begin(), ite = end(); itx += ifirst_false;
    while (itx != ite && *itx) { ++itx; ++(ifirst_false); }
    if (!is_in(ifirst_false)) return ifirst_false;
    else { ifirst_false = ilast_false = size()-1; return size_type(-1); }
  }

  bit_vector::size_type bit_vector::last_true(void) const {
    const_iterator itb = begin(), itx = itb; itx += ilast_true;
    while (itx != itb && !*itx) { --itx; --(ilast_true); }
    if (is_in(ilast_true)) return ilast_true;
    else return size_type(-1);
  }
  
  bit_vector::size_type bit_vector::last_false(void) const {
    const_iterator itb = begin(), itx = itb; itx += ilast_false;
    while (itx != itb && *itx) { --itx; --(ilast_false); }
    return ilast_false;
  }
  
  bit_vector &bit_vector::setminus(const bit_vector& b) {
    for (bv_visitor i(b); !i.finished(); ++i) del(i); 
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
    ifirst_true = std::max(ifirst_true, bv.ifirst_true);
    ilast_true = std::min(ilast_true, bv.ilast_true);
    if (ifirst_true > ilast_true) clear();
    else {
      ilast_false = std::min(size()-1, std::max(ilast_false,bv.ilast_false));
      ifirst_false = std::min(ifirst_false, bv.ifirst_false);
    }
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
    if (nb)
      { add(i+nb-1); std::fill(this->begin()+i, this->begin()+(i+nb), true); }
  }

  void bit_vector::sup(size_type i, size_type nb) {
    if (nb)
      { del(i+nb-1); std::fill(this->begin()+i, this->begin()+(i+nb), false); }
  }

  void bit_vector::del(size_type i, size_type nb) {
    if (nb)
      { del(i+nb-1); std::fill(this->begin()+i, this->begin()+(i+nb), false); }
  }

  bool bit_vector::contains(const dal::bit_vector& other) const {
    for (dal::bv_visitor i(other); !i.finished(); ++i)
      if (!this->is_in(i)) return false;
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

  bool bv_visitor::operator++() {
    while (1) {
      size_type ind_b = (ind&(~WD_MASK));
      while (v) {
	++ind; v >>= 1;
	if (v&1) return true;
      }
      ind = ind_b + WD_BIT;
      if (ind >= ilast) return false; 
      v = *(++it);
      if (v&1) return true;
    }
  }

  const_bv_iterator<bv_iterable> bv_iterable::begin() const{
    return const_bv_iterator<bv_iterable>(this, v_.first_true());
  }
  const_bv_iterator<bv_iterable> bv_iterable::end() const{
    return const_bv_iterator<bv_iterable>(this, v_.last_true() + 1);
  }

  const_bv_iterator<bv_iterable_c> bv_iterable_c::begin() const{
    return const_bv_iterator<bv_iterable_c>(this, v_.first_true());
  }
  const_bv_iterator<bv_iterable_c> bv_iterable_c::end() const{
    return const_bv_iterator<bv_iterable_c>(this, v_.last_true() + 1);
  }
}

