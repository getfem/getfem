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


#include <dal_bit_vector.h>

namespace dal
{

  bit_reference& bit_reference::operator = (bool x)
  {
    if (x) { if (!(*p&mask)) { *p |= mask; bv->change_for_true(ind); } } 
    else  { if (*p & mask) { *p &= ~mask; bv->change_for_false(ind); } }
    return *this;
  }

  bit_iterator::bit_iterator(bit_vector &b, size_type i) : p(b, i / WD_BIT)
  { ind = i; bv = &b, mask = bit_support(1) << (i & WD_MASK); }

  bit_iterator& bit_iterator::operator+=(difference_type i)
  {
    ind+=i; mask = bit_support(1) << (ind & WD_MASK); 
    p=bit_container::iterator(*bv, ind/WD_BIT);
    return *this;
  }

  bit_const_iterator::bit_const_iterator(const bit_vector &b, size_type i)
    : p(b, i / WD_BIT)
  { ind = i; bv = &b, mask = bit_support(1) << (i & WD_MASK); }

  bit_const_iterator& bit_const_iterator::operator+=(difference_type i)
  {
    ind += i; mask = bit_support(1) << (ind & WD_MASK); 
    p = bit_container::const_iterator(*bv, ind/WD_BIT);
    return *this;
  }

  void bit_vector::fill_false(size_type i1, size_type i2)
  {
    size_type f = i1 / WD_BIT, r = i1 & (WD_BIT-1), l = i2 / WD_BIT;
    bit_support *q = & ((*((bit_container *)(this)))[l]);
    
    if (r != 0) f++; l++;
    if (f < l)
    {
      //    q = & ((*((bit_container *)(this)))[f]); 
      std::fill(bit_container::begin()+f, bit_container::begin()+l, 0);
    }

    // *q &= ((bit_support(-1)) >> (WD_BIT - r));

    // cout << "fill false " << i1 << " : " << i2 << " : " << *q << " : " << ((bit_support(-1)) >> (WD_BIT - r)) << endl; getchar();

    ilast_false = i2;
  }

  void bit_vector::swap(bit_vector &da)
  {
    ((bit_container *)(this))->swap(da);
    std::swap(ifirst_true, da.ifirst_true);
    std::swap(ifirst_false, da.ifirst_false);
    std::swap(ilast_true, da.ilast_true);
    std::swap(ilast_false, da.ilast_false);
    std::swap(icard, da.icard);
    std::swap(icard_valid, da.icard_valid);
  }
  
  void bit_vector::swap(size_type i1, size_type i2)
  {
    if (i1 != i2) 
    {
      reference r1 = (*this)[i1], r2 = (*this)[i2];
      bool tmp = r1; r1 = r2; r2 = tmp;
    }
  }

  void bit_vector::clear(void)
  {
    icard = 0; icard_valid = true;
    ifirst_false = ilast_false = ifirst_true = ilast_true = 0;
    fill_false(0,0); 
  }

  bit_vector::size_type bit_vector::card(void) const
  {
    if (!icard_valid)
    {
      size_type *pcard = (size_type *)&icard;
      bool *pvalid = (bool *)&icard_valid;
      const_iterator itb = begin(), ite = end();
      *pcard = 0;
      while (itb != ite) { if (*itb) ++(*pcard); ++itb; }
      *pvalid = true;
    }
    return icard;
  }

  bit_vector::size_type bit_vector::first_true(void) const
  {
    size_type *p = (size_type *)&ifirst_true;
    const_iterator itx = begin(), ite = end(); itx += *p;
    while (itx != ite && !*itx ) { ++itx; ++(*p); }
    return *p;
  }
  
  bit_vector::size_type bit_vector::first_false(void) const
  {
    size_type *p = (size_type *)&ifirst_false;
    const_iterator itx = begin(), ite = end(); itx += *p;
    while (itx != ite && *itx) { ++itx; ++(*p); }
    return *p;
  }

  bit_vector::size_type bit_vector::last_true(void) const
  {
    size_type *p = (size_type *)&ilast_true;
    const_iterator itb = begin(), itx = itb; itx += *p;
    while (itx != itb && !*itx) { --itx; --(*p); }
    return *p;
  }
  
  bit_vector::size_type bit_vector::last_false(void) const
  {
    size_type *p = (size_type *)&ilast_false;
    const_iterator itb = begin(), itx = itb; itx += *p;
    while (itx != itb && *itx) { --itx; --(*p); }
    return *p;
  }
  
  bit_vector &bit_vector::operator |=(const bit_vector &bv)
  {
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
  
  bit_vector &bit_vector::operator &=(const bit_vector &bv)
  {
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

  bool bit_vector::operator ==(const bit_vector &bv) const
  {
    bit_container::const_iterator it1b = bit_container::begin();
    bit_container::const_iterator it1e = bit_container::end();
    bit_container::const_iterator it2b = bv.bit_container::begin();
    bit_container::const_iterator it2e = bv.bit_container::end();
    
    while (it1b!=it1e && it2b!=it2e) if (*it1b++ != *it2b++) return false;
    while (it1b != it1e) if (*it1b++ != 0) return false;
    while (it2b != it2e) if (*it2b++ != 0) return false;
    return true;
  }


  void bit_vector::add(size_type i, size_type j)
  { 
    if (j < i) std::swap(i,j); add(j);
    std::fill(this->begin()+i, this->begin()+j+1, true);
  }

  void bit_vector::sup(size_type i, size_type j)
  { 
    if (j < i) std::swap(i,j); sup(j);
    std::fill(this->begin()+i, this->begin()+j+1, false);
  }

  STD_NEEDED ostream &operator <<(STD_NEEDED ostream &o, const bit_vector &s)
  {
    bit_vector u = s; int i;
    o << "ensemble d'entier";
    for (i << u; i != int(-1); i << u) o << " : " << i;
    o << " ";
    return o;
  }
}

