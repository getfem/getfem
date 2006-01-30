// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_tas.h : tas structure on an array.
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

/**@file dal_tas.h
   @brief Heap implementation.
*/
#ifndef DAL_TAS_H__
#define DAL_TAS_H__

#include <dal_basic.h>
#include <dal_bit_vector.h>

namespace dal
{

  template<class T, unsigned char pks = 5> class dynamic_tas;

  template<class T, unsigned char pks = 5> struct dnt_iterator
  {
    typedef T             value_type;
    typedef value_type*   pointer;
    typedef value_type&   reference;
    typedef size_t        size_type;
    typedef ptrdiff_t     difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    typename dynamic_array<T,pks>::iterator id;
    bit_vector::iterator ib;
    size_type lt;
    
    dnt_iterator(void) {}
    dnt_iterator(dynamic_tas<T,pks> &da, bit_vector &bv, size_type ii)
      : id(da, ii), ib(bv, ii) { lt = da.index().last_true(); }
    
    inline size_type index(void) const { return id.index(); }
    
    dnt_iterator &operator ++();
    dnt_iterator &operator --()
    { while (!*(--ib)) --id; --id; return *this; }
    dnt_iterator operator ++(int)
    {  dnt_iterator tmp = *this; ++(*this); return tmp; }
    dnt_iterator operator --(int)
    { dnt_iterator tmp = *this; --(*this); return tmp; }

    difference_type operator -(const dnt_iterator &i) const
    { return id - i.id; }

	
    reference operator *() const { return (*id); }
    pointer operator->() const { return &(operator*()); }

    bool operator ==(const dnt_iterator &i) const { return i.id==id;}
    bool operator !=(const dnt_iterator &i) const { return i.id!=id;}
    bool operator < (const dnt_iterator &i) const { return id <i.id;}
  };
  
  template<class T, unsigned char pks> dnt_iterator<T, pks> &
  dnt_iterator<T, pks>::operator ++()
  { ++ib; ++id; while(id.in <= lt && !*ib) {++ib; ++id; } return *this; }

  template<class T, unsigned char pks = 5> struct dnt_const_iterator {
    typedef T                   value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    
    typename dynamic_array<T,pks>::const_iterator id;
    bit_vector::const_iterator ib;
    size_type lt; 
    
    dnt_const_iterator(void) {}
    dnt_const_iterator(const dynamic_tas<T,pks> &da, size_type ii)
      : id(da, ii), ib(da.index(), ii) { lt = da.index().last_true(); }
    dnt_const_iterator(const dnt_iterator<T,pks> &it)
      : id(it.id), ib(it.ib), lt(it.lt) { }
    
    inline size_type index(void) const { return id.index(); }

    dnt_const_iterator &operator ++()
    { ++ib; ++id; while(id.in <= lt && !*ib) {++ib; ++id; } return *this; }
    dnt_const_iterator &operator --()
    { while (!*(--ib)) --id; --id; return *this; }
    dnt_const_iterator operator ++(int)
    {  dnt_const_iterator tmp = *this; ++(*this); return tmp; }
    dnt_const_iterator operator --(int)
    { dnt_const_iterator tmp = *this; --(*this); return tmp; }

    difference_type operator -(const dnt_const_iterator &i) const
    { return id - i.id; }
	
    reference operator *() const { return (*id); }
    pointer operator->() const { return &(operator*()); }
    
    bool operator ==(const dnt_const_iterator &i) const
    { return i.id == id;}
    bool operator !=(const dnt_const_iterator &i) const
    { return i.id != id;}
    bool operator < (const dnt_const_iterator &i) const
    { return id < i.id;}
  };

  template<class T, unsigned char pks> class dynamic_tas
    : public dynamic_array<T, pks> {
  protected :
    bit_vector ind;
    
  public :
    typedef typename dynamic_array<T, pks>::iterator iterator;
    typedef typename dynamic_array<T, pks>::const_iterator const_iterator;
    typedef dnt_iterator<T, pks> tas_iterator;
    typedef dnt_const_iterator<T, pks> const_tas_iterator;
    typedef typename dynamic_array<T, pks>::size_type size_type;
    
    size_type memsize(void) const
    {	return dynamic_array<T, pks>::memsize() + ind.memsize(); }
    size_type size(void) const
    { return (ind.card() == 0) ? 0 : (ind.last_true() + 1); }
    size_type ind_first(void) const
    { return (ind.card() == 0) ? 0 : ind.first_true(); }
    size_type ind_last(void) const
    { return (ind.card() == 0) ? 0 : ind.last_true(); }
    size_type card(void) const { return ind.card(); }
    
    tas_iterator tas_begin(void)
    { return tas_iterator(*this, ind, ind_first()); }
    const_tas_iterator tas_begin(void) const
    { return const_tas_iterator(*this, ind_first()); }
    tas_iterator tas_end(void) { return tas_iterator(*this, ind, size()); }
    const_tas_iterator tas_end(void) const
    { return const_tas_iterator(*this, size()); }
    
    const bit_vector &index(void) const { return ind; }
    bool index_valid(size_type i) const { return ind[i]; }
    bool empty(void) const { return (ind.card() == 0); }
    
    void swap(size_type i, size_type j);
    void compact(void);
    size_type add(const T &e)
    { size_type n=ind.first_false(); ind[n]=true; (*this)[n]=e;  return n; }
    void add_to_index(size_type i, const T &e) { ind[i]=true; (*this)[i]=e; }
    void sup(size_type n) { ind[n] = false; }
    void clear(void) { dynamic_array<T,pks>::clear(); ind.clear(); }
  };

  template<class T, unsigned char pks>
    void dynamic_tas<T, pks>::swap(size_type i, size_type j) {
    bool ti = ind[i], tj = ind[j]; ind.swap(i,j);
    if (!ti &&  tj) (*this)[i] = (*this)[j];
    if (ti  && !tj) (*this)[j] = (*this)[i];
    if (ti  &&  tj) std::swap((*this)[i], (*this)[j]);
  }

  template<class T, unsigned char pks>
    void dynamic_tas<T, pks>::compact(void) {
    if (!empty())
      while (ind.last_true() >= ind.card())
	swap(ind.first_false(), ind.last_true());
  }
}

#endif /* DAL_TAS_H__ */
