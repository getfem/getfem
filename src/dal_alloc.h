/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_alloc.h : allocation in an array.                        */
/*     									   */
/*                                                                         */
/* Date : March 06, 1997.                                                  */
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


#ifndef __DAL_ALLOC_H
#define __DAL_ALLOC_H

#include <dal_tree_sorted.h>


namespace dal
{

  struct _fr_sp { size_t ind, size; };
  struct _less1_fr_sp : public std::binary_function<_fr_sp, _fr_sp, int>
  {
    int operator()(const _fr_sp &f1, const _fr_sp &f2) const
      { return (f1.ind < f2.ind) ? -1 : ((f1.ind > f2.ind) ? 1 : 0); }
  };
  struct _less2_fr_sp : public std::binary_function<_fr_sp, _fr_sp, int>
  {
    int operator()(const _fr_sp &f1, const _fr_sp &f2) const
    { return (f1.size < f2.size) ? -1 : ((f1.size > f2.size) ? 1 : 0); }
  };
  

  template<class T,  unsigned char pks = 5> class dynamic_alloc
    : public dynamic_array<T, pks>
  {
    public :

      typedef typename dynamic_array<T, pks>::size_type size_type;
                                                    /* gcc needs this ... */

    protected :

      typedef dynamic_tree_sorted<_fr_sp, _less1_fr_sp, 6> fsptab_t;

      fsptab_t fr_tab;
      dynamic_tree_sorted_index<_fr_sp, fsptab_t, _less2_fr_sp, 6> ind_fr_tab;

      void init(void)
      { 
	_fr_sp fsp; fsp.size = size_type(-2); fsp.ind = 0;
	size_type i = fr_tab.add(fsp); ind_fr_tab.add(i);
      }
    
    public :

      dynamic_alloc(void) : ind_fr_tab(fr_tab) { init(); }
      dynamic_alloc(const dynamic_alloc<T, pks> &da)
           : dynamic_array<T, pks>(da), ind_fr_tab(da.ind_fr_tab),
             fr_tab(da.fr_tab)
      { ind_fr_tab.change_tab(fr_tab); }
      dynamic_alloc<T, pks> &operator =(const dynamic_alloc<T, pks> &da)
      {
	dynamic_array<T, pks>::operator =(da);
	fr_tab = da.fr_tab;
	ind_fr_tab = da.ind_fr_tab;
	ind_fr_tab.change_tab(fr_tab);
	return *this;
      }
      
      size_type alloc(size_type); 
      void free(size_type, size_type);
      void clear(void);
  };
  
  template<class T, unsigned char pks>
    void dynamic_alloc<T,pks>::clear(void)
  { dynamic_array<T,pks>::clear(); fr_tab.clear(); ind_fr_tab.clear(); init();}

  template<class T, unsigned char pks>
    dynamic_alloc<T,pks>::size_type dynamic_alloc<T,pks>::alloc(size_type size)
  {
    size_type res = ST_NIL;
    if (size > 0)
    {
      _fr_sp fsp; fsp.size = size; 
      size_type i = ind_fr_tab.search_ge(fsp);
      if (i != ST_NIL)
      {
	res = fr_tab[i].ind;
	assert(size <= fr_tab[i].size); // internal test
	if (size <  fr_tab[i].size)
	{
	  ind_fr_tab.sup(i); fr_tab[i].ind += size;
	  fr_tab[i].size -= size; ind_fr_tab.add(i);
	}
	else
	{ ind_fr_tab.sup(i); fr_tab.sup(i); }
      }
      else
      { assert(false); } // internal test
    }
    return res;
  }

  
  template<class T, unsigned char pks>
    void dynamic_alloc<T,pks>::free(size_type l, size_type size)
  {
    if (size > 0 && l != ST_NIL)
    {
      _fr_sp fsp; fsp.size = size; fsp.ind = l;
      size_type i = fr_tab.add(fsp);
      fsptab_t::const_sorted_iterator it1(fr_tab);
      fr_tab.find_sorted_iterator(i, it1);
      fsptab_t::const_sorted_iterator it2 = it1;
      size_type i1 = (++it1).index(), i2 = (--it2).index();
      if (i1 != ST_NIL && (*it1).ind <= l + size)
      { 
	fr_tab[i].size = (*it1).ind + (*it1).size - l;
	ind_fr_tab.sup(i1); fr_tab.sup(i1);
      }
      if (i2 != ST_NIL && (*it2).ind + (*it2).size >= l)
      {
	fr_tab[i].size = l + size - (*it2).ind;
	fr_tab[i].ind = (*it2).ind;
	ind_fr_tab.sup(i2); fr_tab.sup(i2);
      }
      ind_fr_tab.add(i);
    }
  }

}

#endif  /* __DAL_ALLOC_H */
