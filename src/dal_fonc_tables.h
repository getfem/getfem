/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_fonc_table.h : deals fonctionalities tables, reservation */
/*            search ...                                                   */
/*                                                                         */
/* Date : August 28, 2001                                                  */
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

#ifndef __DAL_FONC_TABLES_H
#define __DAL_FONC_TABLES_H

#include <dal_tree_sorted.h>

namespace dal
{
  /* ********************************************************************* */
  /* Assumptions :                                                         */
  /*   1 - class LIGHT as a default comparator <.                          */
  /*   2 - class DESC as a constructor with a class LIGHT.                 */
  /* ********************************************************************* */
  
  template<class LIGHT, class DESC> class FONC_TABLE {
  public :
    
    typedef DESC * pDESC;
    typedef typename dynamic_tree_sorted<LIGHT>::size_type size_type;
    typedef dynamic_array<pDESC, 2> desc_table_type;
    
  protected :
    
    dynamic_tree_sorted<LIGHT> _light_table;
    desc_table_type desc_table;
    
  public :
    
    size_type search(const LIGHT &l) const { return _light_table.search(l); }
    
    pDESC add(const LIGHT &l) {
      size_type i = _light_table.search(l);
      if (i == size_type(-1))
	{ i = _light_table.add(l); desc_table[i] = new DESC(l); }
      return desc_table[i];
    }
    void sup(const LIGHT &l) {
      size_type i = _light_table.search(l);
      if (i != size_type(-1))
	{ _light_table.sup(i); delete desc_table[i]; desc_table[i] = 0;}
    }
    const desc_table_type &table(void) { return desc_table; }
    const dynamic_tree_sorted<LIGHT> &light_table(void)
      { return _light_table; }
    const bit_vector &index(void) { return _light_table.index(); }
    ~FONC_TABLE(void) { 
      dal::bit_vector nn = _light_table.index();
      size_type i;
      for (i << nn; i != size_type(-1); i << nn) delete desc_table[i];
    }
  };
  
}

#endif /* __DAL_FONC_TABLES_H */
