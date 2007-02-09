// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================


#include <getfemint_matelemtype.h>
#include <getfem/dal_tree_sorted.h>

namespace getfemint
{
  static dal::dynamic_tree_sorted<getfem::pmat_elem_type> *matelemtype_tab;

  static inline void init_tab(void) // because of problem with initialization
  {                                 // in dynamic libraries.
    static bool initialized = false;
    if (!initialized)
    { 
      initialized = true;
      matelemtype_tab = new dal::dynamic_tree_sorted<getfem::pmat_elem_type>();
    }
  }

  id_type ind_matelemtype(getfem::pmat_elem_type p)
  { init_tab(); return matelemtype_tab->add_norepeat(p); }
  
  getfem::pmat_elem_type addr_matelemtype(id_type i)
  { init_tab(); return (*matelemtype_tab)[i]; }

  bool exists_matelemtype(id_type i)
  { init_tab(); return matelemtype_tab->index()[i]; }


}  /* end of namespace getfemint.                                          */
