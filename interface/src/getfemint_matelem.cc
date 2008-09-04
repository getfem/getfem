// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2001-2008 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================


#include <getfemint_matelem.h>
#include <getfem/dal_tree_sorted.h>

namespace getfemint
{
  static dal::dynamic_tree_sorted<getfem::pmat_elem_computation> *matelem_tab;

  static inline void init_tab(void) // because of problem with initialization
  {                                 // in dynamic libraries.
    static bool initialized = false;
    if (!initialized)
    { 
      initialized = true;
      matelem_tab
	= new dal::dynamic_tree_sorted<getfem::pmat_elem_computation>();
    }
  }

  id_type ind_matelem(getfem::pmat_elem_computation p)
  { init_tab(); return id_type(matelem_tab->add_norepeat(p)); }
  
  getfem::pmat_elem_computation addr_matelem(id_type i)
  { init_tab(); return (*matelem_tab)[i]; }

  bool exists_matelem(id_type i)
  { init_tab(); return matelem_tab->index()[i]; }


}  /* end of namespace getfemint.                                          */
