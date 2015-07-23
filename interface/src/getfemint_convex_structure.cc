/*===========================================================================

 Copyright (C) 2001-2015 Yves Renard.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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


#include <getfemint_convex_structure.h>
#include <getfem/dal_tree_sorted.h>

namespace getfemint
{
  static dal::dynamic_tree_sorted<bgeot::pconvex_structure> *convex_structure_tab;

  static inline void init_tab(void) // because of problem with initialization
  {                                 // in dynamic libraries.
    static bool initialized = false;
    if (!initialized)
    { 
      initialized = true;
      convex_structure_tab = new dal::dynamic_tree_sorted<bgeot::pconvex_structure>();
    }
  }

  id_type ind_convex_structure(bgeot::pconvex_structure p)
  { init_tab(); return id_type(convex_structure_tab->add_norepeat(p)); }
  
  bgeot::pconvex_structure addr_convex_structure(id_type i)
  { init_tab(); return (*convex_structure_tab)[i]; }

  bool exists_convex_structure(id_type i)
  { init_tab(); return convex_structure_tab->index()[i]; }


}  /* end of namespace getfemint.                                          */
