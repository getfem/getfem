/*===========================================================================
 
 Copyright (C) 2009-2015 Luis Saavedra.
 
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
// $Id$
#include <getfemint_global_function.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_global_function::getfemint_global_function(getfem::abstract_xy_function *pabs) {
    assert(workspace == 0);
    pgf  = pabs;
    ikey = getfem_object::internal_key_type(pgf);
  }

  getfemint_global_function::~getfemint_global_function() {
    if (!is_static()) delete pgf;
    pgf = NULL;
  }

  getfemint_global_function*
  getfemint_global_function::get_from(getfem::abstract_xy_function *pabs, int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(pabs));
    getfemint_global_function *gpgf = NULL;
    if (!o) {
      gpgf = new getfemint_global_function(pabs);
      gpgf->set_flags(flags);
      getfemint::workspace().push_object(gpgf);
    } else gpgf = dynamic_cast<getfemint_global_function*>(o);
    assert(gpgf);
    return gpgf;
  }
}
