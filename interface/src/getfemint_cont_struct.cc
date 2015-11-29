/*===========================================================================

 Copyright (C) 2012-2015 Tomas Ligursky, Yves Renard, Konstantinos Poulios.

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
#include <getfemint_cont_struct.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_cont_struct::getfemint_cont_struct
  (getfem::cont_struct_getfem_model *s_) {
    assert(workspace == 0);
    s = s_;
    ikey = getfem_object::internal_key_type(s);
  }

  getfemint_cont_struct::~getfemint_cont_struct() {
    if (!is_static()) delete s;
    s = 0;
  }

  getfemint_cont_struct*
  getfemint_cont_struct::get_from(getfem::cont_struct_getfem_model *s_, int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(s_));
    getfemint_cont_struct *gs = 0;
    if (!o) {
      gs = new getfemint_cont_struct(s_);
      gs->set_flags(flags);
      getfemint::workspace().push_object(gs);
    } else gs = dynamic_cast<getfemint_cont_struct*>(o);
    assert(gs);
    return gs;
  }
}
