/*===========================================================================

 Copyright (C) 2001-2020 Yves Renard.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <getfemint.h>
#include <getfemint_workspace.h>

using namespace getfemint;

/*@GFDOC
  Delete an existing getfem object from memory (mesh, mesh_fem, etc.).

  SEE ALSO:
    gf_workspace, gf_mesh, gf_mesh_fem.
 @*/

void gf_delete(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.remaining() < 1)
    THROW_BADARG("Wrong number of input arguments, should be at least 1.");
  if (!out.narg_in_range(0,0))
    THROW_BADARG("No output argument needed.");
  
  while (in.remaining()) {
    id_type id;
    /*@FUNC ('.list', I[, J, K,...])
      
      I should be a descriptor given by gf_mesh(),
      gf_mesh_im(), gf_slice() etc.
      
      Note that if another object uses I, then object I will be deleted only
      when both have been asked for deletion.
      
      Only objects listed in the output of gf_workspace('stats') can be
      deleted (for example gf_fem objects cannot be destroyed).
      
      You may also use gf_workspace('clear all') to erase everything at
      once.
      @*/
    
    if (in.front().is_object_id()) {
      id_type cid; in.pop().to_object_id(&id,&cid);
      /*if (is_static_object(id, cid)) {
	THROW_BADARG("sorry, this object of type " << 
		     name_of_getfemint_class_id(cid) << 
		     " is static, i.e. it can't be deleted");
      }
      */
    } else if (in.front().is_integer()) {
      id = in.pop().to_integer();
    }
    if (getfemint::workspace().object(id)) {
      /* throw an exeption if object not found */
      getfemint::workspace().delete_object(id);
    } else {
      GFI_WARNING("ouuups strange");
    }
  }
}

