/*===========================================================================

 Copyright (C) 2005-2020 Julien Pommier.

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
// $Id$
#include <getfemint.h>
#include <getfem/getfem_mesh_level_set.h>
#include <getfemint_workspace.h>

using namespace getfemint;

/*@GFDOC
  General constructor for mesh_levelset objects. The role of this object is
  to provide a mesh cut by a certain number of level_set. This object is
  used to build conformal integration method (object mim and enriched finite
  element methods (Xfem)).
  @*/

void gf_mesh_levelset(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (check_cmd("MeshLevelSet", "MeshLevelSet", in, out, 1, 1, 0, 1)) {
    /*@INIT MLS = ('.mesh', @tmesh m)
      Build a new @tmls object from a @tmesh and returns its handle. @*/
    getfem::mesh *mm = extract_mesh_object(in.pop());
    id_type id = store_mesh_levelset_object
      (std::make_shared<getfem::mesh_level_set>(*mm));

    workspace().set_dependence(id, mm);
    out.pop().from_object_id(id, MESH_LEVELSET_CLASS_ID);
  }
}
