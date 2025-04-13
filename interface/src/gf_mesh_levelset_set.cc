/*===========================================================================

 Copyright (C) 2006-2020 Julien Pommier.

 This file is a part of GetFEM

 GetFEM is free software;  you can  redistribute it  and/or modify it under
 the  terms  of the  GNU  Lesser General Public License as published by the
 Free Software Foundation;  either version 3  of  the License,  or (at your
 option) any  later  version  along with  the GCC Runtime Library Exception
 either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/
// $Id$
#include <getfem/getfem_mesh_level_set.h>
#include <getfemint.h>
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>

using namespace getfemint;

/*@GFDOC
  General function for modification of @tmesh_levelset objects.
@*/


void gf_mesh_levelset_set(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  getfem::mesh_level_set &mls = *(to_mesh_levelset_object(in.pop()));
  std::string init_cmd        = in.pop().to_string();
  std::string cmd             = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "add", in, out, 1, 1, 0, 0)) {
    /*@SET ('add', @tls ls)
      Add a link to the @tls `ls`.

      Only a reference is kept, no copy is done. In order to indicate
      that the linked @tmesh is cut by a @tls one has to call this
      method, where `ls` is an @tls object. An arbitrary number of
      @tls can be added.

      **WARNING**

      The @tmesh of `ls` and the linked @tmesh must be the same.@*/
    getfem::level_set *gls = to_levelset_object(in.pop());
    if (&mls.linked_mesh() != &gls->get_mesh_fem().linked_mesh())
      THROW_BADARG("The meshes of the levelset and the mesh_levelset "
                   "are not the same!");
    mls.add_level_set(*gls);
    workspace().set_dependence(&mls, gls);
  } else if (check_cmd(cmd, "sup", in, out, 1, 1, 0, 0)) {
    /*@SET ('sup', @tls ls)
      Remove a link to the @tls `ls`.@*/
    getfem::level_set *gls = to_levelset_object(in.pop());
    mls.sup_level_set(*gls);
    workspace().sup_dependence(&mls, gls);
  } else if (check_cmd(cmd, "adapt", in, out, 0, 0, 0, 0)) {
    /*@SET ('adapt')
      Do all the work (cut the convexes with the levelsets).

      To initialice the @tmls object or to actualize it when the
      value of any levelset function is modified, one has to call
      this method.@*/
    mls.adapt();
  } else
    bad_cmd(init_cmd);

}
