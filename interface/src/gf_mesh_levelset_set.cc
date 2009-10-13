// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Julien Pommier.
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
// $Id$
#include <getfemint.h>
#include <getfemint_mesh_levelset.h>
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_mesh_levelset_get(MLS, ...)
    General function for modification of MESHLEVELSET objects.

  $Id$
MLABCOM*/

void gf_mesh_levelset_set(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mesh_levelset *gmls = in.pop().to_getfemint_mesh_levelset(true);
  getfem::mesh_level_set &mls = gmls->mesh_levelset();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "add", in, out, 1, 1, 0, 0)) {
    /*@SET MESHLEVELSET:SET('add',@tls ls)
    Add a link to the @tls `ls`.

    Only a reference is kept, no copy is done. In order to indicate
    that the linked @tmesh is cut by a @tls one has to call this
    method, where `ls` is an @tls object. An arbitrary number of
    @tls can be added.

    **WARNING**<par>

    The @tmesh of `ls` and the linked @tmesh must be the same.@*/
    getfemint_levelset *gls = in.pop().to_getfemint_levelset();
    if (&mls.linked_mesh() != &gls->levelset().get_mesh_fem().linked_mesh())
      THROW_BADARG("The meshes of the levelset and the mesh_levelset "
                   "are not the same!");
    mls.add_level_set(gls->levelset());
    workspace().set_dependance(gmls, gls);
  } else if (check_cmd(cmd, "sup", in, out, 1, 1, 0, 0)) {
    /*@SET MESHLEVELSET:SET('sup', @tls ls)
    Remove a link to the @tls `ls`.@*/
    getfemint_levelset *gls = in.pop().to_getfemint_levelset();
    mls.sup_level_set(gls->levelset());
    workspace().sup_dependance(gmls, gls);
  } else if (check_cmd(cmd, "adapt", in, out, 0, 0, 0, 0)) {
    /*@SET MESHLEVELSET:SET('adapt')
    Do all the work (cut the convexes with the levelsets).

    To initialice the @tmls object or to actualize it when the
    value of any levelset function is modified, one has to call
    this method.@*/
    mls.adapt();
  } else bad_cmd(cmd);
}
