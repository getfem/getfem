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

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_mesh_levelset_get(MLS, ...)
    General function for querying information about MESHLEVELSET objects.

  @GET M = MESHLEVELSET:GET('cut_mesh')
  @GET LM = MESHLEVELSET:GET('linked_mesh')
  @GET nbls = MESHLEVELSET:GET('nb_ls')
  @GET LS = MESHLEVELSET:GET('levelsets')
  @GET CVIDs = MESHLEVELSET:GET('crack_tip_convexes')
  @GET SIZE = MESHLEVELSET:GET('memsize')

  $Id$
MLABCOM*/

void gf_mesh_levelset_get(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mesh_levelset *gmls = in.pop().to_getfemint_mesh_levelset();
  getfem::mesh_level_set &mls = gmls->mesh_levelset();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "cut_mesh", in, out, 0, 0, 0, 1)) {
    /*@GET M = MESHLEVELSET:GET('cut_mesh')
    Return a @tmesh cut by the linked @tls's.@*/
    getfemint_mesh *mm = getfemint_mesh::get_from(new getfem::mesh);
    mls.global_cut_mesh(mm->mesh());
    out.pop().from_object_id(mm->get_id(), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "linked_mesh", in, out, 0, 0, 0, 1)) {
    /*@GET LM = MESHLEVELSET:GET('linked_mesh')
    Return a reference to the linked @tmesh.@*/
    getfemint_mesh *mm = getfemint_mesh::get_from(&mls.linked_mesh());
    out.pop().from_object_id(mm->get_id(), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "nb_ls", in, out, 0, 0, 0, 1)) {
    /*@GET nbls = MESHLEVELSET:GET('nb_ls')
    Return the number of linked @tls's.@*/
    out.pop().from_integer(int(mls.nb_level_sets()));
  } else if (check_cmd(cmd, "levelsets", in, out, 0, 0, 0, 1)) {
    /*@GET LS = MESHLEVELSET:GET('levelsets')
    Return a list of references to the linked @tls's.@*/
    std::vector<id_type> ids;
    for (unsigned i=0; i < mls.nb_level_sets(); ++i) {
      getfemint_levelset *gls =
        getfemint_levelset::get_from(&(*(mls.get_level_set(i))));
      ids.push_back(gls->get_id());
    }
    out.pop().from_object_id(ids, LEVELSET_CLASS_ID);
  } else if (check_cmd(cmd, "crack_tip_convexes", in, out, 0, 0, 0, 1)) {
    /*@GET CVIDs = MESHLEVELSET:GET('crack_tip_convexes')
    Return the list of convex #id's of the linked @tmesh on
    which have a tip of any linked @tls's.@*/
    out.pop().from_bit_vector(mls.crack_tip_convexes());
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET SIZE = MESHLEVELSET:GET('memsize')
    Return the amount of memory (in bytes) used by the @tmls.@*/
    out.pop().from_integer(int(mls.memsize()));
  } else bad_cmd(cmd);
}
