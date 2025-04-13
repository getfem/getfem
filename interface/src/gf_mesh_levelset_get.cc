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
#include <getfemint_workspace.h>
#include <getfemint_levelset.h>

using namespace getfemint;

/*@GFDOC
  General function for querying information about @tmesh_levelset objects.
@*/


void gf_mesh_levelset_get(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  getfem::mesh_level_set &mls = *(to_mesh_levelset_object(in.pop()));
  std::string init_cmd        = in.pop().to_string();
  std::string cmd             = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "cut_mesh", in, out, 0, 0, 0, 1)) {
    /*@GET M = ('cut_mesh')
      Return a @tmesh cut by the linked @tls's.@*/
    auto mm = std::make_shared<getfem::mesh>();
    mls.global_cut_mesh(*mm);
    out.pop().from_object_id(store_mesh_object(mm), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "linked_mesh", in, out, 0, 0, 0, 1)) {
    /*@GET LM = ('linked_mesh')
      Return a reference to the linked @tmesh.@*/
    id_type id = workspace().object(&mls.linked_mesh());
    if (id == id_type(-1)) THROW_INTERNAL_ERROR;
    out.pop().from_object_id(id, MESH_CLASS_ID);
  } else if (check_cmd(cmd, "nb_ls", in, out, 0, 0, 0, 1)) {
    /*@GET nbls = ('nb_ls')
      Return the number of linked @tls's.@*/
    out.pop().from_integer(int(mls.nb_level_sets()));
  } else if (check_cmd(cmd, "levelsets", in, out, 0, 0, 0, 1)) {
    /*@GET LS = ('levelsets')
      Return a list of references to the linked @tls's.@*/
    std::vector<id_type> ids;
    for (unsigned i=0; i < mls.nb_level_sets(); ++i) {
      id_type id = workspace().object((const void *)(mls.get_level_set(i)));
      GMM_ASSERT1(id != id_type(-1), "Unknown levelset !");
      ids.push_back(id);
    }
    out.pop().from_object_id(ids, LEVELSET_CLASS_ID);
  } else if (check_cmd(cmd, "crack_tip_convexes", in, out, 0, 0, 0, 1)) {
    /*@GET CVIDs = ('crack_tip_convexes')
      Return the list of convex #id's of the linked @tmesh on
      which have a tip of any linked @tls's.@*/
    out.pop().from_bit_vector(mls.crack_tip_convexes());
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET SIZE = ('memsize')
      Return the amount of memory (in bytes) used by the @tmls.@*/
    out.pop().from_integer(int(mls.memsize()));
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET s = ('char')
      Output a (unique) string representation of the @tmlsn.

      This can be used to perform comparisons between two
      different @tmls objects.
      This function is to be completed.@*/
    GMM_ASSERT1(false, "Sorry, function to be done");
    // std::string s = ...;
    // out.pop().from_string(s.c_str());
  } else if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      displays a short summary for a @tmls object.@*/
    infomsg() << "gfMeshLevelSet object in dimension "
              << int(mls.linked_mesh().dim()) << " with "
              << mls.linked_mesh().nb_points() << " points and "
              << mls.linked_mesh().convex_index().card() << " elements\n";
  } else
    bad_cmd(init_cmd);

}
