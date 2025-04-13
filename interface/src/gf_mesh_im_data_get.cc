/*===========================================================================

 Copyright (C) 2014-2025 Konstantinos Poulios.

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

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_im_data.h>


using namespace getfemint;


/*@GFDOC
  General function extracting information from mesh_im_data objects.
@*/

void gf_mesh_im_data_get(getfemint::mexargs_in& in,
                         getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  getfem::im_data *mimd = to_meshimdata_object(in.pop());
  std::string init_cmd  = in.pop().to_string();
  std::string cmd       = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "region", in, out, 0, 0, 0, 1)) {
    /*@GET ('region')
      Output the region that the @tmimd is restricted to.@*/
    out.pop().from_integer(int(mimd->filtered_region()));
  } else if (check_cmd(cmd, "nbpts", in, out, 0, 0, 0, 1)) {
    /*@GET ('nbpts')
      Output the number of integration points (filtered in the considered
      region).@*/
    out.pop().from_integer(int(mimd->nb_filtered_index()));
  } else if (check_cmd(cmd, "nb tensor elements", in, out, 0, 0, 0, 1)) {
    /*@GET ('nb tensor elements')
      Output the size of the stored data (per integration point).@*/
     if (mimd->tensor_size().size())
       out.pop().from_integer(int(mimd->nb_tensor_elem()));
  } else if (check_cmd(cmd, "tensor size", in, out, 0, 0, 0, 1)) {
    /*@GET ('tensor size')
      Output the dimensions of the stored data (per integration point).@*/
    if (mimd->tensor_size().size()) {
      iarray oidx = out.pop().create_iarray_h(unsigned(mimd->tensor_size().size()));
      std::copy(mimd->tensor_size().begin(),
                mimd->tensor_size().end(), &oidx[0]);
    }
  } else if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      displays a short summary for a @tmimd object.@*/
    infomsg() << "gfMeshImData object containing data of size "
              << mimd->tensor_size() << " on a mesh in dimension "
              << int(mimd->linked_mesh_im().linked_mesh().dim())
              << " with " << mimd->linked_mesh_im().linked_mesh().nb_points()
              << " points and " << mimd->linked_mesh_im().linked_mesh()
                                                         .convex_index().card()
              << " elements\n";
  } else if (check_cmd(cmd, "linked mesh", in, out, 0, 0, 0, 1)) {
    /*@GET m = ('linked mesh')
    Returns a reference to the @tmesh object linked to `mim`.@*/
    id_type id = workspace().object(&mimd->linked_mesh_im().linked_mesh());
    if (id == id_type(-1)) {
      auto pst = workspace().hidden_object
        (workspace().object(&mimd->linked_mesh_im()),
         &mimd->linked_mesh_im().linked_mesh());
      if (!pst.get()) THROW_INTERNAL_ERROR;
      std::shared_ptr<getfem::mesh> pm =
        std::const_pointer_cast<getfem::mesh>
        (std::dynamic_pointer_cast<const getfem::mesh>(pst));
      id = store_mesh_object(pm);
    }
    out.pop().from_object_id(id, MESH_CLASS_ID);
  } else
    bad_cmd(init_cmd);

}
