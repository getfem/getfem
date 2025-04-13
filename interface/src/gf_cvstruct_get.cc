/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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
#include <getfem/bgeot_convex_structure.h>

using namespace getfemint;

/*@GFDOC
  General function for querying information about convex_structure objects.

  The convex structures are internal structures of GetFEM. They do not
  contain points positions. These structures are recursive, since the faces
  of a convex structures are convex structures.
@*/

void gf_cvstruct_get(getfemint::mexargs_in& in,
                     getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG( "Wrong number of input arguments");

  bgeot::pconvex_structure cs = to_cvstruct_object(in.pop());
  std::string init_cmd        = in.pop().to_string();
  std::string cmd             = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "nbpts", in, out, 0, 0, 0, 1)) {
    /*@RDATTR n = ('nbpts')
      Get the number of points of the convex structure.@*/
    out.pop().from_scalar(cs->nb_points());
  } else if (check_cmd(cmd, "dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR d = ('dim')
      Get the dimension of the convex structure.@*/
    out.pop().from_scalar(cs->dim());
  } else if (check_cmd(cmd, "basic_structure", in, out, 0, 0, 0, 1)) {
    /*@RDATTR cs = ('basic structure')
      Get the simplest convex structure.

      For example, the 'basic structure' of the 6-node triangle, is the
      canonical 3-noded triangle.@*/
    out.pop().from_object_id(store_cvstruct_object(bgeot::basic_structure(cs)),
                             CVSTRUCT_CLASS_ID);
  } else if (check_cmd(cmd, "face", in, out, 1, 1, 0, 1)) {
    /*@RDATTR cs = ('face', @int F)
      Return the convex structure of the face `F`.@*/
    short_type f = in.pop().to_face_number(cs->nb_faces());
    out.pop().from_object_id(store_cvstruct_object(cs->faces_structure()[f]),
                             CVSTRUCT_CLASS_ID);
  } else if (check_cmd(cmd, "facepts", in, out, 1, 1, 0, 1)) {
    /*@GET I = ('facepts', @int F)
      Return the list of point indices for the face `F`.@*/
    short_type f = short_type(in.pop().to_face_number(cs->nb_faces()));
    iarray w = out.pop().create_iarray_h(cs->nb_points_of_face(f));
    for (size_type i=0; i < w.size(); ++i)
      w[i] = cs->ind_points_of_face(f)[i]+config::base_index();
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET s = ('char')
      Output a string description of the @tcvstruct.@*/
    THROW_ERROR("No output format for a convex structure. To be done.");
  } else if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      displays a short summary for a @tcvstruct object.@*/
    infomsg() << "gfCvStruct (convex structure) in dimension "
              << int(cs->dim()) << " with " << cs->nb_points() << "points.\n";
  } else
    bad_cmd(init_cmd);
}
