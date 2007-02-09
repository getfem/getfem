// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfemint.h>
#include <getfemint_convex_structure.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_cvstruct_get(cs, ...)
    General function for querying information about convex_structure objects.

    @RDATTR CVSTRUCT:GET('nbpts')
    @RDATTR CVSTRUCT:GET('dim')
    @RDATTR CVSTRUCT:GET('basic structure')
    @RDATTR CVSTRUCT:GET('face')
    @GET CVSTRUCT:GET('facepts')
MLABCOM*/

void gf_cvstruct_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  bgeot::pconvex_structure cs = in.pop().to_convex_structure();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "nbpts", in, out, 0, 0, 0, 1)) {
    /*@RDATTR CVSTRUCT:GET('nbpts')
      Get the number of points of the convex structure.@*/
    out.pop().from_scalar(cs->nb_points());
  } else if (check_cmd(cmd, "dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR CVSTRUCT:GET('dim')
      Get the dimension of the convex structure.@*/
    out.pop().from_scalar(cs->dim());
  } else if (check_cmd(cmd, "basic_structure", in, out, 0, 0, 0, 1)) {
    /*@RDATTR  cs=CVSTRUCT:GET('basic structure')
      Get the simplest structure.
      
      For example, the 'basic structure' of the 6-node triangle, is the
      canonical 3-noded triangle.@*/
    out.pop().from_object_id(getfemint::ind_convex_structure(cs->basic_structure()),
			     CVSTRUCT_CLASS_ID);
  } else if (check_cmd(cmd, "face", in, out, 1, 1, 0, 1)) {
    /*@RDATTR  CVSTRUCT:GET('face', @int F)
      Return the structure of the face F.@*/
    size_type f = in.pop().to_integer(1,cs->nb_faces()) - config::base_index();
    out.pop().from_object_id(getfemint::ind_convex_structure(cs->faces_structure()[f]),
			     CVSTRUCT_CLASS_ID);    
  } else if (check_cmd(cmd, "facepts", in, out, 1, 1, 0, 1)) {
    /*@GET I=CVSTRUCT:GET('facepts', @int F)
      Return the list of point indices for the face F.@*/
    short_type f = short_type(in.pop().to_integer(1,cs->nb_faces()) - config::base_index());
    iarray w = out.pop().create_iarray_h(cs->nb_points_of_face(f));
    for (size_type i=0; i < w.size(); ++i) w[i] = cs->ind_points_of_face(f)[i]+config::base_index();
  } else bad_cmd(cmd);
}


