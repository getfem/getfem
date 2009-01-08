// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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

#include <getfemint.h>
#include <getfemint_pgt.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_geotrans_get(GT, ...)
    General function for querying information about geometric transformations
    objects.

    @RDATTR GEOTRANS:GET('dim')
    @RDATTR GEOTRANS:GET('is_linear')
    @RDATTR GEOTRANS:GET('nbpts')
    @GET GEOTRANS:GET('pts')
    @GET GEOTRANS:GET('normals')
    @GET GEOTRANS:GET('transform')
    @GET GEOTRANS:GET('char')
MLABCOM*/

void gf_geotrans_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  bgeot::pgeometric_trans pgt = in.pop().to_pgt();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR d = GEOTRANS:GET('dim')
    Get the dimension of the @tgt.

    This is the dimension of the source space, i.e. the dimension of
    the reference convex.@*/
    out.pop().from_scalar(pgt->dim());
  } else if (check_cmd(cmd, "is_linear", in, out, 0, 0, 0, 1)) {
    /*@RDATTR b = GEOTRANS:GET('is_linear')
    Return 0 if the @tgt is not linear.@*/
    out.pop().from_scalar(pgt->is_linear() ? 1. : 0.);
  } else if (check_cmd(cmd, "nbpts", in, out, 0, 0, 0, 1)) {
    /*@RDATTR n = GEOTRANS:GET('nbpts')
    Return the number of points of the @tgt.@*/
    out.pop().from_scalar(double(pgt->nb_points()));
  } else if (check_cmd(cmd, "pts", in, out, 0, 0, 0, 1)) {
    /*@GET P = GEOTRANS:GET('pts')
    Return the reference convex points of the @tgt.

    The points are stored in the columns of the output matrix.@*/
    out.pop().from_vector_container(pgt->convex_ref()->points());
  } else if (check_cmd(cmd, "normals", in, out, 0, 0, 0, 1)) {
    /*@GET N = GEOTRANS:GET('normals')
    Get the normals for each face of the reference convex of the @tgt.

    The normals are stored in the columns of the output matrix.@*/
    out.pop().from_vector_container(pgt->normals());
#if 0
  } else if (check_cmd(cmd, "poly_str", in, out, 0, 0, 0, 1)) {
    /*@GET GEOTRANS:GET('poly_str')
    Return the @tgt expressed as polynomials

    The result is expressed as a @MATLAB{cell array}@PYTHON{tuple} of
    strings.@*/
    std::vector<std::string> s(pgt->poly_vector().size());
    for (size_type i=0; i < s.size(); ++i) {
      std::stringstream ss; ss << pgt->poly_vector()[i];
      s[i] = ss.str();
    }
    out.pop().from_string_container(s);
#endif
  } else if (check_cmd(cmd, "transform", in, out, 2, 2, 0, 1)) {
    /*@GET Pt = GEOTRANS:GET('transform',@mat G, @mat Pr)
    Apply the @tgt to a set of points.

    `G` is the set of vertices of the real convex, `Pr` is the set
    of points (in the reference convex) that are to be transformed.
    The corresponding set of points in the real convex is returned.@*/
    getfem::base_matrix G = in.pop().to_darray(-1, -1).row_col_to_bm();
    darray pts = in.pop().to_darray(pgt->dim(), -1);
    darray w = out.pop().create_darray(unsigned(G.nrows()), pts.getn());
    for (unsigned i=0; i < pts.getn(); ++i) {
      getfem::base_node P = pgt->transform(pts.col_to_bn(i), G);
      for (size_type k=0; i < P.size(); ++k) w(k,i) = P[k];
    }
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET s = GEOTRANS:GET('char')
    Output a (unique) string representation of the @tgt.

    This can be used to perform comparisons between two
    different @tgt objects. @*/
    std::string s = bgeot::name_of_geometric_trans(pgt);
    out.pop().from_string(s.c_str());
  } else bad_cmd(cmd);
}
