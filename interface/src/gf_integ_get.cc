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
#include <getfemint_integ.h>

using namespace getfemint;

static void check_not_exact(getfem::pintegration_method im) {
  if (im->type() != getfem::IM_APPROX) THROW_ERROR("this has no meaning for exact integration methods");
}
/*MLABCOM
  FUNCTION I = gf_integ_get(F, ...)
    General function for querying information about FEM integration method objects.

    @RDATTR INTEG:GET('is_exact')
    @RDATTR INTEG:GET('dim')
    @RDATTR INTEG:GET('nbpts')
    @GET INTEG:GET('pts')
    @GET INTEG:GET('coeffs')
    @GET INTEG:GET('face_pts')
    @GET INTEG:GET('face_coeffs')
    @GET INTEG:GET('char')
MLABCOM*/

void gf_integ_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfem::pintegration_method im = in.pop().to_integration_method();
  getfem::papprox_integration pai = im->type() == getfem::IM_APPROX ? im->approx_method() : 0;
  size_type imdim = 0;
  if (im->type() == getfem::IM_EXACT) imdim = im->exact_method()->dim();
  else if (im->type() == getfem::IM_APPROX) imdim = im->approx_method()->dim();

  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "is_exact", in, out, 0, 0, 0, 1)) {
    /*@RDATTR INTEG:GET('is_exact')
      Return 0 if the integration is an approximate one.
      @*/
    out.pop().from_scalar(im->type() != getfem::IM_APPROX ? 1. : 0.);
  } else if (check_cmd(cmd, "dim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR INTEG:GET('dim')
      Return the dimension of the ref. convex of the method.
      @*/
    out.pop().from_scalar(imdim);
  } else if (check_cmd(cmd, "nbpts", in, out, 0, 0, 0, 1)) {
    /*@RDATTR INTEG:GET('nbpts')
      Return the total number of integration points. 
      
      Count the points for the volume integration, and points for surface integration on each face of the reference convex. Raises an error for exact integration methods.
      @*/
   check_not_exact(im);
    iarray w = out.pop().create_iarray_h(1+pai->structure()->nb_faces());
    w[0] = pai->nb_points_on_convex();
    for (size_type i=0; i < pai->structure()->nb_faces(); ++i) 
      w[i+1] = pai->nb_points_on_face(i);
  } else if (check_cmd(cmd, "pts", in, out, 0, 0, 0, 1)) {
    /*@GET INTEG:GET('pts')
      Return the list of integration points (only for approximate methods).
      @*/
    check_not_exact(im);
    out.pop().from_vector_container(pai->integration_points());
  } else if (check_cmd(cmd, "face_pts", in, out, 1, 1, 0, 1)) {
   /*@GET INTEG:GET('face_pts',F)
     Return the list of integration points for a face.
     @*/
    check_not_exact(im); 
    size_type nbf = pai->structure()->nb_faces();
    size_type f = in.pop().to_face_number(nbf);
    darray w = out.pop().create_darray(imdim, pai->nb_points_on_face(f));
    for (size_type j=0; j < pai->nb_points_on_face(f); ++j)
      for (size_type i=0; i < imdim; ++i)
	w(i,j)=pai->point_on_face(f,j)[i];
  } else if (check_cmd(cmd, "coeffs", in, out, 0, 0, 0, 1)) {
    /*@GET INTEG:GET('coeffs')
      Returns the coefficients associated to each integration point.
      @*/
    check_not_exact(im);
    out.pop().from_dcvector(im->approx_method()->integration_coefficients());
  } else if (check_cmd(cmd, "face_coeffs", in, out, 1, 1, 0, 1)) {
    /*@GET INTEG:GET('face_coeffs',F)
      Returns the coefficients associated to each integration of a face.
      @*/
    check_not_exact(im); 
    size_type f = in.pop().to_face_number(pai->structure()->nb_faces());
    darray w = out.pop().create_darray_h(pai->nb_points_on_face(f));
    for (size_type j=0; j < pai->nb_points_on_face(f); ++j)
      w[j]=pai->coeff_on_face(f,j);
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET INTEG:GET('char')
      Ouput a (unique) string representation of the integration method. 

      This can be used to  comparisons between two different @tinteg objects.
      @*/    
    std::string s = getfem::name_of_int_method(im);    
    out.pop().from_string(s.c_str());
  } else bad_cmd(cmd);
}

