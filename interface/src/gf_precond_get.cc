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

#include <getfemint_precond.h>

using namespace getfemint;

template <typename T> static void 
mult_or_tmult(gprecond<T>& precond, mexargs_in& in, mexargs_out& out, bool tmult) {
  garray<T> v = in.pop().to_garray(T());
  garray<T> w = out.pop().create_array(v.getm(), v.getn(), T());
  gmm::mult_or_transposed_mult(precond, v, w, tmult);
}


/*MLABCOM
  FUNCTION F=gf_precond_get(...)

  @GET PRECOND:GET('mult')
  @GET PRECOND:GET('tmult')
  @GET PRECOND:GET('type')
  @GET PRECOND:GET('size')
  @GET PRECOND:GET('is_complex')
  @GET PRECOND:GET('info')
MLABCOM*/

void gf_precond_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_precond *precond = in.pop().to_precond();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "mult", in, out, 1, 1, 0, 1)) {
    /*@GET PRECOND:GET('mult', @vec V)
      Apply the preconditioner to the supplied vector.
      @*/
    if (!precond->is_complex()) mult_or_tmult(precond->precond(scalar_type()), in, out, false);
    else                        mult_or_tmult(precond->precond(complex_type()), in, out, false);
  } else if (check_cmd(cmd, "tmult", in, out, 1, 1, 0, 1)) {
    /*@GET PRECOND:GET('tmult', @vec V)
      Apply the transposed preconditioner to the supplied vector.
      @*/
    if (!precond->is_complex()) mult_or_tmult(precond->precond(scalar_type()), in, out, true);
    else                        mult_or_tmult(precond->precond(complex_type()), in, out, true);
  } else if (check_cmd(cmd, "type", in, out, 0, 0, 0, 1)) {
    /*@GET PRECOND:GET('type')
      Return a string describing the type of the preconditioner ('ilu', 'ildlt',..).
      @*/
    out.pop().from_string(precond->bprecond().name());
  } else if (check_cmd(cmd, "size", in, out, 0, 0, 0, 1)) {
    /*@GET PRECOND:GET('size')
      Return the dimensions of the preconditioner.
      @*/
    iarray sz = out.pop().create_iarray_h(2);
    sz[0] = precond->bprecond().nrows();
    sz[1] = precond->bprecond().ncols();
  } else if (check_cmd(cmd, "is_complex", in, out, 0, 0, 0, 1)) {
    /*@GET PRECOND:GET('is_complex')
      Return 1 if the preconditioner stores complex values.
      @*/
    out.pop().from_integer(precond->is_complex());
  } else if (check_cmd(cmd, "info", in, out, 0, 1)) {
    /*@GET PRECOND:GET('info')
      Return a short informative string about the preconditioner.
      @*/
    std::stringstream ss;
    ss << precond->bprecond().nrows() << "x" << precond->bprecond().ncols() << " " 
       << (precond->is_complex() ? "COMPLEX" : "REAL") 
       << " " << precond->bprecond().name() << " [" << precond->memsize() << " bytes]";
    out.pop().from_string(ss.str().c_str());
  } else bad_cmd(cmd);
}
