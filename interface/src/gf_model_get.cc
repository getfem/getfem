// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard.
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

/**\file gf_model_get.cc
   \brief getfemint_model getter.
*/

#include <getfemint.h>
#include <getfemint_models.h>

using namespace getfemint;


#define RETURN_SPARSE(realmeth, cplxmeth)				\
  if (!md->is_complex()) {						\
    gf_real_sparse_by_col M(gmm::mat_nrows(md->model().realmeth()),	\
			    gmm::mat_ncols(md->model().realmeth()));	\
    gmm::copy(md->model().realmeth(), M);				\
    out.pop().from_sparse(M);						\
  } else {								\
    gf_cplx_sparse_by_col M(gmm::mat_nrows(md->model().cplxmeth()),	\
			    gmm::mat_ncols(md->model().cplxmeth()));	\
    gmm::copy(md->model().cplxmeth(), M);				\
    out.pop().from_sparse(M);						\
  }

#define RETURN_VECTOR(realmeth, cplxmeth)			\
  if (!md->is_complex()) {					\
    out.pop().from_dcvector(md->model().realmeth);		\
  } else {							\
    out.pop().from_dcvector(md->model().cplxmeth);		\
  }

/*MLABCOM

  FUNCTION M = gf_model_get(cmd, [, args])
  Get information from a model object.

  @RDATTR MODEL:GET('is_complex')
  @GET MODEL:GET('tangent_matrix')
  @GET MODEL:GET('state')
  @GET MODEL:GET('rhs')
  @GET MODEL:GET('memsize')
  @GET MODEL:GET('varlist')
  @GET MODEL:GET('variable')

MLABCOM*/

void gf_model_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {

  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_model *md  = in.pop().to_getfemint_model();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "is_complex", in, out, 0, 0, 0, 1)) {
    /*@RDATTR b = MODEL:GET('is_complex')
    Return 0 is the model is real, 1 if it is complex.@*/
    out.pop().from_integer(md->is_complex());
  } else if (check_cmd(cmd, "tangent_matrix", in, out, 0, 0, 0, 1)) {
    /*@GET T = MODEL:GET('tangent_matrix')
    Return the tangent matrix stored in the model .@*/
    RETURN_SPARSE(real_tangent_matrix, complex_tangent_matrix);
  } else if (check_cmd(cmd, "state", in, out, 0, 0, 0, 1)) {
    /*@GET MODEL:GET('state')
    Return the vector of unknowns, which contains the solution after a solve.@*/
    RETURN_VECTOR(real_state(), complex_state());
  } else if (check_cmd(cmd, "rhs", in, out, 0, 0, 0, 1)) {
    /*@GET MODEL:GET('rhs')
    Return the right hand side of the tangent problem.@*/
    RETURN_VECTOR(real_rhs(), complex_rhs());
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET z = MODEL:GET('memsize')
    Return a rough approximation of the amount of memory (in bytes) used by
    the model.@*/
    out.pop().from_integer(int(md->memsize()));
  } else if (check_cmd(cmd, "varlist", in, out, 0, 0, 0, 0)) {
    /*@GET MODEL:GET('varlist')
    print to the output the list of variables and constants of the model.@*/
    md->model().varlist(infomsg());
  } else if (check_cmd(cmd, "variable", in, out, 1, 2, 0, 0)) {
    /*@GET V = MODEL:GET('variable', @str name[, @int niter])
    Gives the value of a variable or data.@*/
    std::string name = in.pop().to_string();
    size_type niter = 0;
    if (in.remaining())
      niter = in.pop().to_integer(0,10) - config::base_index();
    RETURN_VECTOR(real_variable(name, niter), complex_variable(name, niter));
  } else bad_cmd(cmd);
}
