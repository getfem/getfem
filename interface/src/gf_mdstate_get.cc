/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Y. Renard, J. Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**\file gf_mdstate_get.cc
   \brief getfemint_mdstate getter.
*/

#include <getfemint.h>
#include <getfemint_mdstate.h>

using namespace getfemint;


#define RETURN_SPARSE(WHAT)						\
  if (!md->is_complex()) {						\
    gf_real_sparse_by_col M(gmm::mat_nrows(md->real_mdstate().WHAT()),	\
			    gmm::mat_ncols(md->real_mdstate().WHAT())); \
    gmm::copy(md->real_mdstate().WHAT(), M);				\
    out.pop().from_sparse(M);						\
  } else {								\
    gf_cplx_sparse_by_col M(gmm::mat_nrows(md->cplx_mdstate().WHAT()),	\
			    gmm::mat_ncols(md->cplx_mdstate().WHAT())); \
    gmm::copy(md->cplx_mdstate().WHAT(), M);				\
    out.pop().from_sparse(M);						\
  }

#define RETURN_VECTOR(WHAT)				\
  if (!md->is_complex()) {				\
    out.pop().from_dcvector(md->real_mdstate().WHAT());	\
  } else {						\
    out.pop().from_dcvector(md->cplx_mdstate().WHAT());	\
  }

/*@GFDOC
  Get information from a model state object.
@*/

void gf_mdstate_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mdstate *md  = in.pop().to_getfemint_mdstate();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "is_complex", in, out, 0, 0, 0, 1)) {
    /*@RDATTR b = ('is_complex')
    Return 0 is the model state is real, 1 if it is complex.@*/
    out.pop().from_integer(md->is_complex());
  } else if (check_cmd(cmd, "tangent_matrix", in, out, 0, 0, 0, 1)) {
    /*@GET T = ('tangent_matrix')
    Return the tangent matrix stored in the model state.@*/
    RETURN_SPARSE(tangent_matrix);
  } else if (check_cmd(cmd, "constraints_matrix", in, out, 0, 0, 0, 1)) {
    /*@GET C = ('constraints_matrix')
    Return the constraints matrix stored in the model state.@*/
    RETURN_SPARSE(constraints_matrix);
  } else if (check_cmd(cmd, "reduced_tangent_matrix", in, out, 0, 0, 0, 1)) {
    /*@GET A = ('reduced_tangent_matrix')
    Return the reduced tangent matrix (i.e. the tangent matrix after
    elimination of the constraints).@*/
    RETURN_SPARSE(reduced_tangent_matrix);
  } else if (check_cmd(cmd, "constraints_nullspace", in, out, 0, 0, 0, 1)) {
    /*@GET ('constraints_nullspace')
    Return the nullspace of the constraints matrix.@*/
    RETURN_SPARSE(constraints_nullspace);
  } else if (check_cmd(cmd, "state", in, out, 0, 0, 0, 1)) {
    /*@GET ('state')
    Return the vector of unknowns, which contains the solution after MDBRICK:GET('solve').@*/
    RETURN_VECTOR(state);
  } else if (check_cmd(cmd, "residual", in, out, 0, 0, 0, 1)) {
    /*@GET ('residual')
    Return the residual.@*/
    RETURN_VECTOR(residual);
  } else if (check_cmd(cmd, "reduced_residual", in, out, 0, 0, 0, 1)) {
    /*@GET ('reduced_residual')
    Return the residual on the reduced system.@*/
    RETURN_VECTOR(reduced_residual);
  } else if (check_cmd(cmd, "unreduce", in, out, 1, 1, 0, 1)) {
    /*@GET ('unreduce', @vec U)
    Reinsert the constraint eliminated from the system.@*/
    if (!md->is_complex()) {
      size_type nred = gmm::vect_size(md->real_mdstate().reduced_residual());
      size_type n    = gmm::vect_size(md->real_mdstate().residual());
      darray Ured = in.pop().to_darray(int(nred));
      darray U    = out.pop().create_darray_v(unsigned(n));
      md->real_mdstate().unreduced_solution(Ured, U);
    } else {
      size_type nred = gmm::vect_size(md->cplx_mdstate().reduced_residual());
      size_type n    = gmm::vect_size(md->cplx_mdstate().residual());
      carray Ured = in.pop().to_carray(int(nred));
      carray U    = out.pop().create_carray_v(unsigned(n));
      md->cplx_mdstate().unreduced_solution(Ured, U);
    }
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET z = ('memsize')
    Return the amount of memory (in bytes) used by the model state.@*/
    out.pop().from_integer(int(md->memsize()));
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET s = ('char')
      Output a (unique) string representation of the @tmdstate.

      This can be used to perform comparisons between two
      different @tmdstate objects.
      This function is to be completed.
      @*/
    GMM_ASSERT1(false, "Sorry, function to be done");
    // std::string s = ...;
    // out.pop().from_string(s.c_str());
  } else if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      displays a short summary for a @tmdstate.@*/
    infomsg() << "gfMdState object\n";
  } else bad_cmd(cmd);
}
