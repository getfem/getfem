/*===========================================================================

 Copyright (C) 2025-2025 Konstantinos Poulios.

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
#include <getfemint_gmumps.h>
#include <getfemint_gsparse.h>

using namespace getfemint;

/*@GFDOC
  General function for modifying mumps_context objects
*/

void gf_mumps_context_set(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  gmumps *pctx         = to_mumps_context_object(in.pop());
  std::string init_cmd = in.pop().to_string();
  std::string cmd      = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "matrix", in, out, 1, 1, 0, 0)) {
    /*@SET ('matrix', @mat A)
      Set the matrix @mat A for the @tmct object.@*/
    std::shared_ptr<gsparse> pmat = in.pop().to_sparse();
    if (pctx->is_complex()) {
      if (!pmat->is_complex()) THROW_ERROR("Complex number matrix expected");
      if (pmat->storage() == gsparse::CSCMAT)
        pctx->set_matrix_c(pmat->cplx_csc());
      else if (pmat->storage() == gsparse::WSCMAT)
        pctx->set_matrix_c(pmat->cplx_wsc());
    } else {
      if (pmat->is_complex())  THROW_ERROR("Real number matrix expected");
      if (pmat->storage() == gsparse::CSCMAT)
        pctx->set_matrix_r(pmat->real_csc());
      else if (pmat->storage() == gsparse::WSCMAT)
        pctx->set_matrix_r(pmat->real_wsc());
    }
  } else if (check_cmd(cmd, "vector", in, out, 1, 1, 0, 0)) {
    /*@SET ('vector', @vec b)
      Set the right hand side @vec b for the @tmct object.@*/
    if (pctx->is_complex()) pctx->set_vector_c(in.pop().to_carray());
    else                    pctx->set_vector_r(in.pop().to_darray());
  } else if (check_cmd(cmd, "ICNTL", in, out, 2, 2, 0, 0)) {
    /*@SET ('ICNTL', @int ind, @int val)
      Set the integer option at 1-based index `ind` in the ICNTL vector
      of the @tmct object to value `val`.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 60)
      THROW_BADARG("Invalid ICNTL parameter index" << ind);
    pctx->ICNTL(ind) = in.pop().to_integer();
  } else if (check_cmd(cmd, "CNTL", in, out, 2, 2, 0, 0)) {
    /*@SET ('CNTL', @int ind, @scalar val)
      Set the scalar option at 1-based index `ind` in the CNTL vector
      of the @tmct object to value `val`.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 15)
      THROW_BADARG("Invalid CNTL parameter index" << ind);
    pctx->CNTL(ind) = in.pop().to_scalar();
  } else if (check_cmd(cmd, "analyze", in, out, 0, 0, 0, 0)) {
    /*@SET ('analyze')
      Run the MUMPS analysis job for the @tmct object.@*/
    pctx->analyze();
  } else if (check_cmd(cmd, "factorize", in, out, 0, 0, 0, 0)) {
    /*@SET ('factorize')
      Run the MUMPS factorization job for the @tmct object.@*/
    pctx->factorize();
  } else if (check_cmd(cmd, "solve", in, out, 0, 0, 0, 1)) {
    /*@SET SOL = ('solve')
      Run the MUMPS solve job (only) for the @tmct object.

      The analysis and factorization jobs need to be run first
      before calling this function.

      It returns the solution vector (on all processes if MPI is used).@*/
    pctx->solve();
    pctx->mpi_broadcast(); // make the solution available on all processes
    if (out.remaining()) {
      if (pctx->is_complex()) out.pop().from_dcvector(pctx->vector_c());
      else                    out.pop().from_dcvector(pctx->vector_r());
    }
  } else if (check_cmd(cmd, "analyze and factorize", in, out, 0, 0, 0, 0)) {
    /*@SET ('analyze and factorize')
      Run the MUMPS analysis and factorization jobs for the @tmct object.@*/
    pctx->analyze_and_factorize();
  } else if (check_cmd(cmd, "factorize and solve", in, out, 0, 0, 0, 1)) {
    /*@SET SOL = ('factorize and solve')
      Run the MUMPS factorization and solve jobs for the @tmct object.

      The analysis job needs to be run first before calling this function.

      It returns the solution vector (on all processes if MPI is used).@*/
    pctx->factorize_and_solve();
    pctx->mpi_broadcast(); // make the solution available on all processes
    if (out.remaining()) {
      if (pctx->is_complex()) out.pop().from_dcvector(pctx->vector_c());
      else                    out.pop().from_dcvector(pctx->vector_r());
    }
  } else if (check_cmd(cmd, "analyze factorize and solve", in, out, 0, 0, 0, 1)) {
    /*@SET SOL = ('analyze factorize and solve')
      Run the MUMPS analysis, factorization and solve jobs for the @tmct
      object.

      It returns the solution vector (on all processes if MPI is used).@*/
    pctx->analyze_factorize_and_solve();
    pctx->mpi_broadcast(); // make the solution available on all processes
    if (out.remaining()) {
      if (pctx->is_complex()) out.pop().from_dcvector(pctx->vector_c());
      else                    out.pop().from_dcvector(pctx->vector_r());
    }
   } else
     bad_cmd(init_cmd);
}
