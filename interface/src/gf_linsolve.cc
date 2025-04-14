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

#include <getfemint_gsparse.h>
#include <getfemint_gprecond.h>
#include <getfemint_misc.h>
#include <gmm/gmm_iter_solvers.h>
#include <gmm/gmm_superlu_interface.h>
#include <gmm/gmm_MUMPS_interface.h>

using namespace getfemint;

typedef enum { GMM_GMRES, GMM_CG, GMM_BICGSTAB /*, GMM_QMR*/ } iterative_gmm_solver_type;

template <typename T> static void
iterative_gmm_solver(iterative_gmm_solver_type stype, gsparse &gsp,
                    getfemint::mexargs_in& in,
                     getfemint::mexargs_out& out, T) {
  garray<T> b = in.pop().to_garray(int(gsp.nrows()), T());
  garray<T> x = out.pop().create_array_v(int(gsp.nrows()), T());

  int restart = 50;
  if (in.remaining() && stype == GMM_GMRES) restart = in.pop().to_integer(1,1000000);

  gprecond<T> id_precond;
  gprecond<T> *precond = &id_precond;
  if (in.remaining())
    precond = dynamic_cast<gprecond<T> *>(to_precond_object(in.pop()));
  precond->set_dimensions(gsp.nrows(), gsp.ncols());

  getfemint::interruptible_iteration iter(1e-16);

  while (in.remaining() && in.front().is_string()) {
    std::string opt = in.pop().to_string();
    if (cmd_strmatch(opt, "noisy")) {
      iter.set_noisy(1);
    } else if (cmd_strmatch(opt, "very noisy")) {
      iter.set_noisy(3);
    } else if (cmd_strmatch(opt, "res")) {
      if (!in.remaining()) THROW_BADARG("missing value after '" << opt << "'");
      iter.set_resmax(in.pop().to_scalar(0., 1e100));
    } else if (cmd_strmatch(opt, "maxiter")) {
      if (!in.remaining()) THROW_BADARG("missing value after '" << opt << "'");
      iter.set_maxiter(in.pop().to_integer(1, INT_MAX));
    }
  }

  if (in.remaining()) THROW_BADARG("too much arguments");

  gsp.to_csc();
  switch (stype) {
    case GMM_GMRES: gmm::gmres(gsp.csc(T()), x, b, *precond, restart, iter); break;
    case GMM_CG: gmm::cg(gsp.csc(T()), x, b, *precond, iter); break;
    case GMM_BICGSTAB: gmm::bicgstab(gsp.csc(T()), x, b, *precond, iter); break;
      //case GMM_QMR: gmm::qmr(gsp.csc(T()), x, b, *precond, iter); break;
  }
}

void iterative_gmm_solver(iterative_gmm_solver_type stype,
                         getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
  gsparse &gsp = *pgsp;
  if (!gsp.is_complex() && in.front().is_complex())
    THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
  if (gsp.is_complex()) iterative_gmm_solver(stype, gsp, in, out, complex_type());
  else                  iterative_gmm_solver(stype, gsp, in, out, scalar_type());
}

#if defined(GMM_USES_SUPERLU)
template <typename T> static void
superlu_solver(gsparse &gsp,
               getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  garray<T> b = in.pop().to_garray(int(gsp.nrows()), T());
  garray<T> x = out.pop().create_array(b.getm(), b.getn(), T());
  double rcond;
  gsp.to_csc();
  gmm::SuperLU_solve(gsp.csc(T()),x,b,rcond,1);
  if (out.remaining())
    out.pop().from_scalar(rcond ? 1./rcond : 0.);
}
#endif

#if defined(GMM_USES_MUMPS)
template <typename T> static void
mumps_solver(gsparse &gsp,
             getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  int nrows = int(gsp.nrows());
  garray<T> b = in.pop().to_garray(T());
  int nrhs = 1; // just 1 rhs by default
  if (b.getn() != 1 && b.getm() != 1) { // multiple rhs
    in.last_popped().check_dimensions(b, nrows, -1);
    nrhs = b.getn();
  } else // check whether b is in fact a vector of proper size
    in.last_popped().check_dimensions(b, nrows);
  bool sym = false;
  if (in.remaining() && in.front().is_string()) {
    std::string opt = in.pop().to_string();
    if (cmd_strmatch(opt, "sym"))
      sym = true;
    else
      THROW_BADARG("unknown linsolve option: " << opt);
  }
  garray<T> x = out.pop().create_array(b.getm(), b.getn(), T());
  gsp.to_csc();
# if GETFEM_PARA_LEVEL > 1
  double t_ref = MPI_Wtime();
  gmm::MUMPS_distributed_matrix_solve(gsp.csc(T()), x, b, sym);
  if (getfem::MPI_IS_MASTER()) cout << "MUMPS solve time " << MPI_Wtime()-t_ref << endl;
# else
  gmm::MUMPS_solve(gsp.csc(T()), x, b, sym);
# endif
}
#endif


/*@GFDOC
  Various linear system solvers.
@*/


void gf_linsolve(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {

  if (in.narg() < 1) THROW_BADARG("Wrong number of input arguments");

  std::string init_cmd = in.pop().to_string();
  std::string cmd      = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "gmres", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ('gmres', @tsp M, @vec b[, @int restart][, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
      Solve `M.X = b` with the generalized minimum residuals method.

      Optionally using `P` as preconditioner. The default value of the
      restart parameter is 50.@*/
    iterative_gmm_solver(GMM_GMRES, in, out);
  } else if (check_cmd(cmd, "cg", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ('cg', @tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
      Solve `M.X = b` with the conjugated gradient method.

      Optionally using `P` as preconditioner.@*/
    iterative_gmm_solver(GMM_CG, in, out);
  } else if (check_cmd(cmd, "bicgstab", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ('bicgstab', @tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
      Solve `M.X = b` with the bi-conjugated gradient stabilized method.

      Optionally using `P` as a preconditioner.@*/
    iterative_gmm_solver(GMM_BICGSTAB, in, out);
#if defined(GMM_USES_SUPERLU)
  } else if (check_cmd(cmd, "lu", in, out, 2, 2, 0, 2)) {
    /*@FUNC @CELL{U, cond} = ('lu', @tsp M, @vec b)
      Alias for ::LINSOLVE('superlu',...)@*/
    std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
    gsparse &gsp = *pgsp;
    if (!gsp.is_complex() && in.front().is_complex())
      THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
    if (gsp.is_complex()) superlu_solver(gsp, in, out, complex_type());
    else                  superlu_solver(gsp, in, out, scalar_type());
  } else if (check_cmd(cmd, "superlu", in, out, 2, 2, 0, 2)) {
    /*@FUNC @CELL{U, cond} = ('superlu', @tsp M, @vec b)
      Solve `M.U = b` apply the SuperLU solver (sparse LU factorization).

      The condition number estimate `cond` is returned with the solution `U`.@*/
    std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
    gsparse &gsp = *pgsp;
    if (!gsp.is_complex() && in.front().is_complex())
      THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
    if (gsp.is_complex()) superlu_solver(gsp, in, out, complex_type());
    else                  superlu_solver(gsp, in, out, scalar_type());
#endif

#if defined(GMM_USES_MUMPS)
  } else if (check_cmd(cmd, "mumps", in, out, 2, 3, 0, 1)) {
    /*@FUNC @CELL{U, cond} = ('mumps', @tsp M, @vec b, ... ['sym'])
      Solve `M.U = b` using the MUMPS solver.

      The right hand side `b` can optionally by a matrix with several columns
      in order to solve multiple right hand sides at once.

      If the option `sym` is provided, the symmetric version of the MUMPS
      solver will be used.@*/
    std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
    gsparse &gsp = *pgsp;
    if (!gsp.is_complex() && in.front().is_complex())
      THROW_BADARG("please use a real right hand side, or convert "
                   "the sparse matrix into a complex one");
    if (gsp.is_complex()) mumps_solver(gsp, in, out, complex_type());
    else                  mumps_solver(gsp, in, out, scalar_type());
#endif
  } else
    bad_cmd(init_cmd);

}
