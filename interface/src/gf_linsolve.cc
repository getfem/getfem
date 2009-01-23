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

#include <getfemint_gsparse.h>
#include <getfemint_precond.h>
#include <getfemint.h>
#include <gmm/gmm_iter_solvers.h>
#include <getfemint_misc.h>
#include <getfem/getfem_superlu.h>

using namespace getfemint;

typedef enum { GMM_GMRES, GMM_CG, GMM_BICGSTAB /*, GMM_QMR*/ } iterative_gmm_solver_type;

template <typename T> static void
iterative_gmm_solver(iterative_gmm_solver_type stype, gsparse &gsp,
		    getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  garray<T> b = in.pop().to_garray(int(gsp.nrows()), T());
  garray<T> x = out.pop().create_array_v(int(gsp.nrows()), T());

  int restart = 50;
  if (in.remaining() && stype == GMM_GMRES) restart = in.pop().to_integer(1,1000000);

  gprecond<T> id_precond;
  gprecond<T> *precond = &id_precond;
  if (in.remaining()) precond = &in.pop().to_precond()->precond(T());
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
  dal::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
  gsparse &gsp = *pgsp;
  if (!gsp.is_complex() && in.front().is_complex())
    THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
  if (gsp.is_complex()) iterative_gmm_solver(stype, gsp, in, out, complex_type());
  else                  iterative_gmm_solver(stype, gsp, in, out, scalar_type());
}

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

/*MLABCOM
  FUNCTION F=gf_linsolve(args)

  @FUNC ::LINSOLVE('gmres')
  @FUNC ::LINSOLVE('cg')
  @FUNC ::LINSOLVE('bicgstab')
  @FUNC ::LINSOLVE('lu')
  @FUNC ::LINSOLVE('superlu')
MLABCOM*/

void gf_linsolve(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }

  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "gmres", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ::LINSOLVE('gmres',@tsp M, @vec b[, @int restart][, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the generalized minimum residuals method.

    Optionally using `P` as preconditioner. The default value of the
    restart parameter is 50.@*/
    iterative_gmm_solver(GMM_GMRES, in, out);
  } else if (check_cmd(cmd, "cg", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ::LINSOLVE('cg',@tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the conjugated gradient method.

    Optionally using `P` as preconditioner.@*/
    iterative_gmm_solver(GMM_CG, in, out);
  } else if (check_cmd(cmd, "bicgstab", in, out, 2, 30, 0, 1)) {
    /*@FUNC X = ::LINSOLVE('bicgstab',@tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the bi-conjugated gradient stabilized method.

    Optionally using `P` as a preconditioner.@*/
    iterative_gmm_solver(GMM_BICGSTAB, in, out);
    /*} else if (check_cmd(cmd, "qmr", in, out, 0, 1, 0, 1)) {
      iterative_gmm_solver(GMM_QMR, in, out);*/
  } else if (check_cmd(cmd, "lu", in, out, 2, 2, 0, 1) ||
	     check_cmd(cmd, "superlu", in, out, 2, 2, 0, 1)) {
    /*@FUNC @CELL{U, cond} = ::LINSOLVE('lu',@tsp M, @vec b)
    Alias for ::LINSOLVE('superlu',...)@*/
    /*@FUNC @CELL{U, cond} = ::LINSOLVE('superlu',@tsp M, @vec b)
    Solve `M.U = b` apply the SuperLU solver (sparse LU factorization).

    The condition number estimate `cond` is returned with the solution `U`.@*/
    dal::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
    gsparse &gsp = *pgsp;
    if (!gsp.is_complex() && in.front().is_complex())
      THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
    if (gsp.is_complex()) superlu_solver(gsp, in, out, complex_type());
    else                  superlu_solver(gsp, in, out, scalar_type());
  } else bad_cmd(cmd);
}
