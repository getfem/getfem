/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfemint_gsparse.h>
#include <getfemint_precond.h>
#include <getfemint.h>
#include <gmm/gmm_iter_solvers.h>
#include <getfemint_misc.h>
#include <getfem/getfem_superlu.h>
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

#if defined(GMM_USES_MUMPS) || defined(HAVE_DMUMPS_C_H)
template <typename T> static void
mumps_solver(gsparse &gsp,
             getfemint::mexargs_in& in, getfemint::mexargs_out& out, T) {
  garray<T> b = in.pop().to_garray(int(gsp.nrows()), T());
  garray<T> x = out.pop().create_array(b.getm(), b.getn(), T());
  gsp.to_csc();
  gmm::MUMPS_solve(gsp.csc(T()),x,b);
}
#endif

/*@GFDOC
  Various linear system solvers.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_linsolve : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out) = 0;
};

typedef std::shared_ptr<sub_gf_linsolve> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_linsolve {				\
      virtual void run(getfemint::mexargs_in& in,			\
                       getfemint::mexargs_out& out)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                        



void gf_linsolve(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@FUNC X = ('gmres', @tsp M, @vec b[, @int restart][, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the generalized minimum residuals method.

    Optionally using `P` as preconditioner. The default value of the
    restart parameter is 50.@*/
    sub_command
      ("gmres", 2, 30, 0, 1,
       iterative_gmm_solver(GMM_GMRES, in, out);
       );


    /*@FUNC X = ('cg', @tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the conjugated gradient method.

    Optionally using `P` as preconditioner.@*/
    sub_command
      ("cg", 2, 30, 0, 1,
       iterative_gmm_solver(GMM_CG, in, out);
       );


    /*@FUNC X = ('bicgstab', @tsp M, @vec b [, @tpre P][,'noisy'][,'res', r][,'maxiter', n])
    Solve `M.X = b` with the bi-conjugated gradient stabilized method.

    Optionally using `P` as a preconditioner.@*/
    sub_command
      ("bicgstab", 2, 30, 0, 1,
       iterative_gmm_solver(GMM_BICGSTAB, in, out);
       );


    /*@FUNC @CELL{U, cond} = ('lu', @tsp M, @vec b)
      Alias for ::LINSOLVE('superlu',...)@*/
    sub_command
      ("lu", 2, 2, 0, 2,
       std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
       gsparse &gsp = *pgsp;
       if (!gsp.is_complex() && in.front().is_complex())
         THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
       if (gsp.is_complex()) superlu_solver(gsp, in, out, complex_type());
       else                  superlu_solver(gsp, in, out, scalar_type());
       );


    /*@FUNC @CELL{U, cond} = ('superlu', @tsp M, @vec b)
    Solve `M.U = b` apply the SuperLU solver (sparse LU factorization).

    The condition number estimate `cond` is returned with the solution `U`.@*/
    sub_command
      ("superlu", 2, 2, 0, 2,
       std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
       gsparse &gsp = *pgsp;
       if (!gsp.is_complex() && in.front().is_complex())
         THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
       if (gsp.is_complex()) superlu_solver(gsp, in, out, complex_type());
       else                  superlu_solver(gsp, in, out, scalar_type());
       );

#if defined(GMM_USES_MUMPS) || defined(HAVE_DMUMPS_C_H)
    /*@FUNC @CELL{U, cond} = ('mumps', @tsp M, @vec b)
    Solve `M.U = b` using the MUMPS solver.@*/
    sub_command
      ("mumps", 2, 2, 0, 1,
       std::shared_ptr<gsparse> pgsp = in.pop().to_sparse();
       gsparse &gsp = *pgsp;
       if (!gsp.is_complex() && in.front().is_complex())
         THROW_BADARG("please use a real right hand side, or convert the sparse matrix to a complex one");
       if (gsp.is_complex()) mumps_solver(gsp, in, out, complex_type());
       else                  mumps_solver(gsp, in, out, scalar_type());
       );
#endif

  }



  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");

  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out);
  }
  else bad_cmd(init_cmd);

}
