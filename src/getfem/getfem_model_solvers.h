/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2015 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**
   @file getfem_model_solvers.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 15, 2004.
   @brief Standard solvers for model bricks
   @see getfem_modeling.h
*/

#ifndef GETFEM_MODEL_SOLVERS_H__
#define GETFEM_MODEL_SOLVERS_H__
#include "getfem_models.h"
#include "gmm/gmm_MUMPS_interface.h"
#include "gmm/gmm_iter.h"
#include "gmm/gmm_iter_solvers.h"
#include "gmm/gmm_dense_qr.h"

//#include "gmm/gmm_inoutput.h"

// metis necessary in this header only due to the old bricks system
#if defined GETFEM_HAVE_METIS && !defined GETFEM_HAVE_METIS_OLD_API
#  include <metis.h>
#endif

namespace getfem {

  template <typename T> struct sort_abs_val_
  { bool operator()(T x, T y) { return (gmm::abs(x) < gmm::abs(y)); } };

  template <typename MAT> void print_eigval(const MAT &M) {
    // just to test a stiffness matrix if needing
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    size_type n = gmm::mat_nrows(M);
    gmm::dense_matrix<T> MM(n, n), Q(n, n);
    std::vector<T> eigval(n);
    gmm::copy(M, MM);
    // gmm::symmetric_qr_algorithm(MM, eigval, Q);
    gmm::implicit_qr_algorithm(MM, eigval, Q);
    std::sort(eigval.begin(), eigval.end(), sort_abs_val_<T>());
    cout << "eival = " << eigval << endl;
//     cout << "vectp : " << gmm::mat_col(Q, n-1) << endl;
//     cout << "vectp : " << gmm::mat_col(Q, n-2) << endl;
//     double emax, emin;
//     cout << "condition number" << condition_number(MM,emax,emin) << endl;
//     cout << "emin = " << emin << endl;
//     cout << "emax = " << emax << endl;
  }


  /* ***************************************************************** */
  /*     Linear solvers definition                                     */
  /* ***************************************************************** */

  template <typename MAT, typename VECT>
  struct abstract_linear_solver {
    virtual void operator ()(const MAT &, VECT &, const VECT &,
                             gmm::iteration &) const  = 0;
    virtual ~abstract_linear_solver() {}
  };

  template <typename MAT, typename VECT>
  struct linear_solver_cg_preconditioned_ildlt
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      gmm::ildlt_precond<MAT> P(M);
      gmm::cg(M, x, b, P, iter);
      if (!iter.converged()) GMM_WARNING2("cg did not converge!");
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_gmres_preconditioned_ilu
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      gmm::ilu_precond<MAT> P(M);
      gmm::gmres(M, x, b, P, 500, iter);
      if (!iter.converged()) GMM_WARNING2("gmres did not converge!");
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_gmres_unpreconditioned
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      gmm::identity_matrix P;
      gmm::gmres(M, x, b, P, 500, iter);
      if (!iter.converged()) GMM_WARNING2("gmres did not converge!");
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_gmres_preconditioned_ilut
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      gmm::ilut_precond<MAT> P(M, 40, 1E-7);
      gmm::gmres(M, x, b, P, 500, iter);
      if (!iter.converged()) GMM_WARNING2("gmres did not converge!");
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_gmres_preconditioned_ilutp
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      gmm::ilutp_precond<MAT> P(M, 20, 1E-7);
      gmm::gmres(M, x, b, P, 500, iter);
      if (!iter.converged()) GMM_WARNING2("gmres did not converge!");
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_superlu
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter)  const {
      double rcond;
      /*gmm::HarwellBoeing_IO::write("test.hb", M);
      std::fstream f("bbb", std::ios::out);
      for (unsigned i=0; i < gmm::vect_size(b); ++i) f << b[i] << "\n";*/
      int info = SuperLU_solve(M, x, b, rcond);
      iter.enforce_converged(info == 0);
      if (iter.get_noisy()) cout << "condition number: " << 1.0/rcond<< endl;
    }
  };


#ifdef GMM_USES_MUMPS
  template <typename MAT, typename VECT>
  struct linear_solver_mumps : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      bool ok = gmm::MUMPS_solve(M, x, b, false);
      iter.enforce_converged(ok);
    }
  };
  template <typename MAT, typename VECT>
  struct linear_solver_mumps_sym : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      bool ok = gmm::MUMPS_solve(M, x, b, true);
      iter.enforce_converged(ok);
    }
  };
#endif

#if GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
  template <typename MAT, typename VECT>
  struct linear_solver_distributed_mumps
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      double tt_ref=MPI_Wtime();
      bool ok = MUMPS_distributed_matrix_solve(M, x, b, false);
      iter.enforce_converged(ok);
      if (MPI_IS_MASTER()) cout<<"UNSYMMETRIC MUMPS time "<< MPI_Wtime() - tt_ref<<endl;
    }
  };

  template <typename MAT, typename VECT>
  struct linear_solver_distributed_mumps_sym
    : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      double tt_ref=MPI_Wtime();
      bool ok = MUMPS_distributed_matrix_solve(M, x, b, true);
      iter.enforce_converged(ok);
      if (MPI_IS_MASTER()) cout<<"SYMMETRIC MUMPS time "<< MPI_Wtime() - tt_ref<<endl;
    }
  };
#endif


  /* ***************************************************************** */
  /*     Newton Line search definition                                 */
  /* ***************************************************************** */

  struct abstract_newton_line_search {
    double conv_alpha, conv_r;
    size_t it, itmax, glob_it;
    //  size_t tot_it;
    virtual void init_search(double r, size_t git, double R0 = 0.0) = 0;
    virtual double next_try(void) = 0;
    virtual bool is_converged(double, double R1 = 0.0) = 0;
    virtual double converged_value(void) {
      // tot_it += it; cout << "tot_it = " << tot_it << endl; it = 0;
      return conv_alpha;
    };
    virtual double converged_residual(void) { return conv_r; };
    virtual ~abstract_newton_line_search() { }
  };


  struct quadratic_newton_line_search : public abstract_newton_line_search {
    double R0_, R1_;

    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    virtual void init_search(double r, size_t git, double R0 = 0.0) {
      GMM_ASSERT1(R0 != 0.0, "You have to specify R0");
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
      R0_ = R0;
    }
    virtual double next_try(void) {
      ++it;
      if (it == 1) return double(1);
      GMM_ASSERT1(R1_ != 0.0, "You have to specify R1");
      double a = R0_ / R1_;
      return (a < 0) ? (a*0.5 + sqrt(a*a*0.25-a)) : a*0.5;
    }
    virtual bool is_converged(double r, double R1 = 0.0) {
      conv_r = r;
      R1_ = R1; return ((gmm::abs(R1_) < gmm::abs(R0_*0.5)) || it >= itmax);
    }
    quadratic_newton_line_search(size_t imax = size_t(-1))  { itmax = imax; }
  };


  struct simplest_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min,
           alpha_threshold_res;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r, double = 0.0) {
      conv_r = r;
      return ((it <= 1 && r < first_res)
              || (r <= first_res * alpha_max_ratio && r <= alpha_threshold_res)
              || (conv_alpha <= alpha_min)
              || it >= itmax);
    }
    simplest_newton_line_search
    (size_t imax = size_t(-1), double a_max_ratio = 6.0/5.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0,
     double a_threshold_res = 1e50)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio), alpha_min(a_min),
        alpha_threshold_res(a_threshold_res)
      { itmax = imax; }
  };

  struct default_newton_line_search : public abstract_newton_line_search {
    // This line search try to detect where is the minimum, dividing the step
    // by a factor two each time.
    //    - it stops if the residual is less than the previous residual
    //      times alpha_min_ratio (= 0.9).
    //    - if the minimal step is reached with a residual greater than
    //      the previous residual times alpha_min_ratio then it decides
    //      between two options :
    //        - return with the step corresponding to the smallest residual
    //        - return with a greater residual corresponding to a residual
    //          less than the previous residual times alpha_max_ratio.
    //      the decision is taken regarding the previous iterations.
    //    - in order to shorten the line search, the process stops when
    //      the residual increases three times consecutively.
    // possible improvment : detect the curvature at the origin
    // (only one evaluation) and take it into account.
    // Fitted to some experiments in contrib/tests_newton

    double alpha, alpha_old, alpha_mult, first_res, alpha_max_ratio;
    double alpha_min_ratio, alpha_min;
    size_type count, count_pat;
    bool max_ratio_reached;
    double alpha_max_ratio_reached, r_max_ratio_reached;
    size_type it_max_ratio_reached;

    virtual void init_search(double r, size_t git, double = 0.0);
    virtual double next_try(void);
    virtual bool is_converged(double r, double = 0.0);
    default_newton_line_search(void) { count_pat = 0; }
  };

  /* the former default_newton_line_search */
  struct basic_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res, alpha_max_ratio;
    double alpha_min, prev_res, alpha_max_augment;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r, double = 0.0) {
      if (glob_it == 0 || (r < first_res / double(2))
          || (conv_alpha <= alpha_min && r < first_res * alpha_max_augment)
          || it >= itmax)
        { conv_r = r; return true; }
      if (it > 1 && r > prev_res && prev_res < alpha_max_ratio * first_res)
        return true;
      prev_res = conv_r = r;
      return false;
    }
    basic_newton_line_search
    (size_t imax = size_t(-1),
     double a_max_ratio = 5.0/3.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0, double a_augm = 2.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio),
        alpha_min(a_min), alpha_max_augment(a_augm) { itmax = imax; }
  };


  struct systematic_newton_line_search : public abstract_newton_line_search {
    double alpha, alpha_mult, first_res;
    double alpha_min, prev_res;
    bool first;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1);
      prev_res = conv_r = first_res = r; it = 0; first = true;
    }
    virtual double next_try(void)
    { double a = alpha; alpha *= alpha_mult; ++it; return a; }
    virtual bool is_converged(double r, double = 0.0) {
      // cout << "a = " << alpha / alpha_mult << " r = " << r << endl;
      if (r < conv_r || first)
        { conv_r = r; conv_alpha = alpha / alpha_mult; first = false; }
      if ((alpha <= alpha_min*alpha_mult) || it >= itmax) return true;
      return false;
    }
    systematic_newton_line_search
    (size_t imax = size_t(-1),
     double a_min = 1.0/10000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_min(a_min)  { itmax = imax; }
  };





  /* ***************************************************************** */
  /*     Newton algorithm.                                             */
  /* ***************************************************************** */

  template <typename PB>
  void classical_Newton(PB &pb, gmm::iteration &iter,
                        const abstract_linear_solver<typename PB::MATRIX,
                        typename PB::VECTOR> &linear_solver) {
    typedef typename gmm::linalg_traits<typename PB::VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.reduce_noisy();
    iter_linsolv0.set_resmax(iter.get_resmax()/20.0);
    iter_linsolv0.set_maxiter(10000); // arbitrary

    pb.compute_residual();

    typename PB::VECTOR dr(gmm::vect_size(pb.residual()));
    typename PB::VECTOR b(gmm::vect_size(pb.residual()));

    scalar_type crit = pb.residual_norm()
      / std::max(1E-25, pb.approx_external_load_norm());
    while (!iter.finished(crit)) {
      gmm::iteration iter_linsolv = iter_linsolv0;
      if (iter.get_noisy() > 1)
        cout << "starting computing tangent matrix" << endl;

      int is_singular = 1;
      while (is_singular) {
        pb.compute_tangent_matrix();
        gmm::clear(dr);
        gmm::copy(gmm::scaled(pb.residual(), pb.scale_residual()), b);
        if (iter.get_noisy() > 1) cout << "starting linear solver" << endl;
        iter_linsolv.init();
        linear_solver(pb.tangent_matrix(), dr, b, iter_linsolv);
        if (!iter_linsolv.converged()) {
          is_singular++;
          if (is_singular <= 4) {
            if (iter.get_noisy())
              cout << "Singular tangent matrix:"
                " perturbation of the state vector." << endl;
            pb.perturbation();
            pb.compute_residual();
          } else {
            if (iter.get_noisy())
              cout << "Singular tangent matrix: perturbation failed, aborting."
                   << endl;
            return;
          }
        }
        else is_singular = 0;
      }

      if (iter.get_noisy() > 1) cout << "linear solver done" << endl;
      R alpha = pb.line_search(dr, iter); //it is assumed that the linesearch
                                          //executes a pb.compute_residual();
      if (iter.get_noisy()) cout << "alpha = " << std::setw(6) << alpha << " ";
      ++iter;
      crit = std::min(pb.residual_norm()
                      / std::max(1E-25, pb.approx_external_load_norm()),
                      gmm::vect_norm1(dr) / std::max(1E-25, pb.state_norm()));
    }
  }


  /* ***************************************************************** */
  /*  Intermediary structure for Newton algorithms with getfem::model. */
  /* ***************************************************************** */

  #define TRACE_SOL 0

  template <typename MAT, typename VEC>
  struct model_pb {

    typedef MAT MATRIX;
    typedef VEC VECTOR;
    typedef typename gmm::linalg_traits<VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;

    model &md;
    abstract_newton_line_search &ls;
    VECTOR stateinit, &state;
    const VECTOR &rhs;
    const MATRIX &K;

    void compute_tangent_matrix(void) {
      md.to_variables(state);
      md.assembly(model::BUILD_MATRIX);
    }

    const MATRIX &tangent_matrix(void) { return K; }

    inline T scale_residual(void) const { return T(1); }

    void compute_residual(void) {
      md.to_variables(state);
      md.assembly(model::BUILD_RHS);
    }

    void perturbation(void) {
      R res = gmm::vect_norm2(state), ampl = std::max(res*R(1E-20), R(1E-50));
      std::vector<R> V(gmm::vect_size(state));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), state);
    }

    const VECTOR &residual(void) { return rhs; }

    R state_norm(void) const
    { return gmm::vect_norm1(state); }

    R approx_external_load_norm(void)
    { return md.approx_external_load(); }

    R residual_norm(void) { // A norm1 seems to be better than a norm2
                            // at least for contact problems.
      return gmm::vect_norm1(rhs);
    }

    R compute_res(bool comp = true) {
      if (comp) compute_residual();
      return residual_norm();
    }


    R line_search(VECTOR &dr, const gmm::iteration &iter) {
      size_type nit = 0;
      gmm::resize(stateinit, md.nb_dof());
      gmm::copy(state, stateinit);
      R alpha(1), res, /* res_init, */ R0;

      /* res_init = */ res = compute_res(false);
      // cout << "first residual = " << residual() << endl << endl;
      R0 = gmm::real(gmm::vect_sp(dr, rhs));

#if TRACE_SOL
      static int trace_number = 0;
      int trace_iter = 0;
      {
        std::stringstream trace_name;
        trace_name << "line_search_state" << std::setfill('0')
                   << std::setw(3) << trace_number << "_000_init";
        gmm::vecsave(trace_name.str(),stateinit);
      }
      trace_number++;
#endif

      ls.init_search(res, iter.get_iteration(), R0);
      do {
        alpha = ls.next_try();
        gmm::add(stateinit, gmm::scaled(dr, alpha), state);
#if TRACE_SOL
        {
          trace_iter++;
          std::stringstream trace_name;
          trace_name  << "line_search_state" << std::setfill('0')
                      << std::setw(3) << trace_number << "_"
                      << std::setfill('0') << std::setw(3) << trace_iter;
          gmm::vecsave(trace_name.str(), state);
        }
#endif
        res = compute_res();
        // cout << "residual = " << residual() << endl << endl;
        R0 = gmm::real(gmm::vect_sp(dr, rhs));

        ++ nit;
      } while (!ls.is_converged(res, R0));

      if (alpha != ls.converged_value()) {
        alpha = ls.converged_value();
        gmm::add(stateinit, gmm::scaled(dr, alpha), state);
        res = ls.converged_residual();
        compute_residual();
      }

      return alpha;
    }

    model_pb(model &m, abstract_newton_line_search &ls_, VECTOR &st,
             const VECTOR &rhs_, const MATRIX &K_)
      : md(m), ls(ls_), state(st), rhs(rhs_), K(K_) {}

  };

  //---------------------------------------------------------------------
  // Default linear solver.
  //---------------------------------------------------------------------

  typedef abstract_linear_solver<model_real_sparse_matrix,
                                 model_real_plain_vector> rmodel_linear_solver;
  typedef std::shared_ptr<rmodel_linear_solver> rmodel_plsolver_type;
  typedef abstract_linear_solver<model_complex_sparse_matrix,
                                 model_complex_plain_vector>
          cmodel_linear_solver;
  typedef std::shared_ptr<cmodel_linear_solver> cmodel_plsolver_type;


  template<typename MATRIX, typename VECTOR>
  std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> >
  default_linear_solver(const model &md) {
    std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> > p;

#if GETFEM_PARA_LEVEL == 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    if (md.is_symmetric())
      p.reset(new linear_solver_mumps_sym<MATRIX, VECTOR>);
    else
      p.reset(new linear_solver_mumps<MATRIX, VECTOR>);
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    if (md.is_symmetric())
      p.reset(new linear_solver_distributed_mumps_sym<MATRIX, VECTOR>);
    else
      p.reset(new linear_solver_distributed_mumps<MATRIX, VECTOR>);
#else
    size_type ndof = md.nb_dof(), max3d = 15000, dim = md.leading_dimension();
# ifdef GMM_USES_MUMPS
    max3d = 250000;
# endif
    if ((ndof<300000 && dim<=2) || (ndof<max3d && dim<=3) || (ndof<1000)) {
# ifdef GMM_USES_MUMPS
      if (md.is_symmetric())
        p.reset(new linear_solver_mumps_sym<MATRIX, VECTOR>);
      else
        p.reset(new linear_solver_mumps<MATRIX, VECTOR>);
# else
      p.reset(new linear_solver_superlu<MATRIX, VECTOR>);
# endif
    }
    else {
      if (md.is_coercive())
        p.reset(new linear_solver_cg_preconditioned_ildlt<MATRIX, VECTOR>);
      else {
        if (dim <= 2)
          p.reset(new
                  linear_solver_gmres_preconditioned_ilut<MATRIX,VECTOR>);
        else
          p.reset(new
                  linear_solver_gmres_preconditioned_ilu<MATRIX,VECTOR>);
      }
    }
#endif
    return p;
  }

  template <typename MATRIX, typename VECTOR>
  std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> >
  select_linear_solver(const model &md, const std::string &name) {
    std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> > p;
    if (bgeot::casecmp(name, "superlu") == 0)
      p.reset(new linear_solver_superlu<MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "mumps") == 0) {
#ifdef GMM_USES_MUMPS
# if GETFEM_PARA_LEVEL <= 1
      p.reset(new linear_solver_mumps<MATRIX, VECTOR>);
# else
      p.reset(new linear_solver_distributed_mumps<MATRIX, VECTOR>);
# endif
#else
      GMM_ASSERT1(false, "Mumps is not interfaced");
#endif
    }
    else if (bgeot::casecmp(name, "cg/ildlt") == 0)
      p.reset(new linear_solver_cg_preconditioned_ildlt<MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilu") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilu<MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilut") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilut<MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilutp") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilutp<MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "auto") == 0)
      p = default_linear_solver<MATRIX, VECTOR>(md);
    else
      GMM_ASSERT1(false, "Unknown linear solver");
    return p;
  }

  inline rmodel_plsolver_type rselect_linear_solver(const model &md,
                                             const std::string &name) {
    return select_linear_solver<model_real_sparse_matrix,
                                model_real_plain_vector>(md, name);
  }

  inline cmodel_plsolver_type cselect_linear_solver(const model &md,
                                             const std::string &name) {
    return select_linear_solver<model_complex_sparse_matrix,
                                model_complex_plain_vector>(md, name);
  }

  //---------------------------------------------------------------------
  // Standard solve.
  //---------------------------------------------------------------------


  /** A default solver for the model brick system.
  Of course it could be not very well suited for a particular
  problem, so it could be copied and adapted to change solvers,
  add a special traitement on the problem, etc ...  This is in
  fact a model for your own solver.

  For small problems, a direct solver is used
  (getfem::SuperLU_solve), for larger problems, a conjugate
  gradient gmm::cg (if the problem is coercive) or a gmm::gmres is
  used (preconditioned with an incomplete factorization).

  When MPI/METIS is enabled, a partition is done via METIS, and a parallel
  solver can be used.

  Note that it is possible to disable some variables
  (see the md.disable_variable(varname) method) in order to
  solve the problem only with respect to a subset of variables (the
  disabled variables are the considered as data) for instance to
  replace the global Newton strategy with a fixed point one.

  @ingroup bricks
  */
  void standard_solve(model &md, gmm::iteration &iter,
                      rmodel_plsolver_type lsolver,
                      abstract_newton_line_search &ls);

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver,
                      abstract_newton_line_search &ls);

  void standard_solve(model &md, gmm::iteration &iter,
                      rmodel_plsolver_type lsolver);

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver);

  void standard_solve(model &md, gmm::iteration &iter);

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODEL_SOLVERS_H__  */
