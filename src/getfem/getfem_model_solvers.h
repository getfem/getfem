/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2020 Yves Renard

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
    typedef MAT MATRIX;
    typedef VECT VECTOR;
    virtual void operator ()(const MAT &, VECT &, const VECT &,
                             gmm::iteration &) const = 0;
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

  template <typename MAT, typename VECT>
  struct linear_solver_dense_lu : public abstract_linear_solver<MAT, VECT> {
    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      typedef typename gmm::linalg_traits<MAT>::value_type T;
      gmm::dense_matrix<T> MM(gmm::mat_nrows(M),gmm::mat_ncols(M));
      gmm::copy(M, MM);
      gmm::lu_solve(MM, x, b);
      iter.enforce_converged(true);
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

  // Dummy linesearch for the newton with step control
  struct newton_search_with_step_control : public abstract_newton_line_search {

    virtual void init_search(double /*r*/, size_t /*git*/, double /*R0*/ = 0.0)
    { GMM_ASSERT1(false, "Not to be used"); }

    virtual double next_try(void)
    { GMM_ASSERT1(false, "Not to be used"); }

    virtual bool is_converged(double /*r*/, double /*R1*/ = 0.0)
    { GMM_ASSERT1(false, "Not to be used"); }

    newton_search_with_step_control() {}
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
              || (conv_alpha <= alpha_min && r < first_res * 1e5)
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
  /*     Newton(-Raphson) algorithm with step control.                 */
  /* ***************************************************************** */
  /* Still solves a problem F(X) = 0 starting at X0 but setting        */
  /* B0 = F(X0) and solving                                            */
  /* F(X) = (1-alpha)B0       (1)                                      */
  /* with alpha growing from 0 to 1.                                   */
  /* A very basic line search is applied.                              */
  /*                                                                   */
  /* Step 0 : set alpha0 = 0, alpha = 1, X0 given and B0 = F(X0).      */
  /* Step 1 : Set Ri = F(Xi) - (1-alpha)B0                             */
  /*          If ||Ri|| < rho, Xi+1 = Xi and go to step 2              */
  /*          Perform Newton step on problem (1)                       */
  /*          If the Newton do not converge (stagnation)               */
  /*            alpha <- max(alpha0+1E-4, (alpha+alpha0)/2)            */
  /*            Loop on step 1 with Xi unchanged                       */
  /* Step 2 : if alpha=1 stop                                          */
  /*          else alpha0 <- alpha,                                    */
  /*              alpha <- min(1,alpha0+2*(alpha-alpha0)),             */
  /*              Go to step 1 with Xi+1                               */
  /*              being the result of Newton iterations of step1.      */
  /*                                                                   */
  /*********************************************************************/

  template <typename PB>
  void Newton_with_step_control(PB &pb, gmm::iteration &iter)
  {
    typedef typename gmm::linalg_traits<typename PB::VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.reduce_noisy();
    iter_linsolv0.set_resmax(iter.get_resmax()/20.0);
    iter_linsolv0.set_maxiter(10000); // arbitrary

    pb.compute_residual();
    R approx_eln = pb.approx_external_load_norm();
    if (approx_eln == R(0)) approx_eln = R(1);

    typename PB::VECTOR b0(gmm::vect_size(pb.residual()));
    gmm::copy(pb.residual(), b0);
    typename PB::VECTOR Xi(gmm::vect_size(pb.residual()));
    gmm::copy(pb.state_vector(), Xi);

    typename PB::VECTOR dr(gmm::vect_size(pb.residual()));

    R alpha0(0), alpha(1), res0(gmm::vect_norm1(b0)), minres(res0);
    // const newton_search_with_step_control *ls
    //  = dynamic_cast<const newton_search_with_step_control *>(&(pb.ls));
    // GMM_ASSERT1(ls, "Internal error");
    size_type nit = 0, stagn = 0;
    R coeff = R(2);

    scalar_type crit = pb.residual_norm() / approx_eln;
    if (iter.finished(crit)) return;
    for(;;) {

      crit = gmm::vect_dist1(pb.residual(), gmm::scaled(b0, R(1)-alpha))
        / approx_eln;
      if (!iter.converged(crit)) {
        gmm::iteration iter_linsolv = iter_linsolv0;

        int is_singular = 1;
        while (is_singular) { // Linear system solve
          gmm::clear(dr);
          pb.add_to_residual(b0, alpha-R(1)); // canceled at next compute_residual
          iter_linsolv.init();
          if (iter.get_noisy() > 1)
            cout << "starting tangent matrix computation" << endl;
          pb.compute_tangent_matrix();
          if (iter.get_noisy() > 1)
            cout << "starting linear solver" << endl;
          pb.linear_solve(dr, iter_linsolv);
          if (!iter_linsolv.converged()) {
            is_singular++;
            if (is_singular <= 4) {
              if (iter.get_noisy())
                cout << "Singular tangent matrix:"
                  " perturbation of the state vector." << endl;
              pb.perturbation();
              pb.compute_residual(); // cancels add_to_residual
            } else {
              if (iter.get_noisy())
                cout << "Singular tangent matrix: perturbation failed, "
                     << "aborting." << endl;
              return;
            }
          }
          else is_singular = 0;
        }
        if (iter.get_noisy() > 1) cout << "linear solver done" << endl;

        gmm::add(dr, pb.state_vector());
        pb.compute_residual(); // cancels add_to_residual
        R res = gmm::vect_dist1(pb.residual(), gmm::scaled(b0, R(1)-alpha));
        R dec = R(1)/R(2), coeff2 = coeff * R(1.5);

        while (dec > R(1E-5) && res >= res0 * coeff) {
          gmm::add(gmm::scaled(dr, -dec), pb.state_vector());
          pb.compute_residual();
          R res2 = gmm::vect_dist1(pb.residual(), gmm::scaled(b0, R(1)-alpha));
          if (res2 < res*R(0.95) || res2 >= res0 * coeff2) {
            dec /= R(2); res = res2; coeff2 *= R(1.5);
          } else {
            gmm::add(gmm::scaled(dr, dec), pb.state_vector());
            break;
          }
        }
        dec *= R(2);

        nit++;
        coeff = std::max(R(1.05), coeff*R(0.93));
        bool near_end = (iter.get_iteration() > iter.get_maxiter()/2);
        bool cut = (alpha < R(1)) && near_end;
        if ((res > minres && nit > 4) || cut) {
          stagn++;

          if ((stagn > 10 && alpha > alpha0 + R(5E-2)) || cut) {
            alpha = (alpha + alpha0) / R(2);
            if (near_end) alpha = R(1);
            gmm::copy(Xi, pb.state_vector());
            pb.compute_residual();
            res = gmm::vect_dist1(pb.residual(), gmm::scaled(b0, R(1)-alpha));
            nit = 0;
            stagn = 0; coeff = R(2);
          }
        }
        if (res < minres || (alpha == R(1) &&  nit < 10)) stagn = 0;
        res0 = res;
        if (nit < 5) minres = res0; else minres = std::min(minres, res0);

        if (iter.get_noisy())
          cout << "step control [" << std::setw(8) << alpha0 << ","
               << std::setw(8) << alpha << "," << std::setw(10) << dec <<  "]";
        ++iter;
        // crit = std::min(res / approx_eln,
        //                gmm::vect_norm1(dr) / std::max(1E-25,pb.state_norm()));
        crit = res / approx_eln;
      }

      if (iter.finished(crit)) {
        if (iter.converged() && alpha < R(1)) {
          R a = alpha;
          alpha = std::min(R(1), alpha*R(3) - alpha0*R(2));
          alpha0 = a;
          gmm::copy(pb.state_vector(), Xi);
          nit = 0; stagn = 0; coeff = R(2);
        } else return;
      }
    }
  }



  /* ***************************************************************** */
  /*     Classical Newton(-Raphson) algorithm.                         */
  /* ***************************************************************** */

  template <typename PB>
  void classical_Newton(PB &pb, gmm::iteration &iter)
  {
    typedef typename gmm::linalg_traits<typename PB::VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.reduce_noisy();
    iter_linsolv0.set_resmax(iter.get_resmax()/20.0);
    iter_linsolv0.set_maxiter(10000); // arbitrary

    pb.compute_residual();
    R approx_eln = pb.approx_external_load_norm();
    if (approx_eln == R(0)) approx_eln = R(1);

    typename PB::VECTOR dr(gmm::vect_size(pb.residual()));

    scalar_type crit = pb.residual_norm() / approx_eln;
    while (!iter.finished(crit)) {
      gmm::iteration iter_linsolv = iter_linsolv0;

      int is_singular = 1;
      while (is_singular) {
        gmm::clear(dr);
        iter_linsolv.init();
        if (iter.get_noisy() > 1)
          cout << "starting computing tangent matrix" << endl;
        pb.compute_tangent_matrix();
        if (iter.get_noisy() > 1)
          cout << "starting linear solver" << endl;
        pb.linear_solve(dr, iter_linsolv);
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
      crit = std::min(pb.residual_norm() / approx_eln,
                      gmm::vect_norm1(dr) / std::max(1E-25, pb.state_norm()));
    }
  }



  //---------------------------------------------------------------------
  // Default linear solver.
  //---------------------------------------------------------------------

  typedef std::shared_ptr<abstract_linear_solver<model_real_sparse_matrix,
                                                 model_real_plain_vector> >
    rmodel_plsolver_type;
  typedef std::shared_ptr<abstract_linear_solver<model_complex_sparse_matrix,
                                                 model_complex_plain_vector> >
    cmodel_plsolver_type;

  template<typename MATRIX, typename VECTOR>
  std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR>>
  default_linear_solver(const model &md) {

#if GETFEM_PARA_LEVEL == 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    if (md.is_symmetric())
      return std::make_shared<linear_solver_mumps_sym<MATRIX, VECTOR>>();
    else
      return std::make_shared<linear_solver_mumps<MATRIX, VECTOR>>();
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    if (md.is_symmetric())
      return std::make_shared
        <linear_solver_distributed_mumps_sym<MATRIX, VECTOR>>();
      else
        return std::make_shared
          <linear_solver_distributed_mumps<MATRIX, VECTOR>>();
#else
    size_type ndof = md.nb_dof(), max3d = 15000, dim = md.leading_dimension();
# ifdef GMM_USES_MUMPS
    max3d = 250000;
# endif
    if ((ndof<300000 && dim<=2) || (ndof<max3d && dim<=3) || (ndof<1000)) {
# ifdef GMM_USES_MUMPS
      if (md.is_symmetric())
        return std::make_shared<linear_solver_mumps_sym<MATRIX, VECTOR>>();
      else
        return std::make_shared<linear_solver_mumps<MATRIX, VECTOR>>();
# else
      return std::make_shared<linear_solver_superlu<MATRIX, VECTOR>>();
# endif
    }
    else {
      if (md.is_coercive())
        return std::make_shared
          <linear_solver_cg_preconditioned_ildlt<MATRIX, VECTOR>>();
      else {
        if (dim <= 2)
          return std::make_shared
            <linear_solver_gmres_preconditioned_ilut<MATRIX,VECTOR>>();
          else
            return std::make_shared
              <linear_solver_gmres_preconditioned_ilu<MATRIX,VECTOR>>();
      }
    }
#endif
    return std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR>>();
  }

  template <typename MATRIX, typename VECTOR>
  std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR>>
  select_linear_solver(const model &md, const std::string &name) {
    std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR>> p;
    if (bgeot::casecmp(name, "superlu") == 0)
      return std::make_shared<linear_solver_superlu<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "dense_lu") == 0)
      return std::make_shared<linear_solver_dense_lu<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "mumps") == 0) {
#ifdef GMM_USES_MUMPS
# if GETFEM_PARA_LEVEL <= 1
      return std::make_shared<linear_solver_mumps<MATRIX, VECTOR>>();
# else
      return std::make_shared
        <linear_solver_distributed_mumps<MATRIX, VECTOR>>();
# endif
#else
      GMM_ASSERT1(false, "Mumps is not interfaced");
#endif
    }
    else if (bgeot::casecmp(name, "cg/ildlt") == 0)
      return std::make_shared
        <linear_solver_cg_preconditioned_ildlt<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "gmres/ilu") == 0)
      return std::make_shared
        <linear_solver_gmres_preconditioned_ilu<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "gmres/ilut") == 0)
      return std::make_shared
        <linear_solver_gmres_preconditioned_ilut<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "gmres/ilutp") == 0)
      return std::make_shared
        <linear_solver_gmres_preconditioned_ilutp<MATRIX, VECTOR>>();
    else if (bgeot::casecmp(name, "auto") == 0)
      return default_linear_solver<MATRIX, VECTOR>(md);
    else
      GMM_ASSERT1(false, "Unknown linear solver");
    return std::shared_ptr<abstract_linear_solver<MATRIX, VECTOR>>();
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
