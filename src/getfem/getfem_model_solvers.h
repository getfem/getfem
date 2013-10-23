/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
    double alpha, alpha_mult, first_res, alpha_max_ratio, alpha_min;
    virtual void init_search(double r, size_t git, double = 0.0) {
      glob_it = git;
      conv_alpha = alpha = double(1); conv_r = first_res = r; it = 0;
    }
    virtual double next_try(void)
    { conv_alpha = alpha; alpha *= alpha_mult; ++it; return conv_alpha; }
    virtual bool is_converged(double r, double = 0.0) {
      conv_r = r;
      return ((it <= 1 && r < first_res)
	      || (r <= first_res * alpha_max_ratio)
	      || (conv_alpha <= alpha_min)
	      || it >= itmax);
    }
    simplest_newton_line_search
    (size_t imax = size_t(-1), double a_max_ratio = 6.0/5.0,
     double a_min = 1.0/1000.0, double a_mult = 3.0/5.0)
      : alpha_mult(a_mult), alpha_max_ratio(a_max_ratio), alpha_min(a_min)
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
      if (iter.get_noisy()) cout << "alpha = " << alpha << " ";
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
    bool is_reduced;
    std::vector<size_type> &sind;
    gmm::sub_index I;
    abstract_newton_line_search &ls;
    VECTOR stateinit, &state;
    const VECTOR &rhs;
    const MATRIX &K;
    MATRIX Kr;
    VECTOR rhsr;
    bool with_pseudo_potential;

    void compute_tangent_matrix(void) {
      md.to_variables(state);
      md.assembly(model::BUILD_MATRIX);
      if (is_reduced) {
	gmm::resize(Kr, sind.size(), sind.size());
	gmm::copy(gmm::sub_matrix(K, I, I), Kr);
      }
    }

    const MATRIX &tangent_matrix(void) { return (is_reduced ? Kr : K); }
    
    inline T scale_residual(void) const { return T(1); }

    void compute_residual(void) {
      md.to_variables(state);
      md.assembly(model::BUILD_RHS);
      if (is_reduced) {
	gmm::resize(rhsr, sind.size());
	gmm::copy(gmm::sub_vector(rhs, I), rhsr);
      }
    }

    void compute_pseudo_potential(void)
    { md.to_variables(state); md.assembly(model::BUILD_PSEUDO_POTENTIAL); }

    void perturbation(void) {
      R res = gmm::vect_norm2(state), ampl = std::max(res*R(1E-20), R(1E-50));
      std::vector<R> V(gmm::vect_size(state));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), state);
    }

    const VECTOR &residual(void) { return (is_reduced ? rhsr : rhs); }
   
    R state_norm(void) const
    { return gmm::vect_norm1(gmm::sub_vector(state, I)); }

    R approx_external_load_norm(void)
    { return md.approx_external_load(); }

    R residual_norm(void) { // A norm1 seems to be better than a norm2
                            // at least for contact problems.
      return (is_reduced ? gmm::vect_norm1(rhsr) : gmm::vect_norm1(rhs));
    }

    R compute_res(bool comp = true) {
      if (with_pseudo_potential) {
	compute_pseudo_potential();
	return md.pseudo_potential();
      } else {
	if (comp) compute_residual();
	return residual_norm();
      }
    }


    R line_search(VECTOR &dr, const gmm::iteration &iter) {
      size_type nit = 0;
      gmm::resize(stateinit, md.nb_dof());
      gmm::copy(state, stateinit);
      R alpha(1), res, /* res_init, */ R0;
	  
      /* res_init = */ res = compute_res(false);      
      // cout << "first residual = " << residual() << endl << endl;
      R0 = (is_reduced ? gmm::real(gmm::vect_sp(dr, rhsr))
	               : gmm::real(gmm::vect_sp(dr, rhs)));
      
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
	gmm::add(gmm::sub_vector(stateinit, I), gmm::scaled(dr, alpha),
		 gmm::sub_vector(state, I));
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
	R0 = (is_reduced ? gmm::real(gmm::vect_sp(dr, rhsr))
	                 : gmm::real(gmm::vect_sp(dr, rhs)));

	++ nit;
      } while (!ls.is_converged(res, R0));

      if (alpha != ls.converged_value() || with_pseudo_potential) {
	alpha = ls.converged_value();
	gmm::add(gmm::sub_vector(stateinit, I), gmm::scaled(dr, alpha),
		 gmm::sub_vector(state, I));
	res = ls.converged_residual();
	compute_residual();
      }

      return alpha;
    }

    model_pb(model &m, abstract_newton_line_search &ls_, VECTOR &st,
	     const VECTOR &rhs_, const MATRIX &K_, bool reduced_,
	     std::vector<size_type> &sind_,
	     bool with_pseudo_pot = false)
      : md(m), is_reduced(reduced_), sind(sind_), I(sind_), ls(ls_), state(st),
	rhs(rhs_), K(K_), with_pseudo_potential(with_pseudo_pot) {}

  };

  //---------------------------------------------------------------------
  // Default linear solver.
  //---------------------------------------------------------------------

  typedef abstract_linear_solver<model_real_sparse_matrix,
                                 model_real_plain_vector> rmodel_linear_solver;
  typedef dal::shared_ptr<rmodel_linear_solver> rmodel_plsolver_type;
  typedef abstract_linear_solver<model_complex_sparse_matrix,
                                 model_complex_plain_vector>
          cmodel_linear_solver;
  typedef dal::shared_ptr<cmodel_linear_solver> cmodel_plsolver_type;


  template<typename MATRIX, typename VECTOR>
  dal::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> >
  default_linear_solver(const model &md) {
    dal::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> > p;
    
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
  dal::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> >
  select_linear_solver(const model &md, const std::string &name) {
    dal::shared_ptr<abstract_linear_solver<MATRIX, VECTOR> > p;
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
  used (preconditionned with an incomplete factorization).

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
                      abstract_newton_line_search &ls,
                      bool with_pseudo_potential = false);

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver,
                      abstract_newton_line_search &ls,
                      bool with_pseudo_potential = false);

  void standard_solve(model &md, gmm::iteration &iter,
                      rmodel_plsolver_type lsolver,
                      bool with_pseudo_potential = false);

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver,
                      bool with_pseudo_potential = false);

  void standard_solve(model &md, gmm::iteration &iter,
                      bool with_pseudo_potential = false);


}






//---------------------------------------------------------------------
//
// Solvers for the old brick system. Kept for compatibility reasons.
//
//---------------------------------------------------------------------


#include "getfem_modeling.h"

namespace getfem {

#if GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == SCHWARZADD_PARA_SOLVER
  template <typename MODEL_STATE, typename MAT, typename VECT>
  struct linear_solver_para_schwarzadd
    : public abstract_linear_solver<MAT, VECT> {

    typedef typename MODEL_STATE::value_type value_type;

    const mdbrick_abstract<MODEL_STATE> &problem;
    int nblocsubdom; // Number of sub-domains per process

    linear_solver_para_schwarzadd(const mdbrick_abstract<MODEL_STATE> &problem_,
                                  int nblocsubdom_)
      : problem(problem_), nblocsubdom(nblocsubdom_) {}

    void operator ()(const MAT &M, VECT &x, const VECT &b,
                     gmm::iteration &iter) const {
      double tt_ref=MPI_Wtime();

      // Meshes sub-partition.
      std::set<const mesh *> mesh_set;
      for (size_type i = 0; i < problem.nb_mesh_fems(); ++i)
        mesh_set.insert(&(problem.get_mesh_fem(i).linked_mesh()));

      std::vector< std::vector<int> > eparts(mesh_set.size());

      int size, rank, nset = 0;
      size_type ndof = problem.nb_dof();
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      for (std::set<const mesh *>::iterator it = mesh_set.begin();
           it != mesh_set.end(); ++it, ++nset) {

        int ne = int((*it)->get_mpi_region().nb_convex());
        std::vector<int> xadj(ne+1), adjncy, numelt((*it)->convex_index().last_true()+1), numeltinv(ne), npart(ne);

        int j = 0, k = 0;
        bgeot::mesh_structure::ind_set s;
        for (mr_visitor ic((*it)->get_mpi_region()); !ic.finished();++ic,++j)
          { numelt[ic.cv()] = j; numeltinv[j] = ic.cv(); }

        j = 0;
        for (mr_visitor ic((*it)->get_mpi_region()); !ic.finished();++ic,++j) {
          xadj[j] = k;
          (*it)->neighbours_of_convex(ic.cv(), s);
          for (bgeot::mesh_structure::ind_set::iterator iti = s.begin();
               iti != s.end(); ++iti)
            if ((*it)->get_mpi_region().is_in(*iti))
              { adjncy.push_back(numelt[*iti]); ++k; }
        }
        xadj[j] = k;

        int wgtflag = 0, edgecut, numflag = 0, options[5] = {0,0,0,0,0};
        int nbbl = nblocsubdom/size;

        // The mpi region is partitionned into nblocsubdom sub-domains.
        METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), 0, 0, &wgtflag,
                            &numflag, &nbbl, options, &edgecut,
                            &(npart[0]));

        eparts[nset].resize(0);
        eparts[nset].resize((*it)->convex_index().last()+1, size_type(-1));

        for (size_type i = 0; i < size_type(ne); ++i)
          eparts[nset][numeltinv[i]] = npart[i];
      }

      // To be completeted for non-finite element dofs
      // nblocsubdom is number of sub dom per proc
//       std::vector<dal::bit_vector> Bidof(nblocsubdom);
      // nblocsubdom is the global number of  sub dom
      std::vector<dal::bit_vector> Bidof(nblocsubdom/size);
      size_type apos = 0;
      for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
        const mesh_fem &mf = problem.get_mesh_fem(i);
        nset = 0;
        for (std::set<const mesh *>::iterator it = mesh_set.begin();
             it != mesh_set.end(); ++it, ++nset)
          if (*it == &(mf.linked_mesh())) break;
        size_type pos = problem.get_mesh_fem_position(i);
        GMM_ASSERT1(pos == apos, "Multipliers are not taken into account");
        size_type length = mf.nb_dof();
        apos += length;
        for (dal::bv_visitor j(mf.convex_index()); !j.finished(); ++j) {
          size_type k = eparts[nset][j];
          if (k != size_type(-1))
            for (size_type l = 0; l < mf.nb_dof_of_element(j); ++l)
              Bidof[k].add(mf.ind_dof_of_element(j)[l] + pos);
        }
        GMM_ASSERT1(apos == ndof, "Multipliers are not taken into account");
      }
      // nblocsubdom is number of sub dom per proc
//       std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bi(nblocsubdom*size);
      // nblocsubdom is the global number of  sub dom
      std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bi(nblocsubdom);
      // nblocsubdom is number of sub dom per proc
//       for (size_type i = 0; i < size_type(nblocsubdom); ++i) {
      // nblocsubdom is the global number of  sub dom
      for (size_type i = 0; i < size_type(nblocsubdom/size); ++i) {
     // nblocsubdom is number of sub dom per proc
//         gmm::resize(Bi[size*rank + i], ndof, Bidof[i].card());
     // nblocsubdom is the global number of  sub dom
        gmm::resize(Bi[(nblocsubdom/size)*rank + i], ndof, Bidof[i].card());
        size_type k = 0;
        for (dal::bv_visitor j(Bidof[i]); !j.finished(); ++j, ++k)
    // nblocsubdom is number of sub dom per proc
//           Bi[size*rank + i](j, k) = value_type(1);
     // nblocsubdom is the global number of  sub dom
          Bi[(nblocsubdom/size)*rank + i](j, k) = value_type(1);
      }

      gmm::mpi_distributed_matrix<MAT> mpiM(ndof, ndof);
      gmm::copy(M, mpiM.local_matrix());

      additive_schwarz(mpiM, x, b, gmm::identity_matrix(), Bi, iter,
                       gmm::using_superlu(), gmm::using_cg());

      cout<<"temps SCHWARZ ADD "<< MPI_Wtime() - tt_ref<<endl;
    }
  };
#endif

  template <typename MODEL_STATE> struct useful_types {

    TYPEDEF_MODEL_STATE_TYPES;
    typedef abstract_linear_solver<T_MATRIX, VECTOR> lsolver_type;
    typedef dal::shared_ptr<lsolver_type> plsolver_type;
  };


  template <typename MODEL_STATE>
  typename useful_types<MODEL_STATE>::plsolver_type
  default_linear_solver(const mdbrick_abstract<MODEL_STATE> &problem) {
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::vector_type VECTOR;

    typename useful_types<MODEL_STATE>::plsolver_type p;

#if GETFEM_PARA_LEVEL == 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    p.reset(new linear_solver_mumps<T_MATRIX, VECTOR>);
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
    p.reset(new linear_solver_distributed_mumps<T_MATRIX, VECTOR>);
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == SCHWARZADD_PARA_SOLVER
    int size;
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
    size=32;//global number of sub_dom
    p.reset(new linear_solver_para_schwarzadd<MODEL_STATE, T_MATRIX, VECTOR>(problem, size));
#else
    size_type ndof = problem.nb_dof(), max3d = 15000, dim = problem.dim();
# ifdef GMM_USES_MUMPS
    max3d = 100000;
# endif
    if ((ndof<200000 && dim<=2) || (ndof<max3d && dim<=3) || (ndof<1000)) {
# ifdef GMM_USES_MUMPS
      p.reset(new linear_solver_mumps<T_MATRIX, VECTOR>);
# else
      p.reset(new linear_solver_superlu<T_MATRIX, VECTOR>);
# endif
    }
    else {
      if (problem.is_coercive())
        p.reset(new linear_solver_cg_preconditioned_ildlt<T_MATRIX, VECTOR>);
      else if (problem.mixed_variables().card() == 0) {
        if (dim <= 2)
          p.reset(new
                  linear_solver_gmres_preconditioned_ilu<T_MATRIX,VECTOR>);
        else
          p.reset(new
                  linear_solver_gmres_preconditioned_ilu<T_MATRIX,VECTOR>);
      }
      else {
        if (dim <= 2)
          p.reset(new
                  linear_solver_gmres_preconditioned_ilu<T_MATRIX,VECTOR>);
        else
          p.reset(new
                  linear_solver_gmres_preconditioned_ilu<T_MATRIX,VECTOR>);
      }
    }
#endif
    return p;
  }


  template <typename MODEL_STATE>
  typename useful_types<MODEL_STATE>::plsolver_type
  select_linear_solver(const mdbrick_abstract<MODEL_STATE> &problem,
                       const std::string &name) {
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::vector_type VECTOR;

    typename useful_types<MODEL_STATE>::plsolver_type p;

    if (bgeot::casecmp(name, "superlu") == 0)
      p.reset(new linear_solver_superlu<T_MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "mumps") == 0) {
#ifdef GMM_USES_MUMPS
# if GETFEM_PARA_LEVEL <= 1
      p.reset(new linear_solver_mumps<T_MATRIX, VECTOR>);
# else
      p.reset(new linear_solver_distributed_mumps<T_MATRIX, VECTOR>);
# endif
#else
      GMM_ASSERT1(false, "Mumps is not interfaced");
#endif
    }
    else if (bgeot::casecmp(name, "cg/ildlt") == 0)
      p.reset(new linear_solver_cg_preconditioned_ildlt<T_MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilu") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilu<T_MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilut") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilut<T_MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "gmres/ilutp") == 0)
      p.reset(new linear_solver_gmres_preconditioned_ilutp<T_MATRIX, VECTOR>);
    else if (bgeot::casecmp(name, "auto") == 0)
      p = default_linear_solver(problem);
    else
      GMM_ASSERT1(false, "Unknown linear solver");
    return p;
  }


  /* ***************************************************************** */
  /*   Intermediary structure for Newton algorithms with MODEL_STATE.  */
  /* ***************************************************************** */

  template <typename MODEL_STATE>
  struct model_problem {

    TYPEDEF_MODEL_STATE_TYPES;
    typedef T_MATRIX MATRIX;
    typedef typename gmm::linalg_traits<VECTOR>::value_type T;


    MODEL_STATE &MS;
    mdbrick_abstract<MODEL_STATE> &pb;
    abstract_newton_line_search &ls;
    VECTOR stateinit, d;

    void compute_tangent_matrix(void) {
      pb.compute_tangent_matrix(MS);
      if (pb.nb_constraints() > 0) {
        pb.compute_residual(MS);
        MS.compute_reduced_system();
      }
    }

    void perturbation(void) {
      R res = gmm::vect_norm2(MS.state()), ampl = res / R(1000);
      if (res == R(0)) ampl = 1E-30;
      VECTOR V(gmm::vect_size(MS.state()));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), MS.state());
    }

    inline T scale_residual(void) const { return T(-1); }

    const T_MATRIX &tangent_matrix(void)
    { return MS.reduced_tangent_matrix(); }

    void compute_residual(void) {
      pb.compute_residual(MS);
      if (pb.nb_constraints() > 0) MS.compute_reduced_residual();
    }

    const VECTOR &residual(void) { return MS.reduced_residual(); }

    R residual_norm(void) { return MS.reduced_residual_norm(); }
    R approx_external_load_norm(void) { return R(1);} // Not taken into account
    R state_norm(void) const { return R(0); } // Not taken into account

    R line_search(VECTOR &dr, const gmm::iteration &iter) {
      gmm::resize(d, pb.nb_dof());
      gmm::resize(stateinit, pb.nb_dof());
      gmm::copy(MS.state(), stateinit);
      MS.unreduced_solution(dr, d);
      R alpha(1), res, R0;
      R0 = gmm::real(gmm::vect_sp(dr, residual()));

      ls.init_search(MS.reduced_residual_norm(), iter.get_iteration(), R0);
      do {
        alpha = ls.next_try();
        gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
        compute_residual();
        res = MS.reduced_residual_norm();
        R0 = gmm::real(gmm::vect_sp(dr, residual()));
      } while (!ls.is_converged(res));

      if (alpha != ls.converged_value()) {
        alpha = ls.converged_value();
        gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
        res = ls.converged_residual();
        compute_residual();
      }
      return alpha;
    }

    model_problem(MODEL_STATE &MS_, mdbrick_abstract<MODEL_STATE> &pb_,
                  abstract_newton_line_search &ls_)
      : MS(MS_), pb(pb_), ls(ls_) {}

  };


  /* ***************************************************************** */
  /*     Standard solve.                                               */
  /* ***************************************************************** */

  /** A default solver for the old model brick system.

  Of course it could be not very well suited for a particular
  problem, so it could be copied and adapted to change solvers,
  add a special traitement on the problem, etc ...  This is in
  fact a model for your own solver.

  For small problems, a direct solver is used
  (getfem::SuperLU_solve), for larger problems, a conjugate
  gradient gmm::cg (if the problem is coercive) or a gmm::gmres is
  used (preconditionned with an incomplete factorization).

  When MPI/METIS is enabled, a partition is done via METIS, and a
  parallel solver can be used.

  @ingroup bricks
  */
  template <typename MODEL_STATE> void
  standard_solve
  (MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
   gmm::iteration &iter,
   typename useful_types<MODEL_STATE>::plsolver_type lsolver,
   abstract_newton_line_search &ls) {

    TYPEDEF_MODEL_STATE_TYPES;
    model_problem<MODEL_STATE> mdpb(MS, problem, ls);

    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before

    if (problem.is_linear()) {
      mdpb.compute_tangent_matrix();
      mdpb.compute_residual();
      VECTOR dr(gmm::vect_size(mdpb.residual())), d(problem.nb_dof());
      VECTOR b(gmm::vect_size(dr));
      gmm::copy(gmm::scaled(mdpb.residual(), value_type(-1)), b);
      // cout << "tg matrix = " << mdpb.tangent_matrix() << endl;
      // print_eigval(mdpb.tangent_matrix());
      (*lsolver)(mdpb.tangent_matrix(), dr, b, iter);
      MS.unreduced_solution(dr, d);
      gmm::add(d, MS.state());
    }
    else
      classical_Newton(mdpb, iter, *lsolver);
  }

  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
                 gmm::iteration &iter,
                 typename useful_types<MODEL_STATE>::plsolver_type lsolver) {
    default_newton_line_search ls;
    standard_solve(MS, problem, iter, lsolver, ls);
  }


  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
                 gmm::iteration &iter) {
    getfem::default_newton_line_search ls;
    standard_solve(MS, problem, iter, default_linear_solver(problem), ls);
  }






}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODEL_SOLVERS_H__  */
