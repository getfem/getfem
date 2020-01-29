/*===========================================================================

 Copyright (C) 2009-2020 Yves Renard

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

===========================================================================*/

#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm_inoutput.h"
#include <iomanip>

namespace getfem {

  #define TRACE_SOL 0

  /* ***************************************************************** */
  /*  Intermediary structure for Newton algorithms with getfem::model. */
  /* ***************************************************************** */

  template <typename PLSOLVER>
  class pb_base {
  public:
    typedef typename PLSOLVER::element_type::MATRIX MATRIX;
    typedef typename PLSOLVER::element_type::VECTOR VECTOR;
    typedef typename gmm::linalg_traits<VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;

  protected:
    PLSOLVER linear_solver;
    const MATRIX &K;
    VECTOR &rhs;
    VECTOR state;

  public:
    virtual const VECTOR &residual() const { return rhs; }
    // A norm1 seems to be better than a norm2 (at least for contact problems).
    virtual R residual_norm() { return gmm::vect_norm1(residual()); }
    virtual const VECTOR &state_vector() const { return state; }
    virtual VECTOR &state_vector() { return state; }
    virtual R state_norm() const { return gmm::vect_norm1(state_vector()); }

    virtual void perturbation() {
      R res = gmm::vect_norm2(state), ampl = std::max(res*R(1E-20), R(1E-50));
      std::vector<R> V(gmm::vect_size(state));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), state);
    }

    virtual void add_to_residual(VECTOR &extra_rhs, R mult=1.) {
      if (mult == R(1))       gmm::add(extra_rhs, rhs);
      else if (mult != R(0))  gmm::add(gmm::scaled(extra_rhs, mult), rhs);
    }

    virtual void linear_solve(VECTOR &dr, gmm::iteration &iter) {
      (*linear_solver)(K, dr, rhs, iter);
    }

    pb_base(PLSOLVER linsolv, const MATRIX &K_, VECTOR &rhs_)
      : linear_solver(linsolv), K(K_), rhs(rhs_), state(gmm::vect_size(rhs_))
    {}
  };

  /* ***************************************************************** */
  /*     Linear model problem.                                         */
  /* ***************************************************************** */
  template <typename PLSOLVER>
  class lin_model_pb : pb_base<PLSOLVER> {
    model &md;
  public:
    using typename pb_base<PLSOLVER>::VECTOR;
    using typename pb_base<PLSOLVER>::R;
    using pb_base<PLSOLVER>::state_vector;
    using pb_base<PLSOLVER>::linear_solve;
    void compute_all() { md.assembly(model::BUILD_ALL); }
    lin_model_pb(model &, PLSOLVER) = delete;
  };

  template <>
  lin_model_pb<rmodel_plsolver_type>::lin_model_pb
    (model &md_, rmodel_plsolver_type linsolv)
    : pb_base<rmodel_plsolver_type>
      (linsolv, md_.real_tangent_matrix(), md_.set_real_rhs()),
      md(md_)
  { md.from_variables(state_vector()); }
  template <>
  lin_model_pb<cmodel_plsolver_type>::lin_model_pb
    (model &md_, cmodel_plsolver_type linsolv)
    : pb_base<cmodel_plsolver_type>
      (linsolv, md_.complex_tangent_matrix(), md_.set_complex_rhs()),
      md(md_)
  { md.from_variables(state_vector()); }


  /* ***************************************************************** */
  /*     Non-linear model problem.                                     */
  /* ***************************************************************** */
  template <typename PLSOLVER>
  class nonlin_model_pb : protected pb_base<PLSOLVER> {
  public:
    using typename pb_base<PLSOLVER>::VECTOR;
    using typename pb_base<PLSOLVER>::R;
  protected:
    model &md;
    abstract_newton_line_search &ls;
  private:
    VECTOR stateinit; // just a temp used in line_search, could also be mutable
  public:
    using pb_base<PLSOLVER>::residual;
    using pb_base<PLSOLVER>::residual_norm;
    using pb_base<PLSOLVER>::state_vector;
    using pb_base<PLSOLVER>::state_norm;
    using pb_base<PLSOLVER>::add_to_residual;
    using pb_base<PLSOLVER>::perturbation;
    using pb_base<PLSOLVER>::linear_solve;

    virtual R approx_external_load_norm() { return md.approx_external_load(); }

    virtual void compute_tangent_matrix() {
      md.to_variables(state_vector());
      md.assembly(model::BUILD_MATRIX);
    }

    virtual void compute_residual() {
      md.to_variables(state_vector());
      md.assembly(model::BUILD_RHS);
    }

    virtual R line_search(VECTOR &dr, const gmm::iteration &iter) {
      size_type nit = 0;
      gmm::resize(stateinit, gmm::vect_size(state_vector()));
      gmm::copy(state_vector(), stateinit);
      R alpha(1), res, /* res_init, */ R0;

      /* res_init = */ res = residual_norm();
      // cout << "first residual = " << residual() << endl << endl;
      R0 = gmm::real(gmm::vect_sp(dr, residual()));
      ls.init_search(res, iter.get_iteration(), R0);
      do {
        alpha = ls.next_try();
        gmm::add(stateinit, gmm::scaled(dr, alpha), state_vector());

        compute_residual();
        res = residual_norm();
        // cout << "residual = " << residual() << endl << endl;
        R0 = gmm::real(gmm::vect_sp(dr, residual()));
        ++ nit;
      } while (!ls.is_converged(res, R0));

      if (alpha != ls.converged_value()) {
        alpha = ls.converged_value();
        gmm::add(stateinit, gmm::scaled(dr, alpha), state_vector());
        res = ls.converged_residual();
        compute_residual();
      }

      return alpha;
    }

    nonlin_model_pb(model &md_, abstract_newton_line_search &ls_,
                    PLSOLVER linear_solver_) = delete;
  };

  template <>
  nonlin_model_pb<rmodel_plsolver_type>::nonlin_model_pb
    (model &md_, abstract_newton_line_search &ls_, rmodel_plsolver_type linsolv)
    : pb_base<rmodel_plsolver_type>
      (linsolv, md_.real_tangent_matrix(), md_.set_real_rhs()),
      md(md_), ls(ls_)
  { md.from_variables(state_vector()); }
  template <>
  nonlin_model_pb<cmodel_plsolver_type>::nonlin_model_pb
    (model &md_, abstract_newton_line_search &ls_, cmodel_plsolver_type linsolv)
    : pb_base<cmodel_plsolver_type>
      (linsolv, md_.complex_tangent_matrix(), md_.set_complex_rhs()),
      md(md_), ls(ls_)
  { md.from_variables(state_vector()); }



  /* ***************************************************************** */
  /*     Non-linear model problem with internal variable condensation. */
  /* ***************************************************************** */
  template <typename PLSOLVER>
  class nonlin_condensed_model_pb : public nonlin_model_pb<PLSOLVER> {
  public:
    using typename pb_base<PLSOLVER>::MATRIX;
    using typename pb_base<PLSOLVER>::VECTOR;
    using typename pb_base<PLSOLVER>::R;

  private:
    using nonlin_model_pb<PLSOLVER>::md;
    bool condensed_vars;
    const MATRIX &internal_K;
    VECTOR &full_rhs;

  public:
    virtual const VECTOR &residual() const { return full_rhs; }
    // A norm1 seems to be better than a norm2 (at least for contact problems).
    virtual R residual_norm() { return gmm::vect_norm1(full_rhs); }

    using pb_base<PLSOLVER>::state_vector;
    using pb_base<PLSOLVER>::state_norm;

    virtual void add_to_residual(VECTOR &extra_rhs, R mult=1.) {
      if (mult == R(1))      gmm::add(extra_rhs, full_rhs);
      else if (mult != R(0)) gmm::add(gmm::scaled(extra_rhs, mult), full_rhs);
    }

    using pb_base<PLSOLVER>::perturbation;

    virtual void linear_solve(VECTOR &dr, gmm::iteration &iter) {
      VECTOR dr0(md.nb_primary_dof(), R(0));
      pb_base<PLSOLVER>::linear_solve(dr0, iter);
      gmm::sub_interval I_prim(0, md.nb_primary_dof()),
                        I_intern(md.nb_primary_dof(), md.nb_internal_dof());
      gmm::copy(dr0, gmm::sub_vector(dr, I_prim));
      // recover solution of condensed variables: intern_K*prim_sol + intern_sol
      gmm::mult(internal_K,
                gmm::scaled(gmm::sub_vector(dr, I_prim), scalar_type(-1)),
                md.internal_solution(),
                gmm::sub_vector(dr, I_intern));
    }

    virtual void compute_tangent_matrix() {
      md.to_variables(state_vector(), condensed_vars);
      md.assembly(condensed_vars ? model::BUILD_MATRIX_CONDENSED
                                 : model::BUILD_MATRIX);
    }

    virtual void compute_residual() {
      md.to_variables(state_vector(), condensed_vars);
      md.assembly(condensed_vars ? model::BUILD_RHS_WITH_INTERNAL // --> full_rhs
                                 : model::BUILD_RHS);             // ---> rhs (= full_rhs)
    }

    nonlin_condensed_model_pb(model &md_, abstract_newton_line_search &ls_,
                              PLSOLVER linear_solver_) = delete;
  };

  template <>
  nonlin_condensed_model_pb<rmodel_plsolver_type>::nonlin_condensed_model_pb
    (model &md_, abstract_newton_line_search &ls_, rmodel_plsolver_type linsolv)
    : nonlin_model_pb<rmodel_plsolver_type>(md_, ls_, linsolv),
      condensed_vars(md_.has_internal_variables()),
      internal_K(md_.real_tangent_matrix(condensed_vars)),
      full_rhs(md_.set_real_rhs(condensed_vars))
  {
    gmm::resize(state_vector(), md.nb_dof(condensed_vars));
    md.from_variables(state_vector(), condensed_vars);
  }

  template <>
  nonlin_condensed_model_pb<cmodel_plsolver_type>::nonlin_condensed_model_pb
    (model &md_, abstract_newton_line_search &ls_, cmodel_plsolver_type linsolv)
    : nonlin_model_pb<cmodel_plsolver_type>(md_, ls_, linsolv),
      condensed_vars(false), internal_K(md_.complex_tangent_matrix()),
      full_rhs(md_.set_complex_rhs())
  {
    GMM_ASSERT1(false, "No support for internal variable condensation in "
                       "complex valued models.");
  }



  /* ***************************************************************** */
  /*   Linear solvers.                                                 */
  /* ***************************************************************** */
  static rmodel_plsolver_type rdefault_linear_solver(const model &md) {
    return default_linear_solver<model_real_sparse_matrix,
                                 model_real_plain_vector>(md);
  }

  static cmodel_plsolver_type cdefault_linear_solver(const model &md) {
    return default_linear_solver<model_complex_sparse_matrix,
                                 model_complex_plain_vector>(md);
  }

  void default_newton_line_search::init_search(double r, size_t git, double) {
    alpha_min_ratio = 0.9;
    alpha_min = 1e-10;
    alpha_max_ratio = 10.0;
    alpha_mult = 0.25;
    itmax = size_type(-1);
    glob_it = git; if (git <= 1) count_pat = 0;
    conv_alpha = alpha = alpha_old = 1.;
    conv_r = first_res = r; it = 0;
    count = 0;
    max_ratio_reached = false;
  }

  double default_newton_line_search::next_try() {
    alpha_old = alpha; ++it;
    // alpha *= 0.5;
    if (alpha >= 0.4) alpha *= 0.5; else alpha *= alpha_mult;
    return alpha_old;
  }

  bool default_newton_line_search::is_converged(double r, double) {
    // cout << "r = " << r << " alpha = " << alpha_old << " count_pat = " << count_pat << endl;
    if (!max_ratio_reached && r < first_res * alpha_max_ratio) {
      alpha_max_ratio_reached = alpha_old; r_max_ratio_reached = r;
      it_max_ratio_reached = it; max_ratio_reached = true;
    }
    if (max_ratio_reached &&
        r < r_max_ratio_reached * 0.5 &&
        r > first_res * 1.1 && it <= it_max_ratio_reached+1) {
      alpha_max_ratio_reached = alpha_old; r_max_ratio_reached = r;
      it_max_ratio_reached = it;
    }
    if (count == 0 || r < conv_r)
      { conv_r = r; conv_alpha = alpha_old; count = 1; }
    if (conv_r < first_res) ++count;

    if (r < first_res *  alpha_min_ratio)
      { count_pat = 0; return true; }
    if (count>=5 || (alpha < alpha_min && max_ratio_reached) || alpha<1e-15) {
      if (conv_r < first_res * 0.99) count_pat = 0;
      if (/*gmm::random() * 50. < -log(conv_alpha)-4.0 ||*/ count_pat >= 3)
        { conv_r=r_max_ratio_reached; conv_alpha=alpha_max_ratio_reached; }
      if (conv_r >= first_res * 0.999) count_pat++;
      return true;
    }
    return false;
  }


  /* ***************************************************************** */
  /*     Computation of initial values of velocity/acceleration for    */
  /*     time integration schemes.                                     */
  /* ***************************************************************** */

  template <typename PLSOLVER>
  void compute_init_values(model &md, gmm::iteration &iter,
                           PLSOLVER lsolver,
                           abstract_newton_line_search &ls) {

    typename PLSOLVER::element_type::VECTOR state(md.nb_dof());
    md.from_variables(state);
    md.cancel_init_step();
    md.set_time_integration(2);
    scalar_type dt = md.get_time_step();
    scalar_type ddt = md.get_init_time_step();
    scalar_type t = md.get_time();

    // Solve for ddt
    md.set_time_step(ddt);
    gmm::iteration iter1 = iter;
    standard_solve(md, iter1, lsolver, ls);
    md.copy_init_time_derivative();

    // Restore the model state
    md.set_time_step(dt);
    md.set_time(t);
    md.to_variables(state);
    md.set_time_integration(1);
  }

  /* ***************************************************************** */
  /*     Standard solve.                                               */
  /* ***************************************************************** */

  template <typename PLSOLVER>
  void standard_solve(model &md, gmm::iteration &iter,
                      PLSOLVER lsolver,
                      abstract_newton_line_search &ls) {

    int time_integration = md.is_time_integration();
    if (time_integration) {
      if (time_integration == 1 && md.is_init_step()) {
        compute_init_values(md, iter, lsolver, ls);
        return;
      }
      md.set_time(md.get_time() + md.get_time_step());
      md.call_init_affine_dependent_variables(time_integration);
    }

    if (md.is_linear()) {
      lin_model_pb<PLSOLVER> mdpb(md, lsolver);
      mdpb.compute_all();
      mdpb.linear_solve(mdpb.state_vector(), iter);
      md.to_variables(mdpb.state_vector()); // copy the state vector into the model variables
    } else {
      std::unique_ptr<nonlin_model_pb<PLSOLVER>> mdpb;
      if (md.has_internal_variables())
        mdpb = std::make_unique<nonlin_condensed_model_pb<PLSOLVER>>(md, ls, lsolver);
      else
        mdpb = std::make_unique<nonlin_model_pb<PLSOLVER>>(md, ls, lsolver);
      if (dynamic_cast<newton_search_with_step_control *>(&ls))
        Newton_with_step_control(*mdpb, iter);
      else
        classical_Newton(*mdpb, iter);
      md.to_variables(mdpb->state_vector()); // copy the state vector into the model variables
    }
  }

  void standard_solve(model &md, gmm::iteration &iter,
                      rmodel_plsolver_type lsolver,
                      abstract_newton_line_search &ls) {
    standard_solve<rmodel_plsolver_type>(md, iter, lsolver, ls);
  }

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver,
                      abstract_newton_line_search &ls) {
    standard_solve<cmodel_plsolver_type>(md, iter, lsolver, ls);
  }


  void standard_solve(model &md, gmm::iteration &iter,
                      rmodel_plsolver_type lsolver) {
    default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls);
  }

  void standard_solve(model &md, gmm::iteration &iter,
                      cmodel_plsolver_type lsolver) {
    newton_search_with_step_control ls;
    // default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls);
  }

  void standard_solve(model &md, gmm::iteration &iter) {
    newton_search_with_step_control ls;
    // default_newton_line_search ls;
    if (md.is_complex())
      standard_solve(md, iter, cdefault_linear_solver(md), ls);
    else
      standard_solve(md, iter, rdefault_linear_solver(md), ls);
  }



}  /* end of namespace getfem.                                             */

