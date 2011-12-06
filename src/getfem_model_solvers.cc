// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2011 Yves Renard
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

#include "getfem/getfem_model_solvers.h"


namespace getfem {


  static rmodel_plsolver_type rdefault_linear_solver(const model &md) {
    return default_linear_solver<model_real_sparse_matrix,
                                 model_real_plain_vector>(md);
  } 

  static cmodel_plsolver_type cdefault_linear_solver(const model &md) {
    return default_linear_solver<model_complex_sparse_matrix,
                                 model_complex_plain_vector>(md);
  }

  /* ***************************************************************** */
  /*     Intermediary structure for Newton algorithms.                 */
  /* ***************************************************************** */

  template <typename MAT, typename VEC> 
  struct model_pb {

    typedef MAT MATRIX;
    typedef VEC VECTOR;
    typedef typename gmm::linalg_traits<VECTOR>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;

    model &md;
    gmm::abstract_newton_line_search &ls;
    VECTOR stateinit, &state;
    const VECTOR &rhs;
    const MATRIX &K;
    bool with_pseudo_potential;

    void compute_tangent_matrix(void)
    { md.to_variables(state); md.assembly(model::BUILD_MATRIX); }

    const MATRIX &tangent_matrix(void) { return K; }
    
    inline T scale_residual(void) const { return T(1); }

    void compute_residual(void)
    { md.to_variables(state); md.assembly(model::BUILD_RHS); }

    void compute_pseudo_potential(void)
    { md.to_variables(state); md.assembly(model::BUILD_PSEUDO_POTENTIAL); }

    void perturbation(void) {
      R res = gmm::vect_norm2(state), ampl = res / R(1000);
      if (res == R(0)) ampl = 1E-30;
      std::vector<R> V(gmm::vect_size(state));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), state);
    }

    const VECTOR &residual(void) { return rhs; }

    R residual_norm(void) { return gmm::vect_norm2(rhs); }

    R line_search(VECTOR &dr, const gmm::iteration &iter) {
      gmm::resize(stateinit, md.nb_dof());
      gmm::copy(state, stateinit);
      R alpha(1), res, R0;
      if (with_pseudo_potential) {
	compute_pseudo_potential();
	res = md.pseudo_potential();
      } else {
	res = residual_norm();
      }
      R0 = gmm::real(gmm::vect_sp(dr, rhs));

      ls.init_search(res, iter.get_iteration(), R0);
      do {
	alpha = ls.next_try();
	gmm::add(stateinit, gmm::scaled(dr, alpha), state);
	if (alpha < 1E-10) break;
	if (with_pseudo_potential) {
	  compute_pseudo_potential();
	  res = md.pseudo_potential();
	} else {
	  compute_residual();
	  res = residual_norm();
	}
	R0 = gmm::real(gmm::vect_sp(dr, rhs));
      } while (!ls.is_converged(res, R0));

      if (alpha != ls.converged_value() || with_pseudo_potential) {
	alpha = ls.converged_value();
	gmm::add(stateinit, gmm::scaled(dr, alpha), state);
	res = ls.converged_residual();
	compute_residual();
      }
      return alpha;
    }

    model_pb(model &m, gmm::abstract_newton_line_search &ls_, VECTOR &st,
	     const VECTOR &rhs_, const MATRIX &K_,
	     bool with_pseudo_pot = false)
      : md(m), ls(ls_), state(st), rhs(rhs_), K(K_),
	with_pseudo_potential(with_pseudo_pot) {}

  };

  /* ***************************************************************** */
  /*     Standard solve.                                               */
  /* ***************************************************************** */

  template <typename MATRIX, typename VECTOR, typename PLSOLVER>
  void standard_solve(model &md, gmm::iteration &iter,
		      PLSOLVER lsolver,
		      gmm::abstract_newton_line_search &ls, const MATRIX &K,
		      const VECTOR &rhs, bool with_pseudo_potential = false) {

    VECTOR state(md.nb_dof());
    
    md.from_variables(state); // copy the model variables in the state vector

    if (md.is_linear()) {
      md.assembly(model::BUILD_ALL);
      (*lsolver)(K, state, rhs, iter);
    }
    else {
      model_pb<MATRIX, VECTOR> mdpb(md, ls, state, rhs, K,
				    with_pseudo_potential);
      classical_Newton(mdpb, iter, *lsolver);
    }

    md.to_variables(state); // copy the state vector into the model variables
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      rmodel_plsolver_type lsolver,
		      gmm::abstract_newton_line_search &ls,
		      bool with_pseudo_potential) {
    standard_solve(md, iter, lsolver, ls, md.real_tangent_matrix(),
		   md.real_rhs(), with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      cmodel_plsolver_type lsolver,
		      gmm::abstract_newton_line_search &ls,
		      bool with_pseudo_potential) {
    standard_solve(md, iter, lsolver, ls, md.complex_tangent_matrix(),
		   md.complex_rhs(), with_pseudo_potential);
  }


  void standard_solve(model &md, gmm::iteration &iter,
			     rmodel_plsolver_type lsolver,
		      bool with_pseudo_potential) {
    gmm::default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls, with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
			     cmodel_plsolver_type lsolver,
		      bool with_pseudo_potential) {
    gmm::default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls, with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      bool with_pseudo_potential) {
    gmm::default_newton_line_search ls;
    if (md.is_complex())
      standard_solve(md, iter, cdefault_linear_solver(md), ls,
		     with_pseudo_potential);
    else
      standard_solve(md, iter, rdefault_linear_solver(md), ls,
		     with_pseudo_potential);
  }



}  /* end of namespace getfem.                                             */

