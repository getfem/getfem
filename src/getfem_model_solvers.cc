// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard
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

    void compute_tangent_matrix(void)
    { md.to_variables(state); md.assembly(model::BUILD_MATRIX); }

    const MATRIX &tangent_matrix(void) { return K; }
    
    inline T scale_residual(void) const { return T(1); }

    void compute_residual(void)
    { md.to_variables(state); md.assembly(model::BUILD_RHS); }

    const VECTOR &residual(void) { return rhs; }

    R residual_norm(void) { return gmm::vect_norm2(rhs); }

    R line_search(VECTOR &dr, const gmm::iteration &iter) {
      gmm::resize(stateinit, md.nb_dof());
      gmm::copy(state, stateinit);
      R alpha(1), res;
      
      ls.init_search(gmm::vect_norm2(residual()), iter.get_iteration());
      do {
	alpha = ls.next_try();
	gmm::add(stateinit, gmm::scaled(dr, alpha), state);
	if (alpha < 1E-10) break;
	compute_residual();
	res = residual_norm();
      } while (!ls.is_converged(res));

      if (alpha != ls.converged_value()) {
	alpha = ls.converged_value();
	gmm::add(stateinit, gmm::scaled(dr, alpha), state);
	res = ls.converged_residual();
	compute_residual();
      }
      return alpha;
    }

    model_pb(model &m, gmm::abstract_newton_line_search &ls_, VECTOR &st,
	     const VECTOR &rhs_, const MATRIX &K_)
      : md(m), ls(ls_), state(st), rhs(rhs_), K(K_) {}

  };

  /* ***************************************************************** */
  /*     Standard solve.                                               */
  /* ***************************************************************** */

  template <typename MATRIX, typename VECTOR, typename PLSOLVER>
  void standard_solve(model &md, gmm::iteration &iter,
		      PLSOLVER lsolver,
		      gmm::abstract_newton_line_search &ls, const MATRIX &K,
		      const VECTOR &rhs) {

    VECTOR state(md.nb_dof());
    
    md.from_variables(state); // copy the model variables in the state vector

    if (md.is_linear()) {
      md.assembly(model::BUILD_ALL);
      (*lsolver)(K, state, rhs, iter);
    }
    else {
      model_pb<MATRIX, VECTOR> mdpb(md, ls, state, rhs, K);
      classical_Newton(mdpb, iter, *lsolver);
    }

    md.to_variables(state); // copy the state vector into the model variables
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      rmodel_plsolver_type lsolver,
		      gmm::abstract_newton_line_search &ls) {
    standard_solve(md, iter, lsolver, ls, md.real_tangent_matrix(),
		   md.real_rhs());
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      cmodel_plsolver_type lsolver,
		      gmm::abstract_newton_line_search &ls) {
    standard_solve(md, iter, lsolver, ls, md.complex_tangent_matrix(),
		   md.complex_rhs());
  }


  void standard_solve(model &md, gmm::iteration &iter,
			     rmodel_plsolver_type lsolver) {
    gmm::default_newton_line_search ls(size_t(-1), 5.0/3.0,
				       1.0/1000.0, 3.0/5.0, 1.6);
    standard_solve(md, iter, lsolver, ls);
  }

  void standard_solve(model &md, gmm::iteration &iter,
			     cmodel_plsolver_type lsolver) {
    gmm::default_newton_line_search ls(size_t(-1), 5.0/3.0,
				       1.0/1000.0, 3.0/5.0, 1.6);
    standard_solve(md, iter, lsolver, ls);
  }

  void standard_solve(model &md, gmm::iteration &iter) {
    gmm::default_newton_line_search ls(size_t(-1), 5.0/3.0,
				       1.0/1000.0, 3.0/5.0, 1.6);
    if (md.is_complex())
      standard_solve(md, iter, cdefault_linear_solver(md), ls);
    else
      standard_solve(md, iter, rdefault_linear_solver(md), ls);
  }



}  /* end of namespace getfem.                                             */

