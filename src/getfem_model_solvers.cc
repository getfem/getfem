/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

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

  double default_newton_line_search::next_try(void) {
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
      { count_pat = 0.; return true; }      
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
  /*     Intermediary structure for Newton algorithms.                 */
  /* ***************************************************************** */

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
      R res = gmm::vect_norm2(state), ampl = res * R(1E-20);
      if (res == R(0)) ampl = R(1E-50);
      std::vector<R> V(gmm::vect_size(state));
      gmm::fill_random(V);
      gmm::add(gmm::scaled(V, ampl), state);
    }

    const VECTOR &residual(void) { return (is_reduced ? rhsr : rhs); }

    R residual_norm(void) { // A norm1 seems to be better than a norm2
                            // at least for contact problems.
      if (is_reduced)
	return gmm::vect_norm1(rhsr)/R(gmm::vect_size(rhsr));
      else
	return gmm::vect_norm1(rhs)/R(gmm::vect_size(rhs));
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
      R alpha(1), res, res_init, R0;
	  
      res_init = res = compute_res(false);      
      // cout << "first residual = " << residual() << endl << endl;
      R0 = (is_reduced ? gmm::real(gmm::vect_sp(dr, rhsr))
	               : gmm::real(gmm::vect_sp(dr, rhs)));

      // store the initial residual for the reprojection step.
//       size_type N = gmm::vect_size(residual());
//       VECTOR first_rhs(N); gmm::copy(residual(), first_rhs); // à éliminer


      // Compute the second derivative at alpha = 0 (value 2*a)
      // not very effective ... precision problem ?
      
//       R EPS = 1e-6;
//       gmm::add(gmm::sub_vector(stateinit, I), gmm::scaled(dr, EPS),
//                                               gmm::sub_vector(state, I));
//       R res2 = compute_res();
      
//       R a = (res2 + EPS * res - res) / gmm::sqr(EPS);
//       R delta = gmm::sqr(res)-R(4)*a*res;
//       R alpha_bis = R(1);
//       if (delta > 0) {
// 	R l1 = (res - sqrt(delta)) / (R(2)*a);
// 	R l2 = (res - sqrt(delta)) / (R(2)*a);
// 	if (l1 > 0) alpha_bis = std::min(alpha_bis, l1);
// 	if (l2 > 0) alpha_bis = std::min(alpha_bis, l2);
//       }
//       if (a > 0)
// 	alpha_bis = std::min(alpha_bis, res / (R(2)*a));
//       cout << "alpha_bis = " << alpha_bis << endl;
      
//       gmm::add(gmm::sub_vector(stateinit, I), gmm::scaled(dr, alpha_bis),
//                                               gmm::sub_vector(state, I));

//       R res3 = compute_res();
//       cout << "res3 = " << res3 << endl;

//       cout.precision(16);
//       cout << "res = " << res << " res2 = " << res2 << endl;
//       cout << "a = " << (res2 + EPS * res - res) / gmm::sqr(EPS) << endl;
//       cout << "aa = " << (res2 + EPS * res - res) << endl;

      
      ls.init_search(res, iter.get_iteration(), R0);
      do {
	alpha = ls.next_try();
	gmm::add(gmm::sub_vector(stateinit, I), gmm::scaled(dr, alpha),
		 gmm::sub_vector(state, I));
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

//     if (1) { // reprojection step ... only for non reduced case ...
//     GMM_ASSERT1(!is_reduced, "to be done ...");
// 	// detection des ddls "coupables"
// 	R mean = R(0);
// 	for (size_type i = 0; i < N; ++i) mean += gmm::abs(residual()[i]);
// 	mean /= R(N);
// 	cout << "mean = " << mean << endl << "selected : ";
	
// 	VECTOR new_dr(N);
// 	size_type Ndof = 0;
// 	for (size_type i = 0; i < N; ++i)
// 	  if (gmm::abs(residual()[i])
// 	      > std::max(gmm::abs(first_rhs[i]), mean) * R(4)) {
// 	    cout << i << " ";
// 	    new_dr[i] = dr[i];
// 	    ++Ndof;
// 	  }
// 	cout << endl;
	
// 	if (Ndof > 0) {

// 	  R EPS = 1e-8;
// 	  R THR = 0.8;
// 	  R mini_alpha = alpha * THR;

// 	  cout << "performing post-correction" << endl;
// 	  cout << "res = " << res << " res_init = " << res_init << endl;
	  

// 	  // divided difference
// 	  gmm::add(gmm::scaled(new_dr, -EPS), state);
// 	  R res2 = compute_res();
// 	  cout << "pente = " << (res - res2) / EPS << endl;
// 	  R beta = (res - res2) / EPS;
// 	  mini_alpha = (res2 - THR*std::min(res_init, res))/((beta != R(0)) ? beta : R(1));

// 	  cout << "mini_alpha = " << mini_alpha << endl;
// 	  if (mini_alpha > alpha * THR)
// 	    mini_alpha = alpha * THR;

// 	  if (mini_alpha > 0) {
	  
// 	    gmm::add(gmm::scaled(new_dr, -mini_alpha+EPS), state);
// 	    R res4 = compute_res(), res3;
// 	    cout << "res4 = " << res4 << endl;
// 	    size_type mini_it = 0;
// 	    do {
// 	      res3 = res4;
// 	      mini_alpha /= R(2); mini_it++;
// 	      gmm::add(gmm::scaled(new_dr, mini_alpha), state);
// 	      res4 = compute_res();
// 	      cout << "res4 = " << res4 << endl;
// 	    } while (res4 < res3 && mini_it <= 6);
	    
// 	    if (res4 < res3) { res3 = res4; }
// 	    else gmm::add(gmm::scaled(new_dr, -mini_alpha), state);
	    
// 	    if (res3 > res) { // cancel correction
// 	      gmm::add(stateinit, gmm::scaled(dr, alpha), state);
// 	    }
// 	  }
// 	  else gmm::add(stateinit, gmm::scaled(dr, alpha), state);
// 	  res = compute_res();
// 	}
//       }
      return alpha;
    }

    model_pb(model &m, abstract_newton_line_search &ls_, VECTOR &st,
	     const VECTOR &rhs_, const MATRIX &K_, bool reduced_,
	     std::vector<size_type> &sind_,
	     bool with_pseudo_pot = false)
      : md(m), is_reduced(reduced_), sind(sind_), I(sind_), ls(ls_), state(st),
	rhs(rhs_), K(K_), with_pseudo_potential(with_pseudo_pot) {}

  };

  /* ***************************************************************** */
  /*     Standard solve.                                               */
  /* ***************************************************************** */

  template <typename MATRIX, typename VECTOR, typename PLSOLVER>
  void standard_solve(model &md, gmm::iteration &iter,
		      PLSOLVER lsolver,
		      abstract_newton_line_search &ls, const MATRIX &K,
		      const VECTOR &rhs, bool with_pseudo_potential = false) {

    VECTOR state(md.nb_dof());
    std::vector<size_type> sind;
    bool is_reduced;
    
    md.from_variables(state); // copy the model variables in the state vector

    is_reduced = md.build_reduced_index(sind);  // sub index of the dof to 
                              //  be solved in case of disabled variables

    if (md.is_linear()) {
      md.assembly(model::BUILD_ALL);
      
      if (is_reduced) {
	gmm::sub_index I(sind);
	MATRIX Kr(sind.size(), sind.size());
	VECTOR rhsr(sind.size()), stater(sind.size());
	gmm::copy(gmm::sub_matrix(K, I, I), Kr);
	gmm::copy(gmm::sub_vector(state, I), stater);
	gmm::copy(gmm::sub_vector(rhs, I), rhsr);
	(*lsolver)(Kr, stater, rhsr, iter);
	gmm::copy(stater, gmm::sub_vector(state, I));
      } else {
	(*lsolver)(K, state, rhs, iter);
      }
    }
    else {
      model_pb<MATRIX, VECTOR> mdpb(md, ls, state, rhs, K, is_reduced, sind,
				    with_pseudo_potential);
      classical_Newton(mdpb, iter, *lsolver);
    }

    md.to_variables(state); // copy the state vector into the model variables
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      rmodel_plsolver_type lsolver,
		      abstract_newton_line_search &ls,
		      bool with_pseudo_potential) {
    standard_solve(md, iter, lsolver, ls, md.real_tangent_matrix(),
		   md.real_rhs(), with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      cmodel_plsolver_type lsolver,
		      abstract_newton_line_search &ls,
		      bool with_pseudo_potential) {
    standard_solve(md, iter, lsolver, ls, md.complex_tangent_matrix(),
		   md.complex_rhs(), with_pseudo_potential);
  }


  void standard_solve(model &md, gmm::iteration &iter,
			     rmodel_plsolver_type lsolver,
		      bool with_pseudo_potential) {
    default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls, with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
			     cmodel_plsolver_type lsolver,
		      bool with_pseudo_potential) {
    default_newton_line_search ls;
    standard_solve(md, iter, lsolver, ls, with_pseudo_potential);
  }

  void standard_solve(model &md, gmm::iteration &iter,
		      bool with_pseudo_potential) {
    getfem::default_newton_line_search ls;
    if (md.is_complex())
      standard_solve(md, iter, cdefault_linear_solver(md), ls,
		     with_pseudo_potential);
    else
      standard_solve(md, iter, rdefault_linear_solver(md), ls,
		     with_pseudo_potential);
  }



}  /* end of namespace getfem.                                             */

