// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_Coulomb_friction.h : 
//           
// Date    : July 6, 2004.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#ifndef GETFEM_COULOMB_FRICTION_H__
#define GETFEM_COULOMB_FRICTION_H__

#include <getfem_modeling.h>

namespace getfem {

  /* ******************************************************************** */
  /*   	Unilateral contact and Coulomb friction condition brick.          */
  /*    (for conformal small displacement problems)                       */
  /* ******************************************************************** */

# define MDBRICK_COULOMB_FRICTION 434245


  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Coulomb_friction : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;
    typedef typename gmm::sub_vector_type<VECTOR *,
				 gmm::sub_interval>::vector_type SUBVECTOR;


    mdbrick_abstract<MODEL_STATE> &sub_problem;
    size_type num_fem;

    T_MATRIX BN, BT;
    VECTOR gap, threshold, WT, WN, friction_coef, RLN, RLT;
    value_type r, alpha, beta;
    size_type d;

    mesh_fem *mf_u;
    gmm::sub_interval SUBU, SUBN, SUBT;
    
    bool Tresca_version, symmetrized, contact_only, stationary;


    template<typename VEC> static void ball_projection(const VEC &x,
						       value_type radius) {
      value_type a = gmm::vect_norm2(x);
      if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
      else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a); 
    }

    template<class VEC, class VECR>
    static void ball_projection_grad_r(const VEC &x, value_type radius,
				       VECR &g) {
      value_type a = gmm::vect_norm2(x);
      if (radius > 0 && a >= radius)
	gmm::copy(gmm::scaled(x, value_type(1)/a), g);
      else gmm::clear(g);
    }

    template <class VEC, class MAT>
    static void ball_projection_grad(const VEC &x, double radius, MAT &g) {
      if (radius <= value_type(0)) { gmm::clear(g); return; }
      gmm::copy(gmm::identity_matrix(), g);
      value_type a = gmm::vect_norm2(x);
      if (a >= radius) { 
	gmm::scale(g, radius/a);
	// gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
	for (size_type i = 0; i < x.size(); ++i)
	  for (size_type j = 0; j < x.size(); ++j)
	    g(i,j) -= radius*x[i]*x[j] / (a*a*a);
      }
    }

    void precomp(MODEL_STATE &MS, size_type i0) {
      size_type i1 = this->mesh_fem_positions[num_fem];
      gmm::resize(RLN, gmm::mat_nrows(BN));
      gmm::resize(RLT, gmm::mat_nrows(BT));
      SUBU = gmm::sub_interval(i0 + i1, mf_u->nb_dof());
      SUBN = gmm::sub_interval(i0 + sub_problem.nb_dof(), gmm::mat_nrows(BN));
      SUBT = gmm::sub_interval(i0 + sub_problem.nb_dof() + gmm::mat_nrows(BN),
			       gmm::mat_nrows(BT));
      gmm::add(gmm::sub_vector(MS.state(), SUBN), gmm::scaled(gap, r), RLN);
      gmm::mult_add(BN, gmm::scaled(WN, -r*alpha), RLN);
      gmm::mult_add(BN, gmm::scaled(gmm::sub_vector(MS.state(), SUBU),
				    -r*alpha),
		    RLN);
      if (!contact_only) {
	gmm::mult(BT, gmm::scaled(WT, -r*beta),
		    gmm::sub_vector(MS.state(), SUBT), RLT);
	if (!stationary)
	  gmm::mult_add(BT, gmm::scaled(gmm::sub_vector(MS.state(), SUBU),
					-r*beta), RLT);
      }
    }

  public :
    
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      b.add(i0 + sub_problem.nb_dof(), gmm::mat_nrows(BN)+gmm::mat_nrows(BT));
    }
    
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }

    virtual size_type nb_dof(void)
    { return sub_problem.nb_dof() + gmm::mat_nrows(BN)+gmm::mat_nrows(BT); }

    inline size_type nb_contact_nodes(void) const
    { return gmm::mat_nrows(BN); }

    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0=0, bool modified=false) {

      sub_problem.compute_tangent_matrix(MS, i0, j0, modified || symmetrized);
      react(MS, i0, modified);
      precomp(MS, i0);
            
      gmm::copy(gmm::scaled(BN, -alpha),
		gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU));
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBN));
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
	if (RLN[i] >= value_type(0)) {
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(),
				gmm::sub_interval(SUBN.first()+i,1), SUBU));
	  MS.tangent_matrix()(SUBN.first()+i, SUBN.first()+i)=-value_type(1)/r;
	}
      }
    
      if (!contact_only) {
	base_matrix pg(d-1, d-1);
	base_vector vg(d-1);
	
	for (size_type i=0; i < nb_contact_nodes(); ++i) {
	  gmm::sub_interval SUBI(i*(d-1), d-1);
	  gmm::sub_interval SUBJ(SUBT.first()+i*(d-1),(d-1));
	  value_type th = Tresca_version ? threshold[i]
	    : - (MS.state())[SUBN.first()+i] * friction_coef[i];
	  
	  ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
	  if (!stationary)
	    gmm::mult(gmm::scaled(pg, -beta), 
		      gmm::sub_matrix(BT, SUBI,
				   gmm::sub_interval(0, gmm::mat_ncols(BT))),
		      gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBU));
	  
	  if (!Tresca_version) {
	    ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
	    for (size_type k = 0; k < d-1; ++k)
	      MS.tangent_matrix()(SUBT.first()+i*(d-1)+k, SUBN.first()+i)
		= - friction_coef[i] * vg[k] / r;
	  }
	  for (size_type j = 0; j < d-1; ++j) pg(j,j) -= value_type(1);
	  gmm::copy(gmm::scaled(pg,value_type(1)/r), 
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ));
	}
      }
      
      if (symmetrized) {
	T_MATRIX tmp(mf_u->nb_dof(), mf_u->nb_dof());
	
	gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BN));
	gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
						  SUBN, SUBU)), tmp);
	gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
	gmm::resize(tmp, mf_u->nb_dof(), mf_u->nb_dof());
	gmm::mult(gmm::transposed(gmm::scaled(BN,-r*alpha)),
		  gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU), tmp);
	gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
	
	if (!contact_only) {
	  gmm::mult(gmm::transposed(gmm::scaled(BT,-r*beta)),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU), tmp);
	  gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
	  gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BT));
	  gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
						    SUBT, SUBU)), tmp);
	  gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
	}
      } 
      else {
	gmm::copy(gmm::scaled(gmm::transposed(BN), value_type(-1)),
		  gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
	if (!contact_only)
	  gmm::copy(gmm::scaled(gmm::transposed(BT), value_type(-1)), 
		    gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
      }
    }
    
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      precomp(MS, i0);
      value_type c1(1);
    
      for (size_type i=0; i < nb_contact_nodes(); ++i) {
	RLN[i] = std::min(value_type(0), RLN[i]);
	if (!contact_only)
	  ball_projection(gmm::sub_vector(RLT, gmm::sub_interval(i*(d-1),d-1)),
			  Tresca_version ? threshold[i]
			  : -friction_coef[i]*(MS.state())[SUBN.first()+i]);
      }

      if (symmetrized) {
	gmm::mult_add(gmm::transposed(BN), gmm::scaled(RLN, -c1),
		      gmm::sub_vector(MS.residu(), SUBU));
	if (!contact_only)
	  gmm::mult_add(gmm::transposed(BT), gmm::scaled(RLT, -c1),
			gmm::sub_vector(MS.residu(), SUBU));
      } else {
	gmm::mult_add(gmm::transposed(BN),
		      gmm::scaled(gmm::sub_vector(MS.state(), SUBN),-c1),
		      gmm::sub_vector(MS.residu(), SUBU));
	if (!contact_only)
	  gmm::mult_add(gmm::transposed(BT),
			gmm::scaled(gmm::sub_vector(MS.state(), SUBT),-c1),
			gmm::sub_vector(MS.residu(), SUBU));
      }
      
      /* residu on LN */
      gmm::add(gmm::scaled(gmm::sub_vector(MS.state(), SUBN), -c1),
	       RLN, gmm::sub_vector(MS.residu(), SUBN));
      gmm::scale(gmm::sub_vector(MS.residu(), SUBN), c1 / r);

      /* residu on LT */
      if (!contact_only) {
	gmm::add(gmm::scaled(gmm::sub_vector(MS.state(),SUBT), -c1),
		 RLT, gmm::sub_vector(MS.residu(), SUBT));
	gmm::scale(gmm::sub_vector(MS.residu(), SUBT), c1 / r);
      }
    }

    void init(size_type nbc) {
      this->add_sub_brick(sub_problem);
      this->proper_is_linear_ = this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ = symmetrized && contact_only;
      this->update_from_context();
      mf_u = this->mesh_fems[num_fem];
      d = mf_u->linked_mesh().dim();
      r = value_type(1);
      beta = value_type(1);
      alpha = value_type(1);
      gmm::resize(BN, nbc, mf_u->nb_dof());
      gmm::resize(BT, nbc*(d-1), mf_u->nb_dof());
      gmm::resize(gap, nbc); gmm::resize(friction_coef, nbc);
      gmm::resize(threshold, nbc); gmm::resize(WT, mf_u->nb_dof());
      gmm::resize(WN, mf_u->nb_dof());
    }

    void set_stationary(bool b) { stationary = b; }
    void set_beta(value_type a) { beta = a; }
    void set_alpha(value_type a) { alpha = a; }

    void set_r(value_type r_) { r = r_; }
    value_type get_r(void) const { return r; }
    template <class VEC> void set_WT(const VEC &WT_) { gmm::copy(WT_, WT); }
    template <class VEC> void set_WN(const VEC &WN_) { gmm::copy(WN_, WN); }

    SUBVECTOR get_LN(MODEL_STATE &MS) {
      SUBN = gmm::sub_interval(this->first_index() + sub_problem.nb_dof(),
			       gmm::mat_nrows(BN));
      return gmm::sub_vector(MS.state(), SUBN);
    }

    SUBVECTOR get_LT(MODEL_STATE &MS) {
      SUBT = gmm::sub_interval(this->first_index() + sub_problem.nb_dof()
			       + gmm::mat_nrows(BN),  gmm::mat_nrows(BT));
      return gmm::sub_vector(MS.state(), SUBT);
    }

    /* contact and friction */
    template <class MAT, class VEC> mdbrick_Coulomb_friction
    (mdbrick_abstract<MODEL_STATE> &problem, const MAT &BN_, const VEC &gap_,
     scalar_type FC_, const MAT &BT_, size_type num_fem_=0)
      : sub_problem(problem), num_fem(num_fem_) {
      contact_only = false; Tresca_version = symmetrized = stationary = false;
      init(gmm::mat_nrows(BN_));
      gmm::copy(BN_, BN); gmm::copy(BT_, BT); gmm::copy(gap_, gap);
      std::fill(friction_coef.begin(), friction_coef.end(), FC_);
    }

    /* contact only.        */
    template <class MAT, class VEC> mdbrick_Coulomb_friction
    (mdbrick_abstract<MODEL_STATE> &problem, const MAT &BN_, const VEC &gap_,
     size_type num_fem_=0) : sub_problem(problem), num_fem(num_fem_) {
      contact_only = true; Tresca_version = symmetrized = stationary = false;
      init(gmm::mat_nrows(BN_));
      gmm::copy(BN_, BN); gmm::copy(gap_, gap);
    }

  };


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NONLINEAR_ELASTICITY_H__ */
