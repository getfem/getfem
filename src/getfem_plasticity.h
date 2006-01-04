// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_plasticity.h : perfect plasticity problem for isotropic
//           materials
// Date    : June 10, 2004.
// Author  : Marc ODUNLAMI, Rémi DELMAS
//
//========================================================================
//
// Copyright (C) 2002-2005 Yves Renard
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

/**@file getfem_plasticity.h
   @brief Plasticity model brick.
*/

#ifndef GETFEM_PLASTICITY__
#define GETFEM_PLASTICITY__

#include <getfem_modeling.h>

namespace getfem {

  /** Abstract projection of a stress tensor onto a set of admissible
      stress tensors.
  */
  class abstract_constraints_projection  {
  protected : 
    size_type flag_hyp;
    
  public :
      /* if flag_proj=0 the output will be Proj(tau)
       * if flag_proj=1 the output will be gradProj(tau)
       * no others values allowed for flag_proj
       */
      virtual void do_projection(const base_matrix& tau,
				     scalar_type stress_threshold,
				     base_matrix& proj,
				     size_type flag_proj) const = 0;
    abstract_constraints_projection (size_type flag_hyp_ = 0) : flag_hyp(flag_hyp_) {}
    virtual ~abstract_constraints_projection () {}
  };

  /** Von Mises projection */
  class VM_projection : public abstract_constraints_projection   {
    /* used to compute the projection */
    template<typename MAT> void tau_m_Id(const MAT& tau, MAT &taumid) const {
      scalar_type trace = gmm::mat_trace(tau);
      size_type size_of_tau = gmm::mat_nrows(tau);
      gmm::copy(gmm::identity_matrix(),taumid);
      gmm::scale(taumid, trace / scalar_type(size_of_tau));
    }
    
    /* used to compute the projection */
    template<typename MAT> void tau_d(const MAT& tau, MAT &taud) const {
      tau_m_Id(tau, taud);
      gmm::scale(taud, scalar_type(-1));
      gmm::add(tau, taud);
    }

    public :      

      /* on input : tau matrix, on output : the projection of tau */
      virtual void do_projection(const base_matrix& tau,
				     scalar_type stress_threshold,
				     base_matrix& proj,
				     size_type flag_proj)  const {
	
	/* be sure that flag_proj has a correct value */
	if(flag_proj!=0 && flag_proj!=1)
	  DAL_THROW(failure_error,
		    "wrong value for the projection flag, must be 0 or 1 ");
      
	/* be sure that stress_threshold has a correct value */
	if(!(stress_threshold>=0.))
	  DAL_THROW(failure_error, "s is not a positive number "
		    << stress_threshold << ". You need to set "
		    << "s as a positive number");
	
	size_type tausize = gmm::mat_nrows(tau);
	size_type projsize = (flag_proj == 0) ? tausize : gmm::sqr(tausize);
	scalar_type normtaud;

	/* calculate tau_m*Id */
	base_matrix taumId(tausize, tausize);
	tau_m_Id(tau, taumId); 

	// calcul du deviateur de tau, taud
	base_matrix taud(tausize,tausize);
	gmm::add(gmm::scaled(taumId, scalar_type(-1)), tau, taud);

	/* plane constraints */    
	if(flag_hyp==1){  // To be done ...
	  if(tausize/=2)
	    DAL_THROW(failure_error,
	     "wrong value for CALCULATION HYPOTHESIS, must be /=1 SINCE n/=2");
	  // we form the 3D tau tensor considering that tau(3,j)=tau(i,3)=0
	  base_matrix tau_aux(3,3); gmm::clear(tau_aux);
	  gmm::copy(tau,gmm::sub_matrix(tau_aux,gmm::sub_interval(0,2)));
	  // we calculate tau deviator and its norms
	  base_matrix taud_aux(3,3);
	  tau_d(tau_aux, taud_aux);            
	  normtaud=gmm::mat_euclidean_norm(taud_aux);  
	} 
	else normtaud=gmm::mat_euclidean_norm(taud);  
	
    
	/* dimension and initialization of proj matrix or its derivative */
	gmm::resize(proj, projsize, projsize);

	if(normtaud <= stress_threshold) {
	  switch(flag_proj) {
	  case 0: gmm::copy(tau, proj); break;
	  case 1: gmm::copy(gmm::identity_matrix(), proj); break;
	  }
	}
	else {
	  switch(flag_proj) {
	  case 0:
	    gmm::copy(gmm::scaled(taud, stress_threshold/normtaud),proj);
	    gmm::add(taumId,proj);
	    break;
	  case 1:
	    base_matrix Igrad(projsize, projsize);
	    gmm::copy(gmm::identity_matrix(),Igrad); 
	    base_matrix Igrad2(projsize, projsize);
	    
	    // build vector[1 0 0 1  0 0 1...] to be copied in certain
	    // columns of Igrad(*)Igrad
	    base_vector aux(projsize);
	    for(size_type i=0; i < tausize; ++i)
	      aux[i*tausize + i] = scalar_type(1);
	    
	    // Copy in a selection of columns of Igrad(*)Igrad
	    for(size_type i=0; i < tausize; ++i)
	      gmm::copy(aux, gmm::mat_col(Igrad2, i*tausize + i)); 
	    
	    // Compute Id_grad
	    base_matrix Id_grad(projsize, projsize);
	    scalar_type rr = scalar_type(1)/scalar_type(tausize);
	    gmm::copy(gmm::scaled(Igrad2, -rr), Id_grad);
	    gmm::add(Igrad, Id_grad);         
	    
	    
	    // Compute ngrad(*)ngrad
	    base_matrix ngrad2(projsize, projsize);
	    // Compute the normal n
	    base_matrix un(tausize, tausize);
	    gmm::copy(gmm::scaled(taud, 1./normtaud),un);  
	    
	    // Copy of the normal in a column vector in the Fortran order
	    std::copy(un.begin(), un.end(), aux.begin());
	    
	    // Loop on the columns of ngrad(*)ngrad
	    for(size_type j=0; j < projsize; ++j)
	      gmm::copy(gmm::scaled(aux,aux[j]), gmm::mat_col(ngrad2,j));
	    
	    
	    // Final computation of the projection gradient
	    gmm::copy(gmm::identity_matrix(), proj);
	    gmm::add(gmm::scaled(ngrad2, scalar_type(-1)), proj);
	    base_matrix aux2(projsize, projsize);
	    gmm::copy(gmm::scaled(proj, stress_threshold/normtaud),aux2);
	    gmm::mult(aux2,Id_grad,proj);
	    gmm::add(gmm::scaled(Igrad2, rr),proj);
	    break;
	  }    
	}
      }
    VM_projection(size_type flag_hyp_ = 0) : abstract_constraints_projection (flag_hyp_) {}
  };
  
  /** Compute the projection of D*e + sigma_bar_ on a Gauss point. */
  class plasticity_projection : public nonlinear_elem_term {
  protected:
    base_vector params, coeff;
    size_type N;
    const mesh_im &mim;
    const mesh_fem &mf;
    const mesh_fem &mf_data;
    const std::vector<scalar_type> &U;
    const std::vector<scalar_type> &stress_threshold;
    const std::vector<scalar_type> &lambda, &mu;  
    bgeot::multi_index sizes_;
    const abstract_constraints_projection  *t_proj; 
    std::vector<std::vector<scalar_type> > &sigma_bar_;

    // to save the projection
    std::vector<std::vector<scalar_type> > &saved_proj_;
 
    const size_type flag_proj;
    bool fill_sigma_bar;
  
  public:  

    std::vector<std::vector<scalar_type> > &sigma_bar() { return sigma_bar_; }
    scalar_type &sigma_bar(size_type cv, size_type ii, int i, int j)
    { return sigma_bar_[cv][ii*N*N + j*N + i]; }
    std::vector<std::vector<scalar_type> > &saved_proj()
    { return saved_proj_; }
    scalar_type &saved_proj(size_type cv, size_type ii, int i, int j)
    { return saved_proj_[cv][ii*N*N + j*N + i]; }

    // constructor
    plasticity_projection(const mesh_im &mim_,
			  const mesh_fem &mf_,
			  const mesh_fem &mf_data_,
			  const std::vector<scalar_type> &U_, 
			  const std::vector<scalar_type> &stress_threshold_, 
			  const std::vector<scalar_type> &lambda_,
			  const std::vector<scalar_type> &mu_, 
			  const abstract_constraints_projection  *t_proj_,
			  std::vector<std::vector<scalar_type> > &sigma_bar__, 
			  std::vector<std::vector<scalar_type> > &saved_proj__,
			  const size_type flag_proj_,
			  const bool fill_sigma) :
      params(3), N(mf_.linked_mesh().dim()), mim(mim_),
      mf(mf_), mf_data(mf_data_),
      U(U_),  stress_threshold(stress_threshold_),
      lambda(lambda_), mu(mu_), sizes_(N, N, N, N), t_proj(t_proj_),
      sigma_bar_( sigma_bar__),saved_proj_(saved_proj__),
      flag_proj(flag_proj_)  {
      
      fill_sigma_bar = fill_sigma;   /* always false during resolution, */
      /*                      true when called from compute_constraints */

      if (mf.get_qdim() != N)
	DAL_THROW(failure_error, "wrong qdim for the mesh_fem");      
      if (flag_proj==0) sizes_.resize(2);
    
      sigma_bar_.resize(mf.linked_mesh().convex_index().last_true()+1);    
      saved_proj_.resize(mf.linked_mesh().convex_index().last_true()+1);
    }

    const bgeot::multi_index &sizes() const { return sizes_; }

    // compute() method from nonlinear_elem, gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
			 bgeot::base_tensor &t){ 

      size_type cv = ctx.convex_num();

      size_type ii = ctx.ii(); 

      pfem pf = ctx.pf();
      coeff.resize(mf.nb_dof_of_element(cv));
      base_matrix gradU(N, N), sigma(N,N);
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))),
		coeff);
      pf->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      
      scalar_type ltrace_eps;
      ltrace_eps = params[0]*gmm::mat_trace(gradU);

      // if needed, we give sigma_bar[cv] and saved_proj[cv] a size equal
      // to the number of integration points on the convexe. Seems that
      // this is rarely needed. 
      if (sigma_bar_[cv].size() == 0){
	size_type nbgausspt = mim.int_method_of_element(cv)
	  ->approx_method()->nb_points_on_convex();
	sigma_bar_[cv].resize(N*N*nbgausspt);
	gmm::clear(sigma_bar_[cv]);
	saved_proj_[cv].resize(N*N*nbgausspt);
	gmm::clear(saved_proj_[cv]);
      }

      t.adjust_sizes(sizes_);
   
      for (size_type i=0; i < N; ++i) {
	for (size_type j=0; j < N; ++j) {
	  sigma(i,j) = 2*params[1]*(gradU(i,j)+gradU(j,i))/2.
	    + sigma_bar(cv,ii,i,j);
	  if(i==j) sigma(i,i) += ltrace_eps;
	}
      }
    
      base_matrix tau_star(N,N), gradproj(N,N), proj;
      t_proj->do_projection(sigma, params[2], proj, flag_proj);

      // we fill sigma_bar only when called from compute_constraints
      // (ie, when fill_sigma_bar is set)
      if (fill_sigma_bar && flag_proj==0) {

	for (size_type i=0; i < N; ++i)
	  for (size_type j=0; j < N; ++j)
	    saved_proj(cv,ii,i,j) = proj(i,j);

	
	gmm::add(gmm::scaled(gradU, -params[1]), proj);
	gmm::add(gmm::scaled(gmm::transposed(gradU), -params[1]), proj);

	for (size_type i=0; i < N; ++i) {
	  proj(i,i) += ltrace_eps;
	  for (size_type j=0; j < N; ++j)
	    sigma_bar(cv,ii,i,j) = proj(i,j);
	}
      }
      std::copy(proj.begin(),proj.end(), t.begin());
    }
    
    virtual void prepare(fem_interpolation_context& ctx, size_type ) {
      size_type cv = ctx.convex_num();

      coeff.resize(mf_data.nb_dof_of_element(cv)*3);
      for (size_type i = 0; i < mf_data.nb_dof_of_element(cv); ++i) {
	coeff[i * 3] = lambda[mf_data.ind_dof_of_element(cv)[i]];
	coeff[i * 3+1] = mu[mf_data.ind_dof_of_element(cv)[i]];
	coeff[i * 3+2] = stress_threshold[mf_data.ind_dof_of_element(cv)[i]];
      }
      ctx.pf()->interpolation(ctx, coeff, params, 3);
    } 

  };


  /** 
     Right hand side vector for plasticity 
      @ingroup asm
  */
  template<typename VECT> 
  void asm_rhs_for_plasticity(VECT &V, 
			      const mesh_im &mim, 
			      const mesh_fem &mf,
			      const mesh_fem &mfdata,
			      nonlinear_elem_term *plast,
			      const mesh_region &rg = mesh_region::all_convexes())
  {
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("t=comp(NonLin(#1,#2).vGrad(#1));"
				   "e=(t{:,:,:,4,5}+t{:,:,:,5,4})/2;"
				   "V(#1) += e(i,j,:,i,j)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_nonlinear_term(plast);
    assem.push_vec(V);
    assem.assembly(rg);
  }

  /** 
      Left hand side matrix for plasticity
      @ingroup asm
  */
  template<typename MAT,typename VECT> 
  void asm_lhs_for_plasticity(MAT &H, 
			      const mesh_im &mim, 
			      const mesh_fem &mf,
			      const mesh_fem &mfdata,
			      const VECT &LAMBDA, const VECT &MU,
			      nonlinear_elem_term *gradplast,
			      const mesh_region &rg = mesh_region::all_convexes())
  {
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
				   "t=comp(NonLin(#1,#2).vGrad(#1).vGrad(#1).Base(#2));"
				   "e=(t{:,:,:,:,:,6,7,:,9,10,:}+t{:,:,:,:,:,7,6,:,9,10,:}+t{:,:,:,:,:,6,7,:,10,9,:}+t{:,:,:,:,:,7,6,:,10,9,:})/4;"
				   "M(#1,#1)+= sym(2*e(i,j,k,l,:,k,l,:,i,j,m).mu(m)+e(i,j,k,k,:,l,l,:,i,j,m).lambda(m))");
    // comp()  to be optimized !!
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_nonlinear_term(gradplast);
    assem.push_mat(H);
    assem.assembly(rg);
  }


  /* ******************************************************************** */
  /*		Plasticity bricks.                                        */
  /* ******************************************************************** */  
# define MDBRICK_SMALL_DEF_PLASTICITY 556433
  
  /**
     Plasticity brick (small deformations, quasi-static).

     @todo plane strain, plane stress  (cf flag_hyp).
     @ingroup bricks
   */

  template<typename MODEL_STATE = standard_model_state> 
    class mdbrick_plasticity : public mdbrick_abstract<MODEL_STATE> {
      
      TYPEDEF_MODEL_STATE_TYPES;

      const mesh_im &mim;
      const mesh_fem &mf_u;
      mdbrick_parameter<VECTOR> lambda_, mu_;
      mdbrick_parameter<VECTOR> stress_threshold_;

      // flag_hyp=0 : 3D case, or '2D plane'
      // other cases : to implement
      size_type N;
      
      std::vector<std::vector<scalar_type> > sigma_bar;
      std::vector<std::vector<scalar_type> > saved_proj;
      
      const abstract_constraints_projection  &t_proj;

      void proper_update(void) {}

      public:
      /** accessor for the lambda lame coefficient */
      mdbrick_parameter<VECTOR> &lambda(void) { return lambda_; }
      const mdbrick_parameter<VECTOR> &lambda(void) const { return lambda_; }
      /** accessor for the mu lame coefficient */
      mdbrick_parameter<VECTOR> &mu(void) { return mu_; }
      const mdbrick_parameter<VECTOR> &mu(void) const { return mu_; }
      /** accessor for the stresh threshold */
      mdbrick_parameter<VECTOR> &stress_threshold(void) { return stress_threshold_; }
      const mdbrick_parameter<VECTOR> &stress_threshold(void) const { return stress_threshold_; }
      
      SUBVECTOR get_solution(MODEL_STATE &MS) {
	gmm::sub_interval SUBU(this->first_index(), mf_u.nb_dof());
	return gmm::sub_vector(MS.state(), SUBU);
      }
      
      void get_proj(std::vector<std::vector<scalar_type> > &p) {
	gmm::resize(p, gmm::vect_size(saved_proj));
	for (size_type cv=0; cv < gmm::vect_size(saved_proj); ++cv) {
	  gmm::resize(p[cv], gmm::vect_size(saved_proj[cv]));
	  gmm::copy(saved_proj[cv], p[cv]);
	}
      }
      
      virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					     size_type) {
	gmm::sub_interval SUBI(i0, mf_u.nb_dof());      
	T_MATRIX K(mf_u.nb_dof(), mf_u.nb_dof());

	plasticity_projection gradproj(mim, mf_u, lambda_.mf(), MS.state(),
				       stress_threshold_.get(), lambda_.get(), mu_.get(), &t_proj,
				       sigma_bar, saved_proj, 1, false);
	
	/* Calculate the actual matrix */
	DAL_TRACE2("Assembling plasticity tangent matrix");
	asm_lhs_for_plasticity(K, mim, mf_u, lambda_.mf(), lambda_.get(), mu_.get(), &gradproj);
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      }
      
      virtual void do_compute_residual(MODEL_STATE &MS, size_type i0, size_type) {
	gmm::sub_interval SUBI(i0, mf_u.nb_dof());        
	VECTOR K(mf_u.nb_dof());
	plasticity_projection proj(mim, mf_u, lambda_.mf(), MS.state(),
				   stress_threshold_.get(),
				   lambda_.get(), mu_.get(), &t_proj, sigma_bar,
				   saved_proj, 0, false);
	
	/* Calculate the actual vector */
	DAL_TRACE2("Assembling plasticity rhs");
	asm_rhs_for_plasticity(K, mim, mf_u, lambda_.mf(), &proj);
	gmm::copy(K, gmm::sub_vector(MS.residual(), SUBI));
      }
      
      void compute_constraints(MODEL_STATE &MS) {
	VECTOR K(mf_u.nb_dof());
	
	plasticity_projection proj(mim, mf_u, lambda_.mf(), MS.state(),
				   stress_threshold_.get(),
				   lambda_.get(), mu_.get(), &t_proj, sigma_bar,
				   saved_proj, 0, true);
	
	/* Calculate the actual vector */
	DAL_TRACE2("Assembling plasticity rhs");
	asm_rhs_for_plasticity(K, mim, mf_u, lambda_.mf(), &proj);
      }

      /** constructor for a homogeneous material.
	  (non homogeneous lamba, mu and stress threshold can be set afterwards).

          @param lambdai 
	  @param mu the Lame coefficients
	  @param stress_th the stress threshold
	  @param t_proj the projection object (projection on the admissible constraints set).
      */
      mdbrick_plasticity(const mesh_im &mim_, const mesh_fem &mf_u_,
			 value_type lambdai, value_type mui,
			 value_type stress_th,
			 const abstract_constraints_projection &t_proj_) 
	: mim(mim_), mf_u(mf_u_), lambda_("lambda", mf_u_.linked_mesh(), this),
	  mu_("mu", mf_u_.linked_mesh(), this),
	  stress_threshold_("stress_threshold", mf_u_.linked_mesh(), this),
	  t_proj(t_proj_) {
	lambda_.set(lambdai); mu_.set(mui); stress_threshold_.set(stress_th);
	this->add_proper_mesh_im(mim);
	this->add_proper_mesh_fem(mf_u, MDBRICK_SMALL_DEF_PLASTICITY);
	this->proper_is_coercive_ = this->proper_is_linear_ = false;
	this->proper_is_symmetric_ = true;
	N = mf_u.linked_mesh().dim();
	this->force_update();
      }

    };

} /* namespace getfem */

#endif
