/************************************************************************* */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  plasticity.h : perfect plasticity problem for isotropic      */
/*                           materials                                     */
/*                                                                         */
/* Date : June 10, 2004.                                                   */
/* Author : Marc ODUNLAMI, Rémi DELMAS                                     */
/*                                                                         */
/* *********************************************************************** */

/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Julien Pommier.                    */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef GETFEM_PLASTICITY
#define GETFEM_PLASTICITY

/* try to enable the SIGFPE if something evaluates to a Not-a-number of infinity
 * during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

#include <getfem_assembling_tensors.h>
#include <getfem_modeling.h>

namespace getfem {

  /* used to calcute the projection */
  template<typename MAT> MAT tau_m_Id(const MAT& tau){
    scalar_type trace=gmm::mat_trace(tau);
    size_type size_of_tau=gmm::mat_nrows(tau);
    MAT taumId(size_of_tau,size_of_tau);
    gmm::copy(gmm::identity_matrix(),taumId);
    gmm::scale(taumId,trace/size_of_tau);
    return taumId;
  } 

  /* used to calcute the projection */
  template<typename MAT> MAT tau_d(const MAT& tau){
    size_type size_of_tau=gmm::mat_nrows(tau);
    MAT taud(size_of_tau,size_of_tau);
    gmm::copy(gmm::scaled(tau_m_Id(tau),-1.0),taud); 
    gmm::add(tau,taud);
    return taud;
  }

  class type_proj {
    public :
      /* if flag_proj=0 il output proj will be Proj(tau)
       * if flag_proj=1 il output proj will be gradProj(tau)
       * no others values allowed for flag_proj
       */
      virtual void compute_type_proj(const base_matrix& tau,const scalar_type stress_threshold, const scalar_type TOL, base_matrix& proj,const size_type flag_proj, const size_type flag_hyp)  const = 0;
  };

  class VM_projection : public type_proj  {
    public :      

      /* on input : tau matrix, on output : the projection of tau */
      virtual void compute_type_proj(const base_matrix& tau ,const scalar_type stress_threshold,
				     const scalar_type TOL, base_matrix& proj,const size_type flag_proj, 
				     const size_type flag_hyp)  const
      {
      
	/* be sure that flag_proj has a correct value */
	if(flag_proj!=0 && flag_proj!=1) DAL_THROW(dal::failure_error, "wrong value for the projection flag, must be 0 or 1 ");
      
	/* be sure that stress_threshold has a correct value */
	if(!(stress_threshold>=0.))
	  DAL_THROW(dal::failure_error, "s is not a positive number "
		    << stress_threshold << ". You need to set "
		    << "s as a positive number");

	size_type size_of_tau=gmm::mat_nrows(tau);

 
	scalar_type normtaud;

	/* calculate tau_m*Id */
	base_matrix taumId(size_of_tau,size_of_tau);
	gmm::copy(tau_m_Id(tau),taumId); 

	// calcul du deviateur de tau, taud
	base_matrix taud(size_of_tau,size_of_tau);
	gmm::copy(tau_d(tau),taud);

	/* plane constraints */    
	if(flag_hyp==1){
	  if(size_of_tau/=2) DAL_THROW(dal::failure_error, "wrong value for CALCULATION HYPOTHESIS, must be /=1 SINCE n/=2");
	  // we form the 3D tau tensor considering that tau(3,j)=tau(i,3)=0
	  base_matrix tau_aux(3,3); gmm::clear(tau_aux);
	  gmm::copy(tau,gmm::sub_matrix(tau_aux,gmm::sub_interval(0,2)));
	  // we calculate tau deviator and its norms
	  base_matrix taud_aux(3,3);
	  gmm::copy(tau_d(tau_aux),taud_aux);             
	  normtaud=gmm::mat_euclidean_norm(taud_aux);  
	} else {
	  normtaud=gmm::mat_euclidean_norm(taud);  
	}
    
	/* we dimentionate and initialize the projection matrix or its derivative */
	switch(flag_proj) {
	case 0:
	  gmm::resize(proj,size_of_tau,size_of_tau);
	  break;
	case 1:
	  gmm::resize(proj,size_of_tau*size_of_tau, size_of_tau*size_of_tau);
	  gmm::copy(gmm::identity_matrix(),proj);
	  break;
	}
    
	if (-TOL<=normtaud && normtaud<=TOL){
	  switch(flag_proj) {
	  case 0:
	    gmm::copy(tau,proj);
	    break;
	  }
	} else {  
	  if (normtaud-stress_threshold< -TOL){
	    switch(flag_proj) {
	    case 0:
	      gmm::copy(tau,proj);
	      break;
	    }
	  }
	  else if(normtaud-stress_threshold>TOL){
	    switch(flag_proj) {
	    case 0:
	      gmm::copy(gmm::scaled(taud, stress_threshold/normtaud),proj);
	      gmm::add(taumId,proj);
	      break;
	    case 1:
	      base_matrix Igrad(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	      gmm::copy(gmm::identity_matrix(),Igrad); 
	      base_matrix Igrad2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	      gmm::clear(Igrad2);
        
	      // on construit le vecteur [1 0 0 1  0 0 1...] qui va etre copie dans certaines colonnes de Igrad(*)Igrad
	      base_vector aux(size_of_tau*size_of_tau);
	      gmm::clear(aux);
	      // boucle sur le vecteur aux pour placer les 1
	      for(size_type i=0;i<size_of_tau;++i)
		for(size_type j=0;j<size_of_tau;++j)
		  if(i==j) aux[j*size_of_tau+i]=1.;

	      // on copie aux dans les colonnes de Igrad(*)Igrad qu'il faut
	      for(size_type i=0;i<size_of_tau;++i)
		for(size_type j=0;j<size_of_tau;++j)
		  if(i==j)  gmm::copy(aux, gmm::mat_col(Igrad2,j*size_of_tau+i)); 

	      //calcul de Id_grad
	      base_matrix Id_grad(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	      gmm::copy(gmm::scaled(Igrad2,-1./size_of_tau),Id_grad);
	      gmm::add(Igrad,Id_grad);         


	      //calcul de ngrad(*)ngrad
	      base_matrix ngrad2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	      // calcul de la normale n
	      base_matrix un(size_of_tau,size_of_tau);
	      gmm::copy(gmm::scaled(taud,1./normtaud),un);  

	      // on copie la normale dasn un vecteur colonne range dans l'ordre du fortran
	      std::copy(un.begin(),un.end(),aux.begin());
	
	      // boucle sur les colonnes de ngrad(*)ngrad
	      for(size_type j=0;j<size_of_tau*size_of_tau;++j)
		gmm::copy(gmm::scaled(aux,aux[j]), gmm::mat_col(ngrad2,j));
         
 
	      //calcul final du gradient de la projection
	      // NB : proj contient l'identite au depart (proj=Igrad)
	      gmm::add(gmm::scaled(ngrad2,-1.),proj);
	      base_matrix aux2(size_of_tau*size_of_tau,size_of_tau*size_of_tau);
	      gmm::copy(gmm::scaled(proj,stress_threshold/normtaud),aux2);
	      gmm::mult(aux2,Id_grad,proj);
	      gmm::add(gmm::scaled(Igrad2,1./size_of_tau),proj);

	      break;
	    }    
	  } else {

	    switch(flag_proj) {
      
	    case 0:
	      gmm::copy(tau,proj);
	      break;

	    case 1:
	      gmm::copy(gmm::identity_matrix(),proj); 
	      break;
	    }      
	  }
	}
      }
  };

  // calculate the projection of D*e + sigma_bar_ on a point of Gauss
  class plasticity_projection_base  {
  protected:
    size_type N;
    const getfem::mesh_fem &mf;
    const std::vector<scalar_type> &U;
    const scalar_type stress_threshold, TOL;
    std::vector<scalar_type> lambda, mu;  
    bgeot::multi_index sizes_;
    const type_proj *t_proj; 
    std::vector<std::vector<scalar_type> > &sigma_bar_;

    // to save the projection
    std::vector<std::vector<scalar_type> > &saved_proj_;
 
    const size_type flag_proj, flag_hyp;
  
    bool fill_sigma_bar;
  
  public:  

    std::vector<std::vector<scalar_type> > &sigma_bar() { return sigma_bar_; }
    scalar_type &sigma_bar(size_type cv, size_type ii, int i, int j) { return sigma_bar_[cv][ii*N*N + j*N + i]; }
    std::vector<std::vector<scalar_type> > &saved_proj() { return saved_proj_; }
    scalar_type &saved_proj(size_type cv, size_type ii, int i, int j) { return saved_proj_[cv][ii*N*N + j*N + i]; }

    // constructor
    plasticity_projection_base(const getfem::mesh_fem &mf_, const std::vector<scalar_type> &U_, 
			       const scalar_type stress_threshold_, const scalar_type TOL_, 
			       const std::vector<scalar_type> lambda_, const std::vector<scalar_type> mu_, 
			       const type_proj *t_proj_, std::vector<std::vector<scalar_type> > &sigma_bar__, 
			       std::vector<std::vector<scalar_type> > &saved_proj__,
			       const size_type flag_proj_, const size_type flag_hyp_, const bool fill_sigma) :
      mf(mf_), U(U_), stress_threshold(stress_threshold_), TOL(TOL_), lambda(lambda_), mu(mu_),t_proj(t_proj_),
      sigma_bar_( sigma_bar__),saved_proj_(saved_proj__), flag_proj(flag_proj_), flag_hyp(flag_hyp_)  {
        
      fill_sigma_bar = fill_sigma;   /* always false during resolution, true when called from compute_constraints */
    
      N = mf.linked_mesh().dim(); 

      /* non-homogeneous case */
      if (gmm::vect_size(lambda_)==1) {
	lambda.resize(mf.linked_mesh().convex_index().last_true()+1);
	mu.resize(mf.linked_mesh().convex_index().last_true()+1);
	std::fill(lambda.begin(), lambda.end(), lambda_[0]);
	std::fill(mu.begin(), mu.end(), mu_[0]);
      }
      
      if (mf.get_qdim() != N) DAL_THROW(dal::failure_error, "wrong qdim for the mesh_fem");      
    
      // calculate the projection : N*N tensor
      if(flag_proj==0){
	sizes_.resize(2); sizes_[0] = N; sizes_[1] = N; 
      }
      // calculate the derivative of the projection : N*N*N*N tensor
      else if(flag_proj==1){
      	sizes_.resize(4); sizes_[0] = N; sizes_[1] = N; sizes_[2] = N; sizes_[3] = N;
      }
      // forbidden values for flag_proj
      else
       	DAL_THROW(dal::failure_error, "wrong value for the projection flag, must be 0 or 1 ");   
    
      sigma_bar_.resize(mf.linked_mesh().convex_index().last_true()+1);    
      saved_proj_.resize(mf.linked_mesh().convex_index().last_true()+1);
    }

    // sizes() method from nonlinear_elem, gives on output the size of the tensor
    const bgeot::multi_index &sizes() const { return sizes_; }

    // compute() method from nonlinear_elem, gives on output the tensor
    void compute_proj(getfem::fem_interpolation_context& ctx, bgeot::base_tensor &t){ 

      size_type cv = ctx.convex_num();

      size_type ii = ctx.ii(); 

      getfem::pfem pf = ctx.pf();
      base_vector coeff(mf.nb_dof_of_element(cv));
      base_matrix gradU(N, N), sigma(N,N);
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))), coeff);
      pf->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      scalar_type ltrace_eps = lambda[cv]*gmm::mat_trace(gradU);
    
      //if needed, we give sigma_bar[cv] and saved_proj[cv] a size equal to the number of integration points on the convexe. Seems that this is rarely needed. 
      if (sigma_bar_[cv].size() == 0){
	size_type nbgausspt = mf.int_method_of_element(cv)->approx_method()->nb_points_on_convex();
	sigma_bar_[cv].resize(N*N*nbgausspt);
	gmm::clear(sigma_bar_[cv]);
	saved_proj_[cv].resize(N*N*nbgausspt);
	gmm::clear(saved_proj_[cv]);
      }

      t.adjust_sizes(sizes_);
   
      for (size_type i=0; i < N; ++i) {
	for (size_type j=0; j < N; ++j) {
	  if(i==j)
	    sigma(i,j) = 2*mu[cv]*(gradU(i,j)+gradU(j,i))/2. + ltrace_eps + sigma_bar(cv,ii,i,j);
	  else
	    sigma(i,j) = 2*mu[cv]*(gradU(i,j)+gradU(j,i))/2. + sigma_bar(cv,ii,i,j);
	}
      }
    
      base_matrix tau_star(N,N), gradproj(N,N), proj;
 
      t_proj->compute_type_proj(sigma, stress_threshold, TOL, proj, flag_proj, flag_hyp);

      // we fill sigma_bar only when called from compute_constraints (ie, when fill_sigma_bar is set)
      if (fill_sigma_bar && flag_proj==0) {

	for (size_type i=0; i < N; ++i)
	  for (size_type j=0; j < N; ++j)
	    saved_proj(cv,ii,i,j) = proj(i,j);

	gmm::add(gmm::scaled(gradU, -mu[cv]), proj);
	gmm::add(gmm::scaled(gmm::transposed(gradU), -mu[cv]), proj);

	for (size_type i=0; i < N; ++i) proj(i,i) += ltrace_eps;    
	for (size_type i=0; i < N; ++i) {
	  for (size_type j=0; j < N; ++j) {          
	    sigma_bar(cv,ii,i,j) = proj(i,j);
	  }
	}
      }
      std::copy(proj.begin(),proj.end(),t.begin());
    }
  };

  class plasticity_projection : public getfem::nonlinear_elem_term {
    plasticity_projection_base *p;
  public:  
    plasticity_projection(plasticity_projection_base *p_) : p(p_) {}
    virtual const bgeot::multi_index &sizes() const { return p->sizes(); }

    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) { return p->compute_proj(ctx,t); }
  };


  /* ******************************************************************** */
  /*		Plasticity bricks.                                        */
  /* ******************************************************************** */  
  /* TODO :
     - contraintes planes, deformations planes  (cf flag_hyp) 
  */
  
  template<typename MODEL_STATE = standard_model_state> 
    class mdbrick_plasticity : public mdbrick_abstract<MODEL_STATE> {
      
      typedef typename MODEL_STATE::vector_type VECTOR;
      typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
      typedef typename MODEL_STATE::value_type value_type;
      mesh_fem &mf_u;
      mesh_fem &mf_data;
      VECTOR lambda_, mu_;
      bool homogeneous;
      value_type stress_threshold, TOL;

      // flag_hyp=0 : 3D case, or '2D plane'
      // other cases : to implement
      size_type N, flag_hyp;
      
      std::vector<std::vector<scalar_type> > sigma_bar;
      std::vector<std::vector<scalar_type> > saved_proj;
      
      VM_projection VM_1;
            
      public:
      
      virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
      virtual bool is_linear(void) { return false; }
      virtual bool is_coercive(void) { return false; } // means we use LU factorisation instead of conjugate gradient
      virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
      virtual size_type nb_constraints(void) { return 0; }
      
      template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
	gmm::sub_interval SUBI(this->first_index(), nb_dof());
	gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
      }
      
      void get_proj(std::vector<std::vector<scalar_type> > &p) {
	gmm::resize(p,gmm::vect_size(saved_proj));
	for (size_type cv=0;cv<gmm::vect_size(saved_proj);++cv) {
	  gmm::resize(p[cv],gmm::vect_size(saved_proj[cv]));
	  gmm::copy(saved_proj[cv], p[cv]);
	}
      }
      
      virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0, size_type = 0, bool = false) {
	
	gmm::sub_interval SUBI(i0, nb_dof());      
	T_MATRIX K;
	gmm::clear(K);
	gmm::resize(K, nb_dof(), nb_dof());
	VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
	
	
	plasticity_projection_base gradproj_base(mf_u, MS.state(), stress_threshold,
						 TOL, lambda_, mu_, &VM_1,
						 sigma_bar, saved_proj,1,flag_hyp, false);
	plasticity_projection gradproj(&gradproj_base);
	
	for (size_type i=0;i<mf_data.nb_dof();++i) {
	  lambda[i]=lambda_[0];
	  mu[i]=mu_[0];
	}
	
	/* Calculate the actual matrix */
	asm_lhs_for_plasticity(K, mf_u, mf_data, lambda, mu, &gradproj);
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      }
      
      virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0, size_type = 0) {
	gmm::sub_interval SUBI(i0, nb_dof());        
	VECTOR K;
	
	/* Initialise K correctly */
	gmm::clear(K);
	gmm::resize(K, nb_dof());
	plasticity_projection_base proj_base(mf_u, MS.state(), stress_threshold,
					     TOL, lambda_, mu_, &VM_1, sigma_bar,
					     saved_proj,0,flag_hyp, false);
	plasticity_projection proj(&proj_base);
	
	/* Calculate the actual matrix */
	asm_rhs_for_plasticity(K, mf_u, &proj);
	gmm::copy(K, gmm::sub_vector(MS.residu(), SUBI));
      }
      
      virtual void compute_constraints(MODEL_STATE &MS, size_type i0 = 0, size_type = 0) {
	gmm::sub_interval SUBI(i0, nb_dof());        
	VECTOR K;
	
	/* Initialise K correctly */
	gmm::clear(K);
	gmm::resize(K, nb_dof());
	plasticity_projection_base proj_base(mf_u, MS.state(), stress_threshold,
					     TOL, lambda_, mu_, &VM_1, sigma_bar,
					     saved_proj,0,flag_hyp, true);
	plasticity_projection proj(&proj_base);
	
	/* Calculate the actual matrix */
	asm_rhs_for_plasticity(K, mf_u, &proj);
      }


      virtual mesh_fem &main_mesh_fem(void) { return mf_u; }
      
      void set_Lame_coeff(value_type lambdai, value_type mui) {
	homogeneous = true;
	gmm::resize(lambda_, 1); lambda_[0] = lambdai;
	gmm::resize(mu_, 1); mu_[0] = mui;
      }
      
      void set_Lame_coeff(const VECTOR &lambdai, const VECTOR &mui) {
	homogeneous = false;
	gmm::resize(lambda_, mf_data.nb_dof()); gmm::copy(lambdai, lambda_);
	gmm::resize(mu_, mf_data.nb_dof()); gmm::copy(mui, mu_);
      }
      
      // constructor for a homogeneous material (constant lambda and mu)
      mdbrick_plasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
			 value_type lambdai, value_type mui, value_type stress_threshold_,
			 value_type TOL_, size_type flag_hyp_) 
      : mf_u(mf_u_), mf_data(mf_data_) {
	set_Lame_coeff(lambdai, mui);
	N = mf_data.linked_mesh().dim();
	stress_threshold = stress_threshold_;
	TOL = TOL_;
	flag_hyp = flag_hyp_;
	this->add_dependency(mf_u); this->add_dependency(mf_data);
      }
    
      // constructor for a non-homogeneous material
      mdbrick_plasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
			 const VECTOR &lambdai, const VECTOR &mui, value_type stress_treshhold_,
			 value_type TOL_, size_type flag_hyp_)
      : mf_u(mf_u_), mf_data(mf_data_) {
	set_Lame_coeff(lambdai, mui);
	N = mf_data.linked_mesh().dim();
	stress_treshold = stress_treshold_;
	TOL = TOL_;
	flag_hyp = flag_hyp_;
	this->add_dependency(mf_u); this->add_dependency(mf_data);
      }
    };

} /* namespace getfem */

#endif
