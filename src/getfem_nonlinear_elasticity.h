// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2000-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**@file getfem_nonlinear_elasticity.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date July 6, 2004.
   @brief Non-linear elasticty brick.
*/
#ifndef GETFEM_NONLINEAR_ELASTICITY_H__
#define GETFEM_NONLINEAR_ELASTICITY_H__

#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>
#include <getfem_derivatives.h>
#include <getfem_interpolation.h>

namespace getfem {


  int check_symmetry(const base_tensor &t);

  /** Base class for material law. 
      Inherit from this class to define a new law.
   */
  struct abstract_hyperelastic_law {
    mutable int uvflag;
    size_type nb_params_;
    void reset_unvalid_flag(void) const { uvflag = 0; }
    void inc_unvalid_flag(void) const { uvflag++; }
    int get_unvalid_flag(void) const { return uvflag; }
    
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params) const = 0;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params) const = 0;
	// the result of grad_sigma has to be completely symmetric.
    virtual void grad_sigma(const base_matrix &E, base_tensor &result, 
			    const base_vector &params) const = 0;
    size_type nb_params(void) const { return nb_params_; }
    abstract_hyperelastic_law() { nb_params_ = 0; }
    virtual ~abstract_hyperelastic_law() {}
    static void random_E(base_matrix &E);
    void test_derivatives(size_type N, scalar_type h,
			  const base_vector& param) const;
  };

  /** Saint-Venant / Kirchhoff hyperelastic law. 
      
      This is the linear law used in linear elasticity, it is not well 
      suited to large strain.. (the convexes may become flat) 
  */
  struct SaintVenant_Kirchhoff_hyperelastic_law : 
    public abstract_hyperelastic_law {
    /* W = lambda*0.5*trace(E)^2 + mu*tr(E^2) */
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params) const;
    /* sigma = lambda*trace(E) + 2 mu * E */
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params) const;
    SaintVenant_Kirchhoff_hyperelastic_law(void) { nb_params_ = 2; }
  };

  /** Mooney-Rivlin hyperelastic law 
      
      To be used for incompressible problems (with getfem::mdbrick_nonlinear_incomp).
  */
  struct Mooney_Rivlin_hyperelastic_law : public abstract_hyperelastic_law {
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params) const;
    Mooney_Rivlin_hyperelastic_law(void) { nb_params_ = 2; }
  };

  /** Ciarlet-Geymonat hyperelastic law ( @f$ W=~_1i_1(L) + \frac{~}{2}i_2(L) + 8ci_3(L) - \frac{~_1}{2} \textrm{log}~\textrm{det}~C @f$ )
      
      A "good" law.
  */
  struct Ciarlet_Geymonat_hyperelastic_law : public abstract_hyperelastic_law {
    // parameters are lambda=params[0], mu=params[1], gamma'(1)=params[2]
    // The parameters gamma'(1) has to verify gamma'(1) in ]-lambda/2-mu, -mu[
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params) const;
    Ciarlet_Geymonat_hyperelastic_law(void) { nb_params_ = 3; }
  };


  template<typename VECT1, typename VECT2> class elasticity_nonlinear_term 
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf;
    const VECT1 &U;
    const mesh_fem &mf_data;
    const VECT2 &PARAMS;
    size_type N;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff;
    base_matrix E, Sigma, gradU;
    base_tensor tt;
    bgeot::multi_index sizes_;
    int version;
  public:
    elasticity_nonlinear_term(const mesh_fem &mf_, const VECT1 &U_,
			      const mesh_fem &mf_data_, const VECT2 &PARAMS_,
			      const abstract_hyperelastic_law &AHL_,
			      int version_) 
      : mf(mf_), U(U_), mf_data(mf_data_), PARAMS(PARAMS_), 
	N(mf_.get_qdim()), AHL(AHL_), params(AHL_.nb_params()),
	E(N, N), Sigma(N, N), gradU(N, N), tt(N, N, N, N), sizes_(N, N, N, N),
	version(version_)
    { if (version == 1) sizes_.resize(2); }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))),
		coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));
      gmm::add(gmm::identity_matrix(), gradU);


      AHL.sigma(E, Sigma, params);

      if (version == 0) {	  
	AHL.grad_sigma(E, tt, params);	
	for (size_type n = 0; n < N; ++n)
	  for (size_type m = 0; m < N; ++m)
	    for (size_type l = 0; l < N; ++l)
	      for (size_type k = 0; k < N; ++k) {
		scalar_type aux = (k == n) ? Sigma(m,l) : 0.0;
		for (size_type j = 0; j < N; ++j)
		  for (size_type i = 0; i < N; ++i) {
		    aux += gradU(n ,j) * gradU(k, i) * tt(j, m, i, l);
		  }
		t(n, m, k, l) = aux;
	      }
      } else {
        if (gmm::lu_det(gradU) < 0) AHL.inc_unvalid_flag();

	for (size_type i = 0; i < N; ++i)
	  for (size_type j = 0; j < N; ++j) {
	    scalar_type aux(0);
	    for (size_type k = 0; k < N; ++k)
	      aux += gradU(i, k) * Sigma(k, j);
	    t(i,j) = aux;
	  }
      }
    }
    virtual void prepare(fem_interpolation_context& ctx, size_type ) {
      size_type cv = ctx.convex_num();
      size_type nb = AHL.nb_params();

      coeff.resize(mf_data.nb_dof_of_element(cv)*nb);
      for (size_type i = 0; i < mf_data.nb_dof_of_element(cv); ++i)
	for (size_type k = 0; k < nb; ++k)
	  coeff[i * nb + k] = PARAMS[mf_data.ind_dof_of_element(cv)[i]*nb+k];
      ctx.pf()->interpolation(ctx, coeff, params, nb);
    } 
    
  };


  /** 
      Tangent matrix for the non-linear elasticity 
      @ingroup asm
   */
  template<typename MAT, typename VECT1, typename VECT2> 
  void asm_nonlinear_elasticity_tangent_matrix(const MAT &K_, 
					       const mesh_im &mim, 
					       const getfem::mesh_fem &mf,
					       const VECT1 &U,
					       const getfem::mesh_fem &mf_data,
					       const VECT2 &PARAMS,
					       const abstract_hyperelastic_law &AHL,
					       const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &K = const_cast<MAT &>(K_);
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT1, VECT2>
      nterm(mf, U, mf_data, PARAMS, AHL, 0);

    getfem::generic_assembly
      /*assem("t=comp(NonLin(#1,#2).vGrad(#1).vGrad(#1)); "
	    "M(#1,#1)+= sym(t(i,j,k,l,:,i,j,:,k,l))");
      */
      assem("M(#1,#1)+=sym(comp(NonLin(#1,#2)(i,j,k,l).vGrad(#1)(:,i,j).vGrad(#1)(:,k,l)))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.assembly(rg);
  }


  /**
      @ingroup asm
  */
  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_rhs(const VECT1 &R_, 
				    const mesh_im &mim, 
				    const getfem::mesh_fem &mf,
				    const VECT2 &U,
				    const getfem::mesh_fem &mf_data,
				    const VECT3 &PARAMS,
				    const abstract_hyperelastic_law &AHL,
				    const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R = const_cast<VECT1 &>(R_);
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT2, VECT3>
      nterm(mf, U, mf_data, PARAMS, AHL, 1);

    getfem::generic_assembly
      assem("t=comp(NonLin(#1,#2).vGrad(#1)); V(#1) += t(i,j,:,i,j)");
    // comp() to be optimized ?
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }


# define MDBRICK_NONLINEAR_ELASTICITY 821357

  /** Non-linear elasticity brick  ( @f$ \int (I+\nabla u)\hat{\hat{\sigma}}:\nabla v = l(v)  @f$ ).
      
      @f$ \hat{\hat{\sigma}} @f$ is known as the second Piola-Kirchhoff stress tensor, and is given by \f[ \hat{\hat{\sigma}} = -\frac{\partial}{\partial L}W(L) \f],
      with @f$W@f$ the strain energy of the material, and @f$L=\frac{1}{2}\left(\nabla u^t\nabla u + \nabla u^t + \nabla u\right)@f$ is the Green-Lagrange strain tensor.

      This brick handle the computation of the tangent matrix and the
      right hand side for large strain problems, with hyperelastic
      material laws.
  
      @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_nonlinear_elasticity : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    const abstract_hyperelastic_law &AHL;
    const mesh_im &mim;
    const mesh_fem &mf_u;
    mdbrick_parameter<VECTOR> PARAMS_;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      asm_nonlinear_elasticity_tangent_matrix
	(gmm::sub_matrix(MS.tangent_matrix(), SUBI), mim, mf_u,
	 gmm::sub_vector(MS.state(), SUBI), params().mf(), params().get(),
	 AHL, this->mf_u.linked_mesh().get_mpi_region());
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::clear(gmm::sub_vector(MS.residual(), SUBI));
      asm_nonlinear_elasticity_rhs(gmm::sub_vector(MS.residual(), SUBI), mim,
				   mf_u, gmm::sub_vector(MS.state(), SUBI), 
				   params().mf(), params().get(), AHL,
				   this->mf_u.linked_mesh().get_mpi_region());
    }

    mdbrick_parameter<VECTOR> &params() { 
      PARAMS_.reshape(AHL.nb_params()); return PARAMS_; 
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    template <class VECTVM>
    void compute_Von_Mises_or_Tresca(MODEL_STATE &MS, const mesh_fem &mf_vm, VECTVM &VM, bool tresca) {
      unsigned int N = mf_u.linked_mesh().dim(), NP = AHL.nb_params();
      VECTOR GRAD(mf_vm.nb_dof()*N*N), PARAMS(mf_vm.nb_dof()*NP);
      
      interpolation(PARAMS_.mf(), mf_vm, PARAMS_.get(), PARAMS);
      compute_gradient(mf_u, mf_vm, get_solution(MS), GRAD);
      if (gmm::vect_size(VM) != mf_vm.nb_dof())
	DAL_THROW(dimension_error, "The vector has not the good size");
      base_matrix E(N, N), gradphi(N,N), Id(N, N), sigmahathat(N,N),
	aux(N,N), sigma(N,N);
      base_vector p(NP);
      base_vector eig(N);
      gmm::copy(gmm::identity_matrix(), Id);
      for (size_type i = 0; i < mf_vm.nb_dof(); ++i) {
	std::copy(GRAD.begin()+i*N*N, GRAD.begin()+(i+1)*N*N,
		  gradphi.begin());
	gmm::add(Id, gradphi);
	gmm::mult(gmm::transposed(gradphi), gradphi, E);
	gmm::add(gmm::scaled(Id, -scalar_type(1)), E);
	gmm::scale(E, scalar_type(1)/scalar_type(2));
	gmm::copy(gmm::sub_vector(PARAMS, gmm::sub_interval(i*NP,NP)), p);
	AHL.sigma(E, sigmahathat, p);
	gmm::mult(gradphi, sigmahathat, aux);
	gmm::mult(aux, gmm::transposed(gradphi), sigma);
	gmm::scale(sigma, scalar_type(1) / gmm::lu_det(gradphi));

	if (!tresca) {
	  /* von mises: norm(deviator(sigma)) */
	  gmm::add(gmm::scaled(Id, -gmm::mat_trace(sigma) / N), sigma);
	  
	  VM[i] = gmm::mat_euclidean_norm(sigma);
	} else {
	  /* else compute the tresca criterion */
	  gmm::symmetric_qr_algorithm(sigma, eig);
	  std::sort(eig.begin(), eig.end());
	  VM[i] = eig.back() - eig.front();
	}
      }
    }

    mdbrick_nonlinear_elasticity(const abstract_hyperelastic_law &AHL_,
				 const mesh_im &mim_,
				 const mesh_fem &mf_u_,
				 const VECTOR &PARAMS)
      : AHL(AHL_), mim(mim_), mf_u(mf_u_),
	PARAMS_("params", mf_u.linked_mesh(), this) {
      params().set(PARAMS);
      this->add_proper_mesh_fem(mf_u, MDBRICK_NONLINEAR_ELASTICITY);
      this->add_proper_mesh_im(mim);
      this->proper_is_linear_ = false;
      this->proper_is_coercive_ = this->proper_is_symmetric_ = true;
      this->force_update();
    }
  };


  /* ******************************************************************** */
  /*		Mixed nonlinear incompressible condition brick.           */
  /* ******************************************************************** */

# define MDBRICK_NONLINEAR_INCOMP 964552

  template<typename VECT1> class incomp_nonlinear_term 
    : public getfem::nonlinear_elem_term {

    const mesh_fem &mf;
    const VECT1 &U;
    size_type N;
    base_vector coeff;
    base_matrix gradPhi;
    bgeot::multi_index sizes_;
    int version; 

  public:
    incomp_nonlinear_term(const mesh_fem &mf_, const VECT1 &U_,
			      int version_) 
      : mf(mf_), U(U_),
	N(mf_.get_qdim()),
	gradPhi(N, N), sizes_(N, N),
	version(version_)
    { if (version == 1) { sizes_.resize(1); sizes_[0] = 1; } }
    const bgeot::multi_index &sizes() const { return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))),
		coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradPhi, mf.get_qdim());
      gmm::add(gmm::identity_matrix(), gradPhi);
      scalar_type det = gmm::lu_inverse(gradPhi);

      if (version != 1) {
	if (version == 2) det = sqrt(gmm::abs(det));
	for (size_type i = 0; i < N; ++i) 
	  for (size_type j = 0; j < N; ++j) {
	    t(i,j) = - det * gradPhi(j,i);
	  }
      }
      else t[0] = scalar_type(1) - det;
     
    }
  };

  /**@ingroup asm*/
  template<typename MAT1, typename MAT2, typename VECT1, typename VECT2> 
  void asm_nonlinear_incomp_tangent_matrix(const MAT1 &K_, const MAT2 &B_,
					   const mesh_im &mim,
					   const mesh_fem &mf_u,
					   const mesh_fem &mf_p,
					   const VECT1 &U, const VECT2 &P,
					   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT1 &K = const_cast<MAT1 &>(K_);
    MAT2 &B = const_cast<MAT2 &>(B_);
    if (mf_u.get_qdim() != mf_u.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    incomp_nonlinear_term<VECT1> ntermk(mf_u, U, 0);
    incomp_nonlinear_term<VECT1> ntermb(mf_u, U, 2);
    getfem::generic_assembly
      assem("P=data(#2);"
	    "t=comp(NonLin$1(#1).vGrad(#1).Base(#2));"
	    "M$2(#1,#2)+= t(i,j,:,i,j,:);"
 	    /*"w=comp(NonLin$2(#1).vGrad(#1).NonLin$2(#1).vGrad(#1).Base(#2));"
	    "M$1(#1,#1)+= w(j,i,:,j,k, m,k,:,m,i,p).P(p)"
	    "-w(i,j,:,i,j, k,l,:,k,l,p).P(p)"*/
            /*
            "w=comp(vGrad(#1).NonLin$2(#1).vGrad(#1).NonLin$2(#1).Base(#2));"
	    "M$1(#1,#1)+= w(:,j,k, j,i, :,m,i, m,k, p).P(p)"
	                "-w(:,j,i, j,i, :,m,l, m,l, p).P(p)"
            */
            "w1=comp(vGrad(#1)(:,j,k).NonLin$2(#1)(j,i).vGrad(#1)(:,m,i).NonLin$2(#1)(m,k).Base(#2)(p).P(p));"
            "w2=comp(vGrad(#1)(:,j,i).NonLin$2(#1)(j,i).vGrad(#1)(:,m,l).NonLin$2(#1)(m,l).Base(#2)(p).P(p));"
	    "M$1(#1,#1)+= w1-w2"            
	    );

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&ntermk);
    assem.push_nonlinear_term(&ntermb);
    assem.push_mat(K);
    assem.push_mat(B);
    assem.push_data(P);
    assem.assembly(rg);
  }


  /**@ingroup asm
   */
  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_incomp_rhs(const VECT1 &R_U_, const VECT1 &R_P_, 
				const mesh_im &mim,
				const getfem::mesh_fem &mf_u,
				const getfem::mesh_fem &mf_p,
				const VECT2 &U, const VECT3 &P,
				const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R_U = const_cast<VECT1 &>(R_U_);
    VECT1 &R_P = const_cast<VECT1 &>(R_P_);
    if (mf_u.get_qdim() != mf_u.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    incomp_nonlinear_term<VECT2> nterm_tg(mf_u, U, 0);
    incomp_nonlinear_term<VECT2> nterm(mf_u, U, 1);

    getfem::generic_assembly
      assem("P=data(#2); "
	    "t=comp(NonLin$1(#1).vGrad(#1).Base(#2));"
	    "V$1(#1) += t(i,j,:,i,j,k).P(k);"
	    "w=comp(NonLin$2(#1).Base(#2)); V$2(#2) += w(1,:)");
    // assem() to be optimized ?

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&nterm_tg);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R_U);
    assem.push_vec(R_P);
    assem.push_data(P);
    assem.assembly(rg);
  }

  /** Incompressible non-linear elasticity brick.
      @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_nonlinear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;
   
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    const mesh_fem &mf_p;
    size_type num_fem;

    virtual void proper_update(void) {
      this->proper_mixed_variables.clear();
      this->proper_mixed_variables.add(sub_problem.nb_dof(), mf_p.nb_dof());
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_p.nb_dof()); /* P */
      gmm::sub_interval SUBJ(i0+i1, mf_u.nb_dof());           /* U */
      T_MATRIX B(mf_u.nb_dof(), mf_p.nb_dof());
      asm_nonlinear_incomp_tangent_matrix(gmm::sub_matrix(MS.tangent_matrix(),
							  SUBJ, SUBJ), B,
					  *(this->mesh_ims[0]), mf_u, mf_p, 
					  gmm::sub_vector(MS.state(), SUBJ), 
					  gmm::sub_vector(MS.state(), SUBI),
					  mf_u.linked_mesh().get_mpi_region());
      gmm::copy(B, gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      gmm::copy(gmm::transposed(B),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }

    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0, size_type) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, mf_u.nb_dof());
      gmm::clear(gmm::sub_vector(MS.residual(), SUBI));
      asm_nonlinear_incomp_rhs(gmm::sub_vector(MS.residual(), SUBJ),
			       gmm::sub_vector(MS.residual(), SUBI),
			       *(this->mesh_ims[0]), mf_u, mf_p, 
			       gmm::sub_vector(MS.state(), SUBJ),
			       gmm::sub_vector(MS.state(), SUBI),
			       mf_u.linked_mesh().get_mpi_region());
    }

    mdbrick_nonlinear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
			     const mesh_fem &mf_p_, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), num_fem(num_fem_) {
      this->add_proper_mesh_fem(mf_p, MDBRICK_NONLINEAR_INCOMP);
      this->add_sub_brick(sub_problem);
      this->proper_is_linear_ = this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ = true;
      this->force_update();
    }
  };



}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NONLINEAR_ELASTICITY_H__ */
