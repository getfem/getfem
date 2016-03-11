/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2000-2015 Yves Renard, Julien Pommier

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

/**@file getfem_nonlinear_elasticity.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date July 6, 2004.
   @brief Non-linear elasticty and incompressibility bricks.
*/
#ifndef GETFEM_NONLINEAR_ELASTICITY_H__
#define GETFEM_NONLINEAR_ELASTICITY_H__

#include "getfem_models.h"
#include "getfem_assembling_tensors.h"
#include "getfem_derivatives.h"
#include "getfem_interpolation.h"
#include "getfem_generic_assembly.h"
#include "gmm/gmm_inoutput.h"

namespace getfem {


  int check_symmetry(const base_tensor &t);

  class abstract_hyperelastic_law;
  typedef std::shared_ptr<const abstract_hyperelastic_law> phyperelastic_law;

  /** Base class for material law. 
      Inherit from this class to define a new law.
  */
  class abstract_hyperelastic_law {
  public:
    mutable int uvflag;
    size_type nb_params_;
    phyperelastic_law pl; /* optional reference */
    void reset_unvalid_flag(void) const { uvflag = 0; }
    void inc_unvalid_flag(void) const { uvflag++; }
    int get_unvalid_flag(void) const { return uvflag; }
    std::string adapted_tangent_term_assembly_fem_data; // should be filled
    std::string adapted_tangent_term_assembly_cte_data; // to replace the
                                                        // default assembly
    
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params,
				      scalar_type det_trans) const = 0;
	/**True Cauchy stress (for Updated Lagrangian formulation)*/
	virtual void cauchy_updated_lagrangian(const base_matrix& F, 
			const base_matrix &E, 
			base_matrix &cauchy_stress,
			const base_vector &params,
			scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params,
		       scalar_type det_trans) const = 0;
    // the result of grad_sigma has to be completely symmetric.
    virtual void grad_sigma(const base_matrix &E, base_tensor &result, 
		 const base_vector &params, scalar_type det_trans) const = 0;

	/**cauchy-truesdel tangent moduli, used in updated lagrangian*/
	virtual void grad_sigma_updated_lagrangian(const base_matrix& F, 
			const base_matrix& E,
			const base_vector &params,
			scalar_type det_trans,
			base_tensor &grad_sigma_ul)const;

    size_type nb_params(void) const { return nb_params_; }
    abstract_hyperelastic_law() { nb_params_ = 0; }
    virtual ~abstract_hyperelastic_law() {}
    static void random_E(base_matrix &E);
    void test_derivatives(size_type N, scalar_type h,
			  const base_vector& param) const;
  };

  /** Saint-Venant Kirchhoff hyperelastic law. 
      
      This is the linear law used in linear elasticity, it is not well 
      suited to large strain. (the convexes may become flat) 
  */
  struct SaintVenant_Kirchhoff_hyperelastic_law : 
    public abstract_hyperelastic_law {
    /* W = lambda*0.5*trace(E)^2 + mu*tr(E^2) */
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params,
				      scalar_type det_trans) const;
    /* sigma = lambda*trace(E) + 2 mu * E */
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
		      const base_vector &params, scalar_type det_trans) const;
	virtual void grad_sigma_updated_lagrangian(const base_matrix& F, 
				const base_matrix& E,
				const base_vector &params,
				scalar_type det_trans,
				base_tensor &grad_sigma_ul)const;
    SaintVenant_Kirchhoff_hyperelastic_law(void);
  };


  /** Linear law for a membrane (plane stress), orthotropic material
      caracterized by Ex, Ey, vYX and G, with optional orthotropic prestresses.
      due to Jean-Yves Heddebaut <jyhed@alpino.be>
  */
  struct membrane_elastic_law : 
    public abstract_hyperelastic_law {
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params,
				      scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
		       const base_vector &params, scalar_type det_trans) const;
    membrane_elastic_law(void) { nb_params_ = 6; }
  };


  /** Mooney-Rivlin hyperelastic law 

      To be used for compressible and incompressible problems.
      Following combinations are possible:
        not compressible, not neohookean (default): 2 parameters (C1,C2)
        not compressible, neohookean: 1 parameter (C1)
        compressible, not neohookean: 3 parameters (C1,C2,D1)
        compressible, neohookean: 2 parameters (C1,D1)
  */
  struct Mooney_Rivlin_hyperelastic_law : public abstract_hyperelastic_law {
    const bool compressible, neohookean;
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params, scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params, scalar_type det_trans) const;
    Mooney_Rivlin_hyperelastic_law(bool compressible_=false,
                                   bool neohookean_=false);
  };

  /** Neo-Hookean hyperelastic law variants

      To be used for compressible problems.
      For incompressible ones Mooney-Rivlin hyperelastic law can be used.
  */
  struct Neo_Hookean_hyperelastic_law : public abstract_hyperelastic_law {
    const bool bonet;
    virtual scalar_type strain_energy(const base_matrix &E,
                                      const base_vector &params, scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
                       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
                            const base_vector &params, scalar_type det_trans) const;
    Neo_Hookean_hyperelastic_law(bool bonet_=true);
  };

  /** Blatz_Ko hyperelastic law 
      
      To be used for compressible problems.
  */
  struct generalized_Blatz_Ko_hyperelastic_law : public abstract_hyperelastic_law {
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params, scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params, scalar_type det_trans) const;
    generalized_Blatz_Ko_hyperelastic_law(void);
  };


  /** Ciarlet-Geymonat hyperelastic law
   */
  struct Ciarlet_Geymonat_hyperelastic_law : public abstract_hyperelastic_law {
    // parameters are lambda=params[0], mu=params[1], a=params[2]
    // The parameter a has to verify a in ]0, mu/2[
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params, scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params, scalar_type det_trans) const;
    Ciarlet_Geymonat_hyperelastic_law(void) { nb_params_ = 3; }
  };

  /** Plane strain hyperelastic law (takes another law as a parameter)
   */
  struct plane_strain_hyperelastic_law : public abstract_hyperelastic_law {
    virtual scalar_type strain_energy(const base_matrix &E,
				      const base_vector &params, scalar_type det_trans) const;
    virtual void sigma(const base_matrix &E, base_matrix &result,
		       const base_vector &params, scalar_type det_trans) const;
    virtual void grad_sigma(const base_matrix &E, base_tensor &result,
			    const base_vector &params, scalar_type det_trans) const;
    plane_strain_hyperelastic_law(const phyperelastic_law &pl_)
    { pl = pl_; nb_params_ = pl->nb_params(); }
  };




  template<typename VECT1, typename VECT2> class elasticity_nonlinear_term 
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf;
    std::vector<scalar_type> U;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    size_type N;
    size_type NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff;
    base_matrix E, Sigma, gradU;
    base_tensor tt;
    bgeot::multi_index sizes_;
    int version;
  public:
    elasticity_nonlinear_term(const mesh_fem &mf_, const VECT1 &U_,
			      const mesh_fem *mf_data_, const VECT2 &PARAMS_,
			      const abstract_hyperelastic_law &AHL_,
			      int version_) 
      : mf(mf_), U(mf_.nb_basic_dof()), mf_data(mf_data_), PARAMS(PARAMS_), 
	N(mf_.linked_mesh().dim()), NFem(mf_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), E(N, N), Sigma(N, N), gradU(NFem, N),
	tt(N, N, N, N), sizes_(NFem, N, NFem, N),
	version(version_) {
      switch (version) {
      case 0 : break; // tangent term
      case 1 : sizes_.resize(2); break; // rhs
      case 2 : sizes_.resize(1); sizes_[0] = 1; break; // strain energy
      case 3 : sizes_.resize(2); break; // Id + grad(u)
      }

      mf.extend_vector(U_, U);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }

    const bgeot::multi_index &sizes(size_type) const {  return sizes_; }

    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      slice_vector_on_basic_dof_of_element(mf, U, cv, coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      
      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha) += scalar_type(1);

      if (version == 3) {
	for (size_type n = 0; n < NFem; ++n)
	  for (size_type m = 0; m < N; ++m)
	    t(n,m) = gradU(n,m);
	return;
      }

      gmm::mult(gmm::transposed(gradU), gradU, E);
      for (unsigned int alpha = 0; alpha < N; ++alpha)
	E(alpha, alpha) -= scalar_type(1);
      gmm::scale(E, scalar_type(0.5));

      scalar_type det_trans = gmm::lu_det(gradU);

      if (version == 2) {
	t[0] = AHL.strain_energy(E, params, det_trans);
	return;
      }

      AHL.sigma(E, Sigma, params, det_trans);

      if (version == 0) {	  
	AHL.grad_sigma(E, tt, params, det_trans);	
	for (size_type n = 0; n < NFem; ++n)
	  for (size_type m = 0; m < N; ++m)
	    for (size_type l = 0; l < N; ++l)
	      for (size_type k = 0; k < NFem; ++k) {
		scalar_type aux = (k == n) ? Sigma(m,l) : 0.0;
		for (size_type j = 0; j < N; ++j)
		  for (size_type i = 0; i < N; ++i) {
		    aux += gradU(n ,j) * gradU(k, i) * tt(j, m, i, l);
		  }
		t(n, m, k, l) = aux;
	      }
      } else {
        if (det_trans < scalar_type(0)) AHL.inc_unvalid_flag();
	for (size_type i = 0; i < NFem; ++i)
	  for (size_type j = 0; j < N; ++j) {
	    scalar_type aux(0);
	    for (size_type k = 0; k < N; ++k)
	      aux += gradU(i, k) * Sigma(k, j);
	    t(i,j) = aux;
	  }
      }
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type ) {
      if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nb = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nb);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nb; ++k)
	    coeff[i * nb + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nb+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nb));
      }
    } 
    
  };


  /** 
      Tangent matrix for the non-linear elasticity 
      @ingroup asm
  */
  template<typename MAT, typename VECT1, typename VECT2> 
  void asm_nonlinear_elasticity_tangent_matrix
  (const MAT &K_, const mesh_im &mim, const getfem::mesh_fem &mf,
   const VECT1 &U, const getfem::mesh_fem *mf_data, const VECT2 &PARAMS,
   const abstract_hyperelastic_law &AHL,
   const mesh_region &rg = mesh_region::all_convexes()) {
    MAT &K = const_cast<MAT &>(K_);
    GMM_ASSERT1(mf.get_qdim() >= mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT1, VECT2>
      nterm(mf, U, mf_data, PARAMS, AHL, 0);
    elasticity_nonlinear_term<VECT1, VECT2>
      nterm2(mf, U, mf_data, PARAMS, AHL, 3);

    getfem::generic_assembly assem;
    if (mf_data)
      if (AHL.adapted_tangent_term_assembly_fem_data.size() > 0)
	assem.set(AHL.adapted_tangent_term_assembly_cte_data);
      else
	assem.set("M(#1,#1)+=sym(comp(NonLin$1(#1,#2)(i,j,k,l).vGrad(#1)(:,i,j).vGrad(#1)(:,k,l)))");
    else
      if (AHL.adapted_tangent_term_assembly_cte_data.size() > 0)
	assem.set(AHL.adapted_tangent_term_assembly_cte_data);
      else
	assem.set("M(#1,#1)+=sym(comp(NonLin$1(#1)(i,j,k,l).vGrad(#1)(:,i,j).vGrad(#1)(:,k,l)))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    if (mf_data) assem.push_mf(*mf_data);
    assem.push_data(PARAMS);
    assem.push_nonlinear_term(&nterm);
    assem.push_nonlinear_term(&nterm2);
    assem.push_mat(K);
    assem.assembly(rg);
  }


  /**
     @ingroup asm
  */
  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_rhs
  (const VECT1 &R_, const mesh_im &mim, const getfem::mesh_fem &mf,
   const VECT2 &U, const getfem::mesh_fem *mf_data, const VECT3 &PARAMS,
   const abstract_hyperelastic_law &AHL,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R = const_cast<VECT1 &>(R_);
    GMM_ASSERT1(mf.get_qdim() >= mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT2, VECT3>
      nterm(mf, U, mf_data, PARAMS, AHL, 1);

    getfem::generic_assembly assem;
    if (mf_data)
      assem.set("t=comp(NonLin(#1,#2).vGrad(#1)); V(#1) += t(i,j,:,i,j)");
    else
      assem.set("t=comp(NonLin(#1).vGrad(#1)); V(#1) += t(i,j,:,i,j)");
    // comp() to be optimized ?
    assem.push_mi(mim);
    assem.push_mf(mf);
    if (mf_data) assem.push_mf(*mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.assembly(rg);
  }

  // added by Jean-Yves Heddebaut <jyhed@alpino.be>
  int levi_civita(int i,int j,int k);


  /**@ingroup asm
   */
  template<typename VECT2, typename VECT3> 
  scalar_type asm_elastic_strain_energy
  (const mesh_im &mim, const getfem::mesh_fem &mf,
   const VECT2 &U, const getfem::mesh_fem *mf_data, const VECT3 &PARAMS,
   const abstract_hyperelastic_law &AHL,
   const mesh_region &rg = mesh_region::all_convexes()) {
 
    GMM_ASSERT1(mf.get_qdim() >= mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT2, VECT3>
      nterm(mf, U, mf_data, PARAMS, AHL, 2);
    std::vector<scalar_type> V(1);

    getfem::generic_assembly assem;
    if (mf_data)
      assem.set("V() += comp(NonLin(#1,#2))(i)");
    else
      assem.set("V() += comp(NonLin(#1))(i)");
    
    assem.push_mi(mim);
    assem.push_mf(mf);
    if (mf_data) assem.push_mf(*mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(V);
    assem.assembly(rg);

    return V[0];
  }









  /* ******************************************************************** */
  /*		Mixed nonlinear incompressibility assembly procedure      */
  /* ******************************************************************** */

  template<typename VECT1> class incomp_nonlinear_term 
    : public getfem::nonlinear_elem_term {

    const mesh_fem &mf;
    std::vector<scalar_type> U;
    size_type N;
    base_vector coeff;
    base_matrix gradPhi;
    bgeot::multi_index sizes_;
    int version; 

  public:
    incomp_nonlinear_term(const mesh_fem &mf_, const VECT1 &U_,
			  int version_) 
      : mf(mf_), U(mf_.nb_basic_dof()),
	N(mf_.get_qdim()),
	gradPhi(N, N), sizes_(N, N),
	version(version_) {
      if (version == 1) { sizes_.resize(1); sizes_[0] = 1; }
      mf.extend_vector(U_, U);
    }

    const bgeot::multi_index &sizes(size_type) const { return sizes_; }

    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      slice_vector_on_basic_dof_of_element(mf, U, cv, coeff);
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
    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

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
  void asm_nonlinear_incomp_rhs
  (const VECT1 &R_U_, const VECT1 &R_P_, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const getfem::mesh_fem &mf_p,
   const VECT2 &U, const VECT3 &P,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R_U = const_cast<VECT1 &>(R_U_);
    VECT1 &R_P = const_cast<VECT1 &>(R_P_);
    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

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



  //===========================================================================
  //
  //  Bricks
  //
  //===========================================================================


  /** Add a nonlinear (large strain) elasticity term to the model with
      respect to the variable
      `varname` (deprecated brick, use add_finite_strain_elaticity instead).
      Note that the constitutive law is described by `AHL` which
      should not be freed while the model is used. `dataname` described the
      parameters of the constitutive laws. It could be a vector of value
      of length the number of parameter of the constitutive law or a vector
      field described on a finite element method.
  */
  size_type add_nonlinear_elasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const phyperelastic_law &AHL, const std::string &dataname,
   size_type region = size_type(-1));



  void compute_Von_Mises_or_Tresca
  (model &md, const std::string &varname, const phyperelastic_law &AHL,
   const std::string &dataname, const mesh_fem &mf_vm,
   model_real_plain_vector &VM, bool tresca);


  void compute_sigmahathat(model &md,
			   const std::string &varname, 
			   const phyperelastic_law &AHL,
			   const std::string &dataname,
			   const mesh_fem &mf_sigma,
			   model_real_plain_vector &SIGMA);


  /**
     Compute the Von-Mises stress or the Tresca stress of a field
     with respect to the constitutive elasticity law AHL (only valid in 3D).
  */
  template <class VECTVM> void compute_Von_Mises_or_Tresca
  (model &md, const std::string &varname, const phyperelastic_law &AHL,
   const std::string &dataname, const mesh_fem &mf_vm,
   VECTVM &VM, bool tresca) {
    model_real_plain_vector VMM(mf_vm.nb_dof());
    compute_Von_Mises_or_Tresca
      (md, varname, AHL, dataname, mf_vm, VMM, tresca);
    gmm::copy(VMM, VM);
  }
  
  /** Add a nonlinear incompressibility term (for large strain elasticity)
      to the model with respect to the variable
      `varname` (the displacement) and `multname` (the pressure).
  */
  size_type add_nonlinear_incompressibility_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region = size_type(-1));

  //===========================================================================
  //  High-level generic assembly bricks
  //===========================================================================

  /** Add a finite strain elasticity brick
      to the model with respect to the variable
      `varname` (the displacement).
      For 2D meshes, switch automatically to plane strain elasticity.
      High-level generic assembly version.
  */
  size_type add_finite_strain_elasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &lawname, const std::string &params,
   size_type region = size_type(-1));


  /** Add a finite strain incompressibility term (for large strain elasticity)
      to the model with respect to the variable
      `varname` (the displacement) and `multname` (the pressure).
      High-level generic assembly version.
  */
  size_type add_finite_strain_incompressibility_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region = size_type(-1));

  /**
     Interpolate the Von-Mises stress of a field `varname`
     with respect to the nonlinear elasticity constitutive law `lawname`
     with parameters `params` (only valid in 3D).
  */
  void finite_strain_elasticity_Von_Mises
  (model &md, const std::string &varname, const std::string &lawname,
   const std::string &params, const mesh_fem &mf_vm,
   model_real_plain_vector &VM,
   const mesh_region &rg=mesh_region::all_convexes());

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NONLINEAR_ELASTICITY_H__ */
