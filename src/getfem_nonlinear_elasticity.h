/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_nonlinear_elasticity.h                                */
/*     									   */
/* Date : July 6, 2004.                                                    */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*           Julien Pommier, pommier@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2004  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
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


#ifndef GETFEM_NONLINEAR_ELASTICITY_H__
#define GETFEM_NONLINEAR_ELASTICITY_H__

#include <getfem_assembling_tensors.h>
namespace getfem {
 
  struct abstract_hyperelastic_law {
    size_type nb_params_;
    virtual void sigma(const base_matrix &L, base_matrix &result,
		       base_vector &params) const = 0;
    virtual void grad_sigma(const base_matrix &L, base_tensor &result, 
			    base_vector &params) const = 0;
    size_type nb_params(void) const { return nb_params_; }
    abstract_hyperelastic_law() { nb_params_ = 0; }
    virtual ~abstract_hyperelastic_law() {}
  };

  struct Hooke_hyperelastic_law : public abstract_hyperelastic_law {
    /* sigma = lambda*trace(L) + 2 mu * L */
    virtual void sigma(const base_matrix &L, base_matrix &result,
		       base_vector &params) const {
      gmm::copy(gmm::identity_matrix(), result);
      gmm::scale(result, params[0] * gmm::mat_trace(L));
      gmm::add(gmm::scaled(L, 2 * params[1]), result);
    }
    virtual void grad_sigma(const base_matrix &L, base_tensor &result,
			    base_vector &params) const {
      std::fill(result.begin(), result.end(), scalar_type(0));
      size_type N = gmm::mat_nrows(L);
      for (size_type i = 0; i < N; ++i)
	for (size_type l = 0; l < N; ++l) {
	  result(i, i, l, l) = params[0];
	  result(i, l, i, l) += params[1];
	  result(i, l, l, i) += params[1];
	}
    }
    Hooke_hyperelastic_law(void) { nb_params_ = 2; }
  };

  struct Mooney_Rivlin_hyperelastic_law : public abstract_hyperelastic_law {
    virtual void sigma(const base_matrix &L, base_matrix &result,
		       base_vector &params) const {
      scalar_type C1 = params[0], C2 = params[1];
      gmm::copy(gmm::identity_matrix(), result);
      gmm::scale(result, 4.0*C2*(gmm::mat_trace(L)+1.0) + C1*2.0);
      gmm::add(gmm::scaled(L, -4.0*C2), result);
    }
    virtual void grad_sigma(const base_matrix &L, base_tensor &result,
			    base_vector &params) const {
      scalar_type C2 = params[1];
      std::fill(result.begin(), result.end(), scalar_type(0));
      size_type N = gmm::mat_nrows(L);
      for (size_type i = 0; i < N; ++i)
	for (size_type l = 0; l < N; ++l) {
	  result(i, i, l, l) = 4.0 * C2;
	  result(i, l, i, l) -= 1.0;
	  result(i, l, l, i) -= 1.0;
	}
    }
    Mooney_Rivlin_hyperelastic_law(void) { nb_params_ = 2; }
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
    base_matrix L, B, gradU;
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
	L(N, N), B(N, N), gradU(N, N), tt(N, N, N, N), sizes_(N, N, N, N),
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

      gmm::mult(gmm::transposed(gradU), gradU, L);
      gmm::add(gradU, L);
      gmm::add(gmm::transposed(gradU), L);
      gmm::scale(L, scalar_type(0.5));
      gmm::add(gmm::identity_matrix(), gradU);

      AHL.sigma(L, B, params);

      if (version == 0) {	  
	AHL.grad_sigma(L, tt, params);
	
	for (size_type n = 0; n < N; ++n)
	  for (size_type m = 0; m < N; ++m)
	    for (size_type l = 0; l < N; ++l)
	      for (size_type k = 0; k < N; ++k) {
		scalar_type aux = (k == l) ? B(m, l) : 0.0;
		for (size_type j = 0; j < N; ++j)
		  for (size_type i = 0; i < N; ++i) {
		    aux += B(n ,j) * B(k, i) * tt(j, m, i, l);
		  }
		t(n, m, k, l) = aux;
	      }
      } else {
	for (size_type i = 0; i < N; ++i)
	  for (size_type j = 0; j < N; ++j) {
	    scalar_type aux(0);
	    for (size_type k = 0; k < N; ++k)
	      aux += gradU(i, k) * B(k, j);
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

  template<typename MAT, typename VECT1, typename VECT2> 
  void asm_nonlinear_elasticity_tangent_matrix(const MAT &K_, 
					       const getfem::mesh_fem &mf,
					       const VECT1 &U,
					       const getfem::mesh_fem &mf_data,
					       const VECT2 &PARAMS,
				     const abstract_hyperelastic_law &AHL) {
    MAT &K = const_cast<MAT &>(K_);
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT1, VECT2>
      nterm(mf, U, mf_data, PARAMS, AHL, 0);

    getfem::generic_assembly
      assem("t=comp(NonLin(#1,#2).vGrad(#1).vGrad(#1));"
	    "M(#1,#1)+= sym(t(i,j,k,l,:,i,j,:,k,l))");

    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(K);
    assem.volumic_assembly();
  }


  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_rhs(const VECT1 &R_, 
				    const getfem::mesh_fem &mf,
				    const VECT2 &U,
				    const getfem::mesh_fem &mf_data,
				    const VECT3 &PARAMS,
				    const abstract_hyperelastic_law &AHL) {
    VECT1 &R = const_cast<VECT1 &>(R_);
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    elasticity_nonlinear_term<VECT2, VECT3>
      nterm(mf, U, mf_data, PARAMS, AHL, 1);

    getfem::generic_assembly
      assem("t=comp(NonLin(#1,#2).vGrad(#1)); V(#1) += t(i,j,:,i,j)");

    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R);
    assem.volumic_assembly();
  }


  /* ******************************************************************** */
  /*		Nonlinear elasticity brick.                               */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_nonlinear_elasticity : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    const abstract_hyperelastic_law &AHL;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR PARAMS_;
    bool homogeneous;

  public :

    virtual bool is_linear(void) { return false; }
    virtual bool is_coercive(void) { return true; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false) {
      react(MS, i0, modified);
      size_type nb = AHL.nb_params();
      VECTOR PARAMS(mf_data.nb_dof() * nb);
      if (homogeneous) {
	for (size_type i = 0; i < mf_data.nb_dof(); ++i) 
	  gmm::copy(PARAMS_,
		    gmm::sub_vector(PARAMS, gmm::sub_interval(i*nb, nb)));
      }
      else
	gmm::copy(PARAMS_, PARAMS);

      gmm::sub_interval SUBI(i0, nb_dof());
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      asm_nonlinear_elasticity_tangent_matrix
	(gmm::sub_matrix(MS.tangent_matrix(), SUBI), mf_u,
	 gmm::sub_vector(MS.state(), SUBI), mf_data, PARAMS,  AHL);
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0) {
      react(MS, i0, false);
      size_type nb = AHL.nb_params();
      VECTOR PARAMS(mf_data.nb_dof() * nb);
      if (homogeneous) {
	for (size_type i = 0; i < mf_data.nb_dof(); ++i) 
	  gmm::copy(PARAMS_,
		    gmm::sub_vector(PARAMS, gmm::sub_interval(i*nb, nb)));
      }
      else
	gmm::copy(PARAMS_, PARAMS);

      gmm::sub_interval SUBI(i0, nb_dof());
      gmm::clear(gmm::sub_vector(MS.residu(), SUBI));
      asm_nonlinear_elasticity_rhs(gmm::sub_vector(MS.residu(), SUBI), mf_u,
				   gmm::sub_vector(MS.state(), SUBI), 
				   mf_data, PARAMS, AHL);
    }
    virtual mesh_fem &main_mesh_fem(void) { return mf_u; }

    void set_params(const VECTOR &PARAMS) {
      homogeneous = gmm::vect_size(PARAMS) == AHL.nb_params();
      gmm::resize(PARAMS_, homogeneous ? AHL.nb_params()
		  : mf_data.nb_dof() * AHL.nb_params());
      gmm::copy(PARAMS, PARAMS_);
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    mdbrick_nonlinear_elasticity(const abstract_hyperelastic_law &AHL_,
				 mesh_fem &mf_u_, mesh_fem &mf_data_,
				 const VECTOR &PARAMS)
      : AHL(AHL_), mf_u(mf_u_), mf_data(mf_data_) {
      set_params(PARAMS);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }
 
    mdbrick_nonlinear_elasticity(const abstract_hyperelastic_law &AHL_,
				 mesh_fem &mf_u_, mesh_fem &mf_data_,
				 value_type p1, value_type p2)
      : AHL(AHL_), mf_u(mf_u_), mf_data(mf_data_) {
      VECTOR PARAMS(2); PARAMS[0] = p1;  PARAMS[1] = p2; 
      set_params(PARAMS);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }

  };

  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NONLINEAR_ELASTICITY_H__ */
