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
	// the result of grad_sigma has to be completely symmetric.
    virtual void grad_sigma(const base_matrix &L, base_tensor &result, 
			    base_vector &params) const = 0;
    size_type nb_params(void) const { return nb_params_; }
    abstract_hyperelastic_law() { nb_params_ = 0; }
    virtual ~abstract_hyperelastic_law() {}
  };

  int check_symmetry(const base_tensor &t) {
    int flags = 7; size_type N = 3;
    for (size_type n = 0; n < N; ++n)
      for (size_type m = 0; m < N; ++m)
	for (size_type l = 0; l < N; ++l)
	  for (size_type k = 0; k < N; ++k) {
	    if (dal::abs(t(n,m,l,k) - t(l,k,n,m))>1e-10) flags &= (~1); 
	    if (dal::abs(t(n,m,l,k) - t(m,n,l,k))>1e-10) flags &= (~2); 
	    if (dal::abs(t(n,m,l,k) - t(n,m,k,l))>1e-10) flags &= (~4);
	  }
    return flags;
  }

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
      assert(check_symmetry(result) == 7);
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
    base_matrix L, Sigma, gradU;
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
	L(N, N), Sigma(N, N), gradU(N, N), tt(N, N, N, N), sizes_(N, N, N, N),
	version(version_)
    { if (version == 1) sizes_.resize(2); }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))),
		coeff);
      base_matrix gradUt(3,3);
      ctx.pf()->interpolation_grad(ctx, coeff, gradUt, mf.get_qdim());
      gmm::copy(gmm::transposed(gradUt),gradU);
      gmm::mult(gmm::transposed(gradU), gradU, L);
      gmm::add(gradU, L);
      gmm::add(gmm::transposed(gradU), L);
      gmm::scale(L, scalar_type(0.5));
      gmm::add(gmm::identity_matrix(), gradU);


      AHL.sigma(L, Sigma, params);

      //cout << "nonlinear_elem_term::compute(version=" << version << ", cv = " << cv << ") -> \n" << "   L=" << L << "\n   E=" << gradU << "\n   Sigma=" << Sigma << "\n";

      if (version == 0) {	  
	AHL.grad_sigma(L, tt, params);


	/*
	  for (size_type n = 0; n < N; ++n)
	    for (size_type m = 0; m <= n; ++m)
	        for (size_type k = 0; k <= m; ++k)
		      for (size_type l = 0; l <= k; ++l) {
		      // scalar_type aux = (k == l) ? B(m, l) : scalar_type(0);
		      scalar_type aux(0);

		      for (size_type j = 0; j < N; ++j)
		        for (size_type i = 0; i < N; ++i) {
			    aux += gradU(n ,j) * gradU(k, i) * tt(j, m, i, l);

			      }
			      t(n, m, k, l) = t(m, n, k, l) = t(n, m, l, k) = aux;
			      t(m, n, l, k) = t(k, l, m, n) = t(l, k, m, n) = aux;
			      t(k, l, n, m) = t(l, k, n, m) = aux;
			            }

				    for (size_type n = 0; n < N; ++n)
				      for (size_type m = 0; m < N; ++m)
				          for (size_type l = 0; l < N; ++l) {
					        t(n, m, n, l) += B(m, l);
						      // t(n, m, l, m) += B(n, l) * 0.5;
						          }
	*/
	
	for (size_type n = 0; n < N; ++n)
	  for (size_type m = 0; m < N; ++m)
	    for (size_type l = 0; l < N; ++l)
	      for (size_type k = 0; k < N; ++k) {
		//scalar_type aux = (k == l) ? Sigma(m, l) : 0.0;
		//scalar_type aux = (m == l) ? Sigma(k,n) : 0.0;
		scalar_type aux = (k == n) ? Sigma(m,l) : 0.0;
		for (size_type j = 0; j < N; ++j)
		  for (size_type i = 0; i < N; ++i) {
		    aux += gradU(n ,j) * gradU(k, i) * tt(j, m, i, l);
		    //aux += gradU(n ,i) * gradU(j, k) * tt(i,m,j,l);
		  }
		t(n, m, k, l) = aux;
	      }
	/*cout << "sym = " << check_symmetry(t) << "\n";
	  if (check_symmetry(t) != 7) cout << "t=" << t << "\n";*/
      } else {
	scalar_type J = gmm::lu_det(gradU);
	if (J < 0) gmm::scale(gradU, 1e10);

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
      //cout << "nonlinear_elem_term::prepare(cv = " << cv << ") -> params=" << params << "\n";
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
      //cout << "mdbrick_nonlinear_elasticity::compute_tangent_matrix\n";
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
      //cout << "mdbrick_nonlinear_elasticity returns " << MS.tangent_matrix() << "\n";
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
      //cout << "mdbrick_nonlinear_elasticity::compute_residu -> " << gmm::vect_norm2(MS.residu()) << "\n";
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





  /* ******************************************************************** */
  /*		Mixed nonlinear incompressible condition brick.           */
  /* ******************************************************************** */


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
      base_matrix gradPhit(3,3);
      ctx.pf()->interpolation_grad(ctx, coeff, gradPhit, mf.get_qdim());
      gmm::copy(gmm::transposed(gradPhit),gradPhi);

      gmm::add(gmm::identity_matrix(), gradPhi);

      scalar_type det = gmm::lu_inverse(gradPhi);

      if (version != 1) {
	if (version == 2) det = sqrt(dal::abs(det));
	for (size_type i = 0; i < N; ++i) 
	  for (size_type j = 0; j < N; ++j) {
	    t(i,j) = - det * gradPhi(j,i);
	  }
      } else if (version == 1) {
	t[0] = 1 - det;
      }
    }
  };

  template<typename MAT1, typename MAT2, typename VECT1, typename VECT2> 
  void asm_nonlinear_incomp_tangent_matrix(const MAT1 &K_, const MAT2 &B_, 
					   const getfem::mesh_fem &mf_u,
					   const getfem::mesh_fem &mf_p,
					   const VECT1 &U, const VECT2 &P) {
    MAT1 &K = const_cast<MAT1 &>(K_);
    MAT2 &B = const_cast<MAT2 &>(B_);
    if (mf_u.get_qdim() != mf_u.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    incomp_nonlinear_term<VECT1>
      ntermk(mf_u, U, 0);
    incomp_nonlinear_term<VECT1>
      ntermb(mf_u, U, 2);

    getfem::generic_assembly
      assem("P=data(#2);"
	    "t=comp(NonLin(#1).vGrad(#1).Base(#2));"
	    "M$2(#1,#2)+= t(i,j,:,i,j,:);"
	    "w=comp(NonLin$2(#1).vGrad(#1).NonLin$2(#1).vGrad(#1).Base(#2));"
	    "M$1(#1,#1)+= w(i,j,:,j,i, k,l,:,l,k,p).P(p) + w(i,j,:,k,j, k,m,:,i,m,p).P(p)");

    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&ntermk);
    assem.push_nonlinear_term(&ntermb);
    assem.push_mat(K);
    assem.push_mat(B);
    assem.push_data(P);
    assem.volumic_assembly();
  }


  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_incomp_rhs(const VECT1 &R_U_, const VECT1 &R_P_, 
				const getfem::mesh_fem &mf_u,
				const getfem::mesh_fem &mf_p,
				const VECT2 &U, const VECT3 &P) {
    VECT1 &R_U = const_cast<VECT1 &>(R_U_);
    VECT1 &R_P = const_cast<VECT1 &>(R_P_);
    if (mf_u.get_qdim() != mf_u.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    incomp_nonlinear_term<VECT2>
      nterm_tg(mf_u, U, 0);
    incomp_nonlinear_term<VECT2>
      nterm(mf_u, U, 1);

    getfem::generic_assembly
      assem("P=data(#2); "
	    "t=comp(NonLin$1(#1).vGrad(#1).Base(#2)); V$1(#1) += t(i,j,:,i,j,k).P(k);"
	    "w=comp(NonLin$2(#1).Base(#2)); V$2(#2) += w(1,:)");

    assem.push_mf(mf_u);
    assem.push_mf(mf_p);
    assem.push_nonlinear_term(&nterm_tg);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(R_U);
    assem.push_vec(R_P);
    assem.push_data(P);
    assem.volumic_assembly();
  }


  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_nonlinear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_p;

  public :
    
    virtual bool is_linear(void)   { return false; }
    virtual bool is_coercive(void) { return false; }
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      b.add(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
    }
    virtual size_type nb_dof(void) {
      return sub_problem.nb_dof() + mf_p.nb_dof();
    }
    
    virtual size_type nb_constraints(void) {
      return sub_problem.nb_constraints();
    }

    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0, bool modified = false) {
      sub_problem.compute_tangent_matrix(MS, i0, j0, modified);
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_p.nb_dof()); /* P */
      gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());           /* U */

      asm_nonlinear_incomp_tangent_matrix(gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBJ),
					  gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI),
					  main_mesh_fem(), mf_p, 
					  gmm::sub_vector(MS.state(), SUBJ), 
					  gmm::sub_vector(MS.state(), SUBI));
      
      gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI)),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }

    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
     
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());

      asm_nonlinear_incomp_rhs(gmm::sub_vector(MS.residu(), SUBJ),
			       gmm::sub_vector(MS.residu(), SUBI),
			       main_mesh_fem(), mf_p, 
			       gmm::sub_vector(MS.state(), SUBJ),
			       gmm::sub_vector(MS.state(), SUBI));
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    // Constructor which does not define the rhs
    mdbrick_nonlinear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_)
      : sub_problem(problem), mf_p(mf_p_) {
      this->add_dependency(mf_p);
      this->add_dependency(sub_problem.main_mesh_fem());
    }
  };



}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NONLINEAR_ELASTICITY_H__ */
