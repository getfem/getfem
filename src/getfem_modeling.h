/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_modeling.h :  Defines model bricks to build           */
/*                                 complete models.                        */
/*     									   */
/* Date : June 15, 2004.                                                   */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
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
/*                                                                         */
/* Required properties of a model brick :                                  */
/*                                                                         */
/*  A model brick is either a fondamental brick (like linearized           */
/*  elasticity brick, platicity brick ...) or a modifier brick which refer */
/*  to a sub brick.                                                        */
/*  The virtual functions of a brick :                                     */
/*  - nb_dof() : number of total variables including the variables of the  */
/*        sub-problem(s) if any.                                           */
/*  - nb_constraints() : number of linear constraints on the system        */
/*        including the constraints defined in the sub-problem(s) if any.  */
/*  - is_linear()   : true if the problem is linear.                       */
/*  - is_coercive() : true if the problem is symmetric coercive.           */
/*  - compute_tangent_matrix(MS, i0, j0, modified) : the brick has to call */
/*        the compute_tangent_matrix(MS, i0,modified) of sub-problem(s) if */
/*        any and to compute  its own part of the tangentand constraint    */
/*        matrices (i0 and j0 are the shifts in the matrices defined in MS)*/
/*        modified indicates if the part of the tangent matrix dedicated   */
/*        to this brick will be modified by other bricks.                  */
/*  - compute_residu(MS, i0, j0) : the brick has to call the               */
/*        compute_residu(MS, i0, modified) of sub-problem(s) if any and    */
/*        has to compute its own part of the residu of the linear system   */
/*        and of the constraint system (i0 and j0 are the shifts in the    */
/*        residu vectors defined in MS)                                    */
/*  - mixed_variables(bv, i0) : indicates in bv the indices of the         */
/*        variables which are considered as multipliers or mixed variables.*/
/*  - main_mesh_fem() : the principal finite element method. For instance  */
/*         a Dirichlet condition will act on this main mesh_fem. For a     */
/*         mixed method, the main fem will be the primal variable fem.     */
/*                                                                         */
/* Dependencies.                                                           */
/*   A brick depends at least on some mesh_fem structures and has to       */
/* react if some changements occur in these mesh_fem structures.           */
/*                                                                         */
/***************************************************************************/

#ifndef GETFEM_MODELING_H__
#define GETFEM_MODELING_H__

#include <getfem_assembling.h>
#include <gmm_solver_cg.h>
#include <gmm_solver_gmres.h>
#include <gmm_precond_ildlt.h>
#include <gmm_precond_ilu.h>
#include <gmm_precond_ilut.h>
#include <gmm_superlu_interface.h>

namespace getfem {

  /* ******************************************************************** */
  /*		Generic definitions.                                      */
  /* ******************************************************************** */
 
  template<typename MODEL_STATE> class mdbrick_abstract;

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  class model_state {
  public :    
    typedef T_MATRIX tangent_matrix_type;
    typedef C_MATRIX constraints_matrix_type;
    typedef VECTOR vector_type;
    typedef typename gmm::linalg_traits<VECTOR>::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type
      magnitude_type;

  protected :
    T_MATRIX tangent_matrix_;
    C_MATRIX constraints_matrix_;
    VECTOR state_, residu_, constraints_rhs_;
    ctx_ident_type ident_;
    
    T_MATRIX SM;
    gmm::col_matrix< gmm::rsvector<value_type> > NS; /* nullspace of constraints */
    VECTOR reduced_residu_, Ud;
  public :

    const T_MATRIX &tangent_matrix(void) const 
    { return tangent_matrix_; }
    T_MATRIX &tangent_matrix(void) { return tangent_matrix_; }
    const C_MATRIX &constraints_matrix(void) const 
    { return constraints_matrix_; }
    C_MATRIX &constraints_matrix(void) { return constraints_matrix_; }
    const VECTOR &constraints_rhs(void) const  { return constraints_rhs_; }
    VECTOR &constraints_rhs(void)  { return constraints_rhs_; }
    const VECTOR &state(void) const  { return state_; }
    VECTOR &state(void)  { return state_; }
    const VECTOR &residu(void) const  { return residu_; }
    const magnitude_type reduced_residu_norm() const {
      if (gmm::mat_nrows(constraints_matrix())) {
	return sqrt(gmm::vect_norm2_sqr(reduced_residu_) + 
		    gmm::vect_norm2_sqr(Ud));
      } else return gmm::vect_norm2(residu_);
    }
    const VECTOR &reduced_residu() const { 
      return gmm::mat_nrows(constraints_matrix()) == 0 ?
	residu_ : reduced_residu_;
    }
    const T_MATRIX &reduced_tangent_matrix() const {
      return gmm::mat_nrows(constraints_matrix()) == 0 ?
	tangent_matrix_ : SM;
    }
    void unreduced_solution(const VECTOR &U_reduced, VECTOR &U) {
      if (gmm::mat_nrows(constraints_matrix()))
	gmm::mult(NS, U_reduced, Ud, U);
      else gmm::copy(U_reduced, U);
    }
    void compute_reduced_system();
    VECTOR &residu(void) { return residu_; }
    ctx_ident_type ident(void) { return ident_; }
    void touch(void) { ident_ = context_dependencies::new_ident(); }
    void adapt_sizes(mdbrick_abstract<model_state> &problem);
    model_state(void) { touch(); }
    model_state(mdbrick_abstract<model_state> &problem)
    { adapt_sizes(problem); }
  };

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::compute_reduced_system() {
    if (gmm::mat_nrows(constraints_matrix()) == 0) return;
    size_type ndof = gmm::mat_ncols(tangent_matrix());
    gmm::resize(NS, ndof, ndof);
    gmm::resize(Ud, ndof);
    
    size_type nbcols=getfem::Dirichlet_nullspace(constraints_matrix(),
						 NS, constraints_rhs(), Ud);
    gmm::resize(NS, ndof, nbcols);
    gmm::resize(SM, nbcols, nbcols);
    VECTOR RHaux(ndof);
    gmm::mult(tangent_matrix(), Ud, residu(), RHaux);
    gmm::resize(reduced_residu_, nbcols);
    gmm::mult(gmm::transposed(NS), RHaux, reduced_residu_);
    T_MATRIX SMaux(nbcols, ndof);
    gmm::col_matrix< gmm::rsvector<value_type> >
      NST(gmm::mat_ncols(NS), gmm::mat_nrows(NS));
    gmm::copy(gmm::transposed(NS), NST);
    gmm::mult(NST, tangent_matrix(), SMaux);
    gmm::mult(SMaux, NS, SM);
  }

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::
  adapt_sizes(mdbrick_abstract<model_state> &problem) {
    size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();

    if (gmm::mat_nrows(tangent_matrix_) != ndof
	|| gmm::mat_nrows(constraints_matrix_) != nc) {
      gmm::resize(tangent_matrix_, ndof, ndof);
      gmm::resize(constraints_matrix_, nc, ndof);
      gmm::resize(constraints_rhs_, nc);
      gmm::resize(state_, ndof); gmm::clear(state_);
      gmm::resize(residu_, ndof);
      touch();
    }
  } 

  template<typename MODEL_STATE>
  class mdbrick_abstract : public context_dependencies {
  protected :
    bool to_compute, to_transfer;
    size_type MS_i0;
    ctx_ident_type ident_ms;

    // to_be_computed : the context has changed (or it is the first call).
    // to_be_transferred : the structure MODEL_STATE has changed, the 
    //                     tangent matrix or constraints system has to be
    //                     copied again if it is stored.
    void react(MODEL_STATE &MS, size_type i0, bool modified) {
      if (this->context_changed()) to_compute = true;
      to_transfer = to_transfer || modified || (ident_ms != MS.ident())
	|| to_compute;
      ident_ms = MS.ident();
      MS_i0 = i0;
    }

    bool to_be_computed(void) { return to_compute; }
    bool to_be_transferred(void) { return to_transfer; }
    void force_recompute(void) { to_compute = to_transfer = true; }
    void computed(void) { to_compute = false; }
    void transferred(void) { to_transfer = false; }
    size_type first_index(void) { return MS_i0; }

  public :

    virtual size_type nb_dof(void) = 0;
    virtual size_type nb_constraints(void) = 0;
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				  size_type j0=0, bool modified = false) = 0;
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) = 0;
    virtual mesh_fem &main_mesh_fem(void) = 0;
    virtual bool is_linear(void) = 0;
    virtual bool is_coercive(void) = 0;
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) = 0;
    mdbrick_abstract(void) :  to_compute(true), to_transfer(true),
			      MS_i0(0), ident_ms(-1) { }
    virtual ~mdbrick_abstract() {}
  };

  typedef gmm::rsvector<scalar_type> modeling_standard_sparse_vector;
  typedef gmm::col_matrix<modeling_standard_sparse_vector>
                                    modeling_standard_sparse_matrix;
  typedef std::vector<scalar_type> modeling_standard_plain_vector;

  typedef gmm::rsvector<complex_type> modeling_standard_complex_sparse_vector;
  typedef gmm::col_matrix<modeling_standard_complex_sparse_vector>
                                    modeling_standard_complex_sparse_matrix;
  typedef std::vector<complex_type> modeling_standard_complex_plain_vector;

  typedef model_state<modeling_standard_sparse_matrix,
		      modeling_standard_sparse_matrix,
		      modeling_standard_plain_vector > standard_model_state;

  typedef model_state<modeling_standard_complex_sparse_matrix,
		      modeling_standard_complex_sparse_matrix,
		      modeling_standard_complex_plain_vector >
    standard_complex_model_state;

  /* ******************************************************************** */
  /*		general scalar elliptic brick.                            */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_scalar_elliptic : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR coeffs_;
    bool homogeneous, laplacian;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void);

  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return true; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false);
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0);
    virtual mesh_fem &main_mesh_fem(void) { return mf_u; }

    
    void set_coeff(value_type a) {
      homogeneous = true; laplacian = true;
      gmm::resize(coeffs_, 1); coeffs_[0] = a;
      this->force_recompute();
    }

    void set_coeff(const VECTOR &coeffs, bool laplace) {
      laplacian = laplace;
      homogeneous = false;
      int N = mf_u.linked_mesh().dim();
      if (laplacian) {
	if (gmm::vect_size(coeffs) == 1) set_coeff(coeffs[0]);
	else gmm::resize(coeffs_,  mf_data.nb_dof());
      }
      else {
	if (gmm::vect_size(coeffs) == dal::sqr(N)) {
	  gmm::resize(coeffs_, dal::sqr(N));
	  homogeneous = true;
	}
	else
	  gmm::resize(coeffs_, mf_data.nb_dof() * dal::sqr(N));
      }
      gmm::copy(coeffs, coeffs_);
      this->force_recompute();
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    // constructor for the Laplace operator
    mdbrick_scalar_elliptic(mesh_fem &mf_u_, mesh_fem &mf_data_,
       value_type a, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored) {
      set_coeff(a);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }

    // constructor for a non-homogeneous material
    mdbrick_scalar_elliptic(mesh_fem &mf_u_, mesh_fem &mf_data_,
       const VECTOR &coeff, bool laplace, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored) {
      set_coeff(coeff, laplace);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }
 
  };

  template<typename MODEL_STATE>
   void mdbrick_scalar_elliptic<MODEL_STATE>::compute_K(void) {
    gmm::clear(K);
    gmm::resize(K, nb_dof(), nb_dof());
    size_type n = laplacian ? 1 : dal::sqr(mf_u.linked_mesh().dim());
    VECTOR coeffs(n * mf_data.nb_dof());
    if (homogeneous) {
      for (size_type i = 0; i < mf_data.nb_dof(); ++i)
	gmm::copy(coeffs_, gmm::sub_vector(coeffs, gmm::sub_interval(i*n, n)));
    }
    else { gmm::copy(coeffs_, coeffs); }
    if (laplacian)
      asm_stiffness_matrix_for_laplacian(K, mf_u, mf_data, coeffs);
    else
      asm_stiffness_matrix_for_scalar_elliptic(K, mf_u, mf_data, coeffs);
    this->computed();
  }

  template<typename MODEL_STATE>
  void mdbrick_scalar_elliptic<MODEL_STATE>::
  compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
			 size_type, bool modified) {
    if (modified && !matrix_stored) 
      DAL_THROW(failure_error, "The residu will not be consistant. "
		"Use this brick with the stiffness matrix stored option");
    react(MS, i0, modified);
    gmm::sub_interval SUBI(i0, nb_dof());
    if (this->to_be_computed()
	|| (!matrix_stored && this->to_be_transferred()))
      compute_K();
    if (this->to_be_transferred()) { 
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      this->transferred();
    }
    if (!matrix_stored) gmm::clear(K);
  }
  
  template<typename MODEL_STATE>
  void mdbrick_scalar_elliptic<MODEL_STATE>::
  compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
    react(MS, i0, false);
    gmm::sub_interval SUBI(i0, nb_dof());
    if (this->to_be_computed()) { 
      compute_K();
      if (!matrix_stored) {
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	gmm::clear(K);
      }
    }
    if (matrix_stored) {
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    } else {
      gmm::mult(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
		gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }
  }


  /* ******************************************************************** */
  /*		Linearized elasticity bricks.                             */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Hooke_linearized_elasticity
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR lambda_, mu_;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void);

  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return true; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false);
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0);
    virtual mesh_fem &main_mesh_fem(void) { return mf_u; }

    void set_Lame_coeff(value_type lambdai, value_type mui) {
      homogeneous = true;
      gmm::resize(lambda_, 1); lambda_[0] = lambdai;
      gmm::resize(mu_, 1); mu_[0] = mui;
      this->force_recompute();
    }

    void set_Lame_coeff(const VECTOR &lambdai, const VECTOR &mui) {
      homogeneous = false;
      gmm::resize(lambda_, mf_data.nb_dof()); gmm::copy(lambdai, lambda_);
      gmm::resize(mu_, mf_data.nb_dof()); gmm::copy(mui, mu_);
      this->force_recompute();
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_Hooke_linearized_elasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
       value_type lambdai, value_type mui, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored) {
      set_Lame_coeff(lambdai, mui);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }

    // constructor for a non-homogeneous material
    mdbrick_Hooke_linearized_elasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
       const VECTOR &lambdai, const VECTOR &mui, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored) {
      set_Lame_coeff(lambdai, mui);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }
 
  };


  template<typename MODEL_STATE>
   void mdbrick_Hooke_linearized_elasticity<MODEL_STATE>::compute_K(void) {
    gmm::clear(K);
    gmm::resize(K, nb_dof(), nb_dof());
    VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
    if (homogeneous) {
      std::fill(lambda.begin(), lambda.end(), value_type(lambda_[0]));
      std::fill(mu.begin(), mu.end(), value_type(mu_[0]));
    }
    else { gmm::copy(lambda_, lambda); gmm::copy(mu_, mu); }
    asm_stiffness_matrix_for_linear_elasticity(K, mf_u, mf_data, lambda, mu);
    this->computed();
  }

  template<typename MODEL_STATE>
  void mdbrick_Hooke_linearized_elasticity<MODEL_STATE>::
  compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
			 size_type, bool modified) {
    if (modified && !matrix_stored) 
      DAL_THROW(failure_error, "The residu will not be consistant. "
		"Use this brick with the stiffness matrix stored option");
    react(MS, i0, modified);
    gmm::sub_interval SUBI(i0, nb_dof());
    if (this->to_be_computed()
	|| (!matrix_stored && this->to_be_transferred()))
      compute_K();
    if (this->to_be_transferred()) { 
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      this->transferred();
    }
    if (!matrix_stored) gmm::clear(K);
  }
  
  template<typename MODEL_STATE>
  void mdbrick_Hooke_linearized_elasticity<MODEL_STATE>::
  compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
    react(MS, i0, false);
    gmm::sub_interval SUBI(i0, nb_dof());
    if (this->to_be_computed()) { 
      compute_K();
      if (!matrix_stored) {
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	gmm::clear(K);
      }
    }
    if (matrix_stored) {
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    } else {
      gmm::mult(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
		gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }
  }

  /* TODO : arbitrary elasticity tensor */




  /* ******************************************************************** */
  /*		Plasticity bricks.                                        */
  /* ******************************************************************** */

  /* TODO :
       - non-homogeneous case
       - saved_proj useful?
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
    value_type stress_threshold, VM_max, TOL;
    size_type N, flag_hyp;
    
    std::vector<std::vector<scalar_type> > sigma_bar;
    std::vector<std::vector<scalar_type> > saved_proj;

    VM_projection VM_1;

    bool save_constraints;

  public:
    
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual bool is_linear(void) { return false; }
    virtual bool is_coercive(void) { return false; } // means we use LU factorisation instead of conjugate gradient
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    void set_save_constraints(bool b) { save_constraints = b; }
    
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
					       TOL, lambda_[0], mu_[0], &VM_1,
					       sigma_bar, saved_proj,1,flag_hyp, save_constraints);
      plasticity_projection gradproj(&gradproj_base);
      
      for (size_type i=0;i<mf_data.nb_dof();++i) {
	lambda[i]=lambda_[0];
	mu[i]=mu_[0];
      }
	
      /* Calculate the actual matrix */
      asm_lhs_for_plasticity(K, mf_u, mf_data, lambda, mu, &gradproj);
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0) {
      gmm::sub_interval SUBI(i0, nb_dof());        
      VECTOR K;
 
      /* Initialise K correctly */
      gmm::clear(K);
      gmm::resize(K, nb_dof());
      plasticity_projection_base proj_base(mf_u, MS.state(), stress_threshold,
						TOL, lambda_[0], mu_[0], &VM_1, sigma_bar,
						//saved_proj,0,flag_hyp, true);
						saved_proj,0,flag_hyp, save_constraints);
      plasticity_projection proj(&proj_base);
	
      /* Calculate the actual matrix */
      asm_rhs_for_plasticity(K, mf_u, &proj);
      gmm::copy(K, gmm::sub_vector(MS.residu(), SUBI));
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
		       value_type VM_max_, value_type TOL_, size_type flag_hyp_,
		       std::vector<std::vector<scalar_type> > sigma_b)
      : mf_u(mf_u_), mf_data(mf_data_) {
      set_Lame_coeff(lambdai, mui);
      N = mf_data.linked_mesh().dim();
      stress_threshold = stress_threshold_;
      VM_max = VM_max_;
      TOL = TOL_;
      flag_hyp = flag_hyp_;
      this->add_dependency(mf_u); this->add_dependency(mf_data);
      set_save_constraints(false);
    }
    
    // constructor for a non-homogeneous material
    mdbrick_plasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
		       const VECTOR &lambdai, const VECTOR &mui, value_type stress_treshhold_,
		       value_type VM_max_, value_type TOL_, size_type flag_hyp_ ,
		       std::vector<std::vector<scalar_type> > sigma_b)
      : mf_u(mf_u_), mf_data(mf_data_) {
      set_Lame_coeff(lambdai, mui);
      N = mf_data.linked_mesh().dim();
      stress_treshold = stress_treshold_;
      VM_max = VM_max_;
      TOL = TOL_;
      flag_hyp = flag_hyp_;
      this->add_dependency(mf_u); this->add_dependency(mf_data);
      set_save_constraints(false);
    }
  };
  
  /* ******************************************************************** */
  /*		Helmholtz brick.                                          */
  /* ******************************************************************** */
 template<typename MODEL_STATE = standard_complex_model_state>
  class mdbrick_Helmholtz
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR wave_number;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void) {
      gmm::clear(K);
      gmm::resize(K, nb_dof(), nb_dof());
      VECTOR wave_number2(mf_data.nb_dof());
      if (homogeneous)
	std::fill(wave_number2.begin(), wave_number2.end(),
		  value_type(dal::sqr(wave_number[0])));
      else
	for (size_type i=0; i < nb_dof(); ++i)
	  wave_number2[i] = dal::sqr(wave_number[i]);
      
      asm_Helmholtz(K, mf_u, mf_data, wave_number2);
      this->computed();
    }

  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return false; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false) {
      if (modified && !matrix_stored) 
	DAL_THROW(failure_error, "The residu will not be consistant. "
		  "Use this brick with the stiffness matrix stored option");
      react(MS, i0, modified);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed()
	  || (!matrix_stored && this->to_be_transferred()))
	  compute_K();
      if (this->to_be_transferred()) { 
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
	this->transferred();
      }
      if (!matrix_stored) gmm::clear(K);
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0) {
      react(MS, i0, false);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed()) {
	compute_K();
	if (!matrix_stored) {
	  gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	  gmm::clear(K);
	}
      }
      if (matrix_stored) {
	gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		     gmm::sub_vector(MS.residu(), SUBI));
      } else {
	gmm::mult(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
		  gmm::sub_vector(MS.state(), SUBI),
		  gmm::sub_vector(MS.residu(), SUBI));
      }
    }
    virtual mesh_fem &main_mesh_fem(void) { return mf_u; }

    void set_wave_number(complex_type k) {
      homogeneous = true;
      gmm::resize(wave_number, 1); wave_number[0] = k;
      this->force_recompute();
    }

    void set_wave_number(const VECTOR &k) {
      homogeneous = false;
      gmm::resize(wave_number, mf_data.nb_dof()); gmm::copy(k, wave_number);
      this->force_recompute();
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    // constructor for a homogeneous wave number
    mdbrick_Helmholtz(mesh_fem &mf_u_, mesh_fem &mf_data_,
       complex_type k, bool mat_stored = true)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored) {
      set_wave_number(k);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }

    // constructor for a non-homogeneous wave number
    mdbrick_Helmholtz(mesh_fem &mf_u_, mesh_fem &mf_data_,
		      const VECTOR &k, bool mat_stored = true)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored) {
      set_wave_number(k);
      this->add_dependency(mf_u); this->add_dependency(mf_data);
    }
  };


  /* ******************************************************************** */
  /*		Source term brick.                                        */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_source_term : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR B_;
    VECTOR F_;
    size_type boundary, qmult;

    void compute_F(void) {
      qmult = sub_problem.main_mesh_fem().get_qdim() / mf_data.get_qdim();
      if (gmm::vect_size(B_) != mf_data.nb_dof() * qmult) 
	DAL_THROW(failure_error, "The data mesh fem structure has changed, "
		  " You have to change the rhs in that case.");
      gmm::resize(F_, sub_problem.main_mesh_fem().nb_dof());
      gmm::clear(F_);
      
      asm_source_term(F_, sub_problem.main_mesh_fem(),mf_data, B_,boundary);
      this->computed();
    }

    void fixing_dimensions(void) {
      size_type q = mf_data.get_qdim();
      size_type qdim = sub_problem.main_mesh_fem().get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

  public :

    virtual bool is_linear(void) { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) { return sub_problem.is_coercive(); }
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem.mixed_variables(b, i0); }
    virtual size_type nb_dof(void) { return sub_problem.nb_dof(); }
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				     size_type j0 = 0, bool modified = false)
    { sub_problem.compute_tangent_matrix(MS, i0, j0, modified); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      react(MS, i0, false);
      if (this->to_be_computed()) {compute_F();}
      
      gmm::add(gmm::scaled(F_, value_type(-1)), gmm::sub_vector(MS.residu(),
	       gmm::sub_interval(i0, sub_problem.main_mesh_fem().nb_dof())));

    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void set_rhs(const VECTOR &B__)
    { fixing_dimensions(); gmm::copy(B__, B_); }    
    // Constructor defining the rhs
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
		       mesh_fem &mf_data_, const VECTOR &B__,
		       size_type bound = size_type(-1))
      : sub_problem(problem), mf_data(mf_data_), boundary(bound) {
      fixing_dimensions();
      gmm::copy(B__, B_);
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
    }
  };

  /* ******************************************************************** */
  /*		Q.U term (for Fourier-Robin conditions)                   */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_QU_term
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR Q;
    bool homogeneous;
    size_type boundary;
    T_MATRIX K;

    void compute_K(void) {
      mesh_fem &mf_u = sub_problem.main_mesh_fem();
      gmm::clear(K);
      gmm::resize(K, nb_dof(), nb_dof());
      size_type N2 = dal::sqr(mf_u.get_qdim());
      VECTOR vQ(mf_data.nb_dof() * N2);
      if (homogeneous) {
	for (size_type i=0; i < mf_data.nb_dof(); ++i) {
	  for (size_type j=0; j < N2; ++j)
	    vQ[i*N2 + j] = (j % (mf_u.get_qdim()+1)) == 0 ? Q[0] : 0.;
	}
      }
      else gmm::copy(Q, vQ);
      asm_qu_term(K, mf_u, mf_data, vQ, boundary);
      this->computed();
    }

  public :

    virtual bool is_linear(void) { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) { return false; }
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem.mixed_variables(b, i0); }
    virtual size_type nb_dof(void) { return sub_problem.nb_dof(); }    
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0 = 0, bool = false) {
      sub_problem.compute_tangent_matrix(MS, i0, j0, true);
      react(MS, i0, true);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed())
	  compute_K();
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      react(MS, i0, false);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed()) { 
	compute_K();
      }
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residu(), SUBI);
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI), SUBV, SUBV);
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void set_Q(value_type q) {
      homogeneous = true;
      gmm::resize(Q, 1); Q[0] = q;
      this->force_recompute();
    }

    void set_Q(const VECTOR &q) {
      homogeneous = false;
      gmm::resize(Q, mf_data.nb_dof()*dal::sqr(main_mesh_fem().get_qdim())); 
      gmm::copy(q, Q);
      this->force_recompute();
    }
    // Constructor which does not define the rhs
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, value_type q=value_type(1),
		    size_type bound = size_type(-1)) 
      : sub_problem(problem), mf_data(mf_data_), boundary(bound) {
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
      set_Q(value_type(q));
    }

    // Constructor defining the rhs
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, const VECTOR &q,
		    size_type bound = size_type(-1))
      : sub_problem(problem), mf_data(mf_data_), boundary(bound) {
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
      set_Q(q);
    }
  };


  /* ******************************************************************** */
  /*		Mixed linear incompressible condition brick.              */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_linear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_p;
    T_MATRIX B;

    void compute_B() {
      mesh_fem &mf_u = sub_problem.main_mesh_fem();
      size_type nd = mf_u.nb_dof(), ndd = mf_p.nb_dof();
      gmm::resize(B, ndd, nd);
      asm_stokes_B(B, mf_u, mf_p);
      this->computed();
    }

  public :
    
    virtual bool is_linear(void)   { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) { return false; }
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      if (this->context_changed()) {
	this->force_recompute();
	compute_B();
      }
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
      react(MS, i0, modified);
      if (this->to_be_computed()) compute_B();
      
      if (this->to_be_transferred()) {
	gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_p.nb_dof());
	gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());
	gmm::copy(B, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	gmm::copy(gmm::transposed(B),
		  gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
	this->transferred();
      }
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      react(MS, i0, false);
      if (this->to_be_computed()) compute_B();
     
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());
      gmm::mult(B, gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(gmm::transposed(B), gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBJ));
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    // Constructor which does not define the rhs
    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_)
      : sub_problem(problem), mf_p(mf_p_) {
      this->add_dependency(mf_p);
      this->add_dependency(sub_problem.main_mesh_fem());
      compute_B();
    }
  };


  /* ******************************************************************** */
  /*		Dirichlet condition bricks.                               */
  /* ******************************************************************** */
  // TODO : Version with local matrices on the boundary

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Dirichlet : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR B_, H_;
    C_MATRIX G;
    VECTOR CRHS;
    size_type boundary, nb_const;
    dal::bit_vector dof_on_bound;
    bool with_H, with_multipliers;

    void fixing_dimensions(void) {
      size_type q = mf_data.get_qdim();
      size_type qdim = sub_problem.main_mesh_fem().get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      size_type qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

    void compute_constraints(int version) {
      mesh_fem &mf_u = sub_problem.main_mesh_fem();
      size_type Q = mf_u.get_qdim();
      size_type nd = mf_u.nb_dof(), ndd = mf_data.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(nd, nd);
      VECTOR V(nd);

      if (!with_H) {
	gmm::resize(H_, dal::sqr(Q) * ndd);
	for (size_type i=0; i < ndd; ++i)
	  for (size_type q=0; q < Q; ++q)  H_[i*Q*Q+q*Q+q] = value_type(1);
      }
      if (!with_multipliers) version |= ASMDIR_SIMPLIFY;
      gmm::clear(M); gmm::clear(V);
      asm_dirichlet_constraints(M, V, sub_problem.main_mesh_fem(),
				mf_data, H_, B_, boundary, version);

      if (!with_H) gmm::resize(H_, 0);

      R tol=gmm::mat_maxnorm(M)*gmm::default_tol(value_type())*R(100);
      if (version & ASMDIR_BUILDH) gmm::clean(M, tol);
      
      std::vector<size_type> ind(0);
      dof_on_bound = mf_u.dof_on_boundary(boundary);
      dal::bit_vector nn = dof_on_bound;
      // The following filter is not sufficient for an arbitrary matrix field
      // H for the multipliers version. To be ameliorated.
      for (size_type i = nn.take_first(); i != size_type(-1); i << nn)
	if (!with_multipliers || gmm::vect_norm2(gmm::mat_row(M, i)) > tol)
	  ind.push_back(i);
      nb_const = ind.size();
      if (version & ASMDIR_BUILDH) gmm::resize(G, nb_const, nd);
      gmm::sub_index SUBI(ind);
      if (version & ASMDIR_BUILDH) 
	gmm::copy(gmm::sub_matrix(M, SUBI, gmm::sub_interval(0, nd)), G);

      gmm::resize(CRHS, nb_const);
      gmm::copy(gmm::sub_vector(V, SUBI), CRHS);
      this->computed();
    }

  public :

    virtual bool is_linear(void)   { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) 
    { return (!with_multipliers && sub_problem.is_coercive()); }
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      if (with_multipliers) {
	if (this->context_changed()) {
	  fixing_dimensions();
	  this->force_recompute();
	  compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
	}
	b.add(i0 + sub_problem.nb_dof(), nb_const);
      }
    }
    virtual size_type nb_dof(void) {
      if (with_multipliers) {
	if (this->context_changed()) {
	  fixing_dimensions();
	  this->force_recompute();
	  compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
	}
	return sub_problem.nb_dof() + nb_const;
      }
      return sub_problem.nb_dof();
    }
    
    virtual size_type nb_constraints(void) {
      if (with_multipliers) return sub_problem.nb_constraints();
      if (this->context_changed()) {
	fixing_dimensions();
	this->force_recompute();
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      }
      return sub_problem.nb_constraints() + nb_const;
    }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0, bool modified = false) {
      sub_problem.compute_tangent_matrix(MS, i0, j0, modified);
      cout << "dans dirichlet::compute_tangent_matrix()\n";
      react(MS, i0, modified);
       
      if (this->to_be_computed()) {
	fixing_dimensions();
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
	
	cout << "dirichlet : computed\n";
      }
      if (this->to_be_transferred()) {
	cout << "dirrichlet : transfered\n";
	if (with_multipliers) {
	  gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), dof_on_bound.card());
	  gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());
	  gmm::copy(G, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	  gmm::copy(gmm::transposed(G),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
	}
	else {
	  cout << "dirrichlet : !with_mutlipliers\n";
	  
	  size_type nd = sub_problem.main_mesh_fem().nb_dof();
	  size_type ncs = sub_problem.nb_constraints();
	  gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0, nd);
	  gmm::copy(G, gmm::sub_matrix(MS.constraints_matrix(), SUBI, SUBJ));
	}
	this->transferred();
      }
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);

      react(MS, i0, false);
      if (this->to_be_computed()) {
	fixing_dimensions();
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      }
      if (with_multipliers) {

	gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), dof_on_bound.card());
	gmm::sub_interval SUBJ(i0, main_mesh_fem().nb_dof());
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residu(), SUBI));

	gmm::mult_add(gmm::transposed(G), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residu(), SUBJ));
      }
      else {

	size_type nd = sub_problem.main_mesh_fem().nb_dof();
	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0, nd);
	gmm::mult(G, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), 
				 value_type(-1)),
		  CRHS, gmm::sub_vector(MS.constraints_rhs(), SUBI));
      }
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void changing_rhs(const VECTOR &B__) {
      if (this->context_changed()) {
	fixing_dimensions(); gmm::copy(B__, B_);
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      }
      else
	{ gmm::copy(B__, B_); compute_constraints(ASMDIR_BUILDR); }
    }

    // Constructor which does not define the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_data_, size_type bound,
		      bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	with_H(false), with_multipliers(with_mult) {
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
      fixing_dimensions();
      gmm::clear(B_);
      compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
    }

    // Constructor defining the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		     mesh_fem &mf_data_, const VECTOR &B__,
		      size_type bound, bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	with_H(false), with_multipliers(with_mult) {
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
      fixing_dimensions();
      gmm::copy(B__, B_);
      compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
    }
    
  };
  
  /* ******************************************************************** */
  /*		Generic solvers.                                          */
  /* ******************************************************************** */

  // faire une version avec using_cg, using_gmres ... (appelée par celle-ci)
  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
	gmm::iteration &iter) {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename gmm::number_traits<value_type>::magnitude_type mtype;

    // TODO : take iter into account for the Newton. compute a consistent 
    //        max residu.

    size_type ndof = problem.nb_dof();

    bool is_linear = problem.is_linear();
    mtype alpha, alpha_min=mtype(1)/mtype(32), alpha_mult=mtype(3)/mtype(4);
    mtype alpha_max_ratio(1);
    dal::bit_vector mixvar;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.set_maxiter(10000);
    if (!is_linear) { iter_linsolv0.reduce_noisy(); iter_linsolv0.set_resmax(iter.get_resmax()/100000.0); }


    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before
    problem.compute_residu(MS);
    problem.compute_tangent_matrix(MS);

    MS.compute_reduced_system();
    
    mtype act_res = MS.reduced_residu_norm(), act_res_new(0);

    while (is_linear || !iter.finished(act_res)) {
    
      gmm::iteration iter_linsolv = iter_linsolv0;
      VECTOR d(ndof), dr(gmm::vect_size(MS.reduced_residu()));

      if (!(iter.first())) {
	problem.compute_tangent_matrix(MS);
	MS.compute_reduced_system();
      }

      cout << "tangent matrix is "
	   << (gmm::is_symmetric(MS.tangent_matrix()) ? "" : "not ")
	   <<  "symmetric. ";

#ifdef GMM_USES_SUPERLU
	  
      double rcond;
      SuperLU_solve(MS.reduced_tangent_matrix(), dr,
		    gmm::scaled(MS.reduced_residu(), value_type(-1)),
		    rcond);
#else
      if (problem.is_coercive()) {
	gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	gmm::cg(MS.reduced_tangent_matrix(), dr, 
		gmm::scaled(MS.reduced_residu(), value_type(-1)),
		P, iter_linsolv);
	if (!iter_linsolv.converged()) DAL_WARNING(2,"cg did not converge!");
      } else {
	problem.mixed_variables(mixvar);
	if (mixvar.card() == 0) {
	  gmm::ilu_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	  
	  gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		     gmm::scaled(MS.reduced_residu(),  value_type(-1)), P,
		     300, iter_linsolv);
	}
	else {
	  cout << "there is " << mixvar.card() << " mixed variables\n";
	  gmm::ilut_precond<T_MATRIX> P(MS.reduced_tangent_matrix(),100,1E-10);
	  // gmm::identity_matrix P;
	  gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		     gmm::scaled(MS.reduced_residu(),  value_type(-1)),
		     P, 300, iter_linsolv);
	}
	if (!iter_linsolv.converged()) DAL_WARNING(2,"gmres did not converge!");
      }
#endif
      MS.unreduced_solution(dr,d);

      if (is_linear) {
	gmm::add(d, MS.state());
	return;
      }
      else { // line search for the non-linear case.
	VECTOR stateinit(ndof);
	gmm::copy(MS.state(), stateinit);
       
	for (alpha = mtype(1); alpha >= alpha_min; alpha *= alpha_mult) {
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  problem.compute_residu(MS);
	  MS.compute_reduced_system(); // The whole reduced system do not
	  // have to be computed, only the RHS. To be adapted.
	  act_res_new = MS.reduced_residu_norm();
	  if (act_res_new <= act_res * alpha_max_ratio) break;
	}
      }
      act_res = act_res_new; ++iter;

      if (iter.get_noisy()) cout << "alpha = " << alpha << " ";
    }

  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
