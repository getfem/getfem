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
/*  - compute_tangent_matrix(MS, i0, modified) : the brick has to compute  */
/*        its own part of the tangent matrix (i0 if the shift in the       */
/*        tangent matrix defined in MS) and has to call the                */
/*        compute_tangent_matrix(MS, i0,modified) of sub-problem(s) if any.*/
/*        modified indicates if the part of the tangent matrix dedicated   */
/*        to this brick will be modified by other bricks.                  */
/*  - compute_residu(MS, i0) : the brick has to compute its own            */
/*        part of the residu (i0 if the shift in the residu vector defined */
/*        in MS) and has to call the compute_residu(MS, i0, modified) of   */
/*        sub-problem(s) if any.                                           */
/*  - constraints_system(MS, i0, j0, modified) : the brick has to compute  */
/*        its own part of the linear constraints system and has to call    */
/*        the constraints_system(MS, i0, j0,modified) of sub-problem(s)    */
/*        if any. modified indicates if the part of the constraint system  */
/*        dedicated to this brick will be modified by other bricks.        */
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
    typedef typename gmm::number_traits<value_type>::magnitude_type magnitude_type;

  protected :
    T_MATRIX tangent_matrix_;
    C_MATRIX constraints_matrix_;
    VECTOR state_, residu_, constraints_rhs_;
    ctx_ident_type ident_;
    
    T_MATRIX SM;
    gmm::col_matrix< gmm::rsvector<value_type> > NS;
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
    void compute_reduced_system() {
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
    VECTOR &residu(void) { return residu_; }
    ctx_ident_type ident(void) { return ident_; }
    void touch(void) { ident_ = context_dependencies::new_ident(); }
    void adapt_sizes(mdbrick_abstract<model_state> &problem) {
      size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();
      gmm::resize(tangent_matrix_, ndof, ndof);
      gmm::resize(constraints_matrix_, nc, ndof);
      gmm::resize(constraints_rhs_, nc);
      gmm::resize(state_, ndof);
      gmm::resize(residu_, ndof);
      touch();
    } 

    model_state(void) { ident_ = context_dependencies::new_ident(); }
  };

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
		      modeling_standard_complex_plain_vector > standard_complex_model_state;


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

    void compute_K(void) {
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

  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return true; }
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

  /* TODO : arbitrary elasticity tensor */


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
      if (homogeneous) {
	std::fill(wave_number2.begin(), wave_number2.end(), value_type(dal::sqr(wave_number[0])));
      }
      else { 
	for (size_type i=0; i < nb_dof(); ++i)
	  wave_number2[i] = dal::sqr(wave_number[i]);
      }
      
      asm_Helmholtz(K, mf_u, mf_data, wave_number2);
      this->computed();
    }

  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return false; }
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
      if (this->to_be_computed()) compute_F();
      gmm::add(gmm::scaled(F_, value_type(-1)), gmm::sub_vector(MS.residu(),
	       gmm::sub_interval(i0, sub_problem.main_mesh_fem().nb_dof())));
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void changing_rhs(const VECTOR &B__)
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
    virtual size_type nb_dof(void) { return sub_problem.nb_dof(); }    
    virtual size_type nb_constraints(void) { return sub_problem.nb_constraints(); }
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
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }
    virtual mesh_fem &main_mesh_fem(void) { return sub_problem.main_mesh_fem(); }

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
		    mesh_fem &mf_data_, value_type q=value_type(1), size_type bound = size_type(-1)) 
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
      react(MS, i0, modified);
      if (this->to_be_computed()) {
	fixing_dimensions();
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      }
      if (this->to_be_transferred()) {
	if (with_multipliers) {
	  gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), dof_on_bound.card());
	  gmm::sub_interval SUBJ(i0, sub_problem.nb_dof());
	  gmm::copy(G, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	  gmm::copy(gmm::transposed(G),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	}
	else {
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
	gmm::sub_interval SUBJ(i0, sub_problem.nb_dof());
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residu(), SUBI));
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
    //        detect the presence of multipliers before using a preconditioner

    size_type ndof = problem.nb_dof();
    bool is_linear = problem.is_linear();
    mtype alpha, alpha_min=mtype(1)/mtype(16), alpha_mult=mtype(3)/mtype(4);
    mtype alpha_max_ratio(2);

    MS.adapt_sizes(problem);
    if (!is_linear) gmm::fill_random(MS.state()); 
    else gmm::clear(MS.state());
    problem.compute_residu(MS);
    problem.compute_tangent_matrix(MS);
    MS.compute_reduced_system();
    mtype act_res = MS.reduced_residu_norm(), act_res_new(0);
    while(act_res > iter.get_resmax()) {
      
      VECTOR d(ndof), dr(gmm::vect_size(MS.reduced_residu()));

      if (!(iter.first())) problem.compute_tangent_matrix(MS);
      if (problem.is_coercive()) {
	gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	gmm::cg(MS.reduced_tangent_matrix(), dr, 
		gmm::scaled(MS.reduced_residu(), value_type(-1)), P, iter);
      } else {
	gmm::ilu_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		   gmm::scaled(MS.reduced_residu(),  value_type(-1)), P, 300, iter);
      }
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
	  MS.compute_reduced_system();
	  act_res_new = MS.reduced_residu_norm();
	  if (act_res_new <= act_res * alpha_max_ratio) break;
	}
      }
      act_res = act_res_new;
    }
  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
