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
/*        sub-problem(s) if any : to be redefine olny if there is          */
/*        additional dof not coming from proper mesh_fems                  */
/*  - nb_constraints() : number of linear constraints on the system        */
/*        including the constraints defined in the sub-problem(s) if any.  */
/*                     coercive.                                           */
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
/*                                                                         */
/* Dependencies.                                                           */
/*   A brick depends at least on some mesh_fem structures and has to       */
/* react if some changements occur in these mesh_fem structures.           */
/*                                                                         */
/***************************************************************************/
/*                                                                         */
/* Brick idents :                                                          */
/* MDBRICK_SCALAR_ELLIPTIC       174397                                    */
/* MDBRICK_LIN_ISO_ELASTICITY    852327                                    */
/* MDBRICK_HELMHOLTZ             354864                                    */
/* MDBRICK_LINEAR_INCOMP         239898                                    */
/* MDBRICK_NONLINEAR_ELASTICITY  821357                                    */
/* MDBRICK_NONLINEAR_INCOMP      964552                                    */
/* MDBRICK_SMALL_DEF_PLASTICITY  556433                                    */
/* MDBRICK_LINEAR_PLATE          897523                                    */
/* MDBRICK_MIXED_LINEAR_PLATE    213456                                    */
/* MDBRICK_COULOMB_FRICTION      434245                                    */
/*                                                                         */
/***************************************************************************/


#ifndef GETFEM_MODELING_H__
#define GETFEM_MODELING_H__

#include <getfem_assembling_tensors.h>
#include <getfem_assembling.h>
#include <gmm_solver_cg.h>
#include <gmm_solver_gmres.h>
#include <gmm_precond_ildlt.h>
#include <gmm_precond_ilu.h>
#include <gmm_precond_ilut.h>
#include <gmm_precond_ilutp.h>
#include <gmm_superlu_interface.h>
#include <gmm_dense_qr.h>

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
    long ident_;
    
    T_MATRIX SM;
    gmm::col_matrix<gmm::rsvector<value_type> > NS; /* constraints nullspace */
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
    void compute_reduced_residu();
    VECTOR &residu(void) { return residu_; }
    long ident(void) { return ident_; }
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
    
    size_type nbcols=Dirichlet_nullspace(constraints_matrix(),
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
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::compute_reduced_residu() {
    // The call to Dirichlet nullspace should be avoided -> we just need Ud
    if (gmm::mat_nrows(constraints_matrix()) == 0) return;
    size_type ndof = gmm::mat_ncols(tangent_matrix());
    gmm::resize(NS, ndof, ndof);
    size_type nbcols=Dirichlet_nullspace(constraints_matrix(),
					 NS, constraints_rhs(), Ud);
    gmm::resize(NS, ndof, nbcols);
    VECTOR RHaux(ndof);
    gmm::mult(tangent_matrix(), Ud, residu(), RHaux);
    gmm::mult(gmm::transposed(NS), RHaux, reduced_residu_);
  }

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::
  adapt_sizes(mdbrick_abstract<model_state> &problem) {
    size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();

    if (gmm::mat_nrows(tangent_matrix_) != ndof
	|| gmm::mat_nrows(constraints_matrix_) != nc) {
      gmm::clear(state_);
      gmm::clear(residu_);
      gmm::clear(tangent_matrix_);
      gmm::clear(constraints_matrix_);
      gmm::clear(constraints_rhs_);
      gmm::resize(tangent_matrix_, ndof, ndof);
      gmm::resize(constraints_matrix_, nc, ndof);
      gmm::resize(constraints_rhs_, nc);
      gmm::resize(state_, ndof);
      gmm::resize(residu_, ndof);
      touch();
    }
  } 


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

  enum bound_cond_type { MDBRICK_UNDEFINED, MDBRICK_DIRICHLET, MDBRICK_NEUMANN,
			 MDBRICK_SIMPLE_SUPPORT, MDBRICK_CLAMPED_SUPPORT,
			 MDBRICK_FOURIER_ROBIN };

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_abstract : public context_dependencies {

  public :

    struct mesh_fem_info_ {
      size_type brick_ident; // basic model brick using the mesh_fem
      size_type info;        // flags
      // type of boundary conditions
      std::map<size_type, bound_cond_type> boundaries;
      mesh_fem_info_(size_type id, size_type in) : brick_ident(id), info(in) {}
      void add_boundary(size_type b, bound_cond_type bc)
      { boundaries[b] = bc; }
      bound_cond_type boundary_type(size_type b) {
	typename std::map<size_type, bound_cond_type>::const_iterator it;
	it = boundaries.find(b);
	return it != boundaries.end() ? it->second : MDBRICK_UNDEFINED;
      }
    };

  protected :

    struct boundary_cond_info_ {
      size_type num_fem, num_bound;
      bound_cond_type bc;
      boundary_cond_info_(size_type a, size_type b, bound_cond_type d)
	: num_fem(a), num_bound(b), bc(d) {}
    };

    mutable bool to_compute, to_transfer;
    size_type MS_i0;
    long ident_ms;

    std::vector<mdbrick_abstract *> sub_bricks;
    mutable std::vector<mesh_fem *> mesh_fems;
    mutable std::vector<mesh_fem_info_> mesh_fems_info;
    mutable std::vector<size_type> mesh_fem_positions;
    std::vector<mesh_fem *> proper_mesh_fems;
    std::vector<mesh_fem_info_> proper_mesh_fems_info;
    std::vector<boundary_cond_info_> proper_boundary_cond_info;


    bool proper_is_linear_, proper_is_symmetric_, proper_is_coercive_;
    mutable bool is_linear_, is_symmetric_, is_coercive_;
    mutable size_type nb_total_dof;

    void update_from_context(void) const {
      to_compute = true;
      nb_total_dof = 0;
      is_linear_ = proper_is_linear_;
      is_symmetric_ = proper_is_symmetric_;
      is_coercive_ = proper_is_coercive_;
      mesh_fems.resize(0); mesh_fem_positions.resize(0);
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	for (size_type j = 0; j < sub_bricks[i]->mesh_fems.size(); ++j) {
	  mesh_fems.push_back(sub_bricks[i]->mesh_fems[j]);
	  mesh_fems_info.push_back(sub_bricks[i]->mesh_fems_info[j]);
	  mesh_fem_positions.push_back(nb_total_dof 
				       + sub_bricks[i]->mesh_fem_positions[j]);
	  is_linear_ = is_linear_ && sub_bricks[i]->is_linear();
	  is_symmetric_ = is_symmetric_ && sub_bricks[i]->is_symmetric();
	  is_coercive_ = is_coercive_ && sub_bricks[i]->is_coercive();
	}
	nb_total_dof += sub_bricks[i]->nb_dof();
      }
      for (size_type j = 0; j < proper_mesh_fems.size(); ++j) {
	mesh_fems.push_back(proper_mesh_fems[j]);
	mesh_fems_info.push_back(proper_mesh_fems_info[j]);
	mesh_fem_positions.push_back(nb_total_dof);
	nb_total_dof += proper_mesh_fems[j]->nb_dof();
      }
      for (size_type j = 0; j < proper_boundary_cond_info.size(); ++j) {
	mesh_fems_info[proper_boundary_cond_info[j].num_fem]
	  .add_boundary(proper_boundary_cond_info[j].num_bound,
			proper_boundary_cond_info[j].bc);
      }
    }

    void add_sub_brick(mdbrick_abstract &mdb) {
      sub_bricks.push_back(&mdb);
      add_dependency(mdb);
    }

    void add_proper_mesh_fem(mesh_fem &mf, size_type brick_ident,
			     size_type info = 0) {
      mesh_fem_info_ mfi(brick_ident, info);
      proper_mesh_fems.push_back(&mf);
      proper_mesh_fems_info.push_back(mfi);
      add_dependency(mf);
    }

    void add_proper_boundary_info(size_type num_fem, size_type num_bound,
				  bound_cond_type bc) {
      boundary_cond_info_ bci(num_fem, num_bound, bc);
      proper_boundary_cond_info.push_back(bci);
    }

    bound_cond_type boundary_type(size_type num_fem, size_type num_bound)
    { return mesh_fems_info[num_fem].boundary_type(num_bound); }

    // to_be_computed : the context has changed (or it is the first call).
    // to_be_transferred : the structure MODEL_STATE has changed, the 
    //                     tangent matrix or constraints system has to be
    //                     copied again if it is stored.
    void react(MODEL_STATE &MS, size_type i0, bool modified) {
      this->context_check();
      to_transfer = to_transfer || modified || (ident_ms != MS.ident());
      ident_ms = MS.ident();
      MS_i0 = i0;
    }

    bool to_be_computed(void) { return to_compute; }
    bool to_be_transferred(void) { return to_transfer; }
    void force_recompute(void) { to_compute = to_transfer = true; }
    void computed(void) { to_compute = false; to_transfer = true; }
    void transferred(void) { to_transfer = false; }
    size_type first_index(void) { return MS_i0; }

  public :

    mesh_fem_info_ &get_mesh_fem_info(size_type i)
    { return mesh_fems_info[i]; }
    mesh_fem &get_mesh_fem(size_type i) { return *(mesh_fems[i]); }
    size_type nb_mesh_fems(void) { return mesh_fems.size(); }

    dim_type dim(void) { return mesh_fems[0]->linked_mesh().dim(); }
    virtual size_type nb_dof(void) { return nb_total_dof; }
    virtual size_type nb_constraints(void) = 0;
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				  size_type j0=0, bool modified = false) = 0;
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) = 0;
    bool is_linear(void) { return is_linear_; }
    bool is_symmetric(void) { return is_symmetric_; }
    bool is_coercive(void) { return is_coercive_; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) = 0;
    mdbrick_abstract(void) : to_compute(true), to_transfer(true),
			     MS_i0(0), ident_ms(-1)
    { proper_is_linear_ = proper_is_symmetric_ = proper_is_coercive_ = true; }
    virtual ~mdbrick_abstract() {}
  };

  /* ******************************************************************** */
  /*		general scalar elliptic brick.                            */
  /* ******************************************************************** */

# define MDBRICK_SCALAR_ELLIPTIC 174397

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_scalar_elliptic : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::sub_vector_type<VECTOR *,
				 gmm::sub_interval>::vector_type SUBVECTOR;

    gmm::sub_interval SUBU;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR coeffs_;
    bool homogeneous, laplacian;
    bool matrix_stored;
    T_MATRIX K;
    

    void compute_K(void);

  public :

    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false);
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0);
    
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

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      SUBU = gmm::sub_interval (this->first_index(), this->nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_u, MDBRICK_SCALAR_ELLIPTIC);
      this->update_from_context();
    }

    // constructor for the Laplace operator
    mdbrick_scalar_elliptic(mesh_fem &mf_u_, mesh_fem &mf_data_,
       value_type a, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored)
    { set_coeff(a); init_(); }

    // constructor for a non-homogeneous material
    mdbrick_scalar_elliptic(mesh_fem &mf_u_, mesh_fem &mf_data_,
       const VECTOR &coeff, bool laplace, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored)
    { set_coeff(coeff, laplace); init_(); }
  };

  template<typename MODEL_STATE>
   void mdbrick_scalar_elliptic<MODEL_STATE>::compute_K(void) {
    gmm::clear(K);
    gmm::resize(K, this->nb_dof(), this->nb_dof());
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
    gmm::sub_interval SUBI(i0, this->nb_dof());
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
    gmm::sub_interval SUBI(i0, this->nb_dof());
    if (this->to_be_computed()) { 
      compute_K();
      if (!matrix_stored) {
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	gmm::clear(K);
	this->transferred();
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

# define MDBRICK_LIN_ISO_ELASTICITY 852327

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_isotropic_linearized_elasticity
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::sub_vector_type<VECTOR *,
				 gmm::sub_interval>::vector_type SUBVECTOR;

    gmm::sub_interval SUBU;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR lambda_, mu_;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void);

  public :

    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false);
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0);

    const T_MATRIX &stiffness_matrix(void)
    { if (this->to_be_computed()) compute_K(); return K; }

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

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      SUBU = gmm::sub_interval (this->first_index(), this->nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_u, MDBRICK_LIN_ISO_ELASTICITY);
      this->update_from_context();
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_isotropic_linearized_elasticity
    (mesh_fem &mf_u_, mesh_fem &mf_data_,
     value_type lambdai, value_type mui, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored)
    { set_Lame_coeff(lambdai, mui); init_(); }

    // constructor for a non-homogeneous material
    mdbrick_isotropic_linearized_elasticity
    (mesh_fem &mf_u_, mesh_fem &mf_data_,
     const VECTOR &lambdai, const VECTOR &mui, bool mat_stored = false)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored)
    { set_Lame_coeff(lambdai, mui); init_(); }
  };


  template<typename MODEL_STATE>
   void mdbrick_isotropic_linearized_elasticity<MODEL_STATE>::compute_K(void) {
    gmm::clear(K);
    gmm::resize(K, this->nb_dof(), this->nb_dof());
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
  void mdbrick_isotropic_linearized_elasticity<MODEL_STATE>::
  compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
			 size_type, bool modified) {
    if (modified && !matrix_stored) 
      DAL_THROW(failure_error, "The residu will not be consistant. "
		"Use this brick with the stiffness matrix stored option");
    react(MS, i0, modified);
    gmm::sub_interval SUBI(i0, this->nb_dof());
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
  void mdbrick_isotropic_linearized_elasticity<MODEL_STATE>::
  compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
    react(MS, i0, false);
    gmm::sub_interval SUBI(i0, this->nb_dof());
    if (this->to_be_computed()) { 
      compute_K();
      if (!matrix_stored) {
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	gmm::clear(K);
	this->transferred();
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
  /*		Helmholtz brick.                                          */
  /* ******************************************************************** */

# define MDBRICK_HELMHOLTZ 354864

 template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Helmholtz
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::sub_vector_type<VECTOR *,
				 gmm::sub_interval>::vector_type SUBVECTOR;

    gmm::sub_interval SUBU;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR wave_number;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void) {
      gmm::clear(K);
      gmm::resize(K, this->nb_dof(), this->nb_dof());
      VECTOR wave_number2(mf_data.nb_dof());
      if (homogeneous)
	std::fill(wave_number2.begin(), wave_number2.end(),
		  value_type(dal::sqr(wave_number[0])));
      else
	for (size_type i=0; i < this->nb_dof(); ++i)
	  wave_number2[i] = dal::sqr(wave_number[i]);
      
      asm_Helmholtz(K, mf_u, mf_data, wave_number2);
      this->computed();
    }

  public :

    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false) {
      if (modified && !matrix_stored) 
	DAL_THROW(failure_error, "The residu will not be consistant. "
		  "Use this brick with the stiffness matrix stored option");
      react(MS, i0, modified);
      gmm::sub_interval SUBI(i0, this->nb_dof());
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
      gmm::sub_interval SUBI(i0, this->nb_dof());
      if (this->to_be_computed()) {
	compute_K();
	if (!matrix_stored) {
	  gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	  gmm::clear(K);
	  this->transferred();
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

    void set_wave_number(value_type k) {
      homogeneous = true;
      gmm::resize(wave_number, 1); wave_number[0] = k;
      this->force_recompute();
    }

    void set_wave_number(const VECTOR &k) {
      homogeneous = false;
      gmm::resize(wave_number, mf_data.nb_dof()); gmm::copy(k, wave_number);
      this->force_recompute();
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      SUBU = gmm::sub_interval (this->first_index(), this->nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_u, MDBRICK_HELMHOLTZ);
      this->proper_is_coercive_ = false;
      this->update_from_context();
    }

    // constructor for a homogeneous wave number
    mdbrick_Helmholtz(mesh_fem &mf_u_, mesh_fem &mf_data_,
       value_type k, bool mat_stored = true)
      : mf_u(mf_u_), mf_data(mf_data_), matrix_stored(mat_stored)
    { set_wave_number(k); init_(); }

    // constructor for a non-homogeneous wave number
    mdbrick_Helmholtz(mesh_fem &mf_u_, mesh_fem &mf_data_,
		      const VECTOR &k, bool mat_stored = true)
      : mf_u(mf_u_), mf_data(mf_data_),	matrix_stored(mat_stored)
    { set_wave_number(k); init_(); }
    
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
    size_type boundary, qmult, num_fem;
    size_type i1, nbd;

    void compute_F(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      qmult = mf_u.get_qdim() / mf_data.get_qdim();
      if (gmm::vect_size(B_) != mf_data.nb_dof() * qmult) 
	DAL_THROW(failure_error, "The data mesh fem structure has changed, "
		  " You have to change the rhs in that case.");
      gmm::resize(F_, mf_u.nb_dof());
      gmm::clear(F_);
      asm_source_term(F_, mf_u, mf_data, B_, boundary);
      this->computed();
    }

    void fixing_dimensions(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type q = mf_data.get_qdim();
      size_type qdim = mf_u.get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

  public :

    const VECTOR &source_term(void) {
      if (this->to_be_computed()) {compute_F();}
      return F_;
    }

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem.mixed_variables(b, i0); }
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
	       gmm::sub_interval(i0+i1, nbd)));

    }

    void set_rhs(const VECTOR &B__)
    { fixing_dimensions(); gmm::copy(B__, B_); this->force_recompute(); }    
    // Constructor defining the rhs
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			mesh_fem &mf_data_, const VECTOR &B__,
			size_type bound = size_type(-1), size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_) {
      this->add_dependency(mf_data);
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->update_from_context();
      fixing_dimensions();
       gmm::copy(B__, B_);
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
    size_type boundary, num_fem, i1, nbd;
    T_MATRIX K;

    void compute_K(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      gmm::clear(K);
      gmm::resize(K, this->nb_dof(), this->nb_dof());
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

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem.mixed_variables(b, i0); }
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0 = 0, bool = false) {
      sub_problem.compute_tangent_matrix(MS, i0, j0, true);
      react(MS, i0, true);
      if (this->to_be_computed()) compute_K();
      gmm::sub_interval SUBI(i0+i1, nbd);
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      react(MS, i0, false);
      if (this->to_be_computed()) compute_K();
      gmm::sub_interval SUBI(i0+i1, nbd);
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residu(), SUBI);
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI), SUBV, SUBV);
    }

    void set_Q(value_type q) {
      homogeneous = true;
      gmm::resize(Q, 1); Q[0] = q;
      this->force_recompute();
    }

    void set_Q(const VECTOR &q) {
      homogeneous = false;
      gmm::resize(Q, mf_data.nb_dof()*dal::sqr(this->mesh_fems[num_fem]->get_qdim())); 
      gmm::copy(q, Q);
      this->force_recompute();
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      if (boundary != size_type(-1))
	this->add_proper_boundary_info(num_fem, boundary,MDBRICK_FOURIER_ROBIN);
      this->update_from_context();
    }

    // Constructor which homogeneous diagonal Q
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, value_type q=value_type(1),
		    size_type bound = size_type(-1), size_type num_fem_=0) 
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_)
    { set_Q(q); init_(); }

    // Constructor defining an arbitrary Q
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, const VECTOR &q,
		    size_type bound = size_type(-1), size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_)
    { set_Q(q); init_(); }
  };


  /* ******************************************************************** */
  /*		Mixed linear incompressible condition brick.              */
  /* ******************************************************************** */

# define MDBRICK_LINEAR_INCOMP 239898

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_linear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_p, &mf_data;
    T_MATRIX B, M;
    bool penalized, homogeneous;
    VECTOR epsilon_; // penalization coefficient if any.
    size_type num_fem, i1, nbd;

    void compute_B() {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      size_type nd = mf_u.nb_dof(), ndd = mf_p.nb_dof();
      gmm::clear(B); gmm::resize(B, ndd, nd);
      asm_stokes_B(B, mf_u, mf_p);
 
//        gmm::dense_matrix<value_type> MM(ndd, ndd);
//        std::vector<value_type> eigval(ndd);
//        gmm::mult(B, gmm::transposed(B), MM);
//        gmm::symmetric_qr_algorithm(MM, eigval);
//        std::sort(eigval.begin(), eigval.end(), Esort);
//        cout << "eival of BBT = " << eigval << endl;

      if (penalized) {
	VECTOR epsilon(mf_data.nb_dof());
	if (homogeneous) std::fill(epsilon.begin(), epsilon.end(),epsilon_[0]);
	else gmm::copy(epsilon_, epsilon);
	gmm::clear(M); gmm::resize(M, ndd, ndd);
	asm_mass_matrix_param(M, mf_p, mf_data, epsilon);
	gmm::scale(M, value_type(-1));
      }
      this->computed();
    }

  public :
    
    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      this->context_check();
      if (this->to_be_computed()) {
	this->force_recompute();
	compute_B();
      }
      b.add(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
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
	gmm::sub_interval SUBJ(i0+i1, nbd);
	gmm::copy(B, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	gmm::copy(gmm::transposed(B),
		  gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	if (penalized)
	  gmm::copy(M, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
	else
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
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::mult(B, gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(gmm::transposed(B), gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBJ));
    }

     void set_coeff(value_type epsiloni) {
      homogeneous = true;
      gmm::resize(epsilon_, 1); epsilon_[0] = epsiloni;
      this->force_recompute();
    }

    void set_coeff(const VECTOR &epsiloni) {
      homogeneous = false;
      gmm::resize(epsilon_, mf_data.nb_dof()); gmm::copy(epsiloni, epsilon_);
      this->force_recompute();
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_p, MDBRICK_LINEAR_INCOMP);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->update_from_context();
    }

    // Constructor for the incompressibility condition
    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), mf_data(mf_p_),
	penalized(false), num_fem(num_fem_)
    { init_(); compute_B(); }

    // Constructor for the nearly incompressibility condition
    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_, mesh_fem &mf_data_, value_type epsilon,
			  size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), mf_data(mf_data_),
	penalized(true), num_fem(num_fem_)
    { set_coeff(epsilon); init_(); compute_B(); }

    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
			  mesh_fem &mf_p_, mesh_fem &mf_data_,
			  const VECTOR& epsilon, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), mf_data(mf_data_),
	penalized(true), num_fem(num_fem_)
    { set_coeff(epsilon); init_(); compute_B(); }

  };


  /* ******************************************************************** */
  /*		Dirichlet condition bricks.                               */
  /* ******************************************************************** */

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
    size_type boundary, nb_const, num_fem;
    bool with_H, with_multipliers;
    gmm::sub_index SUB_CT;
    size_type i1, nbd;

    void fixing_dimensions(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type q = mf_data.get_qdim();
      size_type qdim = mf_u.get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      size_type qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

    void compute_constraints(int version) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      size_type Q = mf_u.get_qdim();
      size_type nd = mf_u.nb_dof(), ndd = mf_data.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(nd, nd);
      VECTOR V(nd);

      if (!with_H) {
	gmm::resize(H_, dal::sqr(Q) * ndd);
	gmm::clear(H_);
	for (size_type i=0; i < ndd; ++i)
	  for (size_type q=0; q < Q; ++q)  H_[i*Q*Q+q*Q+q] = value_type(1);
      }
      if (!with_multipliers) version |= ASMDIR_SIMPLIFY;
      asm_dirichlet_constraints(M, V, mf_u, mf_data,H_, B_, boundary, version);

      if (!with_H) gmm::resize(H_, 0);

      if (version & ASMDIR_BUILDH) {
	R tol=gmm::mat_maxnorm(M)*gmm::default_tol(value_type())*R(100);
	gmm::clean(M, tol);
	std::vector<size_type> ind_ct;
	dal::bit_vector nn = mf_u.dof_on_boundary(boundary);
	// The following filter is not sufficient for an arbitrary matrix field
	// H for the multipliers version. To be ameliorated.
	for (size_type i = nn.take_first(); i != size_type(-1); i << nn)
	  if (!with_multipliers || gmm::vect_norm2(gmm::mat_row(M, i)) > tol)
	    ind_ct.push_back(i);
	nb_const = ind_ct.size();
	SUB_CT = gmm::sub_index(ind_ct);
	gmm::resize(G, nb_const, nd);
	gmm::copy(gmm::sub_matrix(M, SUB_CT, gmm::sub_interval(0, nd)), G);
      }

      gmm::resize(CRHS, nb_const);
      gmm::copy(gmm::sub_vector(V, SUB_CT), CRHS);
      this->computed();
    }

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      this->context_check();
      if (with_multipliers) {
	if (this->to_be_computed()) {
	  fixing_dimensions();
	  this->force_recompute();
	  compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
	}
	b.add(i0 + sub_problem.nb_dof(), nb_const);
      }
    }
    virtual size_type nb_dof(void) {
      this->context_check();
      if (with_multipliers) {
	if (this->to_be_computed()) {
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
      this->context_check();
      if (this->to_be_computed()) {
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
	  gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), nb_const);
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  gmm::copy(G, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	  gmm::copy(gmm::transposed(G),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	  gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
	}
	else {	  
	  size_type ncs = sub_problem.nb_constraints();
	  gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0+i1, nbd);
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

	gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), nb_const);
	gmm::sub_interval SUBJ(i0+i1, nbd);
	
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residu(), SUBI));

	gmm::mult_add(gmm::transposed(G), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residu(), SUBJ));
      }
      else {

	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0+i1, nbd);
	gmm::mult(G, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), 
				 value_type(-1)),
		  CRHS, gmm::sub_vector(MS.constraints_rhs(), SUBI));
      }
    }

    void set_rhs(const VECTOR &B__) {
      this->context_check();
      if (this->to_be_computed()) {
	fixing_dimensions(); gmm::copy(B__, B_);
	compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      }
      else
	{ gmm::copy(B__, B_); compute_constraints(ASMDIR_BUILDR); }
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = !with_multipliers;
      this->add_proper_boundary_info(num_fem, boundary, MDBRICK_DIRICHLET);
      this->update_from_context();
      fixing_dimensions();
    }

    // Constructor which does not define the rhs (0 rhs in fact)
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_data_, size_type bound,
		      size_type num_fem_=0, bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), with_H(false), with_multipliers(with_mult) {
      init_(); gmm::clear(B_);
      compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
    }

    // Constructor defining the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_data_, const VECTOR &B__,
		      size_type bound, size_type num_fem_=0, 
		      bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), with_H(false), with_multipliers(with_mult) {
      init_(); gmm::copy(B__, B_);
      compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
    }
    
  };

  /* ******************************************************************** */
  /*	                 Constraint brick.                                */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_constraint : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    C_MATRIX G;
    VECTOR CRHS;
    size_type num_fem;
    bool with_multipliers;

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      sub_problem.mixed_variables(b, i0);
      if (with_multipliers) b.add(i0+sub_problem.nb_dof(), gmm::mat_nrows(G));
    }
    virtual size_type nb_dof(void)
    { return sub_problem.nb_dof()+(with_multipliers ? gmm::mat_nrows(G) : 0); }
    
    virtual size_type nb_constraints(void) {
      if (with_multipliers) return sub_problem.nb_constraints();
      return sub_problem.nb_constraints() + gmm::mat_nrows(G);
    }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0, bool modified = false) {

      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

      sub_problem.compute_tangent_matrix(MS, i0, j0, modified);
      react(MS, i0, modified);
       
      if (this->to_be_transferred()) {
	if (with_multipliers) {
	  gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), gmm::mat_nrows(G));
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  gmm::copy(G, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	  gmm::copy(gmm::transposed(G),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	}
	else {	  
	  size_type ncs = sub_problem.nb_constraints();
	  gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(G)), SUBJ(i0+i1, nbd);
	  gmm::copy(G, gmm::sub_matrix(MS.constraints_matrix(), SUBI, SUBJ));
	}
	this->transferred();
      }
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

      sub_problem.compute_residu(MS, i0, j0);
      if (with_multipliers) {

	gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), gmm::mat_nrows(G));
	gmm::sub_interval SUBJ(i0+i1, nbd);
	
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residu(), SUBI));

	gmm::mult_add(gmm::transposed(G), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residu(), SUBJ));
      }
      else {

	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(G)), SUBJ(i0+i1, nbd);
	gmm::mult(G, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), 
				 value_type(-1)),
		  CRHS, gmm::sub_vector(MS.constraints_rhs(), SUBI));
      }
    }

    template <class MAT, class VEC>
    void set_constraints(const MAT &G_, const VEC &RHS) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::resize(G, gmm::mat_nrows(G_), mf_u.nb_dof());
      gmm::resize(CRHS, gmm::mat_nrows(G_));
      gmm::copy(G_, G); gmm::copy(RHS, CRHS);
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = !with_multipliers;
      this->update_from_context();
    }

    template <class MAT, class VEC>
    mdbrick_constraint(mdbrick_abstract<MODEL_STATE> &problem,
		       const MAT &G_, const VEC &RHS,		       
		       size_type num_fem_=0, 
		       bool with_mult = false)
      : sub_problem(problem), num_fem(num_fem_), with_multipliers(with_mult)
    { init_(); set_constraints(G_, RHS); }
    
  };
  
  /* ******************************************************************** */
  /*  dynamic brick : not stabilized, could change a lot in the future    */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data, *mf_u;
    VECTOR RHO_, DF;
    T_MATRIX M_;
    size_type num_fem;
    value_type Mcoef, Kcoef;
    gmm::unsorted_sub_index SUBS;
    std::vector<size_type> ind;
    bool have_subs;

    void compute_M(void) {
      gmm::clear(M_); gmm::resize(M_, mf_u->nb_dof(), mf_u->nb_dof());
      asm_mass_matrix_param(M_, *mf_u, mf_data, RHO_);
      if (have_subs) {
	gmm::sub_interval SUBI(0, mf_u->nb_dof());
	gmm::clear(gmm::sub_matrix(M_, SUBS, SUBI));
	gmm::clear(gmm::sub_matrix(M_, SUBI, SUBS));	
      }
      this->computed();
    }

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem.mixed_variables(b, i0); }
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0=0, bool modified=false) {
      sub_problem.compute_tangent_matrix(MS, i0, j0, true);
      react(MS, i0, modified);
      if (this->to_be_computed()) { compute_M(); }
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1)) gmm::scale(MS.tangent_matrix(), Kcoef);
      gmm::add(gmm::scaled(M_, Mcoef),
	       gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type j0 = 0) {
      sub_problem.compute_residu(MS, i0, j0);
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1))  gmm::scale(MS.residu(), Kcoef);
      gmm::add(gmm::scaled(DF, -value_type(1)),
	       gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(M_, gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Mcoef),
		    gmm::sub_vector(MS.residu(), SUBI));
    }

    void set_dynamic_coeff(value_type a, value_type b) { Mcoef=a; Kcoef=b; }
    template <class VEC> void set_DF(const VEC &DF_) { gmm::copy(DF_, DF); }

    const T_MATRIX &mass_matrix(void) const { return M_; }

    void init(void) {
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->update_from_context();
      mf_u = this->mesh_fems[num_fem];
      gmm::resize(RHO_, mf_data.nb_dof());
      gmm::resize(DF, mf_u->nb_dof());
      Mcoef = Kcoef = value_type(1);
      have_subs = false;
    }

    void no_mass_on_boundary(size_type b) {
      dal::bit_vector bv = mf_u->dof_on_boundary(b);
      for (dal::bv_visitor i(bv); !i.finished(); ++i)
	ind.push_back(i);
      SUBS = gmm::unsorted_sub_index(ind);
      have_subs = true;
      if (mf_u->nb_dof() == gmm::mat_nrows(M_)) {
	gmm::sub_interval SUBI(0, mf_u->nb_dof());
	gmm::clear(gmm::sub_matrix(M_, SUBS, SUBI));
	gmm::clear(gmm::sub_matrix(M_, SUBI, SUBS));
      }
      else compute_M();
    }

    template <class VEC>
    mdbrick_dynamic(mdbrick_abstract<MODEL_STATE> &problem, mesh_fem &mf_data_,
		    const VEC &RHO__, size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), num_fem(num_fem_)
    { init(); gmm::copy(RHO__, RHO_); compute_M(); }

    mdbrick_dynamic(mdbrick_abstract<MODEL_STATE> &problem, mesh_fem &mf_data_,
		    value_type RHO__, size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), num_fem(num_fem_)
    { init(); std::fill(RHO_.begin(), RHO_.end(), RHO__); compute_M(); }

  };

  /* ******************************************************************** */
  /*		Generic solvers.                                          */
  /* ******************************************************************** */

  // standard_solve represent a default solver for the model brick system.
  // Of course it could be not adapted for a particular problem, so it could
  // be copied and adapted to change solvers, add a special traitement on the
  // problem, etc ...
  // This is in fact a model for your own solver.

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
    size_type dim = problem.dim();

    bool is_linear = problem.is_linear();
    // mtype alpha, alpha_min=mtype(1)/mtype(20);
    mtype alpha, alpha_min=mtype(1)/mtype(10000);
    mtype alpha_mult=mtype(3)/mtype(5);
    // mtype alpha_max_ratio=mtype(2)/mtype(2);
    dal::bit_vector mixvar;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.set_maxiter(10000);
    if (!is_linear) { 
      iter_linsolv0.reduce_noisy();
      iter_linsolv0.set_resmax(iter.get_resmax()/100.0);
    }


    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before
    problem.compute_residu(MS);
    problem.compute_tangent_matrix(MS);

    // cout << "CM = " << MS.constraints_matrix() << endl;

    MS.compute_reduced_system();

    // cout << "RTM = " << MS.reduced_tangent_matrix() << endl;
    
    
    mtype act_res = MS.reduced_residu_norm(), act_res_new(0);

    while (is_linear || !iter.finished(act_res)) {
    
      size_type nreddof = gmm::vect_size(MS.reduced_residu());
      gmm::iteration iter_linsolv = iter_linsolv0;
      VECTOR d(ndof), dr(nreddof);

      if (!(iter.first())) {
	problem.compute_tangent_matrix(MS);
	MS.compute_reduced_system();
      }

//       if (iter.get_noisy())
// 	cout << "tangent matrix is "
// 	   << (gmm::is_symmetric(MS.tangent_matrix(),
//             1E-6 * gmm::mat_maxnorm(MS.tangent_matrix())) ? "" : "not ")
// 	   <<  "symmetric. ";

      // cout << "MM = " << MS.reduced_tangent_matrix() << endl;

//       gmm::dense_matrix<value_type> MM(nreddof,nreddof), Q(nreddof,nreddof);
//       std::vector<value_type> eigval(nreddof);
//       gmm::copy(MS.reduced_tangent_matrix(), MM);
//       // gmm::symmetric_qr_algorithm(MM, eigval, Q);
//       gmm::implicit_qr_algorithm(MM, eigval, Q);
//       std::sort(eigval.begin(), eigval.end());
//       cout << "eival = " << eigval << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-1) << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-2) << endl;

//       double emax, emin;
//       cout << "condition number" << condition_number(MM,emax,emin) << endl;
//       cout << "emin = " << emin << endl;
//       cout << "emax = " << emax << endl;

      if (iter.get_noisy()) {
	problem.mixed_variables(mixvar);
	cout << "there is " << mixvar.card() << " mixed variables\n";
      }

      // if (0) {
      if ((ndof < 200000 && dim <= 2) || (ndof < 10000 && dim <= 3)
	  || (ndof < 1000)) {
	
	// cout << "M = " << MS.reduced_tangent_matrix() << endl;
	// cout << "L = " << MS.reduced_residu() << endl;
	double rcond;
	SuperLU_solve(MS.reduced_tangent_matrix(), dr,
		      gmm::scaled(MS.reduced_residu(), value_type(-1)),
		      rcond);
	if (iter.get_noisy()) cout << "condition number: " << 1.0/rcond<< endl;
      }
      else {
	if (problem.is_coercive()) {
	  gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	  gmm::cg(MS.reduced_tangent_matrix(), dr, 
		  gmm::scaled(MS.reduced_residu(), value_type(-1)),
		  P, iter_linsolv);
	  if (!iter_linsolv.converged()) DAL_WARNING(2,"cg did not converge!");
	} else {
	  if (mixvar.card() == 0) {
	    gmm::ilu_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	    
	    gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		       gmm::scaled(MS.reduced_residu(),  value_type(-1)), P,
		       300, iter_linsolv);
	  }
	  else {
	    gmm::ilut_precond<T_MATRIX> P(MS.reduced_tangent_matrix(),
					   20, 1E-10);
	    // gmm::identity_matrix P;
	    gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		       gmm::scaled(MS.reduced_residu(),  value_type(-1)),
		       P, 300, iter_linsolv);
	  }
	  if (!iter_linsolv.converged())
	    DAL_WARNING(2,"gmres did not converge!");
	}
      }
      MS.unreduced_solution(dr,d);

      if (is_linear) {
	gmm::add(d, MS.state());
	return;
      }
      else { // line search for the non-linear case.
	VECTOR stateinit(ndof);
	gmm::copy(MS.state(), stateinit);
       
	if ((iter.get_iteration() % 10) || (iter.get_iteration() == 0))
	  alpha = mtype(1); else alpha = mtype(1)/mtype(2);

	mtype previous_res = act_res;
	for (size_type k = 0; alpha >= alpha_min; alpha *= alpha_mult, ++k) {
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  problem.compute_residu(MS);
	  MS.compute_reduced_residu();
	  // Call to Dirichlet nullspace should be avoided -> we just need Ud
	  act_res_new = MS.reduced_residu_norm();
	  // cout << " : " << act_res_new;
	  if (act_res_new <= act_res / mtype(2)) break;
	  if (k > 0 && act_res_new > previous_res) {
	    alpha /= alpha_mult;
	    gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	    act_res_new = previous_res; break;
	  }
	  
	  previous_res = act_res_new;
	}

	// Something should be done to detect oscillating behaviors ...
	// alpha_max_ratio += (mtype(1)-alpha_max_ratio) / mtype(30);
	// alpha_min *= mtype(1) - mtype(1)/mtype(30);
      }
      act_res = act_res_new; ++iter;
      

      if (iter.get_noisy()) cout << "alpha = " << alpha << " ";
    }

  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
