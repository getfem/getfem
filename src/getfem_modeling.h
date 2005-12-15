// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_modeling.h : Defines model bricks to build
//           complete models.
// Date    : June 15, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
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


/**
   @file getfem_modeling.h
   @brief Model Bricks
   @see model_state
   @see mdbrick_abstract
*/

//==============================================
//
// Brick idents :
// MDBRICK_GENERIC_ELLIPTIC      174397
// MDBRICK_LIN_ISO_ELASTICITY    852327
// MDBRICK_MASS_MATRIX           756543
// MDBRICK_HELMHOLTZ             354864
// MDBRICK_LINEAR_INCOMP         239898
// MDBRICK_NONLINEAR_ELASTICITY  821357
// MDBRICK_NONLINEAR_INCOMP      964552
// MDBRICK_SMALL_DEF_PLASTICITY  556433
// MDBRICK_LINEAR_PLATE          897523
// MDBRICK_MIXED_LINEAR_PLATE    213456
// MDBRICK_COULOMB_FRICTION      434245
// MDBRICK_NAVIER_STOKES         394329
//
//==============================================

#ifndef GETFEM_MODELING_H__
#define GETFEM_MODELING_H__

#include <getfem_assembling.h>
#include <gmm_solver_cg.h>
#include <gmm_solver_gmres.h>
#include <gmm_precond_ildlt.h>
#include <gmm_precond_ilu.h>
#include <gmm_precond_ilut.h>
#include <gmm_precond_ilutp.h>
#include <gmm_superlu_interface.h>
#include <gmm_MUMPS_interface.h>
#include <gmm_dense_qr.h>
#include <gmm_matrix.h>
#include <gmm_solver_Schwarz_additive.h>
#include <set>
#include <dal_backtrace.h>

namespace getfem {

  /**@defgroup bricks Model Bricks
   */

  /* ******************************************************************** */
  /*		Generic definitions.                                      */
  /* ******************************************************************** */
 
  template<typename MODEL_STATE> class mdbrick_abstract;

  /** Model State : contains all storage needed by the bricks (@see mdbrick_abstract)
      @ingroup bricks

      The model state is 
        - a tangent matrix (the linear system that is solved)
	- a constraints matrix (the Dirichlet conditions etc.), and its right hand side.
	- a state (the solution of the linear system)
	- a residu (residu_ = A*x - B)
	- the same information, reduced (i.e. after removal of the constraints).
  */
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

    const gmm::col_matrix<gmm::rsvector<value_type> > &nullspace_matrix(void)
    { return NS; }
    const T_MATRIX &tangent_matrix(void) const { return tangent_matrix_; }
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
    /** Apply the constraints to the linear system */
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
    //cout << "NS = " << NS << endl;
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

  /** Default real sparse vector type for bricks @ingroup bricks */
  typedef gmm::rsvector<scalar_type> modeling_standard_sparse_vector;
  /** Default real sparse matrix type for bricks @ingroup bricks */
  typedef gmm::col_matrix<modeling_standard_sparse_vector>
                                    modeling_standard_sparse_matrix;
  /** Default real dense vector type for bricks @ingroup bricks */
  typedef std::vector<scalar_type> modeling_standard_plain_vector;

  /** Default complex sparse vector type for bricks @ingroup bricks */
  typedef gmm::rsvector<complex_type> modeling_standard_complex_sparse_vector;
  /** Default complex sparse matrix type for bricks @ingroup bricks */
  typedef gmm::col_matrix<modeling_standard_complex_sparse_vector>
                                    modeling_standard_complex_sparse_matrix;
  /** Default complex dense vector type for bricks @ingroup bricks */
  typedef std::vector<complex_type> modeling_standard_complex_plain_vector;

  /** Default real model_state for bricks @ingroup bricks */
  typedef model_state<modeling_standard_sparse_matrix,
		      modeling_standard_sparse_matrix,
		      modeling_standard_plain_vector > standard_model_state;

  /** Default complex model_state for bricks @ingroup bricks */
  typedef model_state<modeling_standard_complex_sparse_matrix,
		      modeling_standard_complex_sparse_matrix,
		      modeling_standard_complex_plain_vector >
    standard_complex_model_state;

  enum bound_cond_type { MDBRICK_UNDEFINED, MDBRICK_DIRICHLET, MDBRICK_NEUMANN,
			 MDBRICK_SIMPLE_SUPPORT, MDBRICK_CLAMPED_SUPPORT,
			 MDBRICK_FOURIER_ROBIN, MDBRICK_NAVIERSTOKESNONREF1
  };

  template <typename MAT> struct T_MAT_TYPE {
    typedef MAT T_MATRIX;
  };

#if GETFEM_PARA_LEVEL > 1
  // encore utile ?
  template <typename MAT>
  struct T_MAT_TYPE<gmm::mpi_distributed_matrix<MAT> > {
    typedef MAT T_MATRIX;
  };
#endif

  class mdbrick_abstract_parameter;

#define TYPEDEF_MODEL_STATE_TYPES                                         \
    typedef typename MODEL_STATE::vector_type VECTOR;                     \
    typedef typename T_MAT_TYPE<typename MODEL_STATE::tangent_matrix_type>::T_MATRIX T_MATRIX; \
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;       \
    typedef typename MODEL_STATE::value_type value_type;                  \
    typedef typename gmm::number_traits<value_type>::magnitude_type R;    \
    typedef typename gmm::sub_vector_type<VECTOR *,                       \
		  gmm::sub_interval>::vector_type SUBVECTOR


  /**
     Common base class for real and complex model bricks.
     @ingroup bricks 

     @see getfem::mdbrick_abstract
  */
  class mdbrick_abstract_common_base : public context_dependencies {
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
	std::map<size_type, bound_cond_type>::const_iterator it;
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

    std::vector<mdbrick_abstract_common_base *> sub_bricks;
    
    /** all proper_* specify data which is specific to this brick:
       'proper_mesh_fems' is the list of mesh_fems used by this brick,
       while 'mesh_fems' is 'proper_mesh_fems' plus the list of
       mesh_fems of all the parent bricks.
    */
    std::vector<mesh_fem *> proper_mesh_fems;
    std::vector<mesh_im *> proper_mesh_ims;
    std::vector<mesh_fem_info_> proper_mesh_fems_info;
    std::vector<boundary_cond_info_> proper_boundary_cond_info;
    /** flags indicating how this brick affect the linearity/coercivity
       etc properties of the problem */
    bool proper_is_linear_, proper_is_symmetric_, proper_is_coercive_;
    /** number of new degrees of freedom introduced by this brick */
    size_type proper_additional_dof;
    /** number of new constraints introduced by this brick */
    size_type proper_nb_constraints;
    /** in the dofs, indicates which ones correspound to mixed variables */
    dal::bit_vector proper_mixed_variables;

    /** below is the "global" information, relating to this brick and
	all its parent bricks */
    mutable std::vector<mesh_fem *> mesh_fems;
    mutable std::vector<mesh_im *> mesh_ims;
    mutable std::vector<mesh_fem_info_> mesh_fems_info;
    mutable std::vector<size_type> mesh_fem_positions;
    mutable bool is_linear_, is_symmetric_, is_coercive_;
    mutable size_type nb_total_dof, nb_total_constraints;
    mutable dal::bit_vector total_mixed_variables;

    /** the brick owns the block starting at index @c MS_i0 in the global
       tangent matrix */
    size_type MS_i0;

    /** Brick parameters */
    typedef std::map<std::string, mdbrick_abstract_parameter *> PARAM_MAP;
    PARAM_MAP parameters;
    friend class mdbrick_abstract_parameter;
    
    void parameters_set_uptodate(void);
    bool parameters_is_any_modified(void) const;

    /** the brick-specific update procedure */
    virtual void proper_update(void) = 0;
    // The following function is just for the const cast on "this".
    inline void proper_update_(void) { proper_update(); }

    /** method inherited from getfem::context_dependencies */
    void update_from_context(void) const;

    void force_update(void)
    { if (!this->context_check()) update_from_context(); }

    void add_sub_brick(mdbrick_abstract_common_base &mdb) {
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

    void add_proper_mesh_im(mesh_im &mim) {
      proper_mesh_ims.push_back(&mim);
      add_dependency(mim);
    }

    void add_proper_boundary_info(size_type num_fem, size_type num_bound,
				  bound_cond_type bc) {
      boundary_cond_info_ bci(num_fem, num_bound, bc);
      proper_boundary_cond_info.push_back(bci);
    }

    bound_cond_type boundary_type(size_type num_fem, size_type num_bound)
    { return mesh_fems_info[num_fem].boundary_type(num_bound); }


    size_type first_index(void) { return MS_i0; }


  public :

    mesh_fem_info_ &get_mesh_fem_info(size_type i)
    { context_check(); return mesh_fems_info[i]; }
    mesh_fem &get_mesh_fem(size_type i)
    { context_check(); return *(mesh_fems[i]); }
    size_type get_mesh_fem_position(size_type i)
    { context_check(); return mesh_fem_positions[i]; }
    size_type nb_mesh_fems(void) { context_check(); return mesh_fems.size(); }

    dim_type dim(void)
    { context_check(); return mesh_fems[0]->linked_mesh().dim(); }
    /** total number of variables including the variables of the
	sub-problem(s) if any */
    size_type nb_dof(void) { context_check(); return nb_total_dof; }

    /** number of linear constraints on the system including the
	constraints defined in the sub-problem(s) if any. */
    size_type nb_constraints(void)
    { context_check(); return nb_total_constraints; };

    /** true if the brick (with its sub-bricks) define a linear problem. */
    bool is_linear(void) const { context_check(); return is_linear_; }
    /** true if the brick (with its sub-bricks) define a symmetric problem. */
    bool is_symmetric(void) const { context_check(); return is_symmetric_; }
    /** true if the brick (with its sub-bricks) define a coercive problem. */
    bool is_coercive(void) const { context_check(); return is_coercive_; }
    /** return the set of mixed variables. */
    const dal::bit_vector &mixed_variables(void) const
    { context_check(); return total_mixed_variables; };

    PARAM_MAP &get_parameters() { return parameters; }

    mdbrick_abstract_common_base(void) : proper_additional_dof(0), proper_nb_constraints(0), 
					 MS_i0(0) { 
      proper_is_linear_ = proper_is_symmetric_ = proper_is_coercive_ = true; 
      nb_total_constraints = nb_total_dof = 1000000000;
      proper_additional_dof = proper_nb_constraints = 0;
    }
    virtual ~mdbrick_abstract_common_base() {}
  };

  /**
     Abstract model brick.
     @ingroup bricks 

     Requirements for a model brick :                                  
     
     A model brick is either a fondamental brick (like linearized
     elasticity brick, platicity brick ...) or a modifier brick which
     refer to a sub brick.

     A new brick should define:
     
     - @c proper_update() , which is called each time the brick should
     update itself.  This function is expected to assign the correct
     values to @c proper_nb_dof (the nb of new dof introduced by this
     brick), @c proper_nb_constraints and @c proper_mixed_variables.

     - @c do_compute_tangent_matrix(MS, i0, j0) . This function should
     compute its own part of the tangent and constraint matrices (@c i0 and
     @c j0 are the shifts in the matrices stored in the model_state MS)
     
     - @c do_compute_residu(MS, i0, j0) . Same as above for the residu
     vector stored in MS.
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_abstract : public mdbrick_abstract_common_base {
  public :
    TYPEDEF_MODEL_STATE_TYPES; // usual set of typedefs
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) = 0;
    /** update (if necessary) the tangent matrix stored.
	@param MS the model state (which contains the tangent matrix)
	@param i0,j0 position at which the tangent matrix is to be written in MS.tangent_matrix
    */
    void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0=0) {
      this->context_check();
      size_type i1 = MS_i0 = i0, j1 = j0;
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	((mdbrick_abstract*)sub_bricks[i])->compute_tangent_matrix(MS, i1, j1);
	i1 += sub_bricks[i]->nb_dof();
	j1 += sub_bricks[i]->nb_constraints();
      }
      do_compute_tangent_matrix(MS, i0, j0);
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type j0) = 0;
    void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
			size_type j0 = 0, bool first = true) {
      this->context_check();
      size_type i1 = MS_i0 = i0, j1 = j0;
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	((mdbrick_abstract*)sub_bricks[i])->compute_residu(MS, i1, j1, false);
	i1 += sub_bricks[i]->nb_dof();
	j1 += sub_bricks[i]->nb_constraints();
      }
      do_compute_residu(MS, i0, j0);


#if GETFEM_PARA_LEVEL > 1
      if (first) {
	std::vector<value_type> resloc(gmm::vect_size(MS.residu()));
    double t_init = MPI_Wtime();

	MPI_Allreduce(&((MS.residu())[0]), &(resloc[0]),
		      gmm::vect_size(MS.residu()), gmm::mpi_type(value_type()),
		      MPI_SUM, MPI_COMM_WORLD);
   cout << "reduce residu  time = " << MPI_Wtime() - t_init << endl;
	gmm::copy(resloc, MS.residu());
      }
#endif
    }
  };


  class mdbrick_abstract_parameter {
  protected:
    mdbrick_abstract_common_base *brick_;
    const mesh_fem *pmf_;
    bgeot::multi_index sizes_;
    bool initialized;
    enum { MODIFIED, UPTODATE } state;

    void update_notify() { initialized = true; state = MODIFIED; }

  public:
    const mesh_fem &mf() const { 
      if (!pmf_) DAL_THROW(dal::failure_error, "no mesh fem assigned");
      return *pmf_; 
    }
    const bgeot::multi_index& fsizes() const { return sizes_; }
    size_type fsize() const {
      size_type sz=1;
      for (unsigned i=0; i < sizes_.size(); ++i) 
	sz *= sizes_[i];
      return sz;
    }
    size_type fdim() const { return sizes_.size(); }
    mdbrick_abstract_parameter(const std::string &name,
			       mdbrick_abstract_common_base *b) {
      brick_ = b; pmf_ = 0; state = MODIFIED; initialized = false;
      b->parameters[name] = this;
    }
    mdbrick_abstract_parameter(const std::string &name,
			       const mesh_fem &mf_,
			       mdbrick_abstract_common_base *b, 
			       size_type N=0, size_type M=0, 
			       size_type P=0, size_type Q=0) {
      brick_ = b; pmf_ = &mf_; b->add_dependency(*pmf_);
      reshape(N,M,P,Q);
      state = MODIFIED; initialized = false;
      b->parameters[name] = this;
    }
    void change_mf(const mesh_fem &mf_) {
      if (&mf_ != pmf_) {
	brick_->add_dependency(mf_); pmf_ = &mf_; state = MODIFIED;
	brick_->change_context();
      }
    }
    void redim(unsigned d) {
      if (sizes_.size() != d) { sizes_.resize(d); sizes_.clear(); }
    }
    virtual void reshape(size_type N=0, size_type M=0, size_type P=0, size_type Q=0) {
      sizes_.clear();
      if (N) { sizes_.push_back(N); 
	if (M) { sizes_.push_back(M);
	  if (P) { sizes_.push_back(P);
	    if (Q) { sizes_.push_back(Q); }
	  }
	}
      }
    }
    virtual void check() const = 0;
    virtual ~mdbrick_abstract_parameter() {}
    mdbrick_abstract_common_base &brick() { return *brick_; }
    bool is_modified() const { return state != UPTODATE; }
    bool is_initialized() const { return initialized; }
    void set_uptodate(void) { state = UPTODATE; }
  };

  template <typename VEC> 
  class mdbrick_parameter : public mdbrick_abstract_parameter {
    typedef typename gmm::linalg_traits<VEC>::value_type T;
    VEC value_;

    template <typename W> void set_diagonal_(const W &w, gmm::linalg_false) {
      size_type n = fdim() == 2 ? fsizes()[0] : 1;
      VEC v(n*n);
      for (unsigned i=0; i < n; ++i) v[i*n+i] = w;
      set(w);
    }
    template <typename W> void set_diagonal_(const W &w, gmm::linalg_true) {
      size_type n = fdim() == 2 ? fsizes()[0] : 1;
      int flag = 1;
      if (gmm::vect_size(w) == n) flag = 0;
      else if (gmm::vect_size(w) != mf().nb_dof()*n) 
	DAL_THROW(dal::failure_error, "inconsistent vector dimension for set_diagonal");
      realloc();
      for (unsigned i=0; i < mf().nb_dof(); ++i) {
	for (unsigned j=0; j < n; ++j) {
	  value_[i*n*n*flag + j*n+j] = w[i*n*flag + j];
	}
      }
      update_notify();
    }
    void set_(const mesh_fem &mf_, const T& v, gmm::linalg_false) { 
      change_mf(mf_); realloc(); std::fill(value_.begin(), value_.end(), v);
      update_notify();
    }
    template<typename VEC2> void set_(const mesh_fem &mf_, const VEC2& v, gmm::linalg_true) {
      change_mf(mf_); realloc();
      size_type n = fsize();
      if (gmm::vect_size(v) == n*mf().nb_dof())
	gmm::copy(v, value_);
      else if (gmm::vect_size(v) == n) 
	for (unsigned i=0; i < mf().nb_dof(); ++i)
	  gmm::copy(v, gmm::sub_vector(value_, gmm::sub_interval(i*n, n)));
      else DAL_THROW(dal::failure_error, "inconsistent param value, expected a "
		     << fsizes() << " x " << mf().nb_dof() 
		     << " field, got a vector with " 
		     << gmm::vect_size(v) << " elements");
      update_notify();
    }


  public:
    mdbrick_parameter(const std::string &name, mdbrick_abstract_common_base *b) :
      mdbrick_abstract_parameter(name, b) { }
    mdbrick_parameter(const std::string &name, const mesh_fem &mf_, mdbrick_abstract_common_base *b,
		      size_type N=0, size_type M=0) :
      mdbrick_abstract_parameter(name, mf_,b,N,M) { }
    mdbrick_parameter(const std::string &name, const getfem_mesh &mesh, mdbrick_abstract_common_base *b,
		      size_type N=0, size_type M=0) :
      mdbrick_abstract_parameter(name, classical_mesh_fem(mesh, 0),b,N,M) { }
    void realloc() { gmm::resize(value_, fsize()*mf().nb_dof()); }
    template <typename W> void set(const mesh_fem &mf_, const W &w) {
      this->set_(mf_, w, typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set(const getfem_mesh &mesh, const W &w) {
      this->set_(classical_mesh_fem(mesh, 0), w, typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set(const W &w) {
      this->set_(mf(), w, typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set_diagonal(const W &w) {
      if ((fdim() != 0 && fdim() != 2) || (fdim() == 2 && (fsizes()[0] != fsizes()[1]))) 
	DAL_THROW(dal::failure_error, "wrong field dimension for set_diagonal");
      this->set_diagonal_(w, typename gmm::is_gmm_interfaced<W>::result());
    }
    const VEC &get() const { check(); return value_; }
    virtual void check() const {
      if (gmm::vect_size(value_) != mf().nb_dof() * fsize())
	DAL_THROW(dal::failure_error, "invalid dimension for brick parameter");
    }
  };

  /* ******************************************************************** */
  /*	       Abstract brick for linear PDE problems.                    */
  /* ******************************************************************** */

  /**  Abstract brick for linear PDE problems (such as Helmholtz, Laplacian, etc).

      Bricks which inherit from this one should not redefine
      proper_update, but they should redefine the proper_update_K ,
      whose role is only to update the content of the matrix K.
  
      @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_abstract_linear_pde
    : public mdbrick_abstract<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

  protected:
    mesh_im &mim;   /** the integration method used for the assembly */
    mesh_fem &mf_u; /** the mesh_fem for the PDE unknown */

    T_MATRIX K; /* stores the stiffness matrix. */
    bool K_uptodate;

    /** As proper_update will be called whenever a change occurs
	(such as a modification of the fem or mesh of mf_u, mim
	etc..), we don't want to redo the assembly of the K matrix
	each time, hence a flag is just set in proper_update, and
	the matrix K is computed only when it is needed.
    */
    virtual void proper_update(void) {
      K_uptodate = false; 
    }

    /** Virtual method whose purpose is to update the content of the stiffness matrix K */
    virtual void proper_update_K(void) = 0;

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }

    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }

    /** provide access to the local stiffness matrix K. */
    const T_MATRIX &get_K(void) {
      this->context_check();
      if (!K_uptodate || this->parameters_is_any_modified()) {
	gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
	gmm::clear(K);
	proper_update_K();
	K_uptodate = true;
	this->parameters_set_uptodate();
      }
      return K; 
    }

    /** provide access to the value of the solution corresponding to
	the local mesh_fem mf_u. */
    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU = gmm::sub_interval(this->first_index(),
						 mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    /** constructor 
	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
     */
    mdbrick_abstract_linear_pde(mesh_im &mim_, mesh_fem &mf_u_, size_type brick_id)
      : mim(mim_), mf_u(mf_u_) {
      this->add_proper_mesh_fem(mf_u, brick_id);
      this->add_proper_mesh_im(mim);
      this->force_update();
    }
  };

  /* ******************************************************************** */
  /*		Linearized elasticity bricks.                             */
  /* ******************************************************************** */

# define MDBRICK_LIN_ISO_ELASTICITY 852327

  /** Linear elasticity brick ( @f$ K = \int \sigma(u):\varepsilon(v) @f$ ).

      @see asm_stiffness_matrix_for_linear_elasticity
      @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_isotropic_linearized_elasticity
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;
    
    mdbrick_parameter<VECTOR> lambda_, mu_; /** the Lame coefficients */

    void proper_update_K(void) {
      if (&lambda_.mf() != &mu_.mf()) 
	DAL_THROW(failure_error, "lambda and mu should share the same mesh_fem");
      DAL_TRACE2("Assembling stiffness matrix for linear elasticity");
      asm_stiffness_matrix_for_linear_elasticity
	(this->K, this->mim, this->mf_u, lambda_.mf(), lambda_.get(), mu_.get(),
	 this->mf_u.linked_mesh().get_mpi_region());
    }

  public :
    /** accessor to the Lame coef lamda */
    mdbrick_parameter<VECTOR> &lambda(void) { return lambda_; }
    const mdbrick_parameter<VECTOR> &lambda(void) const { return lambda_; }
    /** accessor to the Lame coef mu */
    mdbrick_parameter<VECTOR> &mu(void) { return mu_; }
    const mdbrick_parameter<VECTOR> &mu(void) const { return mu_; }

    /** constructor for a homogeneous material (constant lambda and
	mu). The lame coefficient can be later changed (and set to non
	homogeneous values).

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param lambdai the (homogeneous) value of the lame coef lambda.
	@param mui the (homogeneous) value of the lame coef mu.
    */
    mdbrick_isotropic_linearized_elasticity(mesh_im &mim_, mesh_fem &mf_u_,
					    value_type lambdai = 100.0, value_type mui = 40.0)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_LIN_ISO_ELASTICITY),
	lambda_("lambda", mf_u_.linked_mesh(), this), mu_("mu", mf_u_.linked_mesh(), this) {
      lambda_.set(lambdai);
      mu_.set(mui);
    }
  };

  /* ******************************************************************** */
  /*		Mass matrix bricks.                                       */
  /* ******************************************************************** */

# define MDBRICK_MASS_MATRIX 756543

  /**
     Mass-matrix brick ( @f$\int \rho u.v @f$ ). 

     @see asm_mass_matrix_param
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mass_matrix
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> rho_;

    void proper_update_K(void) {
      DAL_TRACE2("Assembling mass matrix for mdbrick_mass_matrix");
      asm_mass_matrix_param(this->K, this->mim, this->mf_u, rho().mf(), rho().get(), 
			    this->mf_u.linked_mesh().get_mpi_region());
    }

  public :
    mdbrick_parameter<VECTOR> &rho() { return rho_; }
    const mdbrick_parameter<VECTOR> &rho() const { return rho_; }

    /** constructor for a homogeneous mass matrix. The density can be
	later changed and set to a non-homogeneous density.

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param rhoi default value for the density rho.
    */
    mdbrick_mass_matrix(mesh_im &mim_, mesh_fem &mf_u_, value_type rhoi=1)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_MASS_MATRIX),
	rho_("rho", mf_u_.linked_mesh(), this) {
      rho_.set(rhoi);
    }
  };

  /* ******************************************************************** */
  /*		Helmholtz brick.                                          */
  /* ******************************************************************** */

# define MDBRICK_HELMHOLTZ 354864

  /** Helmholtz brick ( @f$\int k^2 u.v - \nabla u.\nabla v@f$ )
      @see asm_Helmholtz
      @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Helmholtz
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> wave_number_;

    void proper_update_K(void) {
      	VECTOR w(wave_number().get());
	for (unsigned i=0; i < gmm::vect_size(w); ++i) 
	  w[i] = gmm::sqr(w[i]);
	asm_Helmholtz(this->K, this->mim, this->mf_u, wave_number().mf(), w,
		      this->mf_u.linked_mesh().get_mpi_region());
    }

  public :
    /** accessor for the wave_number */
    mdbrick_parameter<VECTOR> &wave_number() { return wave_number_; }
    const mdbrick_parameter<VECTOR> &wave_number() const { return wave_number_; }

    /** constructor for the Helmholtz problem. The wave_number can be
	later changed and set to a non-constant value over the mesh.

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param k default value for the wave number.
    */
    mdbrick_Helmholtz(mesh_im &mim_, mesh_fem &mf_u_, value_type k=1)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_HELMHOLTZ),
	wave_number_("wave_number", mf_u_.linked_mesh(), this) {
      wave_number_.set(k);
    }
   };


  /* ******************************************************************** */
  /*		general scalar elliptic brick.                            */
  /* ******************************************************************** */

# define MDBRICK_GENERIC_ELLIPTIC 174397
  
  /** General elliptic brick ( @f$ (\alpha \nabla u).\nabla v @f$ ).
      
      @f$\alpha@f$ may be a scalar, a (symmetric define positive)
      @f$N\times N@f$ matrix field, or even a (symmetric definite
      positive) @f$N\timesN\timesN\timesN@f tensor field (where N is
      the dimension of the mesh).
      
      If @c mf_u is a vector mesh_fem (qdim > 1), then the assembly is
      done componentwise.

      @see asm_stiffness_matrix_for_laplacian
      @see asm_stiffness_matrix_for_scalar_elliptic
      @see asm_stiffness_matrix_for_laplacian_componentwise
      @see asm_stiffness_matrix_for_scalar_elliptic_componentwise
      
      @ingroup bricks 
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_generic_elliptic
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    /* coeff_ holds the scalar, NxN or NxNxNxN field */
    mdbrick_parameter<VECTOR> coeff_;

    void proper_update_K(void) {
      if (coeff_.fdim() == 0) {
	if (this->mf_u.get_qdim() > 1)
	  asm_stiffness_matrix_for_laplacian_componentwise(this->K, this->mim, this->mf_u,
							   coeff().mf(), coeff().get(),
							   this->mf_u.linked_mesh().get_mpi_region());
	else
	  asm_stiffness_matrix_for_laplacian(this->K, this->mim, this->mf_u,  coeff().mf(),  coeff().get(),
					     this->mf_u.linked_mesh().get_mpi_region());
      }
      else if (coeff_.fdim() == 2) {
	if (this->mf_u.get_qdim() > 1) 
	  asm_stiffness_matrix_for_scalar_elliptic_componentwise(this->K, this->mim, this->mf_u,
								 coeff().mf(), coeff().get(),
								 this->mf_u.linked_mesh().get_mpi_region());
	else
	  asm_stiffness_matrix_for_scalar_elliptic(this->K, this->mim, this->mf_u,
						   coeff().mf(), coeff().get(),
						   this->mf_u.linked_mesh().get_mpi_region());
      }
      else if (coeff_.fdim() == 4) {
	if (this->mf_u.get_qdim() != this->mf_u.linked_mesh().dim())
	  DAL_THROW(failure_error, "Order 4 tensor coefficient applies only to mesh_fem whose Q dim is equal to the mesh dimension");
	asm_stiffness_matrix_for_vector_elliptic(this->K, this->mim, this->mf_u,
						 coeff().mf(), coeff().get(),
						 this->mf_u.linked_mesh().get_mpi_region());
      }
      else DAL_THROW(failure_error, "Bad format for the coefficient of mdbrick_generic_elliptic");
    }

    /** ensure a consistent dimension for the coeff */
    void reshape_coeff() {
      size_type N = this->mf_u.linked_mesh().dim();
      if (coeff_.fdim() == 2) coeff_.reshape(N,N);
      else if (coeff_.fdim() == 4) coeff_.reshape(N,N,N,N);
    }
  public :

    /** accessor to the coefficient k */
    mdbrick_parameter<VECTOR> &coeff() { reshape_coeff(); return coeff_; }
    const mdbrick_parameter<VECTOR> &coeff() const { return  coeff_; }

    /** Switch between a scalar coefficient, a NxN matrix field (with
	N = mf_u.linked_mesh().dim()), and a NxNxNxN tensor field. */
    void set_coeff_dimension(unsigned d) { coeff_.redim(d); }

    /** Constructor, the default coeff is a scalar equal to one
	(i.e. it gives the Laplace operator).

        The coef can be later changed to a matrix or tensor field with
        set_coeff_dimension(2 or 4) and then coeff().set(...).

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param k the scalar value of the coefficient.
    */
    mdbrick_generic_elliptic(mesh_im &mim_, mesh_fem &mf_u_, value_type k = 1.)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_GENERIC_ELLIPTIC),
	coeff_("coeff", mf_u_.linked_mesh(), this) {
      coeff_.set(k);
    }
  };


  /* ******************************************************************** */
  /*		Source term brick.                                        */
  /* ******************************************************************** */

  /**
     Source term brick ( @f$ F = \int b.v @f$ ).
     
     Update the right hand side of the linear system.

     @see asm_source_term
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_source_term : public mdbrick_abstract<MODEL_STATE>  {

    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> B_;
    VECTOR F_, auxF;
    bool F_uptodate;
    size_type boundary, num_fem, i1, nbd;
    bool have_auxF;

    void proper_update(void) {
      mesh_fem &mf_u = this->get_mesh_fem(num_fem);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();

      gmm::resize(F_, mf_u.nb_dof());
      gmm::clear(F_);
      F_uptodate = false;
    }

  public :

    mdbrick_parameter<VECTOR> &source_term(void) {       
      B_.reshape(this->get_mesh_fem(num_fem).get_qdim()); // ensure that the B shape is always consistant with the mesh_fem
      return B_; 
    }
    const mdbrick_parameter<VECTOR> &source_term(void) const { return B_; }

    /// gives the right hand side of the linear system (does no contain the auxilary part).
    const VECTOR &get_F(void) { 
      this->context_check();
      if (!F_uptodate || this->parameters_is_any_modified()) {
	mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
	F_uptodate = true;
	DAL_TRACE2("Assembling a source term");
	asm_source_term(F_, *(this->mesh_ims[0]), mf_u, B_.mf(), B_.get(),
			mf_u.linked_mesh().get_mpi_sub_region(boundary));
	this->parameters_set_uptodate();
      }
      return F_;
    }

    template <class VECT> void set_auxF(const VECT &V) {
      have_auxF = true;
      gmm::resize(auxF, (this->mesh_fems[num_fem])->nb_dof());
      gmm::copy(V, auxF);
    }

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) { }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::add(gmm::scaled(get_F(), value_type(-1)), gmm::sub_vector(MS.residu(),
	       gmm::sub_interval(i0+i1, nbd)));
      if (have_auxF)
	gmm::add(gmm::scaled(auxF, value_type(-1)),
		 gmm::sub_vector(MS.residu(),
				 gmm::sub_interval(i0+i1, nbd)));
    }

    /** Constructor defining the rhs
	@param problem the sub-problem to which this brick applies.
	@param mf_data_ the mesh_fem on which B_ is defined.
	@param B_ the value of the source term.
	@param bound the mesh boundary number on which the source term is applied 
	(by default, it is a volumic source term as the whole mesh is taken).
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			mesh_fem &mf_data_, const VECTOR &B__ = VECTOR(),
			size_type bound = size_type(-1), size_type num_fem_=0)
      : B_("source_term", mf_data_, this), boundary(bound),
	num_fem(num_fem_), have_auxF(false) {
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->force_update();
      B_.reshape(this->get_mesh_fem(num_fem).get_qdim());
      if (gmm::vect_size(B__)) B_.set(B__);
    }
  };


  /* ******************************************************************** */
  /*		Q.U term (for Fourier-Robin conditions)                   */
  /* ******************************************************************** */

  /**
     Q.U term brick ( @f$ \int (qu).v @f$ ) with @f$ q(x) @f$ a @f$N 
     \times N@f$ matrix field (assuming @c N = @c mf_u.get_qdim();)

     This brick updates the tangent matrix (this is not the case for
     the mdbrick_mass_matrix). Integration is done on a boundary or on
     the whole mesh. Can be used for Fourier-Robin boundary
     conditions.

     @see asm_qu_term
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_QU_term
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;
   
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mdbrick_parameter<VECTOR> Q_;
    size_type boundary, num_fem, i1, nbd;
    bool K_uptodate;
    T_MATRIX K;
    

    virtual void proper_update(void) {
      mesh_fem &mf_u = this->get_mesh_fem(num_fem);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      K_uptodate = false;
    }

  public :
    /** the Q parameter. */
    mdbrick_parameter<VECTOR> &Q() {       
      size_type q = this->get_mesh_fem(num_fem).get_qdim();
      Q_.reshape(q,q); // ensure that the shape of Q is coherent with the mesh_fem
      return Q_; 
    }
    const mdbrick_parameter<VECTOR> &Q() const { return Q_; }

    /** Provide access to the assembled Qu mass matrix */
    const T_MATRIX& get_K() {
      this->context_check();
      if (!K_uptodate || this->parameters_is_any_modified()) {
	const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
	gmm::clear(K);
	gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
	asm_qu_term(K, *(this->mesh_ims[0]), mf_u, Q().mf(), Q().get(), boundary);
	K_uptodate = true;
	this->parameters_set_uptodate();
      }
      return K;
    }

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type /*j0*/ ) {
      gmm::sub_interval SUBI(i0+i1, nbd);
      gmm::add(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type /*j0*/ ) {
      gmm::sub_interval SUBI(i0+i1, nbd);
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residu(), SUBI);
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI), SUBV, SUBV);
    }

    /** Constructor. 
	@param problem the sub-problem (for example an Helmholtz brick).
	@param vQ a scalar value that is the value of the (constant) diagonal of the Q matrix field. 
	(you can change it later to a non constant (or non diagonal) matrix field)
	@param bound the boundary number on which the Qu term is computed.
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    value_type vQ = 0,
		    size_type bound = size_type(-1), size_type num_fem_=0)
      : sub_problem(problem), Q_("Q", this), boundary(bound),
	num_fem(num_fem_) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      if (boundary != size_type(-1))
	this->add_proper_boundary_info(num_fem,boundary,MDBRICK_FOURIER_ROBIN);

      this->force_update();
      Q().change_mf(classical_mesh_fem(this->mesh_fems[num_fem]->linked_mesh(),0));
      Q().set_diagonal(vQ);
    }
  };

  /* ******************************************************************** */
  /*		Mixed linear incompressible condition brick.              */
  /* ******************************************************************** */

# define MDBRICK_LINEAR_INCOMP 239898

  /**
     Mixed linear incompressible condition brick.

     Update the tangent matrix with a pressure term:
       @f[
         T \longrightarrow 
          \begin{array}{ll} T & B \\ B^t & M \end{array}
       @f]
       with @f$ B = - \int p.div u @f$ and @f$ M = \int \epsilon p.q @f$ ( @f$ \epsilon @f$ is a penalization coefficient ).

       For nearly incompressible elasticity, 
       @f[ p = -\lambda \textrm{div}~u @f]
       @f[ \sigma = 2 \mu \varepsilon(u) -p I @f]
     @see asm_stokes_B
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_linear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_p;
    T_MATRIX B, M;
    bool penalized, homogeneous, BM_uptodate;
    mdbrick_parameter<VECTOR> epsilon; // penalization coefficient if any.
    size_type num_fem, i1, nbd;

    void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
      i1 = this->mesh_fem_positions.at(num_fem);
      nbd = mf_u.nb_dof();
      BM_uptodate = false;
    }

    void update_M_and_B(void) {
      this->context_check();
      if (!BM_uptodate || this->parameters_is_any_modified()) {
	mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
	size_type nd = mf_u.nb_dof(), ndd = mf_p.nb_dof();
	gmm::clear(B); gmm::resize(B, ndd, nd);
	asm_stokes_B(B, *(this->mesh_ims.at(0)), mf_u, mf_p,
		     mf_u.linked_mesh().get_mpi_region());
	if (penalized) {
	  gmm::clear(M); gmm::resize(M, ndd, ndd);
	  asm_mass_matrix_param(M, *(this->mesh_ims[0]), mf_p,
				epsilon.mf(), epsilon.get(), 
				mf_u.linked_mesh().get_mpi_region());
	  gmm::scale(M, value_type(-1));
	}
	this->proper_mixed_variables.clear();
	this->proper_mixed_variables.add(sub_problem.nb_dof(), mf_p.nb_dof());

	BM_uptodate = true;
	this->parameters_set_uptodate();
      }
    }

  public :
    /** access to the incompressibility term */
    T_MATRIX &get_B(void) { update_M_and_B(); return B; }
    /** access to the local penalized mass matrix. */
    T_MATRIX &get_M(void) { update_M_and_B(); return M; }

    /** access to the penalization term */
    mdbrick_parameter<VECTOR> &penalization_coeff(void) { return epsilon; }
    const mdbrick_parameter<VECTOR> &penalization_coeff(void) const { return epsilon; }
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::copy(get_B(), gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::copy(gmm::transposed(get_B()),
		gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      if (penalized)
	gmm::copy(get_M(), gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
      else
	gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::mult(get_B(), gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(gmm::transposed(get_B()), gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBJ));
      if (penalized) 
	gmm::mult_add(get_M(), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residu(), SUBI));
    }

    /** extract the pressure part from the model state. */
    SUBVECTOR get_pressure(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + sub_problem.nb_dof(),
			     mf_p.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    /** Switch the penalization on or off. */
    void set_penalized(bool f) { 
      if (penalized != f) { penalized = f; BM_uptodate = false; }
    }

    /** Constructor for the incompressibility condition
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), 
	penalized(false), epsilon("epsilon", mf_p, this), num_fem(num_fem_) {
      this->add_proper_mesh_fem(mf_p, MDBRICK_LINEAR_INCOMP);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->force_update();
    }
  };



  /* ******************************************************************** */
  /*		Dirichlet condition bricks.                               */
  /* ******************************************************************** */

  /**
     Dirichlet condition brick.

     This brick updates the constraints matrix and its right hand side
     (except if Lagrange multipliers are used).

     The general form for a Dirichlet condition is @f[ \int h(x)u(x).v
     = \int r(x).v \forall v@f] where @f$ h(x) @f$ is a @f$ N \times N
     @f$ matrix (by default, the identity matrix), and @f$ r(x) @f$ is
     the right hand side for the Dirichlet condition (0 for
     homogeneous conditions).

     @see asm_dirichlet_constraints
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Dirichlet : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;

    mdbrick_parameter<VECTOR> R_, H_;
    C_MATRIX G;
    VECTOR CRHS;
    size_type boundary, nb_const, num_fem;
    bool with_multipliers;
    gmm::sub_index SUB_CT;
    size_type i1, nbd;
    bool mfdata_set;
    
    mesh_fem &mf_u() { return *(this->mesh_fems[num_fem]); }

    void compute_constraints(unsigned version = 0) {
      /* for the occasional reader: this brick is quite complex as it
	 tries hard to avoid unnecessary recomputations */

      if (H_.is_modified()) { version |= ASMDIR_BUILDH; }
      if (R_.is_modified()) { version |= ASMDIR_BUILDR; }

      if (version == 0) return;

      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u().nb_dof();
      // size_type Q = mf_u.get_qdim();
      size_type nd = mf_u().nb_dof(); //  ndd = mf_data.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(nd, nd);
      VECTOR V(nd);
      
      if (!with_multipliers) version |= ASMDIR_SIMPLIFY;
      
      if (!H_.is_initialized()) {
#if GETFEM_PARA_LEVEL > 1
	if (with_multipliers)
	  asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), rhs().mf(),
				    R_.get(), mf_u().linked_mesh().get_mpi_sub_region(boundary), version);   
	else
	  asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), rhs().mf(),
				    R_.get(), boundary, version);
#else
	DAL_TRACE2("Assembling Dirichlet constraints with no H and version " << version);
	asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), rhs().mf(),
				  R_.get(), boundary, version);    
#endif
      } else {
#if GETFEM_PARA_LEVEL > 1
	if (with_multipliers)
	  asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), H().mf(), rhs().mf(),
				    H_.get(), R_.get(), mf_u().linked_mesh().get_mpi_sub_region(boundary),
				    version);
	else
	  asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), H().mf(), rhs().mf(),
				    H_.get(), R_.get(), boundary, version);
#else
	DAL_TRACE2("Assembling Dirichlet constraints with H and version " << version);
	asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u(), H().mf(), rhs().mf(),
				  H_.get(), R_.get(), boundary, version);
#endif
      }
      
      if (version & ASMDIR_BUILDH) {
	R tol=gmm::mat_maxnorm(M)*gmm::default_tol(value_type())*R(100);
	gmm::clean(M, tol);
	std::vector<size_type> ind_ct;
	dal::bit_vector nn = mf_u().dof_on_set(boundary);
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

      this->parameters_set_uptodate();
    }

    virtual void proper_update(void) {
      if (!mfdata_set) {
	// only done once, when proper_update is called by the constructor
	// (cannot be done before since mf_u() cannot be used before)
	rhs().set(classical_mesh_fem(mf_u().linked_mesh(), 0), 0);
	H().change_mf(classical_mesh_fem(mf_u().linked_mesh(), 0));
	mfdata_set = true;
      }
      /* compute_constraints has to be done here because 'nb_const' must be known.. */
      compute_constraints(ASMDIR_BUILDR | ASMDIR_BUILDH);
      this->proper_mixed_variables.clear();
      this->proper_additional_dof = with_multipliers ? nb_const : 0;
      this->proper_nb_constraints = with_multipliers ? 0 : nb_const;
      if (with_multipliers)
	this->proper_mixed_variables.add(sub_problem.nb_dof(), nb_const);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {
      compute_constraints();
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
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      compute_constraints();
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
    /** change the @f$ r(x) @f$ right hand side.
	@param R a vector of size @c Q*mf_data.nb_dof() .
    */
    mdbrick_parameter<VECTOR> &rhs() { 
      R_.reshape(mf_u().get_qdim());
      return R_; 
    }
    /** Accessor to the @f$ h(x) @f$ matrix field. */
    mdbrick_parameter<VECTOR> &H() { 
      H_.reshape(mf_u().get_qdim(), mf_u().get_qdim());
      return H_; 
    }
    
    /** Return true if the brick is using Lagrange multipliers to enforce the Dirichlet condition */
    bool is_using_multipliers() const { return with_multipliers; }
    /** Switch between lagrange multipliers and direct elimination of Dirichlet variables */
    void use_multipliers(bool v) {
      if (v != with_multipliers) {
	with_multipliers = v; 
	this->proper_is_coercive_ = !with_multipliers;
	this->change_context();
      }
    }

    /** Constructor which does not define the rhs (i.e. which sets an homogeneous Dirichlet condition)
	@param problem the sub problem to which this brick is applied.
	@param bound the boundary number for the dirichlet condition.
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      size_type bound,
		      size_type num_fem_=0)
      : sub_problem(problem), 
	R_("R", this), H_("H", this), 
	boundary(bound), num_fem(num_fem_) {
      this->add_sub_brick(sub_problem);
      with_multipliers = false;
      this->proper_is_coercive_ = true;
      this->add_proper_boundary_info(num_fem, boundary, MDBRICK_DIRICHLET);
      mfdata_set = false;
      this->force_update();
    }
  };


  /* ******************************************************************** */
  /*	                 Constraint brick.                                */
  /* ******************************************************************** */
  /** Insert a constraint @c G*U=R into the problem.
      
      @c G is a @c nc * @c mf_u.nb_dof() constraint matrix
      (@c nc is the number of constraints to be added by this brick).

      This is a more general form of the mdbrick_Dirichlet.
      @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_constraint : public mdbrick_abstract<MODEL_STATE>  {
     
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    C_MATRIX G;
    VECTOR CRHS;
    size_type num_fem;
    bool with_multipliers;

    virtual void proper_update(void) {
       size_type nbconst = gmm::mat_nrows(G);
       this->proper_mixed_variables.clear();
       this->proper_additional_dof = with_multipliers ? nbconst : 0;
       this->proper_nb_constraints = with_multipliers ? 0 : nbconst;
       if (with_multipliers)
	 this->proper_mixed_variables.add(sub_problem.nb_dof(), nbconst);
    }

    template <class MAT, class VEC>
    void set_constraints_(const MAT &G_, const VEC &RHS) {
      gmm::resize(G, gmm::mat_nrows(G_), gmm::mat_ncols(G_));
      gmm::resize(CRHS, gmm::mat_nrows(G_));
      gmm::copy(G_, G); gmm::copy(RHS, CRHS);
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = !with_multipliers;
      this->force_update();
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {

      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

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
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

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

    void use_multipliers(bool v) {
      with_multipliers = v; 
      this->proper_is_coercive_ = !with_multipliers;
      this->change_context();
    }

    template <class MAT, class VEC>
    void set_constraints(const MAT &G_, const VEC &RHS) {
      set_constraints_(G_, RHS);
      this->force_update();
    }

    /**
       @param num_fem_ the mesh_fem number on which this brick is is applied.
     */
    template <class MAT, class VEC>
    mdbrick_constraint(mdbrick_abstract<MODEL_STATE> &problem,
		       const MAT &G_, const VEC &RHS,		       
		       size_type num_fem_=0, 
		       bool with_mult = false)
      : sub_problem(problem), num_fem(num_fem_), with_multipliers(with_mult)
    { set_constraints_(G_, RHS); init_(); }
    
  };
  
  /* ******************************************************************** */
  /*  dynamic brick : not stabilized, could change a lot in the future    */
  /* ******************************************************************** */

  /**
     dynamic brick : not stabilized, could change a lot in the future.
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem *mf_u;
    mdbrick_parameter<VECTOR> RHO_;
    VECTOR DF;
    T_MATRIX M_;
    size_type num_fem;
    value_type Mcoef, Kcoef;
    std::set<size_type> boundary_sup;
    bool M_uptodate;

    virtual void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      M_uptodate = false;
    }

    void proper_update_M(void) {
      DAL_TRACE2("Assembling mass matrix for mdbrick_dynamic");
      asm_mass_matrix_param(M_, *(this->mesh_ims[0]), *mf_u, RHO_.mf(), RHO_.get());

      if (!(boundary_sup.empty())) {
	gmm::unsorted_sub_index SUBS;
	std::vector<size_type> ind;
	std::set<size_type> ind_set;
	for (std::set<size_type>::const_iterator it = boundary_sup.begin();
	     it != boundary_sup.end(); ++it) {
	  dal::bit_vector bv = mf_u->dof_on_set(*it);
	  for (dal::bv_visitor i(bv); !i.finished(); ++i) ind_set.insert(i);
	}
	ind.assign(ind_set.begin(), ind_set.end());
	SUBS = gmm::unsorted_sub_index(ind);
	gmm::sub_interval SUBI(0, mf_u->nb_dof());
	gmm::clear(gmm::sub_matrix(M_, SUBS, SUBI));
	gmm::clear(gmm::sub_matrix(M_, SUBI, SUBS));
      }
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1)) gmm::scale(MS.tangent_matrix(), Kcoef);
      gmm::add(gmm::scaled(get_M(), Mcoef),
	       gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1))  gmm::scale(MS.residu(), Kcoef);
      gmm::add(gmm::scaled(DF, -value_type(1)),
	       gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(get_M(), gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Mcoef),
		    gmm::sub_vector(MS.residu(), SUBI));
    }

    void set_dynamic_coeff(value_type a, value_type b) { Mcoef=a; Kcoef=b; }
    template <class VEC> void set_DF(const VEC &DF_)
    { gmm::resize(DF, gmm::vect_size(DF_)); gmm::copy(DF_, DF); }

    const T_MATRIX &get_M(void) {
      this->context_check();
      if (!M_uptodate || this->parameters_is_any_modified()) {
	gmm::clear(M_);
	gmm::resize(M_, mf_u->nb_dof(), mf_u->nb_dof());
	proper_update_M();
	M_uptodate = true;
	this->parameters_set_uptodate();
      }
      return M_; 
    }

    void init(void) {
      Mcoef = Kcoef = value_type(1);
      this->add_sub_brick(sub_problem);
      this->force_update();
    }

    void no_mass_on_boundary(size_type b) {
      if (boundary_sup.find(b) == boundary_sup.end()) {
	boundary_sup.insert(b);
	this->force_update();
      }
    }

    /**
       @param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_dynamic(mdbrick_abstract<MODEL_STATE> &problem,
		    value_type RHO__, size_type num_fem_=0)
      : sub_problem(problem), RHO_("rho", this), num_fem(num_fem_) {
      init();
      RHO_.set(classical_mesh_fem(mf_u->linked_mesh(), 0), RHO__);
    }
  };

  /* ******************************************************************** */
  /*		Generic solvers.                                          */
  /* ******************************************************************** */


  template <typename T> struct sort_abs_val_
  { bool operator()(T x, T y) { return (gmm::abs(x) < gmm::abs(y)); } };

  /** A default solver for the model brick system.  
      
      Of course it could be not very well suited for a particular
      problem, so it could be copied and adapted to change solvers,
      add a special traitement on the problem, etc ...  This is in
      fact a model for your own solver.

      For small problems, a direct solver is used
      (getfem::SuperLU_solve), for larger problems, a conjugate
      gradient gmm::cg (if the problem is coercive) or a gmm::gmres is
      used (preconditionned with an incomplete factorization).

      When MPI/METIS is enabled, a partition is done via METIS, and the
      gmm::additive_schwarz is used as the parallel solver.

      @ingroup bricks
  */
  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
	gmm::iteration &iter) {

    TYPEDEF_MODEL_STATE_TYPES;

    // TODO : take iter into account for the Newton. compute a consistent 
    //        max residu.

    size_type ndof = problem.nb_dof();

    bool is_linear = problem.is_linear();
    // R alpha, alpha_min=R(1)/R(20);
    R alpha, alpha_min=R(1)/R(10000);
    R alpha_mult=R(3)/R(5);
    R alpha_max_ratio=R(2)/R(1);
    dal::bit_vector mixvar;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.set_maxiter(10000);
    if (!is_linear) { 
      iter_linsolv0.reduce_noisy();
      iter_linsolv0.set_resmax(iter.get_resmax()/100.0);
    }


#if GETFEM_PARA_LEVEL > 1
    double t_init = MPI_Wtime();
#endif
    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before
    problem.compute_residu(MS);
#if GETFEM_PARA_LEVEL > 1
   cout << "comput  residu  time = " << MPI_Wtime() - t_init << endl;
   t_init = MPI_Wtime();
#endif
    problem.compute_tangent_matrix(MS);
#if GETFEM_PARA_LEVEL > 1
    cout << "comput tangent time = " << MPI_Wtime() - t_init << endl;
    // cout << "CM = " << MS.constraints_matrix() << endl;
    // cout << "TM = " << MS.tangent_matrix() << endl;
    t_init = MPI_Wtime();
#endif
    MS.compute_reduced_system();

#if GETFEM_PARA_LEVEL > 1
    cout << "comput reduced system time = " << MPI_Wtime() - t_init << endl;
#endif

#if GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == SCHWARZ_ADD /* use a domain partition ? */
    double t_ref = MPI_Wtime();
    std::set<const getfem_mesh *> mesh_set;
    for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
      mesh_set.insert(&(problem.get_mesh_fem(i).linked_mesh()));
    }

    cout << "You have " << mesh_set.size() << " meshes\n";

    std::vector< std::vector<int> > eparts(mesh_set.size());
    size_type nset = 0;
    int nparts = 64;//(ndof / 1000)+1; // number of sub-domains.

    // Quand il y a plusieurs maillages, on découpe tous les maillages en autant de parties
    // et on regroupe les ddl de chaque partition par numéro de sous-partie.
    for (std::set<const getfem_mesh *>::iterator it = mesh_set.begin();
	 it != mesh_set.end(); ++it, ++nset) {
      int ne = int((*it)->nb_convex());
      int nn = int((*it)->nb_points()), k = 0, etype = 0, numflag = 0;
      int edgecut;
      
      bgeot::pconvex_structure cvs
	= (*it)->structure_of_convex((*it)->convex_index().first())->basic_structure();
      
      if (cvs == bgeot::simplex_structure(2)) { k = 3; etype = 1; }
      else if (cvs == bgeot::simplex_structure(3)) { k = 4; etype = 2; }
      else if (cvs == bgeot::parallelepiped_structure(2)) { k = 4; etype = 4; }
      else if (cvs == bgeot::parallelepiped_structure(3)) { k = 8; etype = 3; }
      else DAL_THROW(failure_error, "This kind of element is not taken into account");

      
      std::vector<int> elmnts(ne*k), npart(nn);
      eparts[nset].resize(ne);
      int j = 0;
      // should be adapted for high order geotrans
      for (dal::bv_visitor i((*it)->convex_index()); !i.finished(); ++i, ++j)
	for (int l = 0; l < k; ++l) elmnts[j*k+l] = (*it)->ind_points_of_convex(i)[l];

      METIS_PartMeshNodal(&ne, &nn, &(elmnts[0]), &etype, &numflag, &nparts, &edgecut,
			  &(eparts[nset][0]), &(npart[0]));
    }

    std::vector<dal::bit_vector> Bidof(nparts);
    size_type apos = 0;
    for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
      const mesh_fem &mf = problem.get_mesh_fem(i);
      nset = 0;
      for (std::set<const getfem_mesh *>::iterator it = mesh_set.begin();
	   it != mesh_set.end(); ++it, ++nset)
	if (*it == &(mf.linked_mesh())) break; 
      size_type pos = problem.get_mesh_fem_position(i);
      if (pos != apos)
	DAL_THROW(failure_error, "Multiplicators are not taken into account");
      size_type length = mf.nb_dof();
      apos += length;
      for (dal::bv_visitor j(mf.convex_index()); !j.finished(); ++j) {
	size_type k = eparts[nset][j];
	for (size_type l = 0; l < mf.nb_dof_of_element(i); ++l)
	  Bidof[k].add(mf.ind_dof_of_element(j)[l] + pos);
      }
      if (apos != ndof)
	DAL_THROW(failure_error, "Multiplicators are not taken into account");
    }

    std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bi(nparts);    
    for (size_type i = 0; i < size_type(nparts); ++i) {
      gmm::resize(Bi[i], ndof, Bidof[i].card());
      size_type k = 0;
      for (dal::bv_visitor j(Bidof[i]); !j.finished(); ++j, ++k)
	Bi[i](j, k) = value_type(1);
    }

    std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bib(nparts);
    gmm::col_matrix< gmm::rsvector<value_type> > Bitemp;
    if (problem.nb_constraints() > 0) {
      for (size_type i = 0; i < size_type(nparts); ++i) {
	gmm::resize(Bib[i], gmm::mat_ncols(MS.nullspace_matrix()),
		    gmm::mat_ncols(Bi[i]));
       	gmm::mult(gmm::transposed(MS.nullspace_matrix()), Bi[i], Bib[i]);
	
	gmm::resize(Bitemp, gmm::mat_nrows(Bib[i]), gmm::mat_ncols(Bib[i]));
	gmm::copy(Bib[i], Bitemp);
	scalar_type EPS = gmm::mat_norminf(Bitemp) * gmm::default_tol(scalar_type());
	size_type k = 0;
	for (size_type j = 0; j < gmm::mat_ncols(Bitemp); ++j, ++k) {
	  // should do Schmidt orthogonalization for most sofisticated cases 
	  if (k != j) Bitemp.swap_col(j, k);
	  if (gmm::vect_norm2(gmm::mat_col(Bitemp, k)) < EPS) --k;
	}
	gmm::resize(Bitemp, gmm::mat_nrows(Bib[i]), k);
	gmm::resize(Bib[i], gmm::mat_nrows(Bib[i]), k);
	gmm::copy(Bitemp, Bib[i]);
	// cout << "Bib[" << i << "] = " <<  Bib[i] << endl;
	// cout << "Bi[" << i << "] = " <<  Bi[i] << endl;
      }
    } else std::swap(Bi, Bib);    
    cout << "METIS time = " << MPI_Wtime() - t_ref << endl;
#endif

    // cout << "RTM = " << MS.reduced_tangent_matrix() << endl;
    
    
    R act_res = MS.reduced_residu_norm(), act_res_new(0);

    while (is_linear || !iter.finished(act_res)) {
      
      size_type nreddof = gmm::vect_size(MS.reduced_residu());
      gmm::iteration iter_linsolv = iter_linsolv0;
      VECTOR d(ndof), dr(nreddof);
      
      if (!(iter.first())) {
	problem.compute_tangent_matrix(MS);
	MS.compute_reduced_system();
#if GETFEM_PARA_LEVEL > 1
        DAL_THROW(failure_error, "oups ...");
#endif
      }
      
//       if (iter.get_noisy())
//      cout << "tangent matrix " << MS.tangent_matrix() << endl;
// 	cout << "tangent matrix is "
// 	   << (gmm::is_symmetric(MS.tangent_matrix(),
//             1E-6 * gmm::mat_maxnorm(MS.tangent_matrix())) ? "" : "not ")
// 	   <<  "symmetric. ";

//      cout << "MM = " << MS.reduced_tangent_matrix() << endl;

//       gmm::dense_matrix<value_type> MM(nreddof,nreddof), Q(nreddof,nreddof);
//       std::vector<value_type> eigval(nreddof);
//       gmm::copy(MS.reduced_tangent_matrix(), MM);
//       // gmm::symmetric_qr_algorithm(MM, eigval, Q);
//       gmm::implicit_qr_algorithm(MM, eigval, Q);
//       std::sort(eigval.begin(), eigval.end(), sort_abs_val_<value_type>());
//       cout << "eival = " << eigval << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-1) << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-2) << endl;
//       double emax, emin;
//       cout << "condition number" << condition_number(MM,emax,emin) << endl;
//       cout << "emin = " << emin << endl;
//       cout << "emax = " << emax << endl;

      if (iter.get_noisy()) {
	cout << "there are " << gmm::mat_nrows(MS.constraints_matrix())
	     << " constraints\n";
	mixvar = problem.mixed_variables();
	cout << "there are " << mixvar.card() << " mixed variables\n";
      }
      size_type dim = problem.dim();

#if GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == SCHWARZ_ADD
    double t_ref,t_final;
    t_ref=MPI_Wtime();
    cout<<"begin Seq AS"<<endl;
      additive_schwarz(MS.reduced_tangent_matrix(), dr,
		       gmm::scaled(MS.reduced_residu(), value_type(-1)),
		       gmm::identity_matrix(), Bib, iter_linsolv, gmm::using_superlu(),
		       gmm::using_cg());
    t_final=MPI_Wtime();
    cout<<"temps Seq AS "<< t_final-t_ref<<endl;
#elif GETFEM_PARA_LEVEL == 1 && GETFEM_PARA_SOLVER == MUMPS
     MUMPS_solve(MS.reduced_tangent_matrix(), dr,
		 gmm::scaled(MS.reduced_residu(), value_type(-1)));
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS
     MUMPS_distributed_matrix_solve(MS.reduced_tangent_matrix(), dr,
			     gmm::scaled(MS.reduced_residu(), value_type(-1)));
#else
    // if (0) {
    if (
#ifdef GMM_USES_MUMPS
      (ndof < 200000 && dim <= 2) || (ndof < 100000 && dim <= 3)
	|| (ndof < 1000)
#else  
      (ndof < 200000 && dim <= 2) || (ndof < 13000 && dim <= 3)
      || (ndof < 1000)
#endif
      ) {
	
      // cout << "M = " << MS.reduced_tangent_matrix() << endl;
      // cout << "L = " << MS.reduced_residu() << endl;
	
      double t = ftool::uclock_sec();
#ifdef GMM_USES_MUMPS
      DAL_TRACE2("Solving with MUMPS\n");
      MUMPS_solve(MS.reduced_tangent_matrix(), dr,
		  gmm::scaled(MS.reduced_residu(), value_type(-1)));
#else
      double rcond;
      SuperLU_solve(MS.reduced_tangent_matrix(), dr,
		    gmm::scaled(MS.reduced_residu(), value_type(-1)),
		    rcond);
      if (iter.get_noisy()) cout << "condition number: " << 1.0/rcond<< endl;
#endif
      cout << "resolution time = " << ftool::uclock_sec() - t << endl;
    }
    else {
      if (problem.is_coercive()) {
	gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	gmm::cg(MS.reduced_tangent_matrix(), dr, 
		gmm::scaled(MS.reduced_residu(), value_type(-1)),
		P, iter_linsolv);
	if (!iter_linsolv.converged()) DAL_WARNING2("cg did not converge!");
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
	  DAL_WARNING2("gmres did not converge!");
      }
    } 
#endif // GETFEM_PARA_LEVEL < 2

    MS.unreduced_solution(dr,d);
    
    if (is_linear) {
      gmm::add(d, MS.state());
      return;
    }
    else { // line search for the non-linear case.
      VECTOR stateinit(ndof);
      gmm::copy(MS.state(), stateinit);
      
      if ((iter.get_iteration() % 10) || (iter.get_iteration() == 0))
	alpha = R(1); else alpha = R(1)/R(2);
      
      R previous_res = act_res;
      for (size_type k = 0; alpha >= alpha_min; alpha *= alpha_mult, ++k) {
	gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	problem.compute_residu(MS);
	MS.compute_reduced_residu();
	// Call to Dirichlet nullspace should be avoided -> we just need Ud
	act_res_new = MS.reduced_residu_norm();
	// cout << " : " << act_res_new;
	if (act_res_new <= act_res / R(2)) break;
	if (k > 0 && act_res_new > previous_res
	    && previous_res < alpha_max_ratio * act_res) {
	  alpha /= alpha_mult;
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  act_res_new = previous_res; break;
	}
	
	previous_res = act_res_new;
      }
      
      // Something should be done to detect oscillating behaviors ...
      // alpha_max_ratio += (R(1)-alpha_max_ratio) / R(30);
      alpha_min *= R(1) - R(1)/R(30);
    }
    act_res = act_res_new; ++iter;
    
    
    if (iter.get_noisy()) cout << "alpha = " << alpha << " ";
    }
    
  }
  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
