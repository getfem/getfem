// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**
   @file getfem_modeling.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 15, 2004.
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
// MDBRICK_BILAPLACIAN           783465
//
//==============================================

#ifndef GETFEM_MODELING_H__
#define GETFEM_MODELING_H__

#include "getfem_assembling.h"
#include "getfem_derivatives.h"
#include "gmm/gmm_solver_cg.h"
#include "gmm/gmm_solver_gmres.h"
#include "gmm/gmm_precond_ildlt.h"
#include "gmm/gmm_precond_ilu.h"
#include "gmm/gmm_precond_ilut.h"
#include "gmm/gmm_precond_ilutp.h"
#include "gmm/gmm_superlu_interface.h"
#include "gmm/gmm_dense_qr.h"
#include "gmm/gmm_matrix.h"
#include "gmm/gmm_solver_Schwarz_additive.h"
#include <set>
#include "dal_backtrace.h"

namespace getfem {

  /**@defgroup bricks Model Bricks
   */

  /* ******************************************************************** */
  /*		Generic definitions.                                      */
  /* ******************************************************************** */
 
  template<typename MODEL_STATE> class mdbrick_abstract;

  /** Model State : contains all storage needed by the bricks
   *  (@see mdbrick_abstract)
   *  @ingroup bricks
   *
   *  The model state is 
   *    - a tangent matrix (the linear system that is solved)
   *    - a constraints matrix (the Dirichlet conditions etc.), and its
   *	  right hand side.
   *	- a state (the solution of the linear system)
   *	- a residual (residual_ = A*x - B)
   *	- the same information, reduced (i.e. after removal of constraints).
   */
  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  class model_state {
  public :    
    typedef T_MATRIX tangent_matrix_type;
    typedef C_MATRIX constraints_matrix_type;
    typedef VECTOR vector_type;
    typedef typename gmm::linalg_traits<VECTOR>::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

  protected :
    T_MATRIX tangent_matrix_;
    C_MATRIX constraints_matrix_;
    VECTOR state_, residual_, constraints_rhs_;
    gmm::uint64_type ident_;
    
    T_MATRIX SM;
    gmm::col_matrix<gmm::rsvector<value_type> > NS; /* constraints nullspace */
    VECTOR reduced_residual_, Ud;

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
    gmm::col_matrix<gmm::rsvector<value_type> > &constraints_nullspace() 
    { return NS; }
    const VECTOR &state(void) const  { return state_; }
    VECTOR &state(void)  { return state_; }
    const VECTOR &residual(void) const  { return residual_; }
    const R reduced_residual_norm() const {
      if (gmm::mat_nrows(constraints_matrix())) {
	return sqrt(gmm::vect_norm2_sqr(reduced_residual_) + 
		    gmm::vect_norm2_sqr(Ud));
      } else return gmm::vect_norm2(residual_);
    }
    const VECTOR &reduced_residual() const { 
      return gmm::mat_nrows(constraints_matrix()) == 0 ?
	residual_ : reduced_residual_;
    }
    const T_MATRIX &reduced_tangent_matrix() const {
      return gmm::mat_nrows(constraints_matrix()) == 0 ?
	tangent_matrix_ : SM;
    }
    template <typename VECTOR1, typename VECTOR2> 
    void unreduced_solution(const VECTOR1 &U_reduced, VECTOR2 &U) {
      if (gmm::mat_nrows(constraints_matrix()))
	gmm::mult(NS, U_reduced, Ud, U);
      else gmm::copy(U_reduced, U);
    }
    /** Apply the constraints to the linear system */
    void compute_reduced_system();
    void compute_reduced_residual();
    VECTOR &residual(void) { return residual_; }
    gmm::uint64_type ident(void) { return ident_; }
    void touch(void) { ident_ = act_counter(); }
    void adapt_sizes(mdbrick_abstract<model_state> &problem);
    model_state(void) { touch(); }
    model_state(mdbrick_abstract<model_state> &problem)
    { adapt_sizes(problem); }
  };

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::compute_reduced_system() {

    if (gmm::mat_nrows(constraints_matrix()) == 0) return;
    GMM_TRACE2("Computing reduced system with respect "
	       "to global constraints");
    size_type ndof = gmm::mat_ncols(tangent_matrix());
    gmm::resize(NS, ndof, ndof);
    gmm::resize(Ud, ndof);
    
    size_type nbcols =
      Dirichlet_nullspace(constraints_matrix(), NS,
			  gmm::scaled(constraints_rhs(), value_type(-1)), Ud);
    gmm::resize(NS, ndof, nbcols);
    gmm::resize(SM, nbcols, nbcols);
    VECTOR RHaux(ndof);
    gmm::mult(tangent_matrix(), Ud, residual(), RHaux);
    gmm::resize(reduced_residual_, nbcols);
    gmm::mult(gmm::transposed(NS), RHaux, reduced_residual_);
    T_MATRIX SMaux(nbcols, ndof);
    gmm::col_matrix< gmm::rsvector<value_type> >
      NST(gmm::mat_ncols(NS), gmm::mat_nrows(NS));
    gmm::copy(gmm::transposed(NS), NST);
    gmm::mult(NST, tangent_matrix(), SMaux);
    gmm::mult(SMaux, NS, SM);
  }

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::compute_reduced_residual() {
    // The call to Dirichlet nullspace should be avoided -> we just need Ud
    if (gmm::mat_nrows(constraints_matrix()) == 0) return;
    size_type ndof = gmm::mat_ncols(tangent_matrix());
    gmm::resize(NS, ndof, ndof);
    gmm::resize(Ud, ndof);
    size_type nbcols =
      Dirichlet_nullspace(constraints_matrix(), NS,
			  gmm::scaled(constraints_rhs(), value_type(-1)), Ud);
    gmm::resize(NS, ndof, nbcols);
    gmm::resize(reduced_residual_, nbcols);
    VECTOR RHaux(ndof);
    gmm::mult(tangent_matrix(), Ud, residual(), RHaux);
    gmm::mult(gmm::transposed(NS), RHaux, reduced_residual_);
  }

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  void model_state<T_MATRIX, C_MATRIX, VECTOR>::
  adapt_sizes(mdbrick_abstract<model_state> &problem) {
    size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();

    problem.context_check();
    ndof = problem.nb_dof(); nc = problem.nb_constraints();

    if (gmm::mat_nrows(tangent_matrix_) != ndof
	|| gmm::mat_nrows(constraints_matrix_) != nc) {
      gmm::clear(state_);
      gmm::clear(residual_);
      gmm::clear(tangent_matrix_);
      gmm::clear(constraints_matrix_);
      gmm::clear(constraints_rhs_);
      gmm::resize(tangent_matrix_, ndof, ndof);
      gmm::resize(constraints_matrix_, nc, ndof);
      gmm::resize(constraints_rhs_, nc);
      gmm::resize(state_, ndof);
      gmm::resize(residual_, ndof);
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

  enum bound_cond_type { MDBRICK_UNDEFINED,
			 MDBRICK_DIRICHLET,
			 MDBRICK_NORMAL_DERIVATIVE_DIRICHLET,
			 MDBRICK_NEUMANN,
			 MDBRICK_NORMAL_DERIVATIVE_NEUMANN,
			 MDBRICK_SIMPLE_SUPPORT,
			 MDBRICK_CLAMPED_SUPPORT,
			 MDBRICK_FOURIER_ROBIN,
			 MDBRICK_NAVIERSTOKESNONREF1 };

  class mdbrick_abstract_parameter;

#define TYPEDEF_MODEL_STATE_TYPES					\
  typedef typename MODEL_STATE::vector_type VECTOR;                     \
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;         \
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;	\
    typedef typename MODEL_STATE::value_type value_type;		\
    typedef typename gmm::number_traits<value_type>::magnitude_type R;	\
    typedef typename gmm::sub_vector_type<VECTOR *,			\
			     gmm::sub_interval>::vector_type SUBVECTOR


  /**
     Common base class for real and complex model bricks.
     @ingroup bricks 

     @see getfem::mdbrick_abstract
  */
  class mdbrick_abstract_common_base : public context_dependencies, public boost::noncopyable {
  public :
    struct mesh_fem_info_ {
      size_type brick_ident; // basic model brick using the mesh_fem
      size_type info;        // flags
      // type of boundary conditions
      std::map<size_type, bound_cond_type> boundaries;
      mesh_fem_info_(size_type id, size_type in) : brick_ident(id), info(in) {}
      void add_boundary(size_type b, bound_cond_type bc) { boundaries[b]=bc; }
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
    std::vector<const mesh_fem *> proper_mesh_fems;
    std::vector<const mesh_im *> proper_mesh_ims;
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
    mutable std::vector<const mesh_fem *> mesh_fems;
    mutable std::vector<const mesh_im *> mesh_ims;
    mutable std::vector<mesh_fem_info_> mesh_fems_info;
    mutable std::vector<size_type> mesh_fem_positions;
    mutable bool is_linear_, is_symmetric_, is_coercive_;
    mutable size_type nb_total_dof, nb_total_constraints;
    mutable dal::bit_vector total_mixed_variables;

    /** the brick owns the block starting at index @c MS_i0 in the global
	tangent matrix */
    size_type MS_i0;

    /** Brick parameters */
  public:
    typedef std::map<std::string, mdbrick_abstract_parameter *> PARAM_MAP;
  protected:
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

    void add_proper_mesh_fem(const mesh_fem &mf, size_type brick_ident,
			     size_type info = 0) {
      mesh_fem_info_ mfi(brick_ident, info);
      proper_mesh_fems.push_back(&mf);
      proper_mesh_fems_info.push_back(mfi);
      add_dependency(mf);
    }

    void add_proper_mesh_im(const mesh_im &mim) {
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

    mesh_fem_info_ &get_mesh_fem_info(size_type i) const
    { context_check(); return mesh_fems_info[i]; }
    const mesh_fem &get_mesh_fem(size_type i) const
    { context_check(); return *(mesh_fems[i]); }
    size_type get_mesh_fem_position(size_type i) const
    { context_check(); return mesh_fem_positions[i]; }
    size_type nb_mesh_fems(void) const { context_check(); return mesh_fems.size(); }

    dim_type dim(void) const
    { context_check(); return mesh_fems[0]->linked_mesh().dim(); }
    /** total number of variables including the variables of the
	sub-problem(s) if any */
    size_type nb_dof(void) const { context_check(); return nb_total_dof; }

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

    mdbrick_abstract_common_base(void)
      : proper_additional_dof(0), proper_nb_constraints(0), MS_i0(0) { 
      proper_is_linear_ = proper_is_symmetric_ = proper_is_coercive_ = true; 
      nb_total_constraints = nb_total_dof = 1000000000L;
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
     
     - @c do_compute_residual(MS, i0, j0) . Same as above for the residu
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
	@param i0,j0 position at which the tangent matrix is to be written
	in MS.tangent_matrix
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
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type j0) = 0;
    void compute_residual(MODEL_STATE &MS, size_type i0 = 0,
			size_type j0 = 0, bool
#if GETFEM_PARA_LEVEL > 1			
			first
#endif
			= true) {
      this->context_check();
      size_type i1 = MS_i0 = i0, j1 = j0;
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	((mdbrick_abstract*)sub_bricks[i])->compute_residual(MS, i1, j1, false);
	i1 += sub_bricks[i]->nb_dof();
	j1 += sub_bricks[i]->nb_constraints();
      }
      do_compute_residual(MS, i0, j0);


#if GETFEM_PARA_LEVEL > 1
      if (first) {
	std::vector<value_type> resloc(gmm::vect_size(MS.residual()));

	// MPI_Barrier(MPI_COMM_WORLD);
	double t_init = MPI_Wtime();

	MPI_Allreduce(&((MS.residual())[0]), &(resloc[0]),
		      gmm::vect_size(MS.residual()), gmm::mpi_type(value_type()),
		      MPI_SUM, MPI_COMM_WORLD);
// 	MPI_Reduce(&((MS.residual())[0]), &(resloc[0]),
// 		      gmm::vect_size(MS.residual()), gmm::mpi_type(value_type()),
// 		   MPI_SUM,0, MPI_COMM_WORLD);
	gmm::copy(resloc, MS.residual());
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
    bool isconstant;
    std::string name_;
    enum { MODIFIED, UPTODATE } state;

    void update_notify() { initialized = true; state = MODIFIED; }

  public:
    const std::string name() const { return name_; }
    void rename(const std::string &new_name) { name_ = new_name; }
    const mesh_fem &mf() const { 
      GMM_ASSERT1(pmf_, "no mesh fem assigned to the parameter " << name());
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
    mdbrick_abstract_parameter(const std::string &name__,
			       mdbrick_abstract_common_base *b) {
      brick_ = b; pmf_ = 0; state = MODIFIED; initialized = false;
      name_ = name__;
      b->parameters[name()] = this;
    }
    mdbrick_abstract_parameter(const std::string &name__,
			       const mesh_fem &mf_,
			       mdbrick_abstract_common_base *b, 
			       size_type N=0, size_type M=0, 
			       size_type P=0, size_type Q=0) {
      brick_ = b; pmf_ = &mf_; name_ = name__; b->add_dependency(*pmf_);
      reshape(N,M,P,Q);
      state = MODIFIED; initialized = false; isconstant = false;
      b->parameters[name()] = this;
    }
    void change_mf(const mesh_fem &mf_) {
      if (&mf_ != pmf_) {
	brick_->add_dependency(mf_); pmf_ = &mf_; state = MODIFIED;
	brick_->change_context();
      }
    }
    void redim(unsigned d) {
      if (sizes_.size() != d) { 
	sizes_.resize(d); 
	for (unsigned i=0; i < d; ++i) sizes_[i]=0; 
      }
    }
    virtual void reshape(size_type N=0, size_type M=0, size_type P=0,
			 size_type Q=0) {
      sizes_.clear();
      if (N) { sizes_.push_back(short_type(N)); 
	if (M) { sizes_.push_back(short_type(M));
	  if (P) { sizes_.push_back(short_type(P));
	    if (Q) { sizes_.push_back(short_type(Q)); }
	  }
	}
      }
    }
    virtual void check() const = 0;
    virtual ~mdbrick_abstract_parameter() {}
    mdbrick_abstract_common_base &brick() { return *brick_; }
    bool is_modified() const { return state != UPTODATE; }
    bool is_initialized() const { return initialized; }
    bool is_constant() const { return isconstant; }
    bool is_using_default_mesh_fem() const { 
      return pmf_ == &classical_mesh_fem(pmf_->linked_mesh(),0);
    }
    void set_uptodate(void) { state = UPTODATE; }
  };

  template <typename VEC> 
  class mdbrick_parameter : public mdbrick_abstract_parameter {
    typedef typename gmm::linalg_traits<VEC>::value_type T;
    mutable VEC value_;

    template <typename W> void set_diagonal_(const W &w, gmm::linalg_false) {
      size_type n = fdim() == 2 ? fsizes()[0] : 1;
      VEC v(n*n);
      for (unsigned i=0; i < n; ++i) v[i*n+i] = w;
      set(v);
    }
    template <typename W> void set_diagonal_(const W &w, gmm::linalg_true) {
      size_type n = fdim() == 2 ? fsizes()[0] : 1;
      int flag = 1;
      if (gmm::vect_size(w) == n) flag = 0;
      else GMM_ASSERT1(gmm::vect_size(w) == mf().nb_dof()*n,
		       "inconsistent vector dimension for set_diagonal");
      realloc();
      for (unsigned i=0; i < mf().nb_dof(); ++i) {
	for (unsigned j=0; j < n; ++j) {
	  value_[i*n*n*flag + j*n+j] = w[i*n*flag + j];
	}
      }
      update_notify();
    }
    void set_(const mesh_fem &mf_, const T& v, gmm::linalg_false) {
      isconstant = true;
      change_mf(mf_); realloc(); std::fill(value_.begin(), value_.end(), v);
      update_notify();
    }
    template<typename VEC2>
    void set_(const mesh_fem &mf_, const VEC2& v, gmm::linalg_true) {
      change_mf(mf_); realloc();
      size_type n = fsize();
      if (gmm::vect_size(v) == n*mf().nb_dof()) {
	gmm::copy(v, value_);
	isconstant = false;
      }
      else if (gmm::vect_size(v) == n) {
	for (size_type i=0; i < mf().nb_dof(); ++i)
	  gmm::copy(v, gmm::sub_vector(value_, gmm::sub_interval(i*n, n)));
	isconstant = true;
      }
      else GMM_ASSERT1(false, "inconsistent param value for '" 
		       << name() << "', expected a "
		       << fsizes() << "x" << mf().nb_dof() 
		       << " field, got a vector with " 
		       << gmm::vect_size(v) << " elements");
      update_notify();
    }

  public:
    mdbrick_parameter(const std::string &name__,
		      mdbrick_abstract_common_base *b) :
      mdbrick_abstract_parameter(name__, b) { }
    mdbrick_parameter(const std::string &name__, const mesh_fem &mf_,
		      mdbrick_abstract_common_base *b,
		      size_type N=0, size_type M=0) :
      mdbrick_abstract_parameter(name__, mf_,b,N,M) { }
    mdbrick_parameter(const std::string &name__, const mesh &m,
		      mdbrick_abstract_common_base *b,
		      size_type N=0, size_type M=0) :
      mdbrick_abstract_parameter(name__, classical_mesh_fem(m, 0),b,N,M) { }
    void realloc() const { gmm::resize(value_, fsize()*mf().nb_dof()); }
    template <typename W> void set(const mesh_fem &mf_, const W &w) {
      this->set_(mf_, w, typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set(const mesh &m, const W &w) {
      this->set_(classical_mesh_fem(m, 0), w,
		 typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set(const W &w) {
      this->set_(mf(), w, typename gmm::is_gmm_interfaced<W>::result());
    }
    template <typename W> void set_diagonal(const W &w) {
      GMM_ASSERT1((fdim() == 0 || fdim() == 2)
		  && (fdim() != 2 || (fsizes()[0] == fsizes()[1])),
		  "wrong field dimension for set_diagonal for param '" 
		  << name() << "'");
      this->set_diagonal_(w, typename gmm::is_gmm_interfaced<W>::result());
    }
    const VEC &get() const { check(); return value_; }
    virtual void check() const {

      bool badsize = gmm::vect_size(value_) != mf().nb_dof() * fsize();

      GMM_ASSERT1(is_initialized(), "Parameter " << name()
		  << " is not initialized");
      GMM_ASSERT1(!badsize || (is_constant() && gmm::vect_size(value_) != 0),
		  "invalid dimension for brick parameter '" << name() << 
		  "', expected an array of size " << 
		  mf().nb_dof()*fsize() << "=" << fsize() << "x" << 
		  mf().nb_dof() << ", got an array of size " << 
		  gmm::vect_size(value_));
      if (badsize) {
	realloc();
	size_type n = fsize();
	std::vector<T> v(n);
	gmm::copy(gmm::sub_vector(value_, gmm::sub_interval(0, n)), v);
	for (size_type i=1; i < mf().nb_dof(); ++i)
	  gmm::copy(v, gmm::sub_vector(value_, gmm::sub_interval(i*n, n)));
      }
    }
  };

  /* ******************************************************************** */
  /*	       Abstract brick for linear PDE problems.                    */
  /* ******************************************************************** */

  /** Abstract brick for linear PDE problems (such as Helmholtz,
      Laplacian, etc).
      
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
    const mesh_im &mim;   /** the integration method used for the assembly */
    const mesh_fem &mf_u; /** the mesh_fem for the PDE unknown */

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

    /** Virtual method whose purpose is to update the content of the
	stiffness matrix K
    */
    virtual void proper_update_K(void) = 0;

  public :

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

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }

    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residual(), SUBI));
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
    mdbrick_abstract_linear_pde(const mesh_im &mim_, const mesh_fem &mf_u_,
				size_type brick_id)
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
      GMM_ASSERT1(&lambda_.mf() == &mu_.mf(),
		  "lambda and mu should share the same mesh_fem");
      GMM_TRACE2("Assembling stiffness matrix for linear elasticity");
      this->context_check();

      // mesh_fem mff(this->mf_u.linked_mesh());
      // mff.set_classical_finite_element(0);
      // mff.set_qdim(3);
      // gmm::resize(this->K, mff.nb_dof(), mff.nb_dof());

      asm_stiffness_matrix_for_linear_elasticity
	(this->K, this->mim, this->mf_u, lambda_.mf(), lambda_.get(), mu_.get(),
	 //      asm_stiffness_matrix_for_linear_elasticity
	 // (this->K, this->mim, this->mf_u, lambda_.mf(), lambda_.get(), mu_.get(),
	 this->mf_u.linked_mesh().get_mpi_region());
      
      this->context_check();
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
    mdbrick_isotropic_linearized_elasticity(const mesh_im &mim_, const mesh_fem &mf_u_,
					    value_type lambdai = 100.0, value_type mui = 40.0)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_LIN_ISO_ELASTICITY),
	lambda_("lambda", mf_u_.linked_mesh(), this), mu_("mu", mf_u_.linked_mesh(), this) {
      lambda_.set(lambdai);
      mu_.set(mui);
    }

    template <class VECTVM>
    void compute_Von_Mises_or_Tresca(MODEL_STATE &MS, 
				     const mesh_fem &mf_vm, 
				     VECTVM &VM, bool tresca) {
      getfem::interpolation_von_mises_or_tresca
	(this->mf_u,mf_vm,get_solution(MS),VM,
	 lambda().mf(),lambda().get(),mu().mf(),mu().get(),tresca);
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
      GMM_TRACE2("Assembling mass matrix for mdbrick_mass_matrix");
      gmm::clear(this->K);
      asm_mass_matrix_param
	(this->K, this->mim, this->mf_u, rho().mf(), rho().get(), 
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
    mdbrick_mass_matrix(const mesh_im &mim_, const mesh_fem &mf_u_,
			value_type rhoi=1)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_,
						 MDBRICK_MASS_MATRIX),
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
      gmm::clear(this->K);
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
    mdbrick_Helmholtz(const mesh_im &mim_, const mesh_fem &mf_u_, value_type k=1)
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
  positive) @f$N\times N\times N\times N@f$ tensor field (where N is
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
      gmm::clear(this->K);
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
	GMM_ASSERT1(this->mf_u.get_qdim() == this->mf_u.linked_mesh().dim(),
		    "Order 4 tensor coefficient applies only to mesh_fem "
		    "whose Q dim is equal to the mesh dimension");
	asm_stiffness_matrix_for_vector_elliptic(this->K, this->mim, this->mf_u,
						 coeff().mf(), coeff().get(),
						 this->mf_u.linked_mesh().get_mpi_region());
      }
      else GMM_ASSERT1(false, "Bad format for the coefficient of mdbrick_generic_elliptic");
    }

    /** ensure a consistent dimension for the coeff */
    void reshape_coeff() {
      size_type N = this->mf_u.linked_mesh().dim();
      if (coeff_.fdim() == 0)      coeff_.reshape();
      else if (coeff_.fdim() == 2) coeff_.reshape(N,N);
      else if (coeff_.fdim() == 4) coeff_.reshape(N,N,N,N);
    }
  public :

    /** accessor to the coefficient k */
    mdbrick_parameter<VECTOR> &coeff() { reshape_coeff(); return coeff_; }
    const mdbrick_parameter<VECTOR> &coeff() const { return  coeff_; }

    /** Switch between a scalar coefficient, a NxN matrix field (with
	N = mf_u.linked_mesh().dim()), and a NxNxNxN tensor field. */
    void set_coeff_dimension(unsigned d) { coeff_.redim(d); reshape_coeff(); }

    /** Constructor, the default coeff is a scalar equal to one
	(i.e. it gives the Laplace operator).

        The coef can be later changed to a matrix or tensor field with
        set_coeff_dimension(2 or 4) and then coeff().set(...).

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param k the scalar value of the coefficient.
    */
    mdbrick_generic_elliptic(const mesh_im &mim_,
			     const mesh_fem &mf_u_, value_type k = 1.)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_GENERIC_ELLIPTIC),
	coeff_("A", mf_u_.linked_mesh(), this) {
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
      const mesh_fem &mf_u = this->get_mesh_fem(num_fem);
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
	const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
	F_uptodate = true;
	GMM_TRACE2("Assembling a source term");
	gmm::clear(F_);
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
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::add(gmm::scaled(get_F(), value_type(-1)),
	       gmm::sub_vector(MS.residual(), gmm::sub_interval(i0+i1, nbd)));
      if (have_auxF)
	gmm::add(gmm::scaled(auxF, value_type(-1)),
		 gmm::sub_vector(MS.residual(),
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
			const mesh_fem &mf_data_, const VECTOR &B__ = VECTOR(),
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

    /** Constructor not defining the rhs
	@param problem the sub-problem to which this brick applies.
	@param bound the mesh boundary number on which the source term
	       is applied 
	(by default, it is a volumic source term as the whole mesh is taken).
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			size_type bound = size_type(-1), size_type num_fem_=0)
      : B_("source_term", this), boundary(bound),
	num_fem(num_fem_), have_auxF(false) {
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->force_update();
      B_.reshape(this->get_mesh_fem(num_fem).get_qdim());
    }
  };

  /* ******************************************************************** */
  /*		Normal Source term brick.                                 */
  /* ******************************************************************** */

  /**
     Normal source term brick ( @f$ F = \int (b n).v @f$ ).
     
     Update the right hand side of the linear system.

     @see asm_source_term
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_normal_source_term : public mdbrick_abstract<MODEL_STATE>  {

    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> B_;
    VECTOR F_;
    bool F_uptodate;
    size_type boundary, num_fem, i1, nbd;

    const mesh_fem &mf_u(void) const { return this->get_mesh_fem(num_fem); }

    void proper_update(void) {
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u().nb_dof();
      gmm::resize(F_, nbd); gmm::clear(F_);
      F_uptodate = false;
    }


  public :

    mdbrick_parameter<VECTOR> &normal_source_term(void) {
      // ensure that the B shape is always consistant with the mesh_fem
      B_.reshape(mf_u().get_qdim(), mf_u().linked_mesh().dim());
      return B_; 
    }
    const mdbrick_parameter<VECTOR> &normal_source_term(void) const
    { return B_; }

    /// gives the right hand side of the linear system.
    const VECTOR &get_F(void) { 
      this->context_check();
      if (!F_uptodate || this->parameters_is_any_modified()) {
	F_uptodate = true;
	GMM_TRACE2("Assembling a source term");
	gmm::clear(F_);
	asm_normal_source_term
	  (F_, *(this->mesh_ims[0]), mf_u(), B_.mf(), B_.get(),
	   mf_u().linked_mesh().get_mpi_sub_region(boundary));
	this->parameters_set_uptodate();
      }
      return F_;
    }

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) { }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::add(gmm::scaled(get_F(), value_type(-1)),
	       gmm::sub_vector(MS.residual(), gmm::sub_interval(i0+i1, nbd)));
    }

    /** Constructor defining the rhs
	@param problem the sub-problem to which this brick applies.
	@param mf_data_ the mesh_fem on which B_ is defined.
	@param B_ the value of the source term ( a Qdim x mesh_dim field )
	@param bound the mesh boundary number on which the source term is applied.
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_normal_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			       const mesh_fem &mf_data_, const VECTOR &B__,
			       size_type bound, size_type num_fem_=0)
      : B_("normal_source_term", mf_data_, this), boundary(bound),
	num_fem(num_fem_) {
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->force_update();
      B_.reshape(mf_u().get_qdim(),mf_u().linked_mesh().dim());
      if (gmm::vect_size(B__)) B_.set(B__);
    }

    /** Constructor not defining the rhs
	@param problem the sub-problem to which this brick applies.
	@param bound the mesh boundary number on which the source term is applied.
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_normal_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			       size_type bound, size_type num_fem_=0)
      : B_("normal_source_term", this), boundary(bound),
	num_fem(num_fem_) {
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      this->force_update();
      B_.reshape(mf_u().get_qdim(),mf_u().linked_mesh().dim());
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
      const mesh_fem &mf_u = this->get_mesh_fem(num_fem);
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
	asm_qu_term(K, *(this->mesh_ims[0]), mf_u, Q().mf(), Q().get(),
		    mf_u.linked_mesh().get_mpi_sub_region(boundary));
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
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type /*j0*/ ) {
      gmm::sub_interval SUBI(i0+i1, nbd);
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residual(), SUBI);
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
     with @f$ B = - \int p.div u @f$ and
     @f$ M = \int \epsilon p.q @f$ ( @f$ \epsilon @f$ is a penalization
     coefficient ).

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
    const mesh_fem &mf_p;
    T_MATRIX B, M;
    bool penalized, homogeneous, BM_uptodate;
    mdbrick_parameter<VECTOR> epsilon; // penalization coefficient if any.
    size_type num_fem, i1, nbd;

    void proper_update(void) {
      const mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
      i1 = this->mesh_fem_positions.at(num_fem);
      nbd = mf_u.nb_dof();
      BM_uptodate = false;
    }

    void update_M_and_B(void) {
      this->context_check();
      if (!BM_uptodate || this->parameters_is_any_modified()) {
	const mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
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
    const mdbrick_parameter<VECTOR> &penalization_coeff(void) const
    { return epsilon; }
    
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
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::mult(get_B(), gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residual(), SUBI));
      gmm::mult_add(gmm::transposed(get_B()),
		    gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residual(), SUBJ));
      if (penalized) 
	gmm::mult_add(get_M(), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residual(), SUBI));
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
			  const mesh_fem &mf_p_, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), 
	penalized(false), epsilon("epsilon", mf_p, this), num_fem(num_fem_) {
      this->add_proper_mesh_fem(mf_p, MDBRICK_LINEAR_INCOMP);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->force_update();
    }
  };


  /* ******************************************************************** */
  /*	                 Constraint brick.                                */
  /* ******************************************************************** */
  /** Insert a constraint @c B*U=R into the problem.
   *   
   *  @c B is a @c nc * @c mf_u.nb_dof() constraint matrix
   *  (@c nc is the number of constraints to be added by this brick).
   *
   *  This brick is in particular a base class of mdbrick_Dirichlet.
   *  It can also be used on its own to add a constraint to make
   *  the linear system well-posed for instance when one of the unknown
   *  is defined modulo a constant (typically a pressure term).
   *  @ingroup bricks
   */

  typedef enum { AUGMENTED_CONSTRAINTS = 0,
		 PENALIZED_CONSTRAINTS = 1,
		 ELIMINATED_CONSTRAINTS = 2 } constraints_type;
  
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_constraint : public mdbrick_abstract<MODEL_STATE>  {
  public :

    TYPEDEF_MODEL_STATE_TYPES;
    typedef gmm::row_matrix<gmm::rsvector<value_type> > local_C_MATRIX;

  protected :

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    local_C_MATRIX  B; // Constraint matrix. Supposed to be of maximal rank.
    T_MATRIX optK, optM; // Optional additional matrices (penalisation or
                         // aumentation terms
    VECTOR CRHS;     // right hand side of the constraints (BU = CRHS)
    R eps;           // Parameter for the PENALIZED_CONSTRAINTS option
    size_type num_fem;
    constraints_type co_how;

    virtual void recompute_B(void) {} // for derived classes
    virtual void recompute_B_sizes(void) {}  // for derived classes
    
    virtual void proper_update(void) {
      recompute_B_sizes();
      size_type nbconst = gmm::mat_nrows(B);
      this->proper_mixed_variables.clear();
      this->proper_additional_dof = 0;
      this->proper_nb_constraints = 0;
      switch (co_how) {
      case AUGMENTED_CONSTRAINTS :
	this->proper_additional_dof = nbconst;
	this->proper_mixed_variables.add(sub_problem.nb_dof(), nbconst);
	break;
      case ELIMINATED_CONSTRAINTS :
	this->proper_nb_constraints = nbconst;
	break;
      default : break;
      }
    }

    template <class MAT, class VEC>
    void set_constraints_(const MAT &B_, const VEC &RHS) {
      gmm::resize(B, gmm::mat_nrows(B_), gmm::mat_ncols(B_));
      gmm::resize(CRHS, gmm::mat_nrows(B_));
      gmm::copy(B_, B); gmm::copy(RHS, CRHS);
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = (co_how != AUGMENTED_CONSTRAINTS);
      this->force_update();
    }

  public :

    // For Tomas Ligursky !!
    const VECTOR &get_CRHS(void) const { return CRHS; }
    size_type first_ind(void)
    { return this->first_index()+sub_problem.nb_dof(); }

    /* provide access to the value of the multipliers only for multiplier
     * option (could be extended ?).
     */
    SUBVECTOR get_mult(MODEL_STATE &MS) {
      GMM_ASSERT1(co_how == AUGMENTED_CONSTRAINTS,
		  "Only for Augmented constraint option");
      gmm::sub_interval SUBM
	= gmm::sub_interval(this->first_index()+sub_problem.nb_dof(),
			    gmm::mat_nrows(B));
      return gmm::sub_vector(MS.state(), SUBM);
    }

    template <class MAT1, class MAT2>
    void set_optional_matrices(const MAT1 &K_, const MAT2 &M_) {
      gmm::resize(optK, gmm::mat_nrows(K_), gmm::mat_ncols(K_));
      gmm::resize(optM, gmm::mat_nrows(M_), gmm::mat_ncols(M_));
      gmm::copy(K_, optK); gmm::copy(M_, optM);
    }

    const local_C_MATRIX &get_B() { recompute_B(); return B; }

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {

      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

      switch (co_how) {
      case AUGMENTED_CONSTRAINTS :
	{
	  gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), gmm::mat_nrows(B));
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  gmm::copy(get_B(), gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	  gmm::copy(gmm::transposed(get_B()),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	  if (gmm::mat_nrows(optK) != 0)
	    gmm::add(optK, gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBJ));
	  if (gmm::mat_nrows(optM) != 0)
	    gmm::copy(optM, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
	  else
	    gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));     
	}
	break;
      case PENALIZED_CONSTRAINTS :
	{
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  local_C_MATRIX BTB(nbd, nbd); // could be stored optionally
	  gmm::mult(gmm::transposed(get_B()), get_B(), BTB); // to be optimized
	  gmm::add(gmm::scaled(BTB, value_type(1) / eps),
		   gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBJ));
	}
	break;
      case ELIMINATED_CONSTRAINTS :
	{
	  size_type ncs = sub_problem.nb_constraints();
	  gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(get_B())),
	    SUBJ(i0+i1, nbd);
	  gmm::copy(get_B(), gmm::sub_matrix(MS.constraints_matrix(),
					     SUBI, SUBJ));
	}
	break;
      }
    }

    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type i1 = this->mesh_fem_positions[num_fem];
      size_type nbd = mf_u.nb_dof();

      switch (co_how) {
      case AUGMENTED_CONSTRAINTS :
	{	
	  gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(),
				 gmm::mat_nrows(get_B()));
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  gmm::mult(get_B(), gmm::sub_vector(MS.state(), SUBJ),
		    gmm::scaled(CRHS, value_type(-1)),
		    gmm::sub_vector(MS.residual(), SUBI));
	  if (gmm::mat_nrows(optM) != 0)
	    gmm::mult_add(optM, gmm::sub_vector(MS.state(), SUBI),
			  gmm::sub_vector(MS.residual(), SUBI));
	
	  gmm::mult_add(gmm::transposed(get_B()),
			gmm::sub_vector(MS.state(), SUBI),
			gmm::sub_vector(MS.residual(), SUBJ));
	  if (gmm::mat_nrows(optK) != 0)
	    gmm::mult_add(optK, gmm::sub_vector(MS.state(), SUBJ),
			  gmm::sub_vector(MS.residual(), SUBJ));
	}
	break;
      case PENALIZED_CONSTRAINTS :
	{
	  gmm::sub_interval SUBJ(i0+i1, nbd);
	  std::vector<value_type> Raux(gmm::mat_nrows(get_B()));
	  gmm::mult(get_B(), gmm::sub_vector(MS.state(), SUBJ),
		    gmm::scaled(CRHS, value_type(-1)), Raux);
	  gmm::mult_add(gmm::transposed(get_B()),
			gmm::scaled(Raux, value_type(1) / eps),
			gmm::sub_vector(MS.residual(), SUBJ));
	}
	break;
      case ELIMINATED_CONSTRAINTS :
	{
	  size_type ncs = sub_problem.nb_constraints();
	  gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(get_B())),
	    SUBJ(i0+i1, nbd);
	  gmm::mult(get_B(), gmm::sub_vector(MS.state(), SUBJ),
		    gmm::scaled(CRHS, value_type(-1)),
		    gmm::sub_vector(MS.constraints_rhs(), SUBI));
	  gmm::copy(get_B(), gmm::sub_matrix(MS.constraints_matrix(),
					     SUBI, SUBJ));
	}
	break;
      }
    }

    /** Set the method to take into account the constraints :
     	AUGMENTED_CONSTRAINTS, PENALIZED_CONSTRAINTS or
     	ELIMINATED_CONSTRAINTS.
	
       Remark: the penalization is often a quick and safe choice,
       however you should be aware that the stop criterion of the
       iterative solvers should be lowered accordingly (i.e. for a
       penalization parameter of 1e9, the target residu should be
       multiplied by 1e-9) , or the iterative method may consider it
       has converged to the solution while it has just converged on
       the subspace of penalized constraints!
    */
    void set_constraints_type(constraints_type v) {
      if (co_how != v) {
	co_how = v; 
	this->proper_is_coercive_ = (co_how != AUGMENTED_CONSTRAINTS);
	this->change_context();
      }
    }

    /** Change the penalization parameter for the PENALIZED_CONSTRAINTS
     *	option (the default value is 1e-9)
     */
    void set_penalization_parameter(R new_eps) { eps = new_eps; }

    template <class MAT, class VEC>
    void set_constraints(const MAT &B_, const VEC &RHS) {
      bool fupdate = (gmm::mat_nrows(B_) != gmm::mat_nrows(B))
	|| (gmm::mat_ncols(B_) != gmm::mat_ncols(B));
      set_constraints_(B_, RHS);
      if (fupdate) this->force_update();
    }

    template <class VEC>
    void set_constraints_rhs(const VEC &RHS)  { gmm::copy(RHS, CRHS); }

    /** Constructor with no constraint (to add the constraints use
     *  set_contraints(B, RHS)). 
     *  @param num_fem_ the mesh_fem number on which this brick is is applied.
     */
    mdbrick_constraint(mdbrick_abstract<MODEL_STATE> &problem,		       
		       size_type num_fem_=0)
      : sub_problem(problem), eps(1e-9), num_fem(num_fem_),
	co_how(AUGMENTED_CONSTRAINTS) { 
      init_(); 
    }

    explicit mdbrick_constraint(mdbrick_constraint<MODEL_STATE>& problem) : 
      sub_problem(problem), eps(1e-9), num_fem(0),
      co_how(AUGMENTED_CONSTRAINTS) { init_(); }

  };


  /* ******************************************************************** */
  /*		Standard Dirichlet condition bricks.                      */
  /* ******************************************************************** */

  /** Standard Dirichlet condition brick.
   *
   *  This brick represent a Dirichlet condition on a part of a boundary.
   *  The general form for a Dirichlet condition is @f[ \int u(x)v(x)
   *  = \int r(x)v(x) \forall v@f] where @f$ r(x) @f$ is
   *  the right hand side for the Dirichlet condition (0 for
   *  homogeneous conditions) and @f$ v @f$ is in a space of multipliers
   *  defined by the trace of mf_mult on the considered part of boundary.
   *  (The default is to take the same finite element method as for
   *  the unknown @f$ u(x) @f$. If this fem is not too complex a more
   *  standard method should be prefered).
   *
   *  For the methods of the object see also the mdbrick_constraint which
   *  is a base class of this brick.
   *
   *  @see asm_dirichlet_constraints
   *  @see mdbrick_constraint
   *  @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Dirichlet : public mdbrick_constraint<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> R_;
    
    size_type boundary;
    bool mfdata_set, B_to_be_computed;
    gmm::sub_index SUB_CT;
    const mesh_fem *mf_mult;
    
    const mesh_fem &mf_u() { return *(this->mesh_fems[this->num_fem]); }
    const mesh_im  &mim() { return *(this->mesh_ims[0]); }

    void compute_constraints(unsigned version) {
      size_type ndu = mf_u().nb_dof(), ndm = mf_mult->nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(ndm, ndu);
      VECTOR V(ndm);
      if (this->co_how != AUGMENTED_CONSTRAINTS) version |= ASMDIR_SIMPLIFY;
      GMM_TRACE2("Assembling Dirichlet constraints, version " << version);
      asm_dirichlet_constraints
	(M, V, mim(), mf_u(), *mf_mult, rhs().mf(), R_.get(),
	 mf_u().linked_mesh().get_mpi_sub_region(boundary), version);    
      if (version & ASMDIR_BUILDH)
	gmm::copy(gmm::sub_matrix(M, SUB_CT, gmm::sub_interval(0, ndu)), 
		  this->B);
      gmm::copy(gmm::sub_vector(V, SUB_CT), this->CRHS);
    }

    virtual void recompute_B_sizes(void) {
      if (!mfdata_set) {
	rhs().set(classical_mesh_fem(mf_u().linked_mesh(), 0), 0);
 	mfdata_set = true;
      }
      size_type nd = mf_u().nb_dof();
      dal::bit_vector dof_on_bound;
      if (mf_mult->is_reduced())
	dof_on_bound.add(0, mf_mult->nb_dof());
      else
	dof_on_bound = mf_mult->basic_dof_on_region(boundary);

      size_type nb_const = dof_on_bound.card();
      std::vector<size_type> ind_ct;
      for (dal::bv_visitor i(dof_on_bound); !i.finished(); ++i)
	ind_ct.push_back(i);
      SUB_CT = gmm::sub_index(ind_ct);
      gmm::resize(this->B, nb_const, nd);
      gmm::resize(this->CRHS, nb_const);
      B_to_be_computed = true;
    }

    virtual void recompute_B(void) {
      unsigned version = 0;
      if (R_.is_modified()) { version = ASMDIR_BUILDR; }
      if (B_to_be_computed) { version = ASMDIR_BUILDR | ASMDIR_BUILDH; }
      if (version) { 
	compute_constraints(version);
	this->parameters_set_uptodate();
	B_to_be_computed = false;
      }
    }

  public :

    /** Change the @f$ r(x) @f$ right hand side.
     *	@param R a vector of size @c Q*mf_data.nb_dof() .
     */
    mdbrick_parameter<VECTOR> &rhs()
    { R_.reshape(mf_u().get_qdim()); return R_; }

    /** Constructor which does not define the rhs (i.e. which sets an
     *	homogeneous Dirichlet condition)
     *	@param problem the sub problem to which this brick is applied.
     *	@param bound the boundary number for the dirichlet condition.
     *  @param mf_mult_ the mesh_fem for the multipliers.
     *	@param num_fem_ the mesh_fem number on which this brick is is applied.
     */
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      size_type bound,
		      const mesh_fem &mf_mult_ = dummy_mesh_fem(), 
		      size_type num_fem_=0)
      : mdbrick_constraint<MODEL_STATE>(problem, num_fem_), 
	R_("R", this), boundary(bound) {
      mf_mult = (&mf_mult_ == &dummy_mesh_fem()) ? &(mf_u()) : &mf_mult_;
      this->add_proper_boundary_info(this->num_fem, boundary, 
				     MDBRICK_DIRICHLET);
      this->add_dependency(*mf_mult);
      mfdata_set = false; B_to_be_computed = true;
      this->force_update();
      GMM_ASSERT1(mf_mult->get_qdim() == mf_u().get_qdim(),
		  "The lagrange multipliers mesh fem "
		  "for the Dirichlet brick should have the same Qdim as "
		  "the main mesh_fem");
    }
  };


  /* ******************************************************************** */
  /*		normal component Dirichlet condition bricks.              */
  /* ******************************************************************** */

  /** normal component Dirichlet condition brick.
   *
   *  This brick represent a Dirichlet condition on a part of a boundary for
   *  the normal component of a vectorial unknown.
   *  The general form for this Dirichlet condition is @f[ \int (u(x).n)v(x)
   *  = \int r(x)v(x) \forall v@f] where @f$ r(x) @f$ is the scalar
   *  right hand side for the Dirichlet condition (0 for
   *  homogeneous conditions) and @f$ v @f$ is in a space of scalar multipliers
   *  defined by the trace of mf_mult on the considered part of boundary.
   *
   *  For the methods of the object see also the mdbrick_constraint which
   *  is a base class of this brick.
   *
   *  @see asm_normal_part_dirichlet_constraints
   *  @see mdbrick_constraint
   *  @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_normal_component_Dirichlet
    : public mdbrick_constraint<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> R_;
    
    size_type boundary;
    bool mfdata_set, B_to_be_computed;
    gmm::sub_index SUB_CT;
    const mesh_fem &mf_mult;
    
    const mesh_fem &mf_u() { return *(this->mesh_fems[this->num_fem]); }
    const mesh_im  &mim() { return *(this->mesh_ims[0]); }

    void compute_constraints(unsigned version) {
      size_type ndu = mf_u().nb_dof(), ndm = mf_mult.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(ndm, ndu);
      VECTOR V(ndm);
      if (this->co_how != AUGMENTED_CONSTRAINTS) version |= ASMDIR_SIMPLIFY;
      GMM_TRACE2("Assembling normal component Dirichlet constraints, version "
		 << version);
      asm_normal_component_dirichlet_constraints
	(M, V, mim(), mf_u(), mf_mult, rhs().mf(), rhs().get(),
	 mf_u().linked_mesh().get_mpi_sub_region(boundary), version);    
      if (version & ASMDIR_BUILDH)
	gmm::copy(gmm::sub_matrix(M, SUB_CT, gmm::sub_interval(0, ndu)), 
		  this->B);
      gmm::copy(gmm::sub_vector(V, SUB_CT), this->CRHS);
    }

    virtual void recompute_B_sizes(void) {
      if (!mfdata_set) {
	rhs().set(classical_mesh_fem(mf_u().linked_mesh(), 0), 0);
 	mfdata_set = true;
      }
      size_type nd = mf_u().nb_dof();
      dal::bit_vector dof_on_bound;
      if (mf_mult.is_reduced())
	dof_on_bound.add(0, nd);
      else
	dof_on_bound = mf_mult.basic_dof_on_region(boundary);
      size_type nb_const = dof_on_bound.card();
      std::vector<size_type> ind_ct;
      for (dal::bv_visitor i(dof_on_bound); !i.finished(); ++i)
	ind_ct.push_back(i);
      SUB_CT = gmm::sub_index(ind_ct);
      gmm::resize(this->B, nb_const, nd);
      gmm::resize(this->CRHS, nb_const);
      B_to_be_computed = true;
    }

    virtual void recompute_B(void) {
      unsigned version = 0;
      if (R_.is_modified()) { version = ASMDIR_BUILDR; }
      if (B_to_be_computed) { version = ASMDIR_BUILDR | ASMDIR_BUILDH; }
      if (version) { 
	compute_constraints(version);
	this->parameters_set_uptodate();
	B_to_be_computed = false;
      }
    }

   /** ensure a consistent dimension for the data */
    void reshape_coeff() {
      size_type N = this->mf_u().linked_mesh().dim();
      switch (R_.fdim()) {
	case 0 : R_.reshape(); break;
	case 1 : R_.reshape(N); break;
	case 2 : R_.reshape(mf_mult.get_qdim(),N); break;
      }
    }

  public :

    /** Change the @f$ r(x) @f$ right hand side.
     *	@param R a vector of size @c mf_data.nb_dof() .
     */
    mdbrick_parameter<VECTOR> &rhs() { reshape_coeff(); return R_; }

    /** Switch between a scalar coefficient, a N vector field or a NxN
	matrix field. */
    void set_coeff_dimension(unsigned d) { R_.redim(d); reshape_coeff(); }

    /** Constructor which does not define the rhs (i.e. which sets an
     *	homogeneous Dirichlet condition)
     *	@param problem the sub problem to which this brick is applied.
     *	@param bound the boundary number for the dirichlet condition.
     *  @param mf_mult_ the mesh_fem for the multipliers.
     *	@param num_fem_ the mesh_fem number on which this brick is is applied.
     */
    mdbrick_normal_component_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
				  size_type bound,
				  const mesh_fem &mf_mult_, 
				  size_type num_fem_=0)
      : mdbrick_constraint<MODEL_STATE>(problem, num_fem_), 
	R_("R", this),
	boundary(bound), mf_mult(mf_mult_) {
      this->add_proper_boundary_info(this->num_fem, boundary, 
				     MDBRICK_DIRICHLET);
      this->add_dependency(mf_mult);
      mfdata_set = false; B_to_be_computed = true;
      this->force_update();
      GMM_ASSERT1((mf_u().get_qdim() % mf_u().linked_mesh().dim()) == 0,
		  "This brick is only working for vectorial elements");
    }
  };


  /* ******************************************************************** */
  /*		Generalized Dirichlet condition bricks.                   */
  /* ******************************************************************** */

  /**
     Generalized Dirichlet condition brick.

     The generalized form for a Dirichlet condition is @f[ \int h(x)u(x).v
     = \int r(x).v \forall v@f] where @f$ h(x) @f$ is a @f$ N \times N
     @f$ matrix (by default, the identity matrix), and @f$ r(x) @f$ is
     the right hand side for the Dirichlet condition (0 for
     homogeneous conditions). This brick if slightly outdated and not very
     well stabilized. Particularly, for an arbitrary @f$ h(x) @f$,
     the multipliers option could not work very well.

     @see asm_generalized_dirichlet_constraints
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_generalized_Dirichlet : public mdbrick_abstract<MODEL_STATE>  {
    
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
    
    const mesh_fem &mf_u() { return *(this->mesh_fems[num_fem]); }

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
	GMM_TRACE2("Assembling Dirichlet constraints with no H and version "
		   << version);
	asm_dirichlet_constraints
	  (M, V, *(this->mesh_ims[0]), mf_u(), mf_u(), rhs().mf(), R_.get(),
	   mf_u().linked_mesh().get_mpi_sub_region(boundary), version);
      } else {
	GMM_TRACE2("Assembling Dirichlet constraints with H and version "
		   << version);
	asm_generalized_dirichlet_constraints
	  (M, V, *(this->mesh_ims[0]), mf_u(), H().mf(), rhs().mf(), H_.get(),
	   R_.get(), mf_u().linked_mesh().get_mpi_sub_region(boundary),
	   version);
      }
      
      if (version & ASMDIR_BUILDH) {
	R tol=gmm::mat_maxnorm(M)*gmm::default_tol(value_type())*R(100);
	gmm::clean(M, tol);
	std::vector<size_type> ind_ct;
	GMM_ASSERT1(!mf_u().is_reduced(), "to be adapted");
	dal::bit_vector nn = mf_u().basic_dof_on_region(boundary);
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
      /* compute_constraints has to be done here because 'nb_const' must
	 be known.. */
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
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      compute_constraints();
      if (with_multipliers) {
	gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), nb_const);
	gmm::sub_interval SUBJ(i0+i1, nbd);
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residual(), SUBI));
	
	gmm::mult_add(gmm::transposed(G), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residual(), SUBJ));
      }
      else {
	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0+i1, nbd);
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.constraints_rhs(), SUBI));
	gmm::copy(G, gmm::sub_matrix(MS.constraints_matrix(), SUBI, SUBJ));
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
    
    /** Return true if the brick is using Lagrange multipliers to enforce
	the Dirichlet condition */
    bool is_using_multipliers() const { return with_multipliers; }
    /** Switch between lagrange multipliers and direct elimination of
	Dirichlet variables */
    void use_multipliers(bool v) {
      if (v != with_multipliers) {
	with_multipliers = v; 
	this->proper_is_coercive_ = !with_multipliers;
	this->change_context();
      }
    }

    /** Constructor which does not define the rhs (i.e. which sets an
	homogeneous Dirichlet condition)
	@param problem the sub problem to which this brick is applied.
	@param bound the boundary number for the dirichlet condition.
	@param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_generalized_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
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
  /*  dynamic brick : not stabilized, could change in future versions.    */
  /* ******************************************************************** */

  /**
     dynamic brick : not stabilized, could change a lot in the future.
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    const mesh_fem *mf_u;
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
      GMM_TRACE2("Assembling mass matrix for mdbrick_dynamic");
      gmm::clear(M_);
      asm_mass_matrix_param(M_, *(this->mesh_ims[0]), *mf_u, RHO_.mf(),
			    RHO_.get());

      if (!(boundary_sup.empty())) {
	GMM_ASSERT1(!mf_u->is_reduced(), "To be adapted");

	gmm::unsorted_sub_index SUBS;
	std::vector<size_type> ind;
	dal::bit_vector ind_set;
	
	for (std::set<size_type>::const_iterator it = boundary_sup.begin();
	     it != boundary_sup.end(); ++it) {
	  ind_set = ind_set | mf_u->basic_dof_on_region(*it);
	}

	VECTOR V(mf_u->nb_dof()), MV(mf_u->nb_dof()); 
	for (size_type i=0; i < V.size(); i += mf_u->get_qdim()) V[i] = 1;
	gmm::mult(M_, V, MV);
	cerr << " VMV = " << gmm::vect_sp(V, MV) << "\n";


	redistribute_mass(ind_set);

	gmm::mult(M_, V, MV);
	cerr << " VMV2 = " << gmm::vect_sp(V, MV) << "\n";


	ind.reserve(ind_set.card());
	for (dal::bv_visitor ii(ind_set); !ii.finished(); ++ii)
	  ind.push_back(ii);
	SUBS = gmm::unsorted_sub_index(ind);
	
	gmm::sub_interval SUBI(0, mf_u->nb_dof());

	/*cerr << "gmm::sub_matrix(M_, SUBS, SUBI)) = " << 
	  gmm::sub_matrix(M_, SUBS, SUBI) << "\n";
	  assert(gmm::nnz(gmm::sub_matrix(M_, SUBI, SUBS)) == 0);*/

	gmm::clear(gmm::sub_matrix(M_, SUBS, SUBI));
	gmm::clear(gmm::sub_matrix(M_, SUBI, SUBS));

	gmm::mult(M_, V, MV);
	cerr << " VMV3 = " << gmm::vect_sp(V, MV) << "\n";

      }
    }

    /* valid only for lagrange FEMs */
    void redistribute_mass(const dal::bit_vector &redistributed_dof) {
      GMM_ASSERT1(mf_u->is_reduced(), "To be adapted");
      size_type N = mf_u->linked_mesh().dim();
      size_type Qdim = mf_u->get_qdim();
      size_type nn = mf_u->nb_dof() / Qdim;

      /* extract the "scalar" mass matrix */
      gmm::csc_matrix<value_type> M0;
      gmm::sub_slice IND0(0, nn, Qdim);
      M0.init_with(gmm::sub_matrix(M_, IND0, IND0));

      /* get the non-zero coef as a vector */
      size_type nz = gmm::nnz(M0);
      gmm::array1D_reference<value_type*> vM0(&M0.pr[0], nz);

      //dal::bit_vector removed; removed.sup(0, nz);
      size_type nbmult = 1+N+N*(N+1)/2;

      std::vector<value_type> F(nz); gmm::copy(vM0, F);
      gmm::dense_matrix<value_type> C(nbmult, nz);
      std::vector<value_type> d(nbmult);

      for (size_type j=0, ii=0; j < gmm::mat_ncols(M0); ++j) {
	for (size_type ir=0; ir < M0.jc[j+1] - M0.jc[j]; ++ir, ++ii) {
	  size_type i=M0.ir[ii];
	  const base_node Pi = mf_u->point_of_basic_dof(i * Qdim);
	  const base_node Pj = mf_u->point_of_basic_dof(j * Qdim);
	  C(0, ii) = 1; // X'MX 
	  d[0] += vM0[ii];
	
	  for (unsigned k=0; k < N; ++k) {
	    C(1+k, ii) = Pi[k]; // X'MYk = 0
	    d[1+k] += vM0[ii]*Pi[k];
	  }
	  
	  for (unsigned k=0,cnt=0; k < N; ++k) {
	    for (unsigned l=k; l < N; ++l, ++cnt) {
	      C(1+N+cnt, ii) = Pi[k] * Pj[l];
	      d[1+N+cnt] += vM0[ii] * Pi[k] * Pj[l];
	    }
	  }

	  if (redistributed_dof.is_in(i*Qdim) ||
	      redistributed_dof.is_in(j*Qdim)) {
	    for (unsigned k=0; k < gmm::mat_nrows(C); ++k) {
	      C(k,ii) = 0; 
	    }
	    F[ii] = 0;
	  }
	}
      }

      /* solve [ I C'][X]   [F]
               [ C 0 ][L] = [d] 
      */
      gmm::dense_matrix<value_type> CCt(nbmult, nbmult);
      std::vector<value_type> L(nbmult), CF(nbmult);
      gmm::mult(C, gmm::transposed(C), CCt);
      gmm::mult(C, F, CF);
      gmm::add(gmm::scaled(d, -1), CF);
      gmm::lu_solve(CCt, L, CF);
      gmm::mult(gmm::transposed(C), gmm::scaled(L, -1), vM0); gmm::add(F, vM0);


      /* force the symmetricity */
      gmm::copy(vM0, F);
      for (size_type j=0, ii=0; j < gmm::mat_ncols(M0); ++j) {
	for (size_type ir=0; ir < M0.jc[j+1] - M0.jc[j]; ++ir, ++ii) {
	  size_type i=M0.ir[ii];
	  F[ii] = (M0(i,j) + M0(j,i))/value_type(2);
	}
      }
      gmm::copy(F, vM0);

      
      //gmm::HarwellBoeing_IO::write("redist_mass.hb", M0);


      /* write back the mass matrix */

      for (unsigned k=0; k < Qdim; ++k) {
	gmm::sub_slice IND(k, nn, Qdim);
	gmm::copy(M0, gmm::sub_matrix(M_, IND, IND));
      }

      /* some sanity checks */
      for (unsigned i=0; i < gmm::mat_nrows(M0); ++i) {
	if (gmm::real(M0(i,i)) < 0) {
	  GMM_WARNING1("negative diagonal terms found in the mass matrix!");
	  break;
	}
      }
      std::vector<value_type> e(nbmult);
      for (size_type j=0, ii=0; j < gmm::mat_ncols(M0); ++j) {
	for (size_type ir=0; ir < M0.jc[j+1] - M0.jc[j]; ++ir, ++ii) {
	  size_type i=M0.ir[ii];
	  if (i > j) continue; // deal with upper triangle only
	  const base_node Pi = mf_u->point_of_basic_dof(i * Qdim);
	  const base_node Pj = mf_u->point_of_basic_dof(j * Qdim);
	  e[0] += vM0[ii];
	
	  for (unsigned k=0; k < N; ++k) {
	    e[1+k] += vM0[ii]*Pi[k];
	  }
	  
	  for (unsigned k=0,cnt=0; k < N; ++k) {
	    for (unsigned l=k; l < N; ++l, ++cnt) {
	      e[1+N+cnt] += vM0[ii] * Pi[k] * Pj[l];
	    }
	  }
	}
      }
      // cerr << "d = " << d << "\ne = " << e << "\n";
    }

  public :

    mdbrick_parameter<VECTOR> &rho() { return RHO_; M_uptodate = true; }
    const mdbrick_parameter<VECTOR> &rho() const { return RHO_; }

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1)) gmm::scale(MS.tangent_matrix(), Kcoef);
      gmm::add(gmm::scaled(get_M(), Mcoef),
	       gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1))  gmm::scale(MS.residual(), Kcoef);
      gmm::add(gmm::scaled(DF, -value_type(1)),
	       gmm::sub_vector(MS.residual(), SUBI));
      gmm::mult_add(get_M(),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Mcoef),
		    gmm::sub_vector(MS.residual(), SUBI));
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
  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
