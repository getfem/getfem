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


/*
  Requirements for a model brick :                                  
                                                                         
  A model brick is either a fondamental brick (like linearized
  elasticity brick, platicity brick ...) or a modifier brick which
  refer to a sub brick.

  A new brick should define:

  - proper_update() , which is called each time the brick should
  update itself.
  This function is expected to assign the correct values to
  'proper_nb_dof' (the nb of new dof introduced by this brick),
  'proper_nb_constraints' and 'proper_mixed_variables'.

  - do_compute_tangent_matrix(MS, i0, j0) . This function should
  compute its own part of the tangent and constraint matrices (i0 and
  j0 are the shifts in the matrices stored in the model_state MS)

  - do_compute_residu(MS, i0, j0) . Same as above for the residu
  vector stored in MS.
*/

/***************************************************************************/
/*                                                                         */
/* Brick idents :                                                          */
/* MDBRICK_SCALAR_ELLIPTIC       174397                                    */
/* MDBRICK_LIN_ISO_ELASTICITY    852327                                    */
/* MDBRICK_MASS_MATRIX           756543                                    */
/* MDBRICK_HELMHOLTZ             354864                                    */
/* MDBRICK_LINEAR_INCOMP         239898                                    */
/* MDBRICK_NONLINEAR_ELASTICITY  821357                                    */
/* MDBRICK_NONLINEAR_INCOMP      964552                                    */
/* MDBRICK_SMALL_DEF_PLASTICITY  556433                                    */
/* MDBRICK_LINEAR_PLATE          897523                                    */
/* MDBRICK_MIXED_LINEAR_PLATE    213456                                    */
/* MDBRICK_COULOMB_FRICTION      434245                                    */
/* MDBRICK_NAVIER_STOKES         394329                                    */
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
#include <gmm_precond_ilutp.h>
#include <gmm_superlu_interface.h>
#include <gmm_dense_qr.h>
#include <gmm_matrix.h>
#include <gmm_solver_Schwarz_additive.h>
#include <set>

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

#define TYPEDEF_MODEL_STATE_TYPES                                         \
    typedef typename MODEL_STATE::vector_type VECTOR;                     \
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;           \
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;       \
    typedef typename MODEL_STATE::value_type value_type;                  \
    typedef typename gmm::number_traits<value_type>::magnitude_type R;    \
    typedef typename gmm::sub_vector_type<VECTOR *,                       \
		  gmm::sub_interval>::vector_type SUBVECTOR


  /*
   *   Abstract model brick
   */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_abstract : public context_dependencies {

  public :
    TYPEDEF_MODEL_STATE_TYPES;

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

    std::vector<mdbrick_abstract *> sub_bricks;
    
    /* all proper_* specify data which is specific to this brick:
       'proper_mesh_fems' is the list of mesh_fems used by this brick,
       while 'mesh_fems' is 'proper_mesh_fems' plus the list of
       mesh_fems of all the parent bricks.
    */
    std::vector<mesh_fem *> proper_mesh_fems;
    std::vector<mesh_im *> proper_mesh_ims;
    std::vector<mesh_fem_info_> proper_mesh_fems_info;
    std::vector<boundary_cond_info_> proper_boundary_cond_info;
    /* flags indicating how this brick affect the linearity/coercivity
       etc properties of the problem */
    bool proper_is_linear_, proper_is_symmetric_, proper_is_coercive_;
    /* number of new degrees of freedom introduced by this brick */
    size_type proper_additional_dof;
    /* number of new constraints introduced by this brick */
    size_type proper_nb_constraints;
    /* in the dofs, indicates which ones correspound to mixed variables */
    dal::bit_vector proper_mixed_variables;

    /* below is the "global" information, relating to this brick and
       all its parent bricks */
    mutable std::vector<mesh_fem *> mesh_fems;
    mutable std::vector<mesh_im *> mesh_ims;
    mutable std::vector<mesh_fem_info_> mesh_fems_info;
    mutable std::vector<size_type> mesh_fem_positions;
    mutable bool is_linear_, is_symmetric_, is_coercive_;
    mutable size_type nb_total_dof, nb_total_constraints;
    mutable dal::bit_vector total_mixed_variables;

    /* the brick owns the block starting at index MS_i0 in the global
       tangent matrix */
    size_type MS_i0;

    /* the brick-specific update procedure */
    virtual void proper_update(void) = 0;
    // The following function is just for the const cast on "this".
    inline void proper_update_(void) { proper_update(); }

    /* method inherited from getfem::context_dependencies */
    void update_from_context(void) const {
      nb_total_dof = 0;
      nb_total_constraints = 0;
      total_mixed_variables.clear();
      is_linear_ = proper_is_linear_;
      is_symmetric_ = proper_is_symmetric_;
      is_coercive_ = proper_is_coercive_;
      mesh_fems.resize(0); mesh_ims.resize(0); mesh_fem_positions.resize(0);
      /* get information from parent bricks */
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	for (size_type j = 0; j < sub_bricks[i]->mesh_fems.size(); ++j) {
	  mesh_fems.push_back(sub_bricks[i]->mesh_fems[j]);
	  mesh_fems_info.push_back(sub_bricks[i]->mesh_fems_info[j]);
	  mesh_fem_positions.push_back(nb_total_dof 
				       + sub_bricks[i]->mesh_fem_positions[j]);
	}
	for (size_type j = 0; j < sub_bricks[i]->mesh_ims.size(); ++j) {
	  mesh_ims.push_back(sub_bricks[i]->mesh_ims[j]);
	}
	is_linear_ = is_linear_ && sub_bricks[i]->is_linear();
	is_symmetric_ = is_symmetric_ && sub_bricks[i]->is_symmetric();
	is_coercive_ = is_coercive_ && sub_bricks[i]->is_coercive();
	if (i == 0)
	  total_mixed_variables |= sub_bricks[i]->total_mixed_variables;
	else {
	  for (dal::bv_visitor k(sub_bricks[i]->total_mixed_variables);
	       !k.finished(); ++k)
	    total_mixed_variables.add(k+nb_total_dof);
	}
	nb_total_dof += sub_bricks[i]->nb_total_dof;
	nb_total_constraints += sub_bricks[i]->nb_total_constraints;
      }
      /* merge with information from this brick */
      for (size_type j = 0; j < proper_mesh_fems.size(); ++j) {
	mesh_fems.push_back(proper_mesh_fems[j]);
	mesh_fems_info.push_back(proper_mesh_fems_info[j]);
	mesh_fem_positions.push_back(nb_total_dof);
	nb_total_dof += proper_mesh_fems[j]->nb_dof();
      }
      for (size_type j = 0; j < proper_mesh_ims.size(); ++j)
	mesh_ims.push_back(proper_mesh_ims[j]);
      for (size_type j = 0; j < proper_boundary_cond_info.size(); ++j) {
	mesh_fems_info[proper_boundary_cond_info[j].num_fem]
	  .add_boundary(proper_boundary_cond_info[j].num_bound,
			proper_boundary_cond_info[j].bc);
      }
      /* call the customizable update procedure */
      const_cast<mdbrick_abstract *>(this)->proper_update_();
      nb_total_dof += proper_additional_dof;
      nb_total_constraints += proper_nb_constraints;
      total_mixed_variables |= proper_mixed_variables;
    }

    void force_update(void)
    { if (!this->context_check()) update_from_context(); }

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
    { return mesh_fems_info[i]; }
    mesh_fem &get_mesh_fem(size_type i) { return *(mesh_fems[i]); }
    size_type get_mesh_fem_position(size_type i)
    { return mesh_fem_positions[i]; }
    size_type nb_mesh_fems(void) { return mesh_fems.size(); }

    dim_type dim(void) { return mesh_fems[0]->linked_mesh().dim(); }
    /** total number of variables including the variables of the
	sub-problem(s) if any */
    size_type nb_dof(void) { return nb_total_dof; }

    /** number of linear constraints on the system including the
	constraints defined in the sub-problem(s) if any. */
    size_type nb_constraints(void) { return nb_total_constraints; };

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) = 0;
    void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0=0) {
      this->context_check();
      size_type i1 = MS_i0 = i0, j1 = j0;
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	sub_bricks[i]->compute_tangent_matrix(MS, i1, j1);
	i1 += sub_bricks[i]->nb_dof();
	j1 += sub_bricks[i]->nb_constraints();
      }
      do_compute_tangent_matrix(MS, i0, j0);
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type j0) = 0;
    void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
			size_type j0 = 0) {
      this->context_check();
      size_type i1 = MS_i0 = i0, j1 = j0;
      for (size_type i = 0; i < sub_bricks.size(); ++i) {
	sub_bricks[i]->compute_residu(MS, i1, j1);
	i1 += sub_bricks[i]->nb_dof();
	j1 += sub_bricks[i]->nb_constraints();
      }
      do_compute_residu(MS, i0, j0);
    }
    bool is_linear(void) const { return is_linear_; }
    bool is_symmetric(void) const { return is_symmetric_; }
    bool is_coercive(void) const { return is_coercive_; }
    const dal::bit_vector &mixed_variables(void) const
    { return total_mixed_variables; };
    mdbrick_abstract(void) : proper_additional_dof(0), proper_nb_constraints(0),
			     MS_i0(0)
    { proper_is_linear_ = proper_is_symmetric_ = proper_is_coercive_ = true; }
    virtual ~mdbrick_abstract() {}
  };

  /* ******************************************************************** */
  /*		Function for extracting matrix structure.                 */
  /* ******************************************************************** */
//   void tangent_matrix_structure_single_mf(MODEL_STATE &MS, mesh_fem &mf_u, 
// 					  size_type i0) {
//     for (size_type i = 0; i < mf_u.nb_dof(); ++i) {
//       bgeot::mesh_convex_ind_ct ct = mf_u.convex_to_dof(i);
//       bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin(), ite = ct.end();
//       for (; it != ite; ++it) {
// 	ref_mesh_dof_ind_ct ctd = mf_u.ind_dof_of_element(*it);
// 	ref_mesh_dof_ind_ct::const_iterator itd = ctd.begin(), itde = ctd.end();
// 	  for (; itd != itde; ++itd) {
// 	    (MS.tangent_matrix())(i0 + i, i0 + *itd) = scalar_type(1);
// 	  }
//       }
//     }
//   }

//   void tangent_matrix_structure_double_mf(MODEL_STATE &MS, mesh_fem &mf_u, mesh_fem &mf_v, 
// 					  size_type i0, size_type i1) {
//     for (size_type i = 0; i < mf_u.nb_dof(); ++i) {
//       bgeot::mesh_convex_ind_ct ct = mf_u.convex_to_dof(i);
//       bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin(), ite = ct.end();
//       for (; it != ite; ++it) {
// 	ref_mesh_dof_ind_ct ctd = mf_v.ind_dof_of_element(*it);
// 	ref_mesh_dof_ind_ct::const_iterator itd = ctd.begin(), itde = ctd.end();
// 	  for (; itd != itde; ++itd) {
// 	    (MS.tangent_matrix())(i0 + i, i1 + *itd) = scalar_type(1);
// 	    (MS.tangent_matrix())(i1 + *itd, i0 + i) = scalar_type(1);
// 	  }
//       }
//     }
//   }

//   void tangent_matrix_structure_single_mf_on_boundary(MODEL_STATE &MS, mesh_fem &mf_u, 
// 						      size_type i0, int boundary) {
//     dal::bit_vector nn = mf_u.dof_on_set(boundary);
//     for (dal::bv_visitor i(nn); !i.finished(); ++i) {
//       bgeot::mesh_convex_ind_ct ct = mf_u.convex_to_dof(i);
//       bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin(), ite = ct.end();
//       for (; it != ite; ++it) {
// 	ref_mesh_dof_ind_ct ctd = mf_u.ind_dof_of_element(*it);
// 	ref_mesh_dof_ind_ct::const_iterator itd = ctd.begin(), itde = ctd.end();
// 	  for (; itd != itde; ++itd) {
// 	    (MS.tangent_matrix())(i0 + i, i0 + *itd) = scalar_type(1);
// 	  }
//       }
//     }
//   }

  /* ******************************************************************** */
  /*		general scalar elliptic brick.                            */
  /* ******************************************************************** */

# define MDBRICK_SCALAR_ELLIPTIC 174397

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_scalar_elliptic : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    mesh_im  &mim;
    mesh_fem &mf_u, &mf_data;
    VECTOR coeffs_;
    bool homogeneous, laplacian;
    T_MATRIX K;

    void proper_update(void) {
      gmm::clear(K);
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      size_type n = laplacian ? 1 : gmm::sqr(mf_u.linked_mesh().dim());
      VECTOR coeffs(n * mf_data.nb_dof());
      if (homogeneous) {
	for (size_type i = 0; i < mf_data.nb_dof(); ++i)
	  gmm::copy(coeffs_,gmm::sub_vector(coeffs,gmm::sub_interval(i*n, n)));
      }
      else { gmm::copy(coeffs_, coeffs); }
      if (laplacian) {
	if (mf_u.get_qdim() > 1)
	  asm_stiffness_matrix_for_laplacian_componentwise(K, mim, mf_u,
							   mf_data, coeffs);
	else
	  asm_stiffness_matrix_for_laplacian(K, mim, mf_u, mf_data, coeffs);
      }
      else {
	if (mf_u.get_qdim() > 1) 
	  asm_stiffness_matrix_for_scalar_elliptic_componentwise(K, mim, mf_u,
							      mf_data,coeffs);
	else
	  asm_stiffness_matrix_for_scalar_elliptic(K, mim, mf_u,
						   mf_data,coeffs);
      }
    }
    
    void set_coeff_(const VECTOR &coeffs, bool laplace) {
      laplacian = laplace;
      homogeneous = false;
      int N = mf_u.linked_mesh().dim();
      if (laplacian) {
	if (gmm::vect_size(coeffs) == 1) set_coeff(coeffs[0]);
	else gmm::resize(coeffs_,  mf_data.nb_dof());
      }
      else {
	if (gmm::vect_size(coeffs) == gmm::sqr(N)) {
	  gmm::resize(coeffs_, gmm::sqr(N));
	  homogeneous = true;
	}
	else
	  gmm::resize(coeffs_, mf_data.nb_dof() * gmm::sqr(N));
      }
      gmm::copy(coeffs, coeffs_);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }

    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI)); 
    }
    
    void set_coeff(value_type a) {
      homogeneous = true; laplacian = true;
      gmm::resize(coeffs_, 1); coeffs_[0] = a;
      this->force_update();
    }

    void set_coeff(const VECTOR &coeffs, bool laplace)
    { set_coeff_(coeffs, laplace); this->force_update(); }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_im(mim);
      this->add_proper_mesh_fem(mf_u, MDBRICK_SCALAR_ELLIPTIC);
      this->update_from_context();
    }

    // constructor for the Laplace operator
    mdbrick_scalar_elliptic(mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
			    value_type a)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_),  homogeneous(true),
	laplacian(true)
    { gmm::resize(coeffs_, 1); coeffs_[0] = a; init_(); }

    // constructor for a non-homogeneous material
    mdbrick_scalar_elliptic(mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
       const VECTOR &coeff, bool laplace) : mim(mim_), mf_u(mf_u_),
					    mf_data(mf_data_)
    { set_coeff_(coeff, laplace); init_(); }
  };

  /* ******************************************************************** */
  /*		Linearized elasticity bricks.                             */
  /* ******************************************************************** */

# define MDBRICK_LIN_ISO_ELASTICITY 852327

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_isotropic_linearized_elasticity
    : public mdbrick_abstract<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mesh_im &mim;
    mesh_fem &mf_u, &mf_data;
    VECTOR lambda_, mu_;
    bool homogeneous;
    T_MATRIX K;

    virtual void proper_update(void) {
      gmm::clear(K);
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
      if (homogeneous) {
	std::fill(lambda.begin(), lambda.end(), value_type(lambda_[0]));
	std::fill(mu.begin(), mu.end(), value_type(mu_[0]));
      }
      else { gmm::copy(lambda_, lambda); gmm::copy(mu_, mu); }
      asm_stiffness_matrix_for_linear_elasticity(K, mim, mf_u, mf_data,
						 lambda, mu
#ifdef GMM_USES_MPI
			     , MS.local_domain(mf_u.linked_mesh())
#endif
			     );
    }

    void set_Lame_coeff_(value_type lambdai, value_type mui) {
      homogeneous = true;
      gmm::resize(lambda_, 1); lambda_[0] = lambdai;
      gmm::resize(mu_, 1); mu_[0] = mui;
    }
    
    void set_Lame_coeff_(const VECTOR &lambdai, const VECTOR &mui) {
      homogeneous = false;
      gmm::resize(lambda_, mf_data.nb_dof()); gmm::copy(lambdai, lambda_);
      gmm::resize(mu_, mf_data.nb_dof()); gmm::copy(mui, mu_);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }

    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }
    
    const T_MATRIX &stiffness_matrix(void)
    { this->context_check(); return K; }

    void set_Lame_coeff(value_type lambdai, value_type mui)
    { set_Lame_coeff_(lambdai, mui); this->force_update(); }
    
    void set_Lame_coeff(const VECTOR &lambdai, const VECTOR &mui)
    { set_Lame_coeff_(lambdai, mui); this->force_update(); }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU = gmm::sub_interval(this->first_index(),
						 mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_u, MDBRICK_LIN_ISO_ELASTICITY);
      this->add_proper_mesh_im(mim);
      this->update_from_context();
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_isotropic_linearized_elasticity
    (mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
     value_type lambdai, value_type mui)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_)
    { set_Lame_coeff_(lambdai, mui); init_(); }

    // constructor for a non-homogeneous material
    mdbrick_isotropic_linearized_elasticity
    (mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
     const VECTOR &lambdai, const VECTOR &mui)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_)
    { set_Lame_coeff_(lambdai, mui); init_(); }
  };


  /* TODO : arbitrary elasticity tensor */

  /* ******************************************************************** */
  /*		Mass matrix bricks.                                       */
  /* ******************************************************************** */

# define MDBRICK_MASS_MATRIX 756543

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mass_matrix
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    mesh_im &mim;
    mesh_fem &mf_u, &mf_data;
    VECTOR rho_;
    bool homogeneous;
    T_MATRIX K;

    virtual void proper_update(void) {
      gmm::clear(K);
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      VECTOR rho(mf_data.nb_dof()), mu(mf_data.nb_dof());
      if (homogeneous) {
	std::fill(rho.begin(), rho.end(), value_type(rho_[0]));
      }
      else { gmm::copy(rho_, rho); }
      asm_mass_matrix_param(K, mim, mf_u, mf_data, rho);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }

    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI)); 
    }

    const T_MATRIX &mass_matrix(void)
    { this->context_check(); return K; }

    void set_rho(value_type rhoi) {
      homogeneous = true;
      gmm::resize(rho_, 1); rho_[0] = rhoi; 
      this->force_update();
    }

    void set_rho(const VECTOR &rhoi) {
      homogeneous = false;
      gmm::resize(rho_, mf_data.nb_dof()); gmm::copy(rhoi, rho_);
      this->force_update();
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_im(mim);
      this->add_proper_mesh_fem(mf_u, MDBRICK_MASS_MATRIX);
      this->update_from_context();
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_mass_matrix
    (mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_, value_type rhoi)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_), homogeneous(true)
    { gmm::resize(rho_, 1); rho_[0] = rhoi; init_(); }

    // constructor for a non-homogeneous material
    mdbrick_mass_matrix
    (mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_, const VECTOR &rhoi)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_), homogeneous(false)
    { gmm::resize(rho_, mf_data.nb_dof()); gmm::copy(rhoi, rho_); init_(); }
  };

  /* ******************************************************************** */
  /*		Helmholtz brick.                                          */
  /* ******************************************************************** */

# define MDBRICK_HELMHOLTZ 354864

 template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Helmholtz
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    mesh_im &mim;
    mesh_fem &mf_u, &mf_data;
    VECTOR wave_number;
    bool homogeneous;
    T_MATRIX K;

    virtual void proper_update(void) {
      gmm::clear(K);
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      VECTOR wave_number2(mf_data.nb_dof());
      if (homogeneous)
	std::fill(wave_number2.begin(), wave_number2.end(),
		  value_type(gmm::sqr(wave_number[0])));
      else
	for (size_type i=0; i < mf_u.nb_dof(); ++i)
	  wave_number2[i] = gmm::sqr(wave_number[i]);
      
      asm_Helmholtz(K, mim, mf_u, mf_data, wave_number2);
    }

    void set_wave_number_(value_type k)
    { homogeneous = true; gmm::resize(wave_number, 1); wave_number[0] = k; }

    void set_wave_number_(const VECTOR &k) {
      homogeneous = false;
      gmm::resize(wave_number, mf_data.nb_dof()); gmm::copy(k, wave_number);
    }


  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, mf_u.nb_dof());
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
    }

    void set_wave_number(value_type k)
    { set_wave_number_(k); this->force_update(); }

    void set_wave_number(const VECTOR &k)
    { set_wave_number_(k); this->force_update(); }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_u.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_proper_mesh_fem(mf_u, MDBRICK_HELMHOLTZ);
      this->add_proper_mesh_im(mim);
      this->proper_is_coercive_ = false;
      this->update_from_context();
    }

    // constructor for a homogeneous wave number
    mdbrick_Helmholtz(mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
		      value_type k) : mim(mim_), mf_u(mf_u_), mf_data(mf_data_)
    { set_wave_number_(k); init_(); }

    // constructor for a non-homogeneous wave number
    mdbrick_Helmholtz(mesh_im &mim_, mesh_fem &mf_u_, mesh_fem &mf_data_,
		      const VECTOR &k)
      : mim(mim_), mf_u(mf_u_), mf_data(mf_data_)
    { set_wave_number_(k); init_(); }
    
  };


  /* ******************************************************************** */
  /*		Source term brick.                                        */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_source_term : public mdbrick_abstract<MODEL_STATE>  {

    TYPEDEF_MODEL_STATE_TYPES;

    mesh_fem &mf_data;
    VECTOR B_, F_, auxF;
    size_type boundary, num_fem, i1, nbd;
    bool have_auxF;

    void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      size_type qmult = mf_u.get_qdim() / mf_data.get_qdim();
      if (gmm::vect_size(B_) != mf_data.nb_dof() * qmult) 
	DAL_THROW(failure_error, "The data mesh fem structure has changed, "
		  " You have to change the rhs in that case.");
      gmm::resize(F_, mf_u.nb_dof());
      gmm::clear(F_);
      asm_source_term(F_, *(this->mesh_ims[0]), mf_u, mf_data, B_, boundary);
    }

  public :

    const VECTOR &source_term(void)
    { this->context_check(); return F_; }

    template <class VECT> void set_auxF(const VECT &V) {
      have_auxF = true;
      gmm::resize(auxF, (this->mesh_fems[num_fem])->nb_dof());
      gmm::copy(V, auxF);
    }

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) { }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::add(gmm::scaled(F_, value_type(-1)), gmm::sub_vector(MS.residu(),
	       gmm::sub_interval(i0+i1, nbd)));
      if (have_auxF)
	gmm::add(gmm::scaled(auxF, value_type(-1)),
		 gmm::sub_vector(MS.residu(),
				 gmm::sub_interval(i0+i1, nbd)));
    }

    void set_rhs(const VECTOR &B__) {
      gmm::resize(B_, gmm::vect_size(B__)); gmm::copy(B__, B_);
      this->force_update();
    }    
    // Constructor defining the rhs
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			mesh_fem &mf_data_, const VECTOR &B__,
			size_type bound = size_type(-1), size_type num_fem_=0)
      : mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), have_auxF(false) {
      this->add_dependency(mf_data);
      this->add_sub_brick(problem);
      if (bound != size_type(-1))
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
      gmm::resize(B_, gmm::vect_size(B__));
      gmm::copy(B__, B_);
      this->update_from_context();
    }
  };

  /* ******************************************************************** */
  /*		Q.U term (for Fourier-Robin conditions)                   */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_QU_term
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;
   
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR Q;
    size_type boundary, num_fem, i1, nbd;
    bool homogeneous;
    T_MATRIX K;

    virtual void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      gmm::clear(K);
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      size_type N2 = gmm::sqr(mf_u.get_qdim());
      VECTOR vQ(mf_data.nb_dof() * N2);
      if (homogeneous) {
	for (size_type i=0; i < mf_data.nb_dof(); ++i) {
	  for (size_type j=0; j < N2; ++j)
	    vQ[i*N2 + j] = (j % (mf_u.get_qdim()+1)) == 0 ? Q[0] : 0.;
	}
      }
      else gmm::copy(Q, vQ);
      asm_qu_term(K, *(this->mesh_ims[0]), mf_u, mf_data, vQ, boundary);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+i1, nbd);
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0+i1, nbd);
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residu(), SUBI);
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI), SUBV, SUBV);
    }

    void set_Q(value_type q) {
      homogeneous = true;
      gmm::resize(Q, 1); Q[0] = q;
      this->force_update();
    }

    void set_Q(const VECTOR &q) {
      homogeneous = false;
      gmm::resize(Q, mf_data.nb_dof()
		  * gmm::sqr(this->mesh_fems[num_fem]->get_qdim())); 
      gmm::copy(q, Q);
      this->force_update();
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      if (boundary != size_type(-1))
	this->add_proper_boundary_info(num_fem,boundary,MDBRICK_FOURIER_ROBIN);
      this->update_from_context();
    }

    // Constructor which homogeneous diagonal Q
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, value_type q=value_type(1),
		    size_type bound = size_type(-1), size_type num_fem_=0) 
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), homogeneous(true)
    { gmm::resize(Q, 1); Q[0] = q; init_(); }

    // Constructor defining an arbitrary Q
    mdbrick_QU_term(mdbrick_abstract<MODEL_STATE> &problem,
		    mesh_fem &mf_data_, const VECTOR &q,
		    size_type bound = size_type(-1), size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), homogeneous(false) {
      gmm::resize(Q, mf_data.nb_dof()
		  * gmm::sqr(this->mesh_fems[num_fem]->get_qdim())); 
      gmm::copy(q, Q); init_();
    }
  };


  /* ******************************************************************** */
  /*		Mixed linear incompressible condition brick.              */
  /* ******************************************************************** */

# define MDBRICK_LINEAR_INCOMP 239898

  // for nearly incompressible elasticity, p = -lambda div u
  // sigma = 2 mu epsilon(u) -p I

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_linear_incomp : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_p, &mf_data;
    T_MATRIX B, M;
    bool penalized, homogeneous;
    VECTOR epsilon_; // penalization coefficient if any.
    size_type num_fem, i1, nbd;

    void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
      i1 = this->mesh_fem_positions.at(num_fem);
      nbd = mf_u.nb_dof();
      size_type nd = mf_u.nb_dof(), ndd = mf_p.nb_dof();
      gmm::clear(B); gmm::resize(B, ndd, nd);
      asm_stokes_B(B, *(this->mesh_ims.at(0)), mf_u, mf_p
#ifdef GMM_USES_MPI
		   , MS.local_domain(mf_u.linked_mesh())
#endif
		   );
      if (penalized) {
	VECTOR epsilon(mf_data.nb_dof());
	if (homogeneous) std::fill(epsilon.begin(), epsilon.end(),epsilon_[0]);
	else gmm::copy(epsilon_, epsilon);
	gmm::clear(M); gmm::resize(M, ndd, ndd);
	asm_mass_matrix_param(M, *(this->mesh_ims[0]), mf_p, mf_data, epsilon
#ifdef GMM_USES_MPI
			      , size_type(-1), MS.local_domain(mf_u.linked_mesh())
#endif
			      );
	gmm::scale(M, value_type(-1));
      }
      this->proper_mixed_variables.clear();
      this->proper_mixed_variables.add(sub_problem.nb_dof(), mf_p.nb_dof());
    }

  public :
    
    const T_MATRIX &get_B(void) const { this->context_check(); return B; }
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::copy(B, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::copy(gmm::transposed(B),
		gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      if (penalized)
	gmm::copy(M, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
      else
	gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), mf_p.nb_dof());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::mult(B, gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(gmm::transposed(B), gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBJ));
    }

    SUBVECTOR get_pressure(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + sub_problem.nb_dof(),
			     mf_p.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

     void set_coeff(value_type epsiloni) {
      homogeneous = true;
      gmm::resize(epsilon_, 1); epsilon_[0] = epsiloni;
      this->force_update();
    }

    void set_coeff(const VECTOR &epsiloni) {
      homogeneous = false;
      gmm::resize(epsilon_, mf_data.nb_dof()); gmm::copy(epsiloni, epsilon_);
      this->force_update();
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
    { init_(); }

    // Constructor for the nearly incompressibility condition
    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_p_, mesh_fem &mf_data_, value_type epsilon,
			  size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), mf_data(mf_data_),
	penalized(true), homogeneous(true), num_fem(num_fem_)
    { gmm::resize(epsilon_, 1); epsilon_[0] = epsilon; init_(); }

    mdbrick_linear_incomp(mdbrick_abstract<MODEL_STATE> &problem,
			  mesh_fem &mf_p_, mesh_fem &mf_data_,
			  const VECTOR& epsilon, size_type num_fem_=0)
      : sub_problem(problem), mf_p(mf_p_), mf_data(mf_data_),
	penalized(true), homogeneous(false), num_fem(num_fem_) {
      gmm::resize(epsilon_, mf_data.nb_dof()); gmm::copy(epsilon, epsilon_);
      init_();
    }

  };


  /* ******************************************************************** */
  /*		Dirichlet condition bricks.                               */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Dirichlet : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR B_, H_;
    C_MATRIX G;
    VECTOR CRHS;
    size_type boundary, nb_const, num_fem;
    bool with_H, with_multipliers;
    gmm::sub_index SUB_CT;
    size_type i1, nbd;

    void compute_constraints(int version) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      // size_type Q = mf_u.get_qdim();
      size_type nd = mf_u.nb_dof(); //  ndd = mf_data.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(nd, nd);
      VECTOR V(nd);

      if (!with_multipliers) version |= ASMDIR_SIMPLIFY;
      if (!with_H) {
	asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u, mf_data,
				  B_, boundary, version);
      } else {
	asm_dirichlet_constraints(M, V, *(this->mesh_ims[0]), mf_u, mf_data,mf_data,
				  H_, B_, boundary, version);
      }

      if (version & ASMDIR_BUILDH) {
	R tol=gmm::mat_maxnorm(M)*gmm::default_tol(value_type())*R(100);
	gmm::clean(M, tol);
	std::vector<size_type> ind_ct;
	dal::bit_vector nn = mf_u.dof_on_set(boundary);
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
    }

    virtual void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      if (gmm::vect_size(B_) == 0)
	gmm::resize(B_, mf_data.nb_dof() * mf_u.get_qdim());
      compute_constraints(ASMDIR_BUILDH + ASMDIR_BUILDR);
      this->proper_mixed_variables.clear();
      this->proper_additional_dof = with_multipliers ? nb_const : 0;
      this->proper_nb_constraints = with_multipliers ? 0 : nb_const;
      if (with_multipliers)
	this->proper_mixed_variables.add(sub_problem.nb_dof(), nb_const);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {
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
      gmm::resize(B_, gmm::vect_size(B__));
      gmm::copy(B__, B_);
      if (!(this->context_check())) compute_constraints(ASMDIR_BUILDR);
    }

    void init_(void) {
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = !with_multipliers;
      this->add_proper_boundary_info(num_fem, boundary, MDBRICK_DIRICHLET);
      this->update_from_context();
    }

    // Constructor which does not define the rhs (0 rhs in fact)
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_data_, size_type bound,
		      size_type num_fem_=0, bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), with_H(false), with_multipliers(with_mult) {
      init_();
    }

    // Constructor defining the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		      mesh_fem &mf_data_, const VECTOR &B__,
		      size_type bound, size_type num_fem_=0, 
		      bool with_mult = false)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	num_fem(num_fem_), with_H(false), with_multipliers(with_mult) {
       gmm::resize(B_, gmm::vect_size(B__)); gmm::copy(B__, B_);  init_(); 
    }
  };

  /* ******************************************************************** */
  /*	                 Constraint brick.                                */
  /* ******************************************************************** */

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
      this->update_from_context();
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

    template <class MAT, class VEC>
    void set_constraints(const MAT &G_, const VEC &RHS) {
      set_constraints_(G_, RHS);
      this->force_update();
    }

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

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data, *mf_u;
    VECTOR RHO_, DF;
    T_MATRIX M_;
    size_type num_fem;
    value_type Mcoef, Kcoef;
    gmm::unsorted_sub_index SUBS;
    std::vector<size_type> ind;
    bool have_subs;

    virtual void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      gmm::clear(M_); gmm::resize(M_, mf_u->nb_dof(), mf_u->nb_dof());
      asm_mass_matrix_param(M_, *(this->mesh_ims[0]), *mf_u, mf_data, RHO_);
      if (have_subs) {
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
      gmm::add(gmm::scaled(M_, Mcoef),
	       gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval
	SUBI(i0+this->mesh_fem_positions[num_fem], mf_u->nb_dof());
      if (Kcoef != value_type(1))  gmm::scale(MS.residu(), Kcoef);
      gmm::add(gmm::scaled(DF, -value_type(1)),
	       gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(M_, gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Mcoef),
		    gmm::sub_vector(MS.residu(), SUBI));
    }

    void set_dynamic_coeff(value_type a, value_type b) { Mcoef=a; Kcoef=b; }
    template <class VEC> void set_DF(const VEC &DF_)
    { gmm::resize(DF, gmm::vect_size(DF_)); gmm::copy(DF_, DF); }

    const T_MATRIX &mass_matrix(void) const
    { this->context_check(); return M_; }

    void init(void) {
      Mcoef = Kcoef = value_type(1);
      have_subs = false;
      this->add_dependency(mf_data);
      this->add_sub_brick(sub_problem);
      this->update_from_context();
    }

    void no_mass_on_boundary(size_type b) {
      dal::bit_vector bv = mf_u->dof_on_set(b);
      for (dal::bv_visitor i(bv); !i.finished(); ++i)
	ind.push_back(i);
      SUBS = gmm::unsorted_sub_index(ind);
      have_subs = true;
      if (!(this->context_check())) {
	gmm::sub_interval SUBI(0, mf_u->nb_dof());
	gmm::clear(gmm::sub_matrix(M_, SUBS, SUBI));
	gmm::clear(gmm::sub_matrix(M_, SUBI, SUBS));
      }
    }

    template <class VEC>
    mdbrick_dynamic(mdbrick_abstract<MODEL_STATE> &problem, mesh_fem &mf_data_,
		    const VEC &RHO__, size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), num_fem(num_fem_)
    { gmm::resize(RHO_, mf_data.nb_dof()); gmm::copy(RHO__, RHO_); init();  }

    mdbrick_dynamic(mdbrick_abstract<MODEL_STATE> &problem, mesh_fem &mf_data_,
		    value_type RHO__, size_type num_fem_=0)
      : sub_problem(problem), mf_data(mf_data_), num_fem(num_fem_) {
      gmm::resize(RHO_, mf_data.nb_dof());
      std::fill(RHO_.begin(), RHO_.end(), RHO__);
      init();
    }
  };

  /* ******************************************************************** */
  /*		Generic solvers.                                          */
  /* ******************************************************************** */

  // standard_solve represent a default solver for the model brick system.
  // Of course it could be not adapted for a particular problem, so it could
  // be copied and adapted to change solvers, add a special traitement on the
  // problem, etc ...
  // This is in fact a model for your own solver.

  template <typename T> struct sort_abs_val_
  { bool operator()(T x, T y) { return (gmm::abs(x) < gmm::abs(y)); } };

  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
	gmm::iteration &iter) {

    TYPEDEF_MODEL_STATE_TYPES;
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename gmm::number_traits<value_type>::magnitude_type mtype;

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

#ifdef GMM_USES_MPI
    double t_init = MPI_Wtime();
#endif
    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before
    problem.compute_residu(MS);
    problem.compute_tangent_matrix(MS);

    // cout << "CM = " << MS.constraints_matrix() << endl;

    MS.compute_reduced_system();
#ifdef GMM_USES_MPI
    cout << "comput tangent residu reduction time = " << MPI_Wtime() - t_init << endl;
#endif
#ifdef GMM_USES_METIS
    
    double t_ref = MPI_Wtime();

    std::set<const getfem_mesh *> mesh_set;
    for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
      mesh_set.insert(&(problem.get_mesh_fem(i).linked_mesh()));
    }

    cout << "You have " << mesh_set.size() << " meshes\n";

    std::vector< std::vector<int> > eparts(mesh_set.size());
    size_type nset = 0;
    int nparts = 8;//(ndof / 1000)+1; // number of sub-domains.

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
      // Adapter la boucle aux transformations d'ordre élevé.
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
	  // Il faudrait faire une orthogonalisation de Schmidt ici pour des cas
	  // plus ellaborés.
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
#ifdef GMM_USES_METIS
        DAL_THROW(failure_error, "oups ...");
#endif
      }
      
//       if (iter.get_noisy())
//      cout << "tangent matrix " << MS.tangent_matrix() << endl;
// 	cout << "tangent matrix is "
// 	   << (gmm::is_symmetric(MS.tangent_matrix(),
//             1E-6 * gmm::mat_maxnorm(MS.tangent_matrix())) ? "" : "not ")
// 	   <<  "symmetric. ";

//       cout << "MM = " << MS.reduced_tangent_matrix() << endl;

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
	cout << "there is " << gmm::mat_nrows(MS.constraints_matrix()) << " constraints\n";
	mixvar = problem.mixed_variables();
	cout << "there is " << mixvar.card() << " mixed variables\n";
      }

#ifdef GMM_USES_METIS
#ifdef GMM_USES_MPI
    double t_ref,t_final;
    t_ref=MPI_Wtime();
    cout<<"begin Seq AS"<<endl;
#endif
      sequential_additive_schwarz(MS.reduced_tangent_matrix(), dr,
				  gmm::scaled(MS.reduced_residu(), value_type(-1)),
				  0, Bib, iter_linsolv, gmm::using_superlu(),
				  gmm::using_cg());
#ifdef GMM_USES_MPI
    t_final=MPI_Wtime();
    cout<<"temps Seq AS "<< t_final-t_ref<<endl;
#endif
#else
    size_type dim = problem.dim();

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
#endif

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
