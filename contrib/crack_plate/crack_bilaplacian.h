// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
//  templates or use macros or inline functions from this file, or you compile
//  this file and link it with other files to produce an executable, this
//  file does not by itself cause the resulting executable to be covered by
//  the GNU General Public License.  This exception does not however
//  invalidate any other reasons why the executable file might be covered by
//  the GNU General Public License.
// 
//===========================================================================

#include "getfem/getfem_config.h"
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"
#include "getfem/getfem_modeling.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_mesh_fem_sum.h"


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector; /* dense vector. */
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


/******** Exact Solution *******************************/


scalar_type sol_u(const base_node &x);
scalar_type sol_lapl_u(const base_node &x);
scalar_type sol_f(const base_node &);
base_small_vector sol_du(const base_node &x);
base_small_vector neumann_val(const base_node &x);


struct bilaplacian_singular_functions : public getfem::global_function, public getfem::context_dependencies {
  size_type l ;            // singular function number
  const getfem::level_set &ls;
  scalar_type nu;          // Poison's coefficient
  scalar_type pos ;        // x-position of the crack-tip
  mutable getfem::mesher_level_set mls0, mls1;
  mutable size_type cv;
  
  void update_mls(size_type cv_) const;
  scalar_type sing_function(scalar_type x, scalar_type y) const;
  /* derivative in the levelset coordinates */
  void sing_function_grad(scalar_type x, scalar_type y, base_small_vector &g) const;
  void sing_function_hess(scalar_type x, scalar_type y, base_matrix &he) const;  
  virtual scalar_type val(const getfem::fem_interpolation_context& c) const;
  virtual void grad(const getfem::fem_interpolation_context& c,
		    base_small_vector &v) const;
  virtual void hess(const getfem::fem_interpolation_context& c, base_matrix &he) const;
  void update_from_context(void) const;
  bilaplacian_singular_functions(size_type l_, const getfem::level_set &ls_, scalar_type nu, scalar_type pos_ );
};

inline getfem::pglobal_function bilaplacian_crack_singular(size_type i, const getfem::level_set &ls, scalar_type nu, scalar_type pos){ 
  return new bilaplacian_singular_functions(i, ls, nu, pos);
}

struct exact_solution {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  exact_solution(getfem::mesh &me) : mf(me) {}
  void init(getfem::level_set &ls);
};

inline std::string name_of_dof(getfem::pdof_description dof) {
  char s[200];
  sprintf(s, "UnknownDof[%p]", (void*)dof);
  for (unsigned d = 0; d < 4; ++d) {
    if (dof == getfem::lagrange_dof(d)) {
      sprintf(s, "Lagrange[%d]", d); goto found;
    }
    if (dof == getfem::normal_derivative_dof(d)) {
      sprintf(s, "D_n[%d]", d); goto found;
    }
    if (dof == getfem::global_dof(d)) {
      sprintf(s, "GlobalDof[%d]", d);
    }
    if (dof == getfem::mean_value_dof(d)) {
      sprintf(s, "MeanValue[%d]", d);
    }
    if (getfem::dof_xfem_index(dof) != 0) {
      sprintf(s, "Xfem[idx:%d]", int(dof_xfem_index(dof)));
    }
    
    for (unsigned r = 0; r < d; ++r) {
      if (dof == getfem::derivative_dof(d, r)) {
	sprintf(s, "D_%c[%d]", "xyzuvw"[r], d); goto found;
      }
      for (unsigned t = 0; t < d; ++t) {
	if (dof == getfem::second_derivative_dof(d, r, t)) {
	  sprintf(s, "D2%c%c[%d]", "xyzuvw"[r], "xyzuvw"[t], d); 
	  goto found;
	}
      }
    }
  }
 found:
  return s;
}

/******** Struct for the Bilaplacian problem *******************************/


struct bilaplacian_crack_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3,
	 MORTAR_BOUNDARY_IN = 40, MORTAR_BOUNDARY_OUT = 41};
  
  getfem::mesh mesh;        /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_pre_mortar, mf_pre_mortar_deriv ; 
  getfem::mesh_fem_level_set mfls_u, mfls_mortar, mfls_mortar_deriv ; 
  getfem::mesh_fem_global_function mf_sing_u ;
  getfem::mesh_fem mf_partition_of_unity ;
  getfem::mesh_fem_product mf_u_product ;
  getfem::mesh_fem_sum mf_u_sum ;
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
 
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_mult_d; /* mesh_fem for the Dirichlet condition on    */
                              /* the normal derivative.                     */
  
  scalar_type residual;     /* max residual for the iterative solvers       */
  getfem::constraints_type dirichlet_version, mortar_version, closing_version;
  
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option, mortar_type ;
  size_type NX;
  
  std::string datafilename;
  bgeot::md_param PARAM;

  bool KL;

  dal::bit_vector pm_convexes; /* convexes inside the enrichment 
				  area when point-wise matching is used.*/
  
  scalar_type epsilon ;      /* half-plate thickness */
  scalar_type nu ;
  scalar_type D ;
  
  exact_solution exact_sol;
  size_type sol_ref ;    /* 0 : solution made of pure singularities
                            1 : infinite plate with central crack subject to 
                                moments (pure K1 mode) -> work in progress*/
  bool solve(plain_vector &U);
  bool solve_moment(plain_vector &U) ;
  void init(void);
  void compute_error(plain_vector &U);
  void compute_H2_error_field(const plain_vector &U) ;
  void move_nodes_close_to_crack(void); // bugged
  void compute_error_beta(plain_vector &U) ;
  void init_mixed_elements(void) ;
  void set_matrix_mortar(sparse_matrix &H) ;
  void compute_sif(const plain_vector &U, scalar_type ring_radius);
  
  bilaplacian_crack_problem(void) : ls(mesh, 1, true), 
				    mls(mesh), mim(mls), mf_pre_u(mesh),  
				    mf_pre_mortar(mesh), mf_pre_mortar_deriv(mesh), mfls_u(mls, mf_pre_u), 
				    mfls_mortar(mls, mf_pre_mortar), mfls_mortar_deriv(mls, mf_pre_mortar_deriv),  
				    mf_sing_u(mesh), mf_partition_of_unity(mesh),
				    mf_u_product(mf_partition_of_unity, mf_sing_u), mf_u_sum(mesh),
				    mf_rhs(mesh), mf_mult(mesh), mf_mult_d(mesh), exact_sol(mesh)   
				    { KL = true; } 
};


