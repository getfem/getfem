/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2012-2012 Yves Renard.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**
 * Large sliding contact.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
 */


#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_contact_and_friction_integral.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include <fstream>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
 * structure for the friction problem
 */
struct contact_problem {

  enum {
    DIRICHLET_BOUNDARY, CONTACT_BOUNDARY1, CONTACT_BOUNDARY2
  };
  getfem::mesh mesh1, mesh2;  /* the mesh */
  getfem::mesh_im  mim1, mim2; /* integration methods.                       */
  getfem::mesh_fem mf_u1, mf_u2; /* main mesh_fems.                          */
  getfem::mesh_fem mf_rhs1, mf_rhs2;   /* mesh_fems for the right hand side. */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type friction_coef; /* friction coefficient.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  size_type N, noisy;

  std::string datafilename, INTEGRATION;
  bgeot::md_param PARAM;

  void solve(void);
  void init(void);
  contact_problem(void) : mim1(mesh1), mim2(mesh2), mf_u1(mesh1), mf_u2(mesh2), mf_rhs1(mesh1), mf_rhs2(mesh2) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void contact_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  INTEGRATION = PARAM.string_value("INTEGRATION","Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  std::string meshname1
    (PARAM.string_value("MESHNAME1", "Name of file for the first mesh"));

  getfem::import_mesh(meshname1, mesh1);
  N = mesh1.dim();
 
  std::string meshname2
    (PARAM.string_value("MESHNAME2", "Name of file for the second mesh"));

  getfem::import_mesh(meshname2, mesh2);
  GMM_ASSERT1(N == mesh2.dim(), "Meshes of incompatible dimensions");;

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");

  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");
 
  noisy = (PARAM.int_value("NOISY", "verbosity of iterative methods") != 0);
  
  mf_u1.set_qdim(dim_type(N));
  mf_u2.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim1.set_integration_method(mesh1.convex_index(), ppi);
  mim2.set_integration_method(mesh2.convex_index(), ppi);
  mf_u1.set_finite_element(mesh1.convex_index(), pf_u);
  mf_u2.set_finite_element(mesh2.convex_index(), pf_u);
 
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_ASSERT1(false, "You are using a non-lagrange FEM. "
		  << "In that case you need to set "
		  << "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs1.set_finite_element(mesh1.convex_index(), pf_u);
    mf_rhs2.set_finite_element(mesh2.convex_index(), pf_u);
  } else {
    mf_rhs1.set_finite_element(mesh1.convex_index(), 
			       getfem::fem_descriptor(data_fem_name));
    mf_rhs2.set_finite_element(mesh2.convex_index(), 
			       getfem::fem_descriptor(data_fem_name));
  }
  
  getfem::mesh_region border_faces1;
  getfem::outer_faces_of_mesh(mesh1, border_faces1);
  for (getfem::mr_visitor i(border_faces1); !i.finished(); ++i)
    mesh1.region(CONTACT_BOUNDARY1).add(i.cv(), i.f());

  getfem::mesh_region border_faces2;
  getfem::outer_faces_of_mesh(mesh2, border_faces2);
  for (getfem::mr_visitor i(border_faces2); !i.finished(); ++i)
    mesh2.region(CONTACT_BOUNDARY2).add(i.cv(), i.f());
}


/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void contact_problem::solve(void) {
  cout << "Number of dof for u: " << mf_u1.nb_dof() + mf_u2.nb_dof() << endl;

  getfem::model model;

  // Main unknowns of the problem.
  model.add_fem_variable("u1", mf_u1);
  model.add_fem_variable("u2", mf_u2);

  // Linearized elasticity bricks.
  model.add_initialized_scalar_data("lambda", lambda);
  model.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim1, "u1", "lambda", "mu");
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim2, "u2", "lambda", "mu");
  
  cout << "Number of dof of the model: " << model.nb_dof() << endl;

  getfem::contact_frame cf(N);
  cf.add_obstacle("x");
  cf.add_boundary(mf_u1, model.real_variable("u1"), CONTACT_BOUNDARY1);
  cf.add_boundary(mf_u2, model.real_variable("u2"), CONTACT_BOUNDARY2);

  getfem::test_contact_frame(cf, mim1, mim2);
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {


  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    contact_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.solve();
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
