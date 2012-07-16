/*===========================================================================
 
 Copyright (C) 2012-2012 Yves Renard.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
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
#include "getfem/getfem_contact_and_friction_common.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include <fstream>

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

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





base_small_vector normal_vect(base_matrix A, base_small_vector n0) {
  base_small_vector n = n0;
  gmm::lu_inverse(A);
  gmm::mult(gmm::transposed(A), n0, n);
  n /= gmm::vect_norm2(n);
  return n;
}



base_small_vector grad_normal_vect(base_matrix A, base_small_vector n0, base_matrix Dir) {

  base_small_vector d = n0, d2 = n0;
  base_small_vector n = normal_vect(A, n0);

  gmm::lu_inverse(A);
  gmm::mult(gmm::transposed(Dir), n, d);
  gmm::mult(gmm::transposed(A), gmm::scaled(d, -1.), d2);
  gmm::add(d2, gmm::scaled(n, -gmm::vect_sp(d2, n)), d);
  return d;
}


void grad_normal_test(void) {

  size_type N = 3;
  base_small_vector n0(N), grad1(N), grad2(N);
  base_matrix dir(N, N), x(N, N);
  scalar_type EPS = 1E-8;
  
  
  for (size_type i = 0; i < 1000; ++i) {
    // The normal vector n is not supposed to be unitary
    gmm::fill_random(x); gmm::fill_random(dir); gmm::fill_random(n0);
    
    
    base_matrix xe = x;
    gmm::add(gmm::scaled(dir, EPS), xe);
    
    base_small_vector n1 = normal_vect(x, n0);
    base_small_vector n2 = normal_vect(xe, n0);
    grad1 = (n2-n1)/EPS;
    grad2 = grad_normal_vect(x, n0, dir);
    
    scalar_type err = gmm::vect_norm2(grad1 - grad2);

    // cout << "grad1 = " << grad1 << "grad2 = " << grad2 << endl;
    // cout << "test " << i << " err = " << err << endl; getchar();
    GMM_ASSERT1(err < 4e-6, "Erroneous gradient of normal vector ");
  }
}




// Test the validity of the two gradients of the De Saxce's projection
void De_Saxce_projection_test(void) {

  size_type N = 3;
  base_small_vector x(N), n(N), dir(N), grad1(N), grad2(N);
  base_matrix g(N, N);
  scalar_type f, EPS = 1E-8;
  
  
  for (size_type i = 0; i < 1000; ++i) {
    // The normal vector n is not supposed to be unitary
    gmm::fill_random(x); gmm::fill_random(dir); gmm::fill_random(n);
    f = gmm::random();
    
    base_small_vector xx = x, xxe = x;
    base_small_vector ne = n + dir*EPS;
    
    getfem::De_Saxce_projection(xx, n, f);
    getfem::De_Saxce_projection(xxe, ne, f);
    grad1 = (xxe-xx)/EPS;
    getfem::De_Saxce_projection_gradn(x, n, f, g);
    gmm::mult(g, dir, grad2);
    
    scalar_type err = gmm::vect_norm2(grad1 - grad2);
    GMM_ASSERT1(err < 1e-6, "Erroneous gradient of De Saxcé projection "
		"with respect to the normal vector"); 
    
    xx = x; xxe = x + dir*EPS;
    getfem::De_Saxce_projection(xx, n, f);
    getfem::De_Saxce_projection(xxe, n, f);
    grad1 = (xxe-xx)/EPS;
    getfem::De_Saxce_projection_grad(x, n, f, g);
    gmm::mult(g, dir, grad2);
    err = gmm::vect_norm2(grad1 - grad2);
    GMM_ASSERT1(err < 1e-6, "Erroneous gradient of De Saxcé projection");
  }
}



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
  getfem::mesh_fem mf_lambda1, mf_lambda2; /* mesh_fems for multipliers.     */
  getfem::mesh_fem mf_rhs1, mf_rhs2;   /* mesh_fems for the right hand side. */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type friction_coef; /* friction coefficient.                        */
  scalar_type R;           /* Augmentation parameter.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  size_type N, noisy;

  std::string datafilename, INTEGRATION;
  bgeot::md_param PARAM;

  void solve(void);
  void init(void);
  contact_problem(void) : mim1(mesh1), mim2(mesh2), mf_u1(mesh1), mf_u2(mesh2), mf_lambda1(mesh1), mf_lambda2(mesh2), mf_rhs1(mesh1), mf_rhs2(mesh2) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void contact_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string MULT_FEM_TYPE  = PARAM.string_value("MULT_FEM_TYPE",
						  "FEM name for multipliers");
  INTEGRATION = PARAM.string_value("INTEGRATION","Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "MULT_FEM_TYPE="  << MULT_FEM_TYPE << "\n";
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
  R = PARAM.real_value("R", "Augmentation parameter");
  if (residual == 0.) residual = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");
 
  noisy = (PARAM.int_value("NOISY", "verbosity of iterative methods") != 0);
  
  mf_u1.set_qdim(dim_type(N));
  mf_u2.set_qdim(dim_type(N));
  mf_lambda1.set_qdim(dim_type(N));
  mf_lambda2.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_lambda = getfem::fem_descriptor(MULT_FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim1.set_integration_method(mesh1.convex_index(), ppi);
  mim2.set_integration_method(mesh2.convex_index(), ppi);
  mf_u1.set_finite_element(mesh1.convex_index(), pf_u);
  mf_u2.set_finite_element(mesh2.convex_index(), pf_u);
  mf_lambda1.set_finite_element(mesh1.convex_index(), pf_lambda);
  mf_lambda2.set_finite_element(mesh2.convex_index(), pf_lambda);

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

  dal::bit_vector dol1 = mf_lambda1.basic_dof_on_region(CONTACT_BOUNDARY1);
  mf_lambda1.reduce_to_basic_dof(dol1);
  dal::bit_vector dol2 = mf_lambda2.basic_dof_on_region(CONTACT_BOUNDARY2);
  mf_lambda2.reduce_to_basic_dof(dol2);
}



void test_tangent_matrix(getfem::model &md) {

  size_type nbdof = md.nb_dof();
  scalar_type EPS = 1E-6;
  scalar_type errmax = scalar_type(0);
  size_type NN = 5;
  scalar_type scale = 0.0001;


  if (md.is_linear()) cout << "Non relevant test" << endl;
  else {
    std::vector<scalar_type> U(nbdof), DIR(nbdof);
    std::vector<scalar_type> D1(nbdof), D2(nbdof), Ddiff(nbdof);
    for (size_type i = 0; i < NN; ++i) {
      gmm::fill_random(U); gmm::scale(U, scale);
      gmm::fill_random(DIR); gmm::scale(DIR, scale);
      md.to_variables(U);
      md.assembly(getfem::model::BUILD_ALL);
      gmm::copy(md.real_rhs(), D2);
      gmm::mult(md.real_tangent_matrix(), DIR, D1);
      gmm::add(gmm::scaled(DIR, EPS), U);
      md.to_variables(U);
      md.assembly(getfem::model::BUILD_RHS);
      gmm::add(gmm::scaled(md.real_rhs(), -scalar_type(1)), D2);
      gmm::scale(D2, scalar_type(1)/EPS);
      scalar_type err = gmm::vect_dist2(D1, D2);
      scalar_type ratio = 0;
      for (size_type j = 0; j < nbdof; ++j) {
	scalar_type r = (D1[j] - D2[j]) / (gmm::abs(D1[j]) + 1E-10);
	ratio = std::max(ratio, r);
	if (r > 0.01) {
	  cout << "WARNING, ratio " << r << " j = " << j << " D1 = "
	       << D1[j] << " D2 = " << D2[j] << " D1[j-1] = " << D1[j-1]
	       << " D1[j+1] = " << D1[j+1] << endl; // getchar();
	}
      }
      

      gmm::add(gmm::scaled(D1, -1.), D2, Ddiff);
      gmm::clean(Ddiff, 1E-9); // à enlever ?
      // cout << "Ddiff = " << Ddiff << endl;
      // cout << "norminf(Ddiff) = " << gmm::vect_norminf(Ddiff) << endl;
      // cout << "D1 = " << gmm::sub_vector(D1, gmm::sub_interval(0, 300)) << endl;
      // cout << "D2 = " << gmm::sub_vector(D2, gmm::sub_interval(0, 300)) << endl;
      // cout << "norm(D1) = " << gmm::vect_norm2(D1)
      //   << " norm(D2) = " << gmm::vect_norm2(D2) << endl;
      // cout << "Max ratio " << ratio << endl;
      // cout << "Error at step " << i << " : " << err << endl; // getchar();
      GMM_ASSERT1(err < 1e-6, "Erroneous tangent matrix");
      errmax = std::max(err, errmax);
    }
  }
}





/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

base_small_vector vol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N); res[N-1] = -1.0;
  return res;
}

void contact_problem::solve(void) {
  cout << "Number of dof for u: " << mf_u1.nb_dof() + mf_u2.nb_dof() << endl;

  getfem::model model;
  std::vector<scalar_type> F;

  bool two_bodies = false;

  // Main unknowns of the problem.
  if (two_bodies) model.add_fem_variable("u1", mf_u1);
  model.add_fem_variable("u2", mf_u2);
  if (two_bodies) model.add_fem_variable("lambda1", mf_lambda1);
  model.add_fem_variable("lambda2", mf_lambda2);
  
  gmm::resize(F, mf_rhs1.nb_dof() * N);
  getfem::interpolation_function(mf_rhs1, F, vol_f);
  model.add_initialized_fem_data("VolumicData1", mf_rhs1, F);
  if (two_bodies)
    getfem::add_source_term_brick(model, mim1, "u1", "VolumicData1");
  gmm::resize(F, mf_rhs2.nb_dof() * N);
  getfem::interpolation_function(mf_rhs2, F, vol_f);
  model.add_initialized_fem_data("VolumicData2", mf_rhs2, F);
  getfem::add_source_term_brick(model, mim2, "u2", "VolumicData2");


  // Pointwise constraints for frictionless problems
  GMM_ASSERT1(N == 2, "Pointwises constraints have to be adapted for 3D");
  std::vector<scalar_type> cpoints(N);
  cpoints[0] = 0.0; cpoints[1] = 0.1;
  std::vector<scalar_type> cunitv(N);
  cunitv[0] = 1.0; cunitv[1] = 0.0;
  model.add_initialized_fixed_size_data("cpoints1", cpoints);
  model.add_initialized_fixed_size_data("cunitv1", cunitv);
  if (two_bodies)
    getfem::add_pointwise_constraints_with_multipliers
      (model, "u1", "cpoints1", "cunitv1");
  cpoints[0] = 0.0; cpoints[1] = 0.0;
  model.add_initialized_fixed_size_data("cpoints2", cpoints);
  model.add_initialized_fixed_size_data("cunitv2", cunitv);
  getfem::add_pointwise_constraints_with_multipliers
    (model, "u2", "cpoints2", "cunitv2");
  
  
  // Linearized elasticity bricks.
  model.add_initialized_scalar_data("lambda", lambda);
  model.add_initialized_scalar_data("mu", mu);
  if (two_bodies)
    getfem::add_isotropic_linearized_elasticity_brick
      (model, mim1, "u1", "lambda", "mu");
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim2, "u2", "lambda", "mu");
  
  cout << "Number of dof of the model: " << model.nb_dof() << endl;

  // Contact brick.
  model.add_initialized_scalar_data("r", R);
  model.add_initialized_scalar_data("f", friction_coef);
  size_type indb = add_integral_large_sliding_contact_brick
    (model, mim2, "u2", "lambda2", "r", "f", CONTACT_BOUNDARY2);

  if (two_bodies) 
    getfem::add_boundary_to_large_sliding_contact_brick
      (model, indb, mim1, "u1", "lambda1", CONTACT_BOUNDARY1);
  
  getfem::add_rigid_obstacle_to_large_sliding_contact_brick(model, indb, "y");


  test_tangent_matrix(model);


  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(model, iter);
  
//   gmm::col_matrix<gmm::wsvector<double> > T(model.nb_dof(), model.nb_dof());
//   gmm::copy(model.real_tangent_matrix(), T);
//   gmm::clean(T, 1E-10);
//   cout << T << endl;

}





  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {


  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  // FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  De_Saxce_projection_test();
  grad_normal_test();

  contact_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.solve();

  return 0; 
}
