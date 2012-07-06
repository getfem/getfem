/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard.
 
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

/**@file heat_equation.cc
   @brief Transient heat equation.
 
   The transient heat equation is solved on a regular mesh of the unit
   square, and is compared to an analytical solution.

   This program is used to check that getfem++ is working. This is
   also a good example of use of Getfem++. This program  does not use the
   model bricks intentionally in order to serve as an exemple of solving
   a pde directly with the assembly procedures.
*/

#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
#include "gmm/gmm.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;
typedef std::vector<scalar_type> plain_vector;

/* Definitions for the exact solution of the Heat_Equation problem,
 *  i.e. Delta(sol_u) + sol_f = 0
 */

base_small_vector sol_K; /* a coefficient */
scalar_type sol_c, sol_t;
/* exact solution */
scalar_type sol_u(const base_node &x) {
  scalar_type k2 = gmm::vect_sp(sol_K, sol_K);
  return (1. - exp(-sol_c*sol_t*k2))*sin(gmm::vect_sp(sol_K, x));
}
/* righ hand side */
scalar_type sol_f(const base_node &x) {
  scalar_type k2 = gmm::vect_sp(sol_K, sol_K);
  return sol_c * k2 * sin(gmm::vect_sp(sol_K, x));
}
/* gradient of the exact solution */
base_small_vector sol_grad(const base_node &x) {
  scalar_type k2 = gmm::vect_sp(sol_K, sol_K);
  return (1. - exp(-sol_c*sol_t*k2))*sol_c*sol_K*cos(gmm::vect_sp(sol_K, x));
}

/*
  structure for the Heat_Equation problem
  (not mandatory, just to gather the variables)
*/
struct heat_equation_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  enum { DIRICHLET_WITH_MULTIPLIERS = 0, DIRICHLET_WITH_PENALIZATION = 1};
  getfem::mesh mesh;        /* the mesh */
  getfem::mesh_im mim;      /* the integration methods. */
  getfem::mesh_fem mf_u;    /* the main mesh_fem, for the Heat_Equation solution */
  getfem::mesh_fem mf_rhs;  /* the mesh_fem for the right hand side(f(x),..) */

  scalar_type residual;        /* max residual for the iterative solvers     */
  size_type N, dirichlet_version;
  scalar_type dirichlet_coefficient; /* Penalization parameter.              */
  plain_vector U;

  scalar_type dt, T, theta;

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(void);
  void init(void);
  void compute_error();
  heat_equation_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void heat_equation_problem::init(void) {
  
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  
  bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  if (N>1) { M(0,1) = PARAM.real_value("INCLINE") * PARAM.real_value("LY"); }

  /* scale the unit mesh to [LX,LY,..] and incline it */
  mesh.transformation(M);

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");
  dt = PARAM.real_value("DT", "Time step");
  T = PARAM.real_value("T", "final time");
  theta = PARAM.real_value("THETA", "Theta method parameter");
  sol_c = PARAM.real_value("C", "Diffusion coefficient");
  residual = PARAM.real_value("RESIDUAL");
  dirichlet_version = PARAM.int_value("DIRICHLET_VERSION",
				      "Type of Dirichlet contion");
  if (dirichlet_version == 1)
    dirichlet_coefficient = PARAM.real_value("DIRICHLET_COEFFICIENT",
					     "Penalization coefficient for "
					     "Dirichlet condition");
  if (residual == 0.) residual = 1e-10;
  sol_K.resize(N);
  for (size_type j = 0; j < N; j++)
    sol_K[j] = ((j & 1) == 0) ? FT : -FT;
  sol_t = 0.;

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    } else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}

bool heat_equation_problem::solve(void) {

  dal::bit_vector transient_bricks;

  getfem::model model;

  // Main unknown of the problem
  model.add_fem_variable("u", mf_u, 2);

  // Laplacian term on u.
  model.add_initialized_scalar_data("c", sol_c);
  transient_bricks.add(getfem::add_generic_elliptic_brick(model, mim,
							  "u", "c"));

  // Volumic source term.
  std::vector<scalar_type> F(mf_rhs.nb_dof());
  getfem::interpolation_function(mf_rhs, F, sol_f);
  model.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(model, mim, "u", "VolumicData");

  // Neumann condition.
  gmm::resize(F, mf_rhs.nb_dof()*N);
  getfem::interpolation_function(mf_rhs, F, sol_grad, NEUMANN_BOUNDARY_NUM);
  // The two version of the data make only a difference for midpoint scheme
  model.add_fem_data("NeumannData", mf_rhs, bgeot::dim_type(N), 2);
  gmm::copy(F, model.set_real_variable("NeumannData", 0));
  gmm::copy(F, model.set_real_variable("NeumannData", 1));
  transient_bricks.add(getfem::add_normal_source_term_brick
		       (model, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM));

  // Dirichlet condition.
  gmm::resize(F, mf_rhs.nb_dof());
  getfem::interpolation_function(mf_rhs, F, sol_u);
  model.add_initialized_fem_data("DirichletData", mf_rhs, F);

  if (dirichlet_version == DIRICHLET_WITH_MULTIPLIERS)
    getfem::add_Dirichlet_condition_with_multipliers
      (model, mim, "u", mf_u,
       DIRICHLET_BOUNDARY_NUM, "DirichletData");
  else
    getfem::add_Dirichlet_condition_with_penalization
      (model, mim, "u", dirichlet_coefficient,
       DIRICHLET_BOUNDARY_NUM, "DirichletData");

  // transient part.
  model.add_initialized_scalar_data("dt", dt);
  getfem::add_basic_d_on_dt_brick(model, mim, "u", "dt");
#if 1
  model.add_initialized_scalar_data("theta", theta);
  getfem::add_theta_method_dispatcher(model, transient_bricks, "theta");
#else
  // You can also test the midpoint scheme (but not really different from
  // crank-Nicholson scheme for linear problems).
  getfem::add_midpoint_dispatcher(model, transient_bricks);
#endif

  gmm::iteration iter(residual, 0, 40000);

  model.first_iter();

  // Null initial value for the temperature. Can be modified.
  gmm::resize(U, mf_u.nb_dof());
  gmm::clear(U);
  gmm::copy(U, model.set_real_variable("u", 1));

  for (scalar_type t = 0; t < T; t += dt) {
    sol_t = t+dt;
    
    gmm::resize(F, mf_rhs.nb_dof()*N);
    getfem::interpolation_function(mf_rhs, F, sol_grad, NEUMANN_BOUNDARY_NUM);
    gmm::copy(F, model.set_real_variable("NeumannData"));

    gmm::resize(F, mf_rhs.nb_dof());
    getfem::interpolation_function(mf_rhs, F, sol_u);
    gmm::copy(F, model.set_real_variable("DirichletData"));

    cout << "solving for t = " << sol_t << endl;
    iter.init();
    getfem::standard_solve(model, iter);
    gmm::copy(model.real_variable("u"), U);
    if (PARAM.int_value("EXPORT_SOLUTION") != 0) {
      char s[100]; sprintf(s, "step%d", int(t/dt)+1);
      gmm::vecsave(datafilename + s + ".U", U);
    }
    model.next_iter();
  }

  return (iter.converged());
}

/* compute the error with respect to the exact solution */
void heat_equation_problem::compute_error() {
  plain_vector V(mf_rhs.nb_basic_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_basic_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_basic_dof(i));
  cout.precision(16);
  cout << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V) << endl
       << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V) << endl
       << "Linfty error = " << gmm::vect_norminf(V) << endl;     
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  heat_equation_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  if (!p.solve()) GMM_ASSERT1(false, "Solve procedure has failed");
  if (p.PARAM.int_value("EXPORT_SOLUTION") != 0)
    p.mf_u.write_to_file(p.datafilename + ".mf", true);
  p.compute_error();

  return 0; 
}
