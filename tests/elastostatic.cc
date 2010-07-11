// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard, Julien Pommier.
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
//===========================================================================

/**
   @file elastostatic.cc
   @brief Linear Elastostatic problem. A dummy linear
   elastotatic problem is solved on a regular mesh, and is compared to
   the analytical solution.

   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.

   @see laplacian.cc
   @see nonlinear_elastostatic.cc
*/

#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_error_estimate.h"
#include "getfem/getfem_import.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

gmm::row_matrix<base_small_vector> sol_K;
static scalar_type sol_lambda, sol_mu, alph = 0.3;
int sol_sing;

base_small_vector sol_u(const base_node &x) {
  int N = x.size(); base_small_vector res(N);
  switch(sol_sing) {
    case 0 :
      for (int i = 0; i < N; ++i)
	res[i] = alph * sin(gmm::vect_sp(sol_K.row(i), x));
      break;
    case 1 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	scalar_type a = gmm::vect_norm2(y);
	for (int i = 0; i < N; ++i) res[i] = a;
	break;
      }
    case 2 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	scalar_type a = gmm::sqrt(gmm::vect_norm2(y));
	for (int i = 0; i < N; ++i) res[i] = a;
	break;
      }
  }
  return res;
}

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  switch (sol_sing) {
    case 0 :
      for (int i = 0; i < N; i++) {
	res[i] = alph * ( sol_mu * gmm::vect_sp(sol_K.row(i), sol_K.row(i)) )
	  * sin(gmm::vect_sp(sol_K.row(i), x));
	for (int j = 0; j < N; j++)
	  res[i] += alph * ( (sol_lambda + sol_mu) * sol_K(j,j) * sol_K(j,i))
	    * sin(gmm::vect_sp(sol_K.row(j), x));
      }
      break;
    case 1 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	scalar_type r = gmm::vect_norm2(y) + 1e-100;
	scalar_type r2 = r*r;
	scalar_type tr(0); tr = std::accumulate(y.begin(), y.end(), tr);

	for (int i = 0; i < N; i++)
	  res[i] = sol_lambda * (y[i]*tr / r2 - 1.0) / r
	    + sol_mu * (y[i]*tr/r2 - N) / r;
      }
      break;
    case 2 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	
	scalar_type r = gmm::vect_norm2(y) + 1e-100;
	scalar_type a = gmm::sqrt(1.0/r); 
	scalar_type b = a*a*a, c = b*b*a; 
	scalar_type tr(0); tr = std::accumulate(y.begin(), y.end(), tr);
	for (int i = 0; i < N; ++i)
	  res[i] -= b * (sol_lambda + sol_mu * (N+1-3.0/2.0)) * 0.5
	    - c * 3.0 * (sol_lambda + sol_mu) * (y[i]*tr) / 4.0;
      }
      break;
  }
  return res;
}

base_matrix sol_sigma(const base_node &x) {
  int N = x.size();
  base_matrix res(N,N);
  switch (sol_sing) {
    case 0 :
      for (int i = 0; i < N; i++)
	for (int j = 0; j <= i; j++) {
	  res(j,i) = res(i,j) = alph * sol_mu *
	    ( sol_K(i,j) * cos(gmm::vect_sp(sol_K.row(i), x))
	      +  sol_K(j,i) * cos(gmm::vect_sp(sol_K.row(j), x))
	      );
	  if (i == j)
	    for (int k = 0; k < N; k++)
	      res(i,j) += alph * sol_lambda * sol_K(k,k)
		* cos(gmm::vect_sp(sol_K.row(k), x));
	}
      break;
    case 1 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	scalar_type r = gmm::vect_norm2(y) + 1e-100;
	scalar_type tr(0); tr = std::accumulate(y.begin(), y.end(), tr);
	for (int i = 0; i < N; i++) {
	  res(i, i) += sol_lambda * tr / r;
	  for (int j = 0; j < N; j++)
	    res(i, j) += sol_mu * (y[i] + y[j]) / r;
	}
      }
 
      break;
    case 2 :
      {
	base_small_vector trans(x.size());
	gmm::fill(trans,  M_PI / 10.0);
	base_node y = x - trans;
	scalar_type r = gmm::vect_norm2(y) + 1e-100;
	scalar_type a = gmm::sqrt(1.0/r); 
	scalar_type b = a*a*a;
	scalar_type tr(0); tr = std::accumulate(y.begin(), y.end(), tr);
	for (int i = 0; i < N; i++) {
	  res(i, i) += sol_lambda * tr * b * 0.5;
	  for (int j = 0; j < N; j++)
	    res(i, j) += sol_mu * b * (y[i] + y[j]) * 0.5;
	}
      }
  }

  return res;
}

/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_im mim;       /* the integration methods.                     */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_mult;  /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure for mixed form     */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */

  scalar_type residual;       /* max residual for iterative solvers          */
  bool mixed_pressure, refine;
  getfem::constraints_type dirichlet_version;

  std::string datafilename;
  bgeot::md_param PARAM;
  
  void select_boundaries(void);
  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  elastostatic_problem(void) : mim(mesh),mf_u(mesh), mf_mult(mesh),
			       mf_rhs(mesh),mf_p(mesh) {}
};

/* Selects the boundaries */

void elastostatic_problem::select_boundaries(void) {
  size_type N = mesh.dim();
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 0.5) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    } else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}


/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void elastostatic_problem::init(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");

  cout << "MESH_FILE=" << MESH_FILE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

#if GETFEM_PARA_LEVEL > 1
  double t_init=MPI_Wtime();
#endif

  size_type NX = PARAM.int_value("NX");
  size_type N = PARAM.int_value("N");
  std::stringstream filename; filename << MESH_FILE;
  if ((MESH_FILE.compare(0,11,"structured:") == 0) && NX > 0) {
    filename << ";NSUBDIV=[" << NX;
    for (size_type i = 1; i < N; ++i) filename << "," << NX;
    filename << "];";
  }
  getfem::import_mesh(filename.str(), mesh);
  
  GMM_ASSERT1(N == mesh.dim(), "The mesh has not the right dimension");

#if GETFEM_PARA_LEVEL > 1
  cout<<"temps creation maillage "<< MPI_Wtime()-t_init<<endl;
#endif

  dirichlet_version
    = getfem::constraints_type(PARAM.int_value("DIRICHLET_VERSION",
					       "Dirichlet version"));
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  gmm::resize(sol_K, N, N);
  for (size_type i = 0; i < N; i++)
    for (size_type j = 0; j < N; j++)
      sol_K(i,j) = (i == j) ? FT : -FT;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  sol_sing = int(PARAM.int_value("SOL_SING", "Optional singular solution"));
  refine = (PARAM.int_value("REFINE", "Optional refinement") != 0);
  sol_lambda = lambda; sol_mu = mu;
  mf_u.set_qdim(dim_type(N));
  mf_mult.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(ppi);
  mf_u.set_finite_element(pf_u);

  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0)
    mf_mult.set_finite_element(pf_u);
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(getfem::fem_descriptor(dirichlet_fem_name));
  }

  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  if (mixed_pressure) {
    std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
    mf_p.set_finite_element(getfem::fem_descriptor(FEM_TYPE_P));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(pf_u);
  } else {
    mf_rhs.set_finite_element(getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  select_boundaries();

#if GETFEM_PARA_LEVEL > 1
  
  
  t_init = MPI_Wtime();
  
  mf_u.nb_dof(); mf_rhs.nb_dof(); mf_mult.nb_dof();

  cout<<"enumerate dof time "<< MPI_Wtime()-t_init<<endl;
#else
  double t_init = gmm::uclock_sec();
  mf_u.nb_dof(); mf_rhs.nb_dof(); mf_mult.nb_dof();
  cout << "enumerate dof time " << gmm::uclock_sec() - t_init << endl;
#endif
}

/* compute the error with respect to the exact solution */
void elastostatic_problem::compute_error(plain_vector &U) {
  size_type N = mesh.dim();
  std::vector<scalar_type> V(mf_rhs.nb_basic_dof()*N);
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_basic_dof(); ++i) {
    gmm::add(gmm::scaled(sol_u(mf_rhs.point_of_basic_dof(i)), -1.0),
	     gmm::sub_vector(V, gmm::sub_interval(i*N, N)));
  }

  

  cout.precision(16);
  mf_rhs.set_qdim(dim_type(N));
  scalar_type l2 = getfem::asm_L2_norm(mim, mf_rhs, V);
  scalar_type h1 = getfem::asm_H1_norm(mim, mf_rhs, V);

  if (getfem::MPI_IS_MASTER())
    cout << "L2 error = " << l2 << endl
	 << "H1 error = " << h1 << endl
	 << "Linfty error = " << gmm::vect_norminf(V) << endl;
  
  getfem::vtk_export exp(datafilename + "_err.vtk",
			 PARAM.int_value("VTK_EXPORT")==1);
  exp.exporting(mf_rhs); 
  exp.write_point_data(mf_rhs, V, "elastostatic_displacement");

  mf_rhs.set_qdim(1);
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/


bool elastostatic_problem::solve(plain_vector &U) {

  size_type N = mesh.dim();

  if (mixed_pressure) cout << "Number of dof for P: " << mf_p.nb_dof() << endl;
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  getfem::model model;

  // Main unknown of the problem.
  model.add_fem_variable("u", mf_u);

  // Linearized elasticity brick.
  model.add_initialized_scalar_data("lambda", mixed_pressure ? 0.0 : lambda);
  model.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim, "u", "lambda", "mu");

  // Linearized incompressibility condition brick.
  if (mixed_pressure) {
    model.add_initialized_scalar_data("incomp_coeff", 1.0/lambda);
    model.add_fem_variable("p", mf_p); // Adding the pressure as a variable
    add_linear_incompressibility
      (model, mim, "u", "p", size_type(-1), "incomp_coeff");
  }

  // Volumic source term.
  std::vector<scalar_type> F(mf_rhs.nb_dof()*N);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  model.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(model, mim, "u", "VolumicData");
    
  // Neumann condition.
  gmm::resize(F, mf_rhs.nb_dof()*N*N);
  getfem::interpolation_function(mf_rhs, F, sol_sigma, NEUMANN_BOUNDARY_NUM);
  model.add_initialized_fem_data("NeumannData", mf_rhs, F);
  getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM);

  // Dirichlet condition.
  gmm::resize(F, mf_rhs.nb_dof()*N);
  getfem::interpolation_function(mf_rhs, F, sol_u);
  model.add_initialized_fem_data("DirichletData", mf_rhs, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (model, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, "DirichletData");

  gmm::iteration iter(residual, 1, 40000);
#if GETFEM_PARA_LEVEL > 1
  double t_init=MPI_Wtime();
#endif
  dal::bit_vector cvref;
  
  do { // solve with optional refinement
    
    cout << "Total number of variables : " << model.nb_dof() << endl;

    // Defining the volumic source term.
    size_type nb_dof_rhs = mf_rhs.nb_dof();
    gmm::resize(F, nb_dof_rhs * N);
    getfem::interpolation_function(mf_rhs, F, sol_f);
    gmm::copy(F, model.set_real_variable("VolumicData"));

    // Defining the Neumann source term.
    gmm::resize(F, nb_dof_rhs * N * N);
    getfem::interpolation_function(mf_rhs, F, sol_sigma, NEUMANN_BOUNDARY_NUM);
    gmm::copy(F, model.set_real_variable("NeumannData"));

    // Defining the Dirichlet condition value.
    gmm::resize(F, nb_dof_rhs * N);
    getfem::interpolation_function(mf_rhs, F, sol_u, DIRICHLET_BOUNDARY_NUM);
    gmm::copy(F, model.set_real_variable("DirichletData"));
    
    iter.init();
    getfem::standard_solve(model, iter);
    gmm::resize(U, mf_u.nb_dof());
    gmm::copy(model.real_variable("u"), U);

    if (refine) {
      plain_vector ERR(mesh.convex_index().last_true()+1);
      getfem::error_estimate(mim, mf_u, U, ERR);
      
      cout << "max = " << gmm::vect_norminf(ERR) << endl;
      // scalar_type threshold = gmm::vect_norminf(ERR) * 0.7;
      scalar_type threshold = 0.0001, min_ = 1e18;
      cvref.clear();
      for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
	if (ERR[i] > threshold) cvref.add(i);
	min_ = std::min(min_, ERR[i]);
      }
      cout << "min = " << min_ << endl;
      cout << "Nb elt to be refined : " << cvref.card() << endl;
     
      mesh.Bank_refine(cvref);
    }

  } while (refine && cvref.card() > 0);

#if GETFEM_PARA_LEVEL > 1
    cout<<"temps standard solve "<< MPI_Wtime()-t_init<<endl;
#endif
 

  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GETFEM_MPI_INIT(argc, argv); // For parallelized version

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //  try {



    elastostatic_problem p;
    p.PARAM.read_command_line(argc, argv);
#if GETFEM_PARA_LEVEL > 1
    double t_ref=MPI_Wtime();
#endif
    p.init();
#if GETFEM_PARA_LEVEL > 1
    cout << "temps init "<< MPI_Wtime()-t_ref << endl;
#endif
    if (getfem::MPI_IS_MASTER())
      p.mesh.write_to_file(p.datafilename + ".mesh");

    plain_vector U;

#if GETFEM_PARA_LEVEL > 1
    t_ref=MPI_Wtime();
    cout<<"begining resol"<<endl;
#endif
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

#if GETFEM_PARA_LEVEL > 1
    cout << "temps Resol "<< MPI_Wtime()-t_ref << endl;
    t_ref = MPI_Wtime();
#endif
    p.compute_error(U);
#if GETFEM_PARA_LEVEL > 1
    cout << "temps error "<< MPI_Wtime()-t_ref << endl;
    t_ref = MPI_Wtime();
#endif

    // if (getfem::MPI_IS_MASTER()) { p.mesh.write_to_file("toto.mesh"); }

    if (p.PARAM.int_value("VTK_EXPORT") && getfem::MPI_IS_MASTER()) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi2 -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m Surface -m Outline\n";
    }

    //  } GMM_STANDARD_CATCH_ERROR;

  GETFEM_MPI_FINALIZE;



  return 0; 
}
