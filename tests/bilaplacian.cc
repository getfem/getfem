// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2006-2006 Yves Renard, Julien Pommier.                    */
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

/**
   @file bilaplacian.cc
   @brief Bilaplacian problem. A dummy
   bilaplacian problem is solved on a regular mesh, and is compared to
   the analytical solution.

   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.

   @see laplacian.cc
*/

#include <getfem_config.h>
#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_orderfourpdes.h>
#include <gmm.h>
#include <getfem_superlu.h>
#include <numeric>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_modeling.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/


scalar_type sol_u(const base_node &x) {
  return sin(std::accumulate(x.begin(), x.end(), 0.0));
}

scalar_type sol_f(const base_node &x) {
  return gmm::sqr(x.size())*sin(std::accumulate(x.begin(), x.end(), 0.0));
}

base_small_vector neumann_val(const base_node &x) {
  base_small_vector res(x.size());
  std::fill(res.begin(), res.end(),
	    sin(std::accumulate(x.begin(), x.end(), 0.0))* x.size());
  return res;
}

/*
  structure for the bilaplacian problem
*/
struct bilaplacian_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_im mim;      /* the integration methods.                     */
  getfem::mesh_fem mf_u;    /* main mesh_fem, for the bilaplacian solution  */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */

  scalar_type residual;    /* max residual for the iterative solvers        */
  int dirichlet_version;

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  bilaplacian_problem(void) : mim(mesh),mf_u(mesh), mf_mult(mesh),
			      mf_rhs(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void bilaplacian_problem::init(void) {
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
  size_type N = pgt->dim();
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

  dirichlet_version=PARAM.int_value("DIRICHLET_VERSION","Dirichlet version");
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);

  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0)
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
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
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) >= 1.0E-7) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    } else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}

/* compute the error with respect to the exact solution */
void bilaplacian_problem::compute_error(plain_vector &U) {
  std::vector<scalar_type> V(mf_rhs.nb_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_dof(i));
  cout.precision(16);
  cout << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V) << endl
       << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V) << endl
       << "Linfty error = " << gmm::vect_norminf(V) << endl;
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/


bool bilaplacian_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u);

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = sol_f(mf_rhs.point_of_dof(i));
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(BIL, mf_rhs, F);

  // Defining the Neumann condition right hand side.
  base_small_vector un(N);

  for (getfem::mr_visitor i(mesh.region(NEUMANN_BOUNDARY_NUM));
       !i.finished(); ++i) {
    size_type cv = i.cv(), f = i.f();
    getfem::pfem pf = mf_rhs.fem_of_element(cv);
    for (size_type l = 0; l< pf->structure(cv)->nb_points_of_face(f); ++l) {
      size_type n = pf->structure(cv)->ind_points_of_face(f)[l];
      un = mesh.normal_of_face_of_convex(cv, f, pf->node_of_dof(cv, n));
      un /= gmm::vect_norm2(un);
      size_type dof = mf_rhs.ind_dof_of_element(cv)[n];
      F[dof] = gmm::vect_sp(neumann_val(mf_rhs.point_of_dof(dof)), un);
    }
  }

  // Neumann condition brick.
  getfem::mdbrick_source_term<>
    NEUMANN(VOL_F, mf_rhs, F,NEUMANN_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = sol_u(mf_rhs.point_of_dof(i));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> final_model(NEUMANN, DIRICHLET_BOUNDARY_NUM,
					  mf_mult);
  final_model.set_constraints_type(getfem::constraints_type(dirichlet_version));
  final_model.rhs().set(mf_rhs, F);

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(BIL.get_solution(MS), U);
  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GETFEM_MPI_INIT(argc, argv); // For parallelized version

  DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    bilaplacian_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");

    p.compute_error(U);

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "bilaplacian_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	   << "mayavi -d " << p.datafilename
	   << ".vtk -f BandedSurfaceMap -m Outline\n";
    }
  }
  DAL_STANDARD_CATCH_ERROR;

  GETFEM_MPI_FINALIZE;

  return 0; 
}
