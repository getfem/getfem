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
#include <getfem_fourth_order.h>
#include <getfem_model_solvers.h>
#include <gmm.h>
#include <getfem_superlu.h>
#include <getfem_derivatives.h>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
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

scalar_type FT = 0.0;

// scalar_type sol_u(const base_node &x)
// { return sin(FT*std::accumulate(x.begin(), x.end(), 0.0)); }

// scalar_type sol_lapl_u(const base_node &x)
// { return -FT*FT*sol_u(x) * x.size(); }

// scalar_type sol_f(const base_node &x)
// { return FT*FT*FT*FT*sol_u(x)*gmm::sqr(x.size()); }

// base_small_vector sol_du(const base_node &x) {
//   base_small_vector res(x.size());
//   std::fill(res.begin(), res.end(),
// 	    FT * cos(std::accumulate(x.begin(), x.end(), 0.0)));
//   return res;
// }

// base_small_vector neumann_val(const base_node &x)
// { return -FT*FT*FT*sol_du(x) * scalar_type(x.size()); }


scalar_type sol_u(const base_node &x)
{ return x[0]*x[1]; }

scalar_type sol_lapl_u(const base_node &)
{ return 0.0; }

scalar_type sol_f(const base_node &x) { return 0.0; }

base_small_vector sol_du(const base_node &x) {
  base_small_vector res(x.size()); res[0] = x[1]; res[1] = x[0]; 
  return res;
}

base_small_vector neumann_val(const base_node &x)
{ return base_small_vector(x.size());  }



/**************************************************************************/
/*  Structure for the bilaplacian problem.                                */
/**************************************************************************/

struct bilaplacian_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3};
  getfem::mesh mesh;        /* the mesh */
  getfem::mesh_im mim;      /* the integration methods.                     */
  getfem::mesh_fem mf_u;    /* main mesh_fem, for the bilaplacian solution  */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */

  scalar_type residual;     /* max residual for the iterative solvers       */
  getfem::constraints_type dirichlet_version;

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
  std::string MESH_FILE = PARAM.string_value("MESH_FILE");
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  size_type N;
  if (!MESH_FILE.empty()) {
    mesh.read_from_file(MESH_FILE);
    MESH_TYPE = bgeot::name_of_geometric_trans(mesh.trans_of_convex(mesh.convex_index().first_true()));
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    N = mesh.dim();
  } else {
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    
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
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
  }

  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  dirichlet_version = getfem::constraints_type(dv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  FT = PARAM.real_value("FT"); if (FT == 0.0) FT = 1.0;


  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
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
    if (0 & gmm::abs(un[N-1] - 1.0) <= 0.5) {
      mesh.region(FORCE_BOUNDARY_NUM).add(i.cv(), i.f());
      mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
    }
    else {
      mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
      if (gmm::abs(un[N-1] + 1.0) <= 0.5)
	mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
      else
	mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
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

  // Defining the prescribed momentum.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = sol_lapl_u(mf_rhs.point_of_dof(i));
  
  // Prescribed momentum on the boundary
  getfem::mdbrick_normal_derivative_source_term<>
    MOMENTUM(VOL_F, mf_rhs, F, MOMENTUM_BOUNDARY_NUM);


  // Defining the Neumann condition right hand side.
  base_small_vector un(N);
  gmm::clear(F);
  for (getfem::mr_visitor i(mesh.region(FORCE_BOUNDARY_NUM));
       !i.finished(); ++i) {
    size_type cv = i.cv(), f = i.f();
    getfem::pfem pf = mf_rhs.fem_of_element(cv);
    for (size_type l = 0; l< pf->structure(cv)->nb_points_of_face(f); ++l) {
      size_type n = pf->structure(cv)->ind_points_of_face(f)[l];
      un = mesh.normal_of_face_of_convex(cv, f, pf->node_of_dof(cv, n));
      un /= gmm::vect_norm2(un);
      size_type dof = mf_rhs.ind_dof_of_element(cv)[n];
      F[dof] = -gmm::vect_sp(neumann_val(mf_rhs.point_of_dof(dof)), un);
    }
  }

  // Neumann condition brick.
  getfem::mdbrick_source_term<>
    NEUMANN(MOMENTUM, mf_rhs, F, FORCE_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = sol_u(mf_rhs.point_of_dof(i));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NEUMANN, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  DIRICHLET.set_constraints_type(dirichlet_version);
  DIRICHLET.rhs().set(mf_rhs, F);

  // Defining the normal derivative Dirichlet condition value.
  gmm::clear(F);
  for (getfem::mr_visitor i(mesh.region(CLAMPED_BOUNDARY_NUM));
       !i.finished(); ++i) {
    size_type cv = i.cv(), f = i.f();
    getfem::pfem pf = mf_rhs.fem_of_element(cv);
    for (size_type l = 0; l< pf->structure(cv)->nb_points_of_face(f); ++l) {
      size_type n = pf->structure(cv)->ind_points_of_face(f)[l];
      un = mesh.normal_of_face_of_convex(cv, f, pf->node_of_dof(cv, n));
      un /= gmm::vect_norm2(un);
      size_type dof = mf_rhs.ind_dof_of_element(cv)[n];
      F[dof] = gmm::vect_sp(sol_du(mf_rhs.point_of_dof(dof)), un) * 0.1;
    }
  }
  
  // Normal derivative Dirichlet condition brick.
  getfem::mdbrick_normal_derivative_Dirichlet<> 
    final_model(DIRICHLET, CLAMPED_BOUNDARY_NUM, mf_mult);
  final_model.set_constraints_type(dirichlet_version);
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
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error, "Solve has failed");

    p.compute_error(U);

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "bilaplacian_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	   << "mayavi -d " << p.datafilename
	   << ".vtk -m BandedSurfaceMap -m Outline -f WarpScalar\n";
    }
  }
  DAL_STANDARD_CATCH_ERROR;

  GETFEM_MPI_FINALIZE;

  return 0; 
}
