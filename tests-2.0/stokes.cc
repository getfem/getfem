// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**
   @file stokes.cc
   @brief Solve the stokes problem (incompressible viscuous fluid).
   
   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.
   
   This file is almost identical to @link elastostatic.cc
   tests/elastostatic.cc@endlink, except than an
   getfem::mdbrick_linear_incomp brick is inserted.
*/

#include <getfem/getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem/getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_model_solvers.h>
#include <gmm/gmm.h>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
 * structure for the stokes problem
 */
struct stokes_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the velocity              */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure                    */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  scalar_type nu;            /* Lamé coefficients.                           */

  scalar_type residual;        /* max residual for the iterative solvers         */

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  stokes_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh),
			 mf_rhs(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void stokes_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
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
  
  /* scale the unit mesh to [LX,LY,..] and incline it */
   bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  mesh.transformation(M);

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

  nu = PARAM.real_value("NU", "Viscosité");
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION); 

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_p.set_finite_element(mesh.convex_index(),
			  getfem::fem_descriptor(FEM_TYPE_P));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_THROW(dal::failure_error, "You are using a non-lagrange FEM "
		<< data_fem_name << ". In that case you need to set "
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
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
    } else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(), it.f());
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

base_small_vector sol_f(const base_small_vector &P) {
  base_small_vector res(P.size());
  res[P.size()-1] = -1.0;
  return res;
}


bool stokes_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u, 0.0, nu);

  // Pressure term
  getfem::mdbrick_linear_incomp<> INCOMP(ELAS, mf_p);
  
  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_f(mf_rhs.point_of_dof(i)),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(INCOMP, mf_rhs, F);
  
  // Defining the Dirichlet condition value.
  gmm::clear(F);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> final_model(VOL_F, DIRICHLET_BOUNDARY_NUM);
  final_model.rhs().set(mf_rhs,F);

  // Generic solve.
  cout << "Number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);
  
  return (iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    stokes_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    // p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) GMM_THROW(dal::failure_error,"Solve has failed");
    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "stokes_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline\n";
    }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
