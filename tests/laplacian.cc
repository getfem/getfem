// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard, Julien Pommier.
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

/**@file laplacian.cc
   @brief Laplacian (Poisson) problem.
 
   The laplace equation is solved on a regular mesh of the unit
   square, and is compared to the analytical solution.

   This program is used to check that getfem++ is working. This is
   also a good example of use of Getfem++. This program  does not use the
   model bricks intentionally in order to serve as an exemple of solving
   a pde directly with the assembly procedures.
*/

#include "getfem/getfem_assembling.h" /* assembly methods (and comp. of norms) */
#include "getfem/getfem_export.h"   /* export functions (save solutions in a file) */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_superlu.h"
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

/* Definitions for the exact solution of the Laplacian problem,
 *  i.e. Delta(sol_u) + sol_f = 0
 */

base_small_vector sol_K; /* a coefficient */
/* exact solution */
scalar_type sol_u(const base_node &x) { return sin(gmm::vect_sp(sol_K, x)); }
/* righ hand side */
scalar_type sol_f(const base_node &x)
{ return gmm::vect_sp(sol_K, sol_K) * sin(gmm::vect_sp(sol_K, x)); }
/* gradient of the exact solution */
base_small_vector sol_grad(const base_node &x)
{ return sol_K * cos(gmm::vect_sp(sol_K, x)); }

/*
  structure for the Laplacian problem
*/
struct laplacian_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;        /* the mesh */
  getfem::mesh_im mim;      /* the integration methods. */
  getfem::mesh_fem mf_u;    /* the main mesh_fem, for the Laplacian solution */
  getfem::mesh_fem mf_rhs;  /* the mesh_fem for the right hand side(f(x),..) */
  getfem::mesh_fem mf_coef; /* the mesh_fem to represent pde coefficients    */

  scalar_type residual;        /* max residual for the iterative solvers */
  size_type est_degree, N;
  bool gen_dirichlet;

  sparse_matrix_type SM;     /* stiffness matrix.                           */
  std::vector<scalar_type> U, B;      /* main unknown, and right hand side  */

  std::vector<scalar_type> Ud; /* reduced sol. for gen. Dirichlet condition. */
  col_sparse_matrix_type NS; /* Dirichlet NullSpace 
			      * (used if gen_dirichlet is true)
			      */
  std::string datafilename;
  bgeot::md_param PARAM;

  void assembly(void);
  bool solve(void);
  void init(void);
  void compute_error();
  laplacian_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh),
			    mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void laplacian_problem::init(void) {
  
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
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  sol_K.resize(N);
  for (size_type j = 0; j < N; j++)
    sol_K[j] = ((j & 1) == 0) ? FT : -FT;

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  est_degree = pf_u->estimated_degree();
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
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  gen_dirichlet = PARAM.int_value("GENERIC_DIRICHLET");

  if (!pf_u->is_lagrange() && !gen_dirichlet)
    GMM_WARNING2("With non lagrange fem prefer the generic "
		 "Dirichlet condition option");

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

void laplacian_problem::assembly(void) {
  size_type nb_dof = mf_u.nb_dof();
  size_type nb_dof_rhs = mf_rhs.nb_dof();

  gmm::resize(B, nb_dof); gmm::clear(B);
  gmm::resize(U, nb_dof); gmm::clear(U); 
  gmm::resize(SM, nb_dof, nb_dof); gmm::clear(SM);
  
  cout << "Number of dof : " << nb_dof << endl;
  cout << "Assembling stiffness matrix" << endl;
  getfem::asm_stiffness_matrix_for_laplacian(SM, mim, mf_u, mf_coef, 
     std::vector<scalar_type>(mf_coef.nb_dof(), 1.0));
  
  cout << "Assembling source term" << endl;
  std::vector<scalar_type> F(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  getfem::asm_source_term(B, mim, mf_u, mf_rhs, F);
  
  cout << "Assembling Neumann condition" << endl;
  gmm::resize(F, nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_grad);
  getfem::asm_normal_source_term(B, mim, mf_u, mf_rhs, F,
				 NEUMANN_BOUNDARY_NUM);

  cout << "take Dirichlet condition into account" << endl;  
  if (!gen_dirichlet) {    
    std::vector<scalar_type> D(nb_dof);
    getfem::interpolation_function(mf_u, D, sol_u);
    getfem::assembling_Dirichlet_condition(SM, B, mf_u, 
					   DIRICHLET_BOUNDARY_NUM, D);
  } else {
    gmm::resize(F, nb_dof_rhs);
    getfem::interpolation_function(mf_rhs, F, sol_u);
    
    gmm::resize(Ud, nb_dof);
    gmm::resize(NS, nb_dof, nb_dof);
    col_sparse_matrix_type H(nb_dof_rhs, nb_dof);
    std::vector<scalar_type> R(nb_dof_rhs);
    std::vector<scalar_type> RHaux(nb_dof);

    /* build H and R such that U mush satisfy H*U = R */
    getfem::asm_dirichlet_constraints(H, R, mim, mf_u, mf_rhs,
      mf_rhs, F, DIRICHLET_BOUNDARY_NUM);

    gmm::clean(H, 1e-12);
//     cout << "H = " << H << endl;
//     cout << "R = " << R << endl;
    int nbcols = getfem::Dirichlet_nullspace(H, NS, R, Ud);
    // cout << "Number of irreductible unknowns : " << nbcols << endl;
    gmm::resize(NS, gmm::mat_ncols(H),nbcols);

    gmm::mult(SM, Ud, gmm::scaled(B, -1.0), RHaux);
    gmm::resize(B, nbcols);
    gmm::resize(U, nbcols);
    gmm::mult(gmm::transposed(NS), gmm::scaled(RHaux, -1.0), B);
    sparse_matrix_type SMaux(nbcols, nb_dof);
    gmm::mult(gmm::transposed(NS), SM, SMaux);
    gmm::resize(SM, nbcols, nbcols);
    /* NSaux = NS, but is stored by rows instead of by columns */
    sparse_matrix_type NSaux(nb_dof, nbcols); gmm::copy(NS, NSaux);
    gmm::mult(SMaux, NSaux, SM);
  }
}

bool laplacian_problem::solve(void) {
  cout << "Compute preconditionner\n";
  gmm::iteration iter(residual, 1, 40000);
  double time = gmm::uclock_sec();
  if (1) {
    gmm::identity_matrix P;
    // gmm::diagonal_precond<sparse_matrix_type> P(SM);
    // gmm::mr_approx_inverse_precond<sparse_matrix_type> P(SM, 10, 10E-17);
    // gmm::ildlt_precond<sparse_matrix_type> P(SM);
    // gmm::ildltt_precond<sparse_matrix_type> P(SM, 50, 1E-9);
    // gmm::ilut_precond<sparse_matrix_type> P(SM, 50, 1E-9);
    // gmm::ilutp_precond<sparse_matrix_type> P(SM, 50, 1E-9);
    // gmm::ilu_precond<sparse_matrix_type> P(SM);
    cout << "Time to compute preconditionner : "
	 << gmm::uclock_sec() - time << " seconds\n";

  
    //gmm::HarwellBoeing_IO::write("SM", SM);

    gmm::cg(SM, U, B, P, iter);
    // gmm::gmres(SM, U, B, P, 50, iter);
  } else {
    double rcond; 
    gmm::SuperLU_solve(SM, U, B, rcond); 
    cout << "cond = " << 1/rcond << "\n";
  }
  
  cout << "Total time to solve : "
       << gmm::uclock_sec() - time << " seconds\n";

  if (gen_dirichlet) {
    std::vector<scalar_type> Uaux(mf_u.nb_dof());
    gmm::mult(NS, U, Ud, Uaux);
    gmm::resize(U, mf_u.nb_dof());
    gmm::copy(Uaux, U);
  }

  return (iter.converged());
}

/* compute the error with respect to the exact solution */
void laplacian_problem::compute_error() {
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
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    laplacian_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.assembly();
    if (!p.solve()) GMM_ASSERT1(false, "Solve procedure has failed");
    p.compute_error();
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
