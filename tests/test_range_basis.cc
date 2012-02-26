/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard, Julien Pommier.
 
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

/**@file test_rang_basis.cc
   @brief A small test for gmm_range_basis function.
*/

#include "gmm/gmm.h"
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_regular_meshes.h"

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
  getfem::mesh_fem mf_mult; /* the mesh_fem to represent pde coefficients    */

  scalar_type residual;        /* max residual for the iterative solvers */
  size_type est_degree, N;
  bool gen_dirichlet;

  sparse_matrix_type SM;     /* stiffness matrix.                           */
  sparse_matrix_type B;     
  
  std::vector<scalar_type> U, L;      /* main unknown, and right hand side  */

  std::vector<scalar_type> Ud; /* reduced sol. for gen. Dirichlet condition. */
  col_sparse_matrix_type NS; /* Dirichlet NullSpace 
			      * (used if gen_dirichlet is true)
			      */
  std::string datafilename;
  bgeot::md_param PARAM;

  void assembly(void);
  void init(void);
  laplacian_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh),
			    mf_mult(mesh) {}
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
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
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
  
  std::string mult_fem_name = PARAM.string_value("MULT_FEM_TYPE",
						 "Mult fem type");
  mf_mult.set_finite_element(mesh.convex_index(),
			      getfem::fem_descriptor(mult_fem_name));

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
  size_type nb_dof_mult = mf_mult.nb_dof();

  gmm::resize(B, nb_dof_mult, nb_dof); gmm::clear(B);
  
  cout << "Number of dof : " << nb_dof << endl;
  cout << "Number of dof mult : " << nb_dof_mult << endl;
  cout << "Assembly of the mass matrix" << endl;
  getfem::asm_mass_matrix(B, mim, mf_mult, mf_u, NEUMANN_BOUNDARY_NUM);
  
  
  std::set<size_type> columns;

  double t = gmm::uclock_sec();
  cout << "depart range basis" << endl;
  gmm::range_basis(gmm::transposed(B), columns);
  cout << "Elaps time for range basis : " << gmm::uclock_sec() - t << endl;
  cout << "Rank of B : " << columns.size() << " null space dimension : "
       << nb_dof-columns.size() << endl;
  
  

  // compute the kernel on u.
  t = gmm::uclock_sec();
  NS.resize(nb_dof, nb_dof);
  cout << "depart Dirichlet_nullspace" << endl;
  plain_vector U0(nb_dof);
  col_sparse_matrix_type BB( nb_dof_mult, nb_dof);
  gmm::copy(B, BB);
  size_type nk
    = getfem::Dirichlet_nullspace(BB, NS, plain_vector(gmm::mat_nrows(B)), U0);
  cout << "Elaps time for Dirichlet_nullspace : " << gmm::uclock_sec() - t
       << endl;		  
  cout << "Null space dimension : " << nk << endl;

  GMM_ASSERT1(nk == nb_dof-columns.size(),
	      "Different results for the dimension of the null space");


  // the same test with complex numbers
  cout << "Repeated the test with complex numbers" << endl;
  gmm::row_matrix< gmm::rsvector< std::complex<double> > >
    B2(nb_dof_mult, nb_dof);
  getfem::asm_mass_matrix(gmm::imag_part(B2), mim, mf_mult, mf_u,
			  NEUMANN_BOUNDARY_NUM);
  // gmm::copy(B, gmm::imag_part(B2));
  t = gmm::uclock_sec();
  cout << "depart range basis" << endl;
  gmm::range_basis(gmm::transposed(B2), columns);
  cout << "Elaps time for range basis : " << gmm::uclock_sec() - t << endl;
  cout << "Rank of B2 : " << columns.size() << " null space dimension : "
       << nb_dof-columns.size() << endl;

  GMM_ASSERT1(nk == nb_dof-columns.size(),
	      "Different results for the dimension of the null space");
}


/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.


  std::set<size_type> columns;
  sparse_matrix_type B(0,10);
  gmm::range_basis(gmm::transposed(B), columns);
  GMM_ASSERT1(columns.size() == 0, "Wrong behavior of range_basis !");

  laplacian_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  p.assembly();
  
  return 0; 
}
