// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Yves Renard, Julien Pommier.                    */
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
 * Linear Elastostatic problem.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_superlu.h>
#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <gmm.h>

#ifdef GMM_USES_MPI
#include <mpi.h>
#include <mpi++.h>
#endif
/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
#ifndef GMM_USES_MPI
typedef getfem::modeling_standard_sparse_matrix Tsparse_matrix;
typedef getfem::modeling_standard_sparse_matrix Csparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;
#else
typedef gmm::mpi_distributed_matrix<getfem::modeling_standard_sparse_matrix> Tsparse_matrix;
typedef getfem::modeling_standard_sparse_matrix Csparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;
#endif

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

gmm::row_matrix<base_small_vector> sol_K;
static scalar_type sol_lambda, sol_mu, alph = 0.3;

base_small_vector sol_u(const base_node &x) {
  int N = x.size(); base_small_vector res(N);
  for (int i = 0; i < N; ++i)
    res[i] = alph * sin(gmm::vect_sp(sol_K.row(i), x));
  return res;
}

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  for (int i = 0; i < N; i++) {
    res[i] = alph * ( sol_mu * gmm::vect_sp(sol_K.row(i), sol_K.row(i)) )
                  * sin(gmm::vect_sp(sol_K.row(i), x));
    for (int j = 0; j < N; j++)
      res[i] += alph * ( (sol_lambda + sol_mu) * sol_K(j,j) * sol_K(j,i))
	          * sin(gmm::vect_sp(sol_K.row(j), x));
  }
  return res;
}

base_matrix sol_sigma(const base_node &x) {
  int N = x.size();
  base_matrix res(N,N);
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
  return res;
}

/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_im mim;       /* the integration methods.                     */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure for mixed form     */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */

  scalar_type residue;        /* max residue for the iterative solvers         */
  bool mixed_pressure;

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_p(mesh),
			       mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void elastostatic_problem::init(void) {
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
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


  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");
  residue = PARAM.real_value("RESIDUE"); if (residue == 0.) residue = 1e-10;
  gmm::resize(sol_K, N, N);
  for (size_type i = 0; i < N; i++)
    for (size_type j = 0; j < N; j++)
      sol_K(i,j) = (i == j) ? FT : -FT;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  sol_lambda = lambda; sol_mu = mu;
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  
  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  if (mixed_pressure) {
    const char *FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
    mf_p.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM_TYPE_P));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
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
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    base_node un = mesh.normal_of_face_of_convex(it->cv, it->f);
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) >= 1.0E-7) { // new Neumann face
      mesh.add_face_to_set(NEUMANN_BOUNDARY_NUM, it->cv, it->f);
    } else {
      mesh.add_face_to_set(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
    }
  }
}

/* compute the error with respect to the exact solution */
void elastostatic_problem::compute_error(plain_vector &U) {
  size_type N = mesh.dim();
  std::vector<scalar_type> V(mf_rhs.nb_dof()*N);
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i) {
    gmm::add(gmm::scaled(sol_u(mf_rhs.point_of_dof(i)), -1.0),
	     gmm::sub_vector(V, gmm::sub_interval(i*N, N)));
  }
  cout.precision(16);
  mf_rhs.set_qdim(N);
  cout << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V) << endl
       << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V) << endl
       << "Linfty error = " << gmm::vect_norminf(V) << endl;
  mf_rhs.set_qdim(1);
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

typedef getfem::model_state<Tsparse_matrix, Csparse_matrix, plain_vector> Model_State;

bool elastostatic_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  if (mixed_pressure) cout << "Number of dof for P: " << mf_p.nb_dof() << endl;
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<Model_State>
    ELAS(mim, mf_u, mf_coef, mixed_pressure ? 0.0 : lambda, mu);

  getfem::mdbrick_linear_incomp<Model_State> INCOMP(ELAS, mf_p, mf_coef, 1.0/lambda);

  getfem::mdbrick_abstract<Model_State> *pINCOMP;
  if (mixed_pressure) pINCOMP = &INCOMP; else pINCOMP = &ELAS;

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_f(mf_rhs.point_of_dof(i)),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<Model_State> VOL_F(*pINCOMP, mf_rhs, F);

  // Defining the Neumann condition right hand side.
  base_small_vector un(N), v(N);
  for (dal::bv_visitor cv(mesh.convexes_in_set(NEUMANN_BOUNDARY_NUM));
       !cv.finished(); ++cv) {
    getfem::pfem pf = mf_rhs.fem_of_element(cv);
    getfem::mesh_cvf_set::face_bitset fb = 
      mesh.faces_of_convex_in_set(NEUMANN_BOUNDARY_NUM, cv);
    for (unsigned f = 0; f < MAX_FACES_PER_CV; ++f) if (fb[f]) {
      for (size_type l = 0; l< pf->structure(cv)->nb_points_of_face(f); ++l) {
	size_type n = pf->structure(cv)->ind_points_of_face(f)[l];
	un = mesh.normal_of_face_of_convex(cv, f, pf->node_of_dof(cv, n));
	un /= gmm::vect_norm2(un);
	size_type dof = mf_rhs.ind_dof_of_element(cv)[n];
	gmm::mult(sol_sigma(mf_rhs.point_of_dof(dof)), un, v);
	gmm::copy(v, gmm::sub_vector(F, gmm::sub_interval(dof*N, N)));
      }
    }
  }

  // Neumann condition brick.
  getfem::mdbrick_source_term<Model_State> NEUMANN(VOL_F, mf_rhs, F,NEUMANN_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_u(mf_rhs.point_of_dof(i)), 
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<Model_State> final_model(NEUMANN, mf_rhs,
					  F, DIRICHLET_BOUNDARY_NUM);

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  Model_State MS(final_model);
  gmm::iteration iter(residue, 1, 40000);
#ifdef GMM_USES_MPI
    double t_init=MPI_Wtime();
#endif
  getfem::standard_solve(MS, final_model, iter);

#ifdef GMM_USES_MPI
    cout<<"temps standard solve "<< MPI_Wtime()-t_init<<endl;
#endif
#ifdef GMM_USES_MPI
    double t_ref=MPI_Wtime();
#endif
  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);
#ifdef GMM_USES_MPI
    cout<<"temps copy elas "<< MPI_Wtime()-t_ref<<endl;
#endif
  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
#ifdef GMM_USES_MPI
  int rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {
#ifdef GMM_USES_MPI
    static double t_resol = 0.0;
    double t_ref,t_final;
#endif

    elastostatic_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u.nb_dof());
#ifdef GMM_USES_MPI
    t_ref=MPI_Wtime();
    cout<<"begining resol"<<endl;
#endif
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");

#ifdef GMM_USES_MPI
    t_final=MPI_Wtime();
    t_resol += t_final-t_ref;
    cout<<"end resol"<<endl;
    cout<<"["<< rank <<"] temps Resol "<< t_final-t_ref << " t_tot = " << t_resol << endl;
#endif
    p.compute_error(U);

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline\n";
    }
  }
  DAL_STANDARD_CATCH_ERROR;
#ifdef GMM_USES_MPI
   MPI_Finalize();
#endif
  return 0; 
}
