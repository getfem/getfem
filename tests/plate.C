/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Michel Salaün.                     */
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
 * Linear Elastostatic plate problem.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_linearized_plates.h>
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <gmm.h>

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
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the elastostatic problem
*/
struct plate_problem {

  enum { SIMPLY_FIXED_BOUNDARY_NUM = 0 };
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_fem mf_ut;
  getfem::mesh_fem mf_u3;
  getfem::mesh_fem mf_theta;
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type epsilon;       /* thickness of the plate.                      */
  scalar_type pressure;
  scalar_type residu;        /* max residu for the iterative solvers         */
  scalar_type LX;
  bool mixed, symmetrized;

  std::string datafilename;
  ftool::md_param PARAM;

  base_small_vector theta_exact(base_node P);
  scalar_type u3_exact(base_node P);

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  plate_problem(void) : mf_ut(mesh), mf_u3(mesh), mf_theta(mesh),
			       mf_rhs(mesh), mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void plate_problem::init(void) {
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE_UT  = PARAM.string_value("FEM_TYPE_UT","FEM name");
  const char *FEM_TYPE_U3  = PARAM.string_value("FEM_TYPE_U3","FEM name");
  const char *FEM_TYPE_THETA = PARAM.string_value("FEM_TYPE_THETA","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  const char *INTEGRATION_CT = PARAM.string_value("INTEGRATION_CT",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE_UT="  << FEM_TYPE_UT << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  cout << "INTEGRATION_CT=" << INTEGRATION_CT << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  if (N != 2)
    DAL_THROW(getfem::failure_error, "For a plate problem, N should be 2");
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  
  LX =  PARAM.real_value("LX");
  bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  mesh.transformation(M);

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;
  mixed = (PARAM.int_value("MIXED", "Mixed version ?") != 0);
  symmetrized = (PARAM.int_value("SYMMETRIZED",
				 "Mixed symmetrized version ?") != 0);

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  pressure = PARAM.real_value("PRESSURE",
			      "pressure on the top surface of the plate.");
  epsilon = PARAM.real_value("EPSILON", "thickness of the plate");
  mf_ut.set_qdim(N);
  mf_theta.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_ut = getfem::fem_descriptor(FEM_TYPE_UT);
  getfem::pfem pf_u3 = getfem::fem_descriptor(FEM_TYPE_U3);
  getfem::pfem pf_theta = getfem::fem_descriptor(FEM_TYPE_THETA);
  
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method ppi_ct = 
    getfem::int_method_descriptor(INTEGRATION_CT);

  mf_ut.set_finite_element(mesh.convex_index(), pf_ut, ppi);
  mf_u3.set_finite_element(mesh.convex_index(), pf_u3, ppi_ct);
  mf_theta.set_finite_element(mesh.convex_index(), pf_theta, ppi);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_ut->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_ut, ppi);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name), ppi);
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0), ppi);

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
    if (dal::abs(un[1]) <= 1.0E-7) { // new Neumann face
      mf_ut.add_boundary_elt(SIMPLY_FIXED_BOUNDARY_NUM, it->cv, it->f);
      mf_theta.add_boundary_elt(SIMPLY_FIXED_BOUNDARY_NUM, it->cv, it->f);
      mf_u3.add_boundary_elt(SIMPLY_FIXED_BOUNDARY_NUM, it->cv, it->f);
    }
  }
}

base_small_vector plate_problem::theta_exact(base_node P) {
  base_small_vector theta(2);
  theta[0] = (-pressure / (32. * mu * epsilon * epsilon * epsilon))
    * (4. * pow(P[0] - LX * .5, 3) - 3 * LX * LX * (P[0] - LX * .5));
  theta[1] = 0;
  return theta;
}

scalar_type plate_problem::u3_exact(base_node P) {
  return (pressure / (32. * mu * epsilon * epsilon * epsilon))
    * P[0] * (P[0] - LX)
    * (dal::sqr(P[0] - LX * .5) -1.25*LX*LX-(mixed ? 0 : 8.*epsilon*epsilon));
}


/* compute the error with respect to the exact solution */
void plate_problem::compute_error(plain_vector &U) {
  cout.precision(16);

  std::vector<scalar_type> V(mf_rhs.nb_dof()*2);
  
  size_type i1 = mf_ut.nb_dof();
  size_type i2 = mf_u3.nb_dof();
  size_type i3 = mf_theta.nb_dof();
  
  getfem::interpolation(mf_ut, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(0, i1)), V);
  mf_rhs.set_qdim(2);
  scalar_type l2 = dal::sqr(getfem::asm_L2_norm(mf_rhs, V));
  scalar_type h1 = dal::sqr(getfem::asm_H1_norm(mf_rhs, V));
  scalar_type linf = gmm::vect_norminf(V);
  mf_rhs.set_qdim(1);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  getfem::interpolation(mf_theta, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(i1+i2, i3)), V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i) {
    gmm::add(gmm::scaled(theta_exact(mf_rhs.point_of_dof(i)), -1.0),
	     gmm::sub_vector(V, gmm::sub_interval(i*2, 2)));
  }
  mf_rhs.set_qdim(2);
  l2 += dal::sqr(getfem::asm_L2_norm(mf_rhs, V));
  h1 += dal::sqr(getfem::asm_H1_norm(mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));
  mf_rhs.set_qdim(1);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  gmm::resize(V, mf_rhs.nb_dof());
  getfem::interpolation(mf_u3, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(i1, i2)), V);

  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= u3_exact(mf_rhs.point_of_dof(i));

  l2 += dal::sqr(getfem::asm_L2_norm(mf_rhs, V));
  h1 += dal::sqr(getfem::asm_H1_norm(mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));

  cout.precision(16);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool plate_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  
  cout << "Number of dof for ut: " << mf_ut.nb_dof() << endl;
  cout << "Number of dof for u3: " << mf_u3.nb_dof() << endl;
  cout << "Number of dof for theta: " << mf_theta.nb_dof() << endl;

  getfem::mdbrick_abstract<> *ELAS;

  // Linearized plate brick.
  getfem::mdbrick_isotropic_linearized_plate<>
    ELAS1(mf_ut, mf_u3, mf_theta, mf_coef, lambda, mu, epsilon);

  getfem::mdbrick_mixed_isotropic_linearized_plate<>
    ELAS2(mf_ut, mf_u3, mf_theta, mf_coef, lambda, mu, epsilon, symmetrized);

  if (mixed) ELAS = &ELAS2; else ELAS = &ELAS1;

  // Defining the surface source term.
  plain_vector F(nb_dof_rhs * 3);
  for (size_type i = 0; i < nb_dof_rhs; ++i) F[3*i+2] = pressure;
  getfem::mdbrick_plate_source_term<> VOL_F(*ELAS, mf_rhs, F);
  
  getfem::mdbrick_plate_simple_support<> SIMPLE
    (VOL_F, mf_rhs, SIMPLY_FIXED_BOUNDARY_NUM, 0, 1);

  getfem::mdbrick_plate_closing<> final_model(SIMPLE, 0, 1);
  

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residu, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);


  gmm::resize(U, mf_u3.nb_dof());
  if (mixed) gmm::copy(ELAS2.get_u3(MS), U);
  else gmm::copy(ELAS1.get_u3(MS), U);
  
  if (PARAM.int_value("VTK_EXPORT")) {
    cout << "export to " << datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(datafilename + ".vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(mf_u3); 
    exp.write_point_data(mf_u3, U, "plate_normal_displacement");
    cout << "export done, you can view the data file with (for example)\n"
       "mayavi -d " << datafilename << ".vtk -f "
       "WarpScalar -m BandedSurfaceMap -m Outline\n";
//     cout << "export done, you can view the data file with (for example)\n"
//       "mayavi -d " << datafilename << ".vtk -f ExtractVectorNorm -f "
//       "WarpVector -m BandedSurfaceMap -m Outline\n";
  }

  // Solution extraction
  gmm::resize(U, ELAS->nb_dof());
  if (mixed) gmm::copy(ELAS2.get_solution(MS), U);
  else gmm::copy(ELAS1.get_solution(MS), U);

  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {    
    plate_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.pressure *= p.epsilon * p.epsilon * p.epsilon;
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U;
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");
    p.compute_error(U);
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
