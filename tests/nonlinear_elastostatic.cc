/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
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
   @file nonlinear_elastostatic.cc
   @brief Nonlinear Elastostatic problem (large strain).

   A rubber bar is submitted to a large torsion.
   
   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "getfem/getfem_superlu.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
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
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure.                   */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type p1, p2, p3;    /* elastic coefficients.                        */

  scalar_type residual;        /* max residual for the iterative solvers         */

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh), mf_rhs(mesh), mf_coef(mesh) {}
};


/* Read parameters from the .param file, build the mesh, set finite element
   and integration methods and selects the boundaries.

   (this is boilerplate code, not very interesting)
 */
void elastostatic_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name for the pressure");
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
  nsubdiv[1] = PARAM.int_value("NY") ? PARAM.int_value("NY") : nsubdiv[0];
  if (N>2) nsubdiv[2] = PARAM.int_value("NZ") ? PARAM.int_value("NZ") : nsubdiv[0];
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
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

  p1 = PARAM.real_value("P1", "First Elastic coefficient");
  p2 = PARAM.real_value("P2", "Second Elastic coefficient");
  p3 = PARAM.real_value("P3", "Third Elastic coefficient");
  
  mf_u.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(ppi);
  mf_u.set_finite_element(pf_u);

  mf_p.set_finite_element(getfem::fem_descriptor(FEM_TYPE_P));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM"
		". In that case you need to set "
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
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { 
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
    } else if (gmm::abs(un[N-1] + 1.0) < 1.0E-7) {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool elastostatic_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  size_type law_num = PARAM.int_value("LAW");
  // Linearized elasticity brick.
  base_vector p(3); p[0] = p1; p[1] = p2; p[2] = p3;

  /* choose the material law */
  getfem::abstract_hyperelastic_law *pl1 = 0, *pl = 0;
  switch (law_num) {
    case 0:
    case 1: pl1 = new getfem::SaintVenant_Kirchhoff_hyperelastic_law(); break;
    case 2: pl1 = new getfem::Ciarlet_Geymonat_hyperelastic_law(); break;
    case 3: pl1 = new getfem::Mooney_Rivlin_hyperelastic_law(); break;
    default: GMM_ASSERT1(false, "no such law");
  }

  if (N == 2) {
    cout << "2D plane strain hyper-elasticity\n";
    pl = new getfem::plane_strain_hyperelastic_law(pl1);
  } else pl = pl1;

  p.resize(pl->nb_params());


  pl->test_derivatives(3, 5e-9, p);


  getfem::model model;

  // Main unknown of the problem (displacement).
  model.add_fem_variable("u", mf_u);

  // Nonlinear elasticity brick
  model.add_initialized_fixed_size_data("params", p);
  getfem::add_nonlinear_elasticity_brick(model, mim, "u", *pl, "params");

  // Incompressibility
  if (law_num == 1 || law_num == 3) {
    model.add_fem_variable("p", mf_p);
    getfem::add_nonlinear_incompressibility_brick(model, mim, "u", "p");
  }

  // Defining the volumic source term.
  base_vector f(N);
  f[0] = PARAM.real_value("FORCEX","Amplitude of the gravity");
  f[1] = PARAM.real_value("FORCEY","Amplitude of the gravity");
  if (N>2)
    f[2] = PARAM.real_value("FORCEZ","Amplitude of the gravity");
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i) {
    gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  }
  // Volumic source term brick.
  model.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(model, mim, "u", "VolumicData");
  // getfem::mdbrick_source_term<> VOL_F(*pINCOMP, mf_rhs, F);

  // Dirichlet condition
  plain_vector F2(nb_dof_rhs * N);
  model.add_initialized_fem_data("DirichletData", mf_rhs, F2);
  if (PARAM.int_value("DIRICHLET_VERSION") == 0)
    getfem::add_Dirichlet_condition_with_multipliers
      (model, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, "DirichletData");
  else
    getfem::add_Dirichlet_condition_with_penalization
      (model, mim, "u", 1E15, DIRICHLET_BOUNDARY_NUM, "DirichletData");


  gmm::iteration iter;


  /* prepare the export routine for OpenDX (all time steps will be exported) 
     (can be viewed with "dx -edit nonlinear_elastostatic.net")
  */
  getfem::dx_export exp(datafilename + ".dx",
			PARAM.int_value("VTK_EXPORT")==1);
  getfem::stored_mesh_slice sl; sl.build(mesh, getfem::slicer_boundary(mesh),8); 
  exp.exporting(sl,true); exp.exporting_mesh_edges();
  //exp.begin_series("deformationsteps");
  exp.write_point_data(mf_u, U, "stepinit"); 
  exp.serie_add_object("deformationsteps");

  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted for reduced mesh_fems");
  
  int nb_step = int(PARAM.int_value("NBSTEP"));
  size_type maxit = PARAM.int_value("MAXITER");
 
  for (int step = 0; step < nb_step; ++step) {
    plain_vector DF(F);

    gmm::copy(gmm::scaled(F, (step+1.)/(scalar_type)nb_step), DF);
    gmm::copy(DF, model.set_real_variable("VolumicData"));

    if (N>2) {
      /* Apply the gradual torsion/extension */
      scalar_type torsion = PARAM.real_value("TORSION","Amplitude of the torsion");
      torsion *= (step+1)/scalar_type(nb_step);
      scalar_type extension = PARAM.real_value("EXTENSION","Amplitude of the extension");
      extension *= (step+1)/scalar_type(nb_step);
      base_node G(N); G[0] = G[1] = 0.5;
      for (size_type i = 0; i < nb_dof_rhs; ++i) {
	const base_node P = mf_rhs.point_of_basic_dof(i) - G;
	scalar_type r = sqrt(P[0]*P[0]+P[1]*P[1]),
	  theta = atan2(P[1],P[0]);    
	F2[i*N+0] = r*cos(theta + (torsion*P[2])) - P[0]; 
	F2[i*N+1] = r*sin(theta + (torsion*P[2])) - P[1]; 
	F2[i*N+2] = extension * P[2];
      }
    }
    /* update the imposed displacement  */
    gmm::copy(F2, model.set_real_variable("DirichletData"));

    cout << "step " << step << ", number of variables : " << model.nb_dof() << endl;
    iter = gmm::iteration(residual, int(PARAM.int_value("NOISY")), maxit ? maxit : 40000);

    /* let the default non-linear solve (Newton) do its job */
    getfem::standard_solve(model, iter);

    pl->reset_unvalid_flag();
    model.assembly(getfem::model::BUILD_RHS);
    if (pl->get_unvalid_flag()) 
      GMM_WARNING1("The solution is not completely valid, the determinant "
		   "of the transformation is negative on "
		   << pl->get_unvalid_flag() << " gauss points");

    gmm::copy(model.real_variable("u"), U);
    //char s[100]; sprintf(s, "step%d", step+1);

    /* append the new displacement to the exported opendx file */
    exp.write_point_data(mf_u, U); //, s);
    exp.serie_add_object("deformationsteps");
  }

  // Solution extraction
  gmm::copy(model.real_variable("u"), U);
  
  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    elastostatic_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.mf_u.write_to_file(p.datafilename + ".mf", true);
    p.mf_rhs.write_to_file(p.datafilename + ".mfd", true);
    plain_vector U(p.mf_u.nb_dof());
    GMM_ASSERT1(p.solve(U), "Solve has failed");
    if (p.PARAM.int_value("VTK_EXPORT")) {
      gmm::vecsave(p.datafilename + ".U", U);
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
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
