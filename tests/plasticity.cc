/*===========================================================================

 Copyright (C) 2000-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
   @file plasticity.cc
   @brief Small deformation plasticity problem.

   This program is used to check that getfem++ is working. 
   This is also a good example of use of GetFEM.
*/

#include "getfem/getfem_assembling.h" 
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_plasticity.h"
#include "getfem/getfem_export.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some GetFEM types that we will be using */

/* special class for small (dim<16) vectors */
using bgeot::base_small_vector; 
/* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_node;  
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
   using the predefined types in Gmm++ */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

template<typename VEC> 
static void vecsave(std::string fname, const VEC& V);

//function to save a vector
template<typename VEC> 
static void vecsave( std::string fname, const VEC& V) {

  std::ofstream f(fname.c_str()); f.precision(16);
  for (size_type i=0; i < V.size(); ++i) 
    f << V[i] << "\n"; 
}
 


//=================================================================
// structure for the elastoplastic problem
//=================================================================

struct elastoplasticity_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, 
	 NEUMANN_BOUNDARY_NUM = 1};

  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im  mim;      /* integration methods. */
  getfem::im_data mim_data;  /* Mim data for the pastic strain. */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the
				elastoplastic displacement */
  getfem::mesh_fem mf_xi;    /* mesh_fem, for the plastic multiplier. */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side 
				(f(x),..)   */
  scalar_type lambda, mu;    /* Lamé coefficients. */

  scalar_type residual; /* max residual for the iterative solvers */
  scalar_type sigma_y;
  size_type flag_hyp;
  std::string datafilename;
  bgeot::md_param PARAM;
  bool do_export;

  bool solve(plain_vector &U);
  void init(void);

  elastoplasticity_problem(void) : mim(mesh), mim_data(mim), mf_u(mesh), 
			     mf_xi(mesh), mf_rhs(mesh) {}

};




/* Read parameters from the .param file, build the mesh, 
   set finite element and integration methods 
   and selects the boundaries.
*/
void elastoplasticity_problem::init(void) {

  std::string MESH_TYPE = 
    PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = 
    PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_XI = 
    PARAM.string_value("FEM_TYPE_XI","FEM name");
  std::string INTEGRATION = 
    PARAM.string_value("INTEGRATION", 
		       "Name of integration method");
  do_export = (PARAM.int_value("EXPORT", "Perform or not the vtk export") != 0);

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  residual = PARAM.real_value("RESIDUAL", "residual");

  //  file to save the mesh
  datafilename = PARAM.string_value("ROOTFILENAME",
				    "Filename for saving");

  /* First step : build the mesh */
  size_type N;
  bgeot::pgeometric_trans pgt = 0; 
  
  if (MESH_TYPE != "load") {
    std::cout << "created getfem mesh"  << "\n"; 
    pgt = bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    nsubdiv[0]=PARAM.int_value
      ("NX", "Number of space steps in x direction ");
    nsubdiv[1]=PARAM.int_value
      ("NY", "Number of space steps in y direction ");

    if(N==3)
      nsubdiv[2]=PARAM.int_value
	("NZ", "Number of space steps in z direction ");
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED")!= 0);
    
    bgeot::base_matrix M(N,N);

    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }

    if (N>1) { 
      M(0,1) = PARAM.real_value("INCLINE") * 
	PARAM.real_value("LY"); 
    }

    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);

  } else {
    std ::cout << "mesh from pdetool"  << "\n"; 
    std::string MESH_FILE 
      = PARAM.string_value("MESH_FILE", "Mesh file name");

    mesh.read_from_file(MESH_FILE);
    
    N = mesh.dim();
    pgt = mesh.trans_of_convex
      (mesh.convex_index().first_true());
  }

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  /*  muT = PARAM.real_value("MUT", "Lamé coefficient muT");
  lambdaB = PARAM.real_value("LAMBDAB", 
			    "Lamé coefficient lambdaB");
  */
  lambda = PARAM.real_value("LAMBDA", 
			     "Lamé coefficient lambda");
  mf_u.set_qdim(bgeot::dim_type(N));


  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mim_data.set_tensor_size(bgeot::multi_index(N,N));
  mf_u.set_finite_element(mesh.convex_index(), pf_u);


  /* set the finite element on the mf_sigma */
  getfem::pfem pf_xi = getfem::fem_descriptor(FEM_TYPE_XI);
  mf_xi.set_finite_element(mesh.convex_index(), pf_xi);


  /* set the finite element on mf_rhs 
     (same as mf_u is DATA_FEM_TYPE is not used in the .param file)*/
  std::string data_fem_name 
    = PARAM.string_value("DATA_FEM_TYPE");
  
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), 
		"You are using a non-lagrange FEM. "
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

  for (getfem::mr_visitor it(border_faces); 
       !it.finished(); ++it) {
    assert(it.is_face());
    base_node un 
      = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);

    if (gmm::abs(un[0] - 1.0) < 1.0E-7)
      mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
    else if (gmm::abs(un[0] + 1.0) < 1.0E-7) 
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
      
  }
 
  // Plasticity part 
  sigma_y = PARAM.real_value("SIGMA_Y", "plasticity yield stress");
  flag_hyp=PARAM.int_value("FLAG_HYP");
}



//==================================================================
// Model.                                                           
//==================================================================

bool elastoplasticity_problem::solve(plain_vector &U) {

  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  getfem::model model;

  gmm::set_traces_level(1);

  // Main unknown of the problem.
  model.add_fem_variable("u", mf_u);
  model.add_fem_data("Previous_u", mf_u);

  model.add_initialized_scalar_data("lambda", lambda);
  model.add_initialized_scalar_data("mu", mu);
  model.add_initialized_scalar_data("sigma_y", sigma_y);

  model.add_fem_data("xi", mf_xi);
  model.add_fem_data("Previous_xi", mf_xi);
  model.add_im_data("Previous_Ep", mim_data);

  /* choose the projection type */
  getfem::pconstraints_projection
    proj = std::make_shared<getfem::VM_projection>(0);

  std::vector<std::string> plastic_variables = {"u", "xi", "Previous_Ep"};
  std::vector<std::string> plastic_data = {"lambda", "mu", "sigma_y"};
  

  add_small_strain_elastoplasticity_brick
    (model, mim, "Prandtl Reuss", getfem::DISPLACEMENT_ONLY,
     plastic_variables, plastic_data);
  
  plain_vector F(nb_dof_rhs * N);
  model.add_initialized_fem_data("NeumannData", mf_rhs, F);
  getfem::add_source_term_brick
    (model, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM);
  
  model.add_initialized_fem_data("DirichletData", mf_rhs, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (model, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, 
     "DirichletData");

  const size_type Nb_t = 19;
  std::vector<scalar_type> T(19);
  T[0] = 0; T[1] = 0.9032; T[2] = 1; T[3] = 1.1; T[4] = 1.3;
  T[5] = 1.5; T[6] = 1.7; T[7] = 1.74; T[8] = 1.7; T[9] = 1.5;
  T[10] = 1.3; T[11] = 1.1; T[12] = 1; T[13] = 0.9032; T[14] = 0.7;
  T[15] = 0.5; T[16] = 0.3; T[17] = 0.1; T[18] = 0;


  getfem::mesh_fem mf_vm(mesh);
  mf_vm.set_classical_discontinuous_finite_element(1);
  getfem::base_vector VM(mf_vm.nb_dof());
  getfem::base_vector plast(mf_vm.nb_dof());
  
  for (size_type nb = 0; nb < Nb_t; ++nb) {
    cout << "=============iteration number : " << nb << "==========" << endl;
   
    scalar_type t = T[nb];
    
    // Defining the Neumann condition right hand side.
    base_small_vector v(N);
    v[N-1] = -PARAM.real_value("FORCE");
    gmm::scale(v,t); 
    
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(v, gmm::sub_vector
		(F, gmm::sub_interval(i*N, N)));
    
    gmm::copy(F, model.set_real_variable("NeumannData"));
    
    // Generic solve.
    cout << "Number of variables : " 
	 << model.nb_dof() << endl;
    
    getfem::newton_search_with_step_control ls;
    // getfem::simplest_newton_line_search ls;
    gmm::iteration iter(residual, 2, 40000);
    getfem::standard_solve(model, iter,
#ifdef GMM_USES_MUMPS
			   getfem::rselect_linear_solver(model, "mumps"), ls);
#else
			   getfem::rselect_linear_solver(model, "superlu"), ls);
#endif
 
    getfem::small_strain_elastoplasticity_next_iter
      (model, mim, "Prandtl Reuss", getfem::DISPLACEMENT_ONLY,
       plastic_variables, plastic_data);
    
    // Get the solution and save it
    gmm::copy(model.real_variable("u"), U);
    
 
    getfem::compute_small_strain_elastoplasticity_Von_Mises
      (model, mim, "Prandtl Reuss", getfem::DISPLACEMENT_ONLY,
       plastic_variables, plastic_data, mf_vm, VM);
    
    std::stringstream fname; fname << datafilename << "_" << nb << ".vtk";

    if (do_export) {
      getfem::vtk_export exp(fname.str());
      exp.exporting(mf_vm);
      exp.write_point_data(mf_vm,VM, "Von Mises stress");
      exp.write_point_data(mf_u, U, "displacement");
    }
    
  }

  if (do_export) {
    cout << "export done, you can view the data file with "
      "(for example)\n"
      "mayavi2 -d " << datafilename << "_1.vtk -f "
      "WarpVector -m Surface -m Outline\n";
  }

  return true;
}

  
//==================================================================
// main program.                                                    
//==================================================================

int main(int argc, char *argv[]) {

  GETFEM_MPI_INIT(argc, argv);
  GMM_SET_EXCEPTION_DEBUG; 
  // Exceptions make a memory fault, to debug.

  FE_ENABLE_EXCEPT;        
  // Enable floating point exception for Nan.
   
  elastoplasticity_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  plain_vector U(p.mf_u.nb_dof());
  if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");
  vecsave(p.datafilename + ".U", U);
  cout << "Resultats dans fichier : "
       <<p.datafilename<<".U \n";
  p.mf_u.write_to_file(p.datafilename + ".meshfem",true);
  scalar_type t[2]={p.mu,p.lambda};
  vecsave(p.datafilename+".coef", 
	  std::vector<scalar_type>(t, t+2));    

  GETFEM_MPI_FINALIZE;
  
  return 0; 
}
