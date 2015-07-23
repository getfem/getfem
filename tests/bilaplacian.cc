/*===========================================================================

 Copyright (C) 2006-2015 Yves Renard, Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
   @file bilaplacian.cc
   @brief Bilaplacian problem. A dummy
   bilaplacian problem is solved on a regular mesh, and is compared to
   the analytical solution.

   This program is used to check that getfem++ is working. This is also 
   a good example of use of GetFEM++.

   @see laplacian.cc
*/

#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm_inoutput.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some GetFEM++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

typedef getfem::model_real_plain_vector  plain_vector;

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

scalar_type FT = 0.0;
scalar_type D  = 1.0;
scalar_type nu = 0.0;
scalar_type pressure  = 1.0 ;

#if 1

scalar_type sol_u(const base_node &x)
{ return sin(FT*std::accumulate(x.begin(), x.end(), 0.0)); }
scalar_type sol_lapl_u(const base_node &x)
{ return -FT*FT*sol_u(x) * x.size(); }
scalar_type sol_f(const base_node &x)
{ return FT*FT*FT*FT*sol_u(x)*gmm::sqr(x.size()); }
base_small_vector sol_du(const base_node &x) {
  base_small_vector res(x.size());
  std::fill(res.begin(), res.end(),
	    FT * cos(FT*std::accumulate(x.begin(), x.end(), 0.0)));
  return res;
}
base_small_vector neumann_val(const base_node &x)
{ return -FT*FT*sol_du(x) * scalar_type(x.size()); }

base_matrix sol_hessian(const base_node &x) {
  base_matrix m(x.size(), x.size());
  std::fill(m.begin(), m.end(), -FT*FT*sol_u(x));
  return m;
}

#else

scalar_type sol_u(const base_node &x) {
 return pressure/(2. * D) * x[0] * x[0] * ( x[0] * x[0] / 12. - x[0] / 3. + 1./2. ); }
 
scalar_type sol_lapl_u(const base_node &) { 
return pressure/(2. * D) * ( x[0] * x[0] - 2. * x[0] + 1. ); }

scalar_type sol_f(const base_node &) { return pressure ; }
base_small_vector sol_du(const base_node &x) {
  base_small_vector res(x.size()); 
  res[0] = pressure/(2. * D) * x[0] * ( x[0] * x[0] / 3. - x[0] + 1. ); 
  res[1] = 0.; 
  return res;
}
base_small_vector neumann_val(const base_node &x)
{ cout << "erreur : KL = 1 doit etre selectionne \n";
  return 0 ; /*base_small_vector(x.size()); */ }

base_matrix sol_hessian(const base_node &x)
{ base_matrix m(x.size(), x.size()); 
  m(0,0) = sol_lapl_u(x);
  return m; }


// scalar_type sol_u(const base_node &x) { return x[0]*x[1]; }
// scalar_type sol_lapl_u(const base_node &) { return 0.0; }
// scalar_type sol_f(const base_node &) { return 0.0; }
// base_small_vector sol_du(const base_node &x) {
//   base_small_vector res(x.size()); res[0] = x[1]; res[1] = x[0]; 
//   return res;
// }
// base_small_vector neumann_val(const base_node &x)
// { return base_small_vector(x.size());  }
// 
// base_matrix sol_hessian(const base_node &x)
// { base_matrix m(x.size(), x.size()); m(1,0) = m(0,1) = 1.0; return m; }

#endif

base_matrix sol_mtensor(const base_node &x) { 
  base_matrix m = sol_hessian(x), mm(x.size(), x.size());
  scalar_type l = sol_lapl_u(x);
  for (size_type i = 0; i < x.size(); ++i) mm(i,i) = l * nu;
  gmm::scale(m, (1-nu));
  gmm::add(mm, m);
  gmm::scale(m, -D);
  return m;
}

base_small_vector sol_bf(const base_node &x)
{ return -D * neumann_val(x); }

/**************************************************************************/
/*  Structure for the bilaplacian problem.                                */
/**************************************************************************/

struct bilaplacian_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3};
  /* Boundary conditions :
    CLAMPED_BOUNDARY_NUM        --> normal derivative of the solution
    SIMPLE_SUPPORT_BOUNDARY_NUM --> value of the solution
    FORCE_BOUNDARY_NUM          --> the "big" condition
    MOMENTUM_BOUNDARY_NUM       --> flexure moment condition  */

  getfem::mesh mesh;        /* the mesh */
  getfem::mesh_im mim;      /* the integration methods.                     */
  getfem::mesh_fem mf_u;    /* main mesh_fem, for the bilaplacian solution  */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */

  scalar_type residual;     /* max residual for the iterative solvers       */
  size_type dirichlet_version;

  std::string datafilename;
  bgeot::md_param PARAM;

  bool KL;
  size_type NX ;
  size_type boundary_ref ;
  /* boundary_ref = 0 corresponds to clamped edge on the 4 edges   
     boundary_ref = 1 corresponds to :    
                                    
                                               free edge
                                              __________
                                             |          |
                                   simple    |          |   simple
                                   support   |          |   support
                                             |          |
                                             |__________|
                                    
                                               clamped
                                               edge
  */
                                 
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
    MESH_TYPE = bgeot::name_of_geometric_trans
      (mesh.trans_of_convex(mesh.convex_index().first_true()));
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    N = mesh.dim();
  } else {
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    
    /* First step : build the mesh */
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    NX = PARAM.int_value("NX", "Nomber of space steps ") ;
    std::fill(nsubdiv.begin(),nsubdiv.end(), NX );
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

  int dv = int(PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version"));
  dirichlet_version = size_type(dv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  FT = PARAM.real_value("FT"); if (FT == 0.0) FT = 1.0;
  pressure = PARAM.real_value("pressure") ; 
  KL = (PARAM.int_value("KL", "Kirchhoff-Love model or not") != 0);
  D = PARAM.real_value("D", "Flexion modulus");
  if (KL) nu = PARAM.real_value("NU", "Poisson ratio");
  boundary_ref = PARAM.int_value("BOUNDARY_REF");
  
  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);

  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0) {
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
    cout << "No DIRICHLET_FEM_TYPE\n";
  } else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_ASSERT1(false, "You are using a non-lagrange FEM. "
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
  
  if (boundary_ref == 0) { // clamped plate
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
	mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
     }
  }
  
  if (boundary_ref == 1) { // mixed boundary conditions
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
       base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
       un /= gmm::vect_norm2(un);
       if (gmm::abs(un[N-1] - 1.0) <= 0.5) {
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
   
   if (boundary_ref == 2) { // clamped plate on the left edge, free boundary elsewhere
     for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
	if (un[0] < - 0.5 ) { 
           mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
 	   mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f()); }
	else
	   mesh.region(FORCE_BOUNDARY_NUM).add(i.cv(), i.f());
           mesh.region(MOMENTUM_BOUNDARY_NUM).add(i.cv(), i.f());
     }
  }
}

/* compute the relative error with respect to the exact solution */
void bilaplacian_problem::compute_error(plain_vector &U) {

  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  std::vector<scalar_type> V(mf_rhs.nb_basic_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_basic_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_basic_dof(i));
  cout.precision(16);
  cout  << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V)  << endl
        << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V)  << endl
        << "H2 error = " << getfem::asm_H2_norm(mim, mf_rhs, V)  << endl
        /*<< "Linfty error = " << gmm::vect_norminf(V)  << endl*/; 
  cout  << "semi-norme H1 = " << getfem::asm_H1_semi_norm(mim, mf_rhs, V)  << endl 
        << "semi-norme H2 = " << getfem::asm_H2_semi_norm(mim, mf_rhs, V)  << endl ;
       
}



/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/


bool bilaplacian_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  
  getfem::model model;
  
  // Main unknown of the problem.
  model.add_fem_variable("u", mf_u);
  
  // Bilaplacian brick.
  model.add_initialized_scalar_data("D", D);
  if (KL) {
    model.add_initialized_scalar_data("nu", nu);
    getfem::add_bilaplacian_brick_KL(model, mim, "u", "D", "nu");
  } else {
    getfem::add_bilaplacian_brick(model, mim, "u", "D");
  }
  
  // Volumic source term brick.
  plain_vector F(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  model.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(model, mim, "u", "VolumicData");
  
  
  // Defining the prescribed momentum.
  if (KL) {
    gmm::resize(F, nb_dof_rhs*N*N);
    getfem::interpolation_function(mf_rhs, F,
                                   sol_mtensor, MOMENTUM_BOUNDARY_NUM);
    gmm::scale(F, -1.0);
  } else {
    gmm::resize(F, nb_dof_rhs);
    getfem::interpolation_function(mf_rhs, F, sol_lapl_u,
                                   MOMENTUM_BOUNDARY_NUM);
  }
  
  // Prescribed momentum on the boundary
  model.add_initialized_fem_data("Momentum", mf_rhs, F);
  getfem::add_normal_derivative_source_term_brick
    (model, mim, "u", "Momentum", MOMENTUM_BOUNDARY_NUM);
  
  // Defining the Neumann condition right hand side.
  plain_vector H(nb_dof_rhs*N*N);
  gmm::resize(F, nb_dof_rhs*N);
  if (KL) {
    getfem::interpolation_function(mf_rhs, F,sol_bf,FORCE_BOUNDARY_NUM);
    getfem::interpolation_function(mf_rhs, H,sol_mtensor,FORCE_BOUNDARY_NUM);
    model.add_initialized_fem_data("M", mf_rhs, F);
    model.add_initialized_fem_data("H", mf_rhs, H);
    add_Kirchoff_Love_Neumann_term_brick
      (model, mim, "u", "H", "M", FORCE_BOUNDARY_NUM);
  }
  else {
    getfem::interpolation_function(mf_rhs, F,neumann_val,FORCE_BOUNDARY_NUM);
    gmm::scale(F, -1.0);
    model.add_initialized_fem_data("Neumanndata", mf_rhs, F);
    add_normal_source_term_brick
      (model, mim, "u", "Neumanndata", FORCE_BOUNDARY_NUM);
  }
  
  // Normal derivative Dirichlet condition brick
  gmm::resize(F, nb_dof_rhs*N);
  gmm::clear(F);
  getfem::interpolation_function(mf_rhs, F, sol_du, CLAMPED_BOUNDARY_NUM);
  model.add_initialized_fem_data("DerDirdata", mf_rhs, F);
  if (dirichlet_version == 0)
    add_normal_derivative_Dirichlet_condition_with_multipliers
      (model, mim, "u", mf_u, CLAMPED_BOUNDARY_NUM, "DerDirdata", false);
  else
    add_normal_derivative_Dirichlet_condition_with_penalization
      (model, mim, "u", PARAM.real_value("EPS_DIRICHLET_PENAL"),
       CLAMPED_BOUNDARY_NUM, "DerDirdata", false);
  
  
  // Dirichlet condition brick.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function
    (mf_rhs, F, sol_u, SIMPLE_SUPPORT_BOUNDARY_NUM);
  model.add_initialized_fem_data("Dirichletdata", mf_rhs, F);
  if (dirichlet_version == 0)
    add_Dirichlet_condition_with_multipliers
      (model, mim, "u", mf_u, SIMPLE_SUPPORT_BOUNDARY_NUM, "Dirichletdata");
  else
    add_Dirichlet_condition_with_penalization
      (model, mim, "u", PARAM.real_value("EPS_DIRICHLET_PENAL"),
       SIMPLE_SUPPORT_BOUNDARY_NUM, "Dirichletdata", &mf_mult);
  
  // Generic solve.
  cout << "Total number of variables : " << model.nb_dof() << endl;
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(model, iter);
  gmm::resize(U, mf_u.nb_dof());
  gmm::copy(model.real_variable("u"), U);
  return (iter.converged());
 
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GETFEM_MPI_INIT(argc, argv); // For parallelized version
  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    bilaplacian_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U(p.mf_u.nb_dof()), V(p.mf_u.nb_dof()) ;
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

    p.compute_error(U);

    if (p.PARAM.int_value("VTK_EXPORT") ) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "bilaplacian_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	   << "mayavi2 -d " << p.datafilename
	   << ".vtk -f WarpScalar -m Surface -m Outline\n";
    }
  }
  GMM_STANDARD_CATCH_ERROR;

  GETFEM_MPI_FINALIZE;

  return 0; 
}
