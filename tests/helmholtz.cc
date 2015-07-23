/*===========================================================================

 Copyright (C) 2002-2015 Yves Renard, Julien Pommier.

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
   @file helmholtz.cc
   @brief Helmholtz problem (Delta(u) + k^2 u = 0)

   Diffraction of a plane wave by a circular obstacle. 

   This program is used to check that getfem++ is working. This is also 
   a good example of use of GetFEM++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some GetFEM++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::complex_type; /* = std::complex<double> */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
   using the predefined types in Gmm++ */
typedef getfem::modeling_standard_complex_sparse_vector sparse_vector;
typedef getfem::modeling_standard_complex_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_complex_plain_vector  plain_vector;

/*
  structure for the Helmholtz problem
*/
struct Helmholtz_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, ROBIN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  complex_type wave_number;

  scalar_type residual;        /* max residual for the iterative solvers     */
  int with_mult;

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  Helmholtz_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh) {}
};

complex_type __wave_number;

complex_type incoming_field(const base_node &P) {
  return complex_type(cos(__wave_number.real()*P[1]+.2),
		      sin(__wave_number.real()*P[1]+.2));
  /*scalar_type s = 0;
  for (size_type i=1; i < P.size(); ++i) s += P[i]*(1.-P[i]);
  s = rand()*3. / RAND_MAX;
  return s;
  */
}

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void Helmholtz_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  size_type Nt = PARAM.int_value("NTHETA", "Nomber of space steps "), 
    Nr=PARAM.int_value("NR", "Nomber of space steps ");
  size_type gt_order = PARAM.int_value("GTDEGREE",
		       "polynomial degree of geometric transformation");
  scalar_type dtheta=2*M_PI*1./scalar_type(Nt);
  scalar_type R0 = PARAM.real_value("R0","R0");
  scalar_type R1 = PARAM.real_value("R1","R1");
  scalar_type dR = (R1-R0)/scalar_type(Nr-1);
  bgeot::pgeometric_trans
    pgt = bgeot::parallelepiped_geotrans(2, short_type(gt_order));
  for (size_type i=0; i < Nt; ++i) {
    for (size_type j=0; j < Nr-1; ++j) {
      std::vector<size_type> ipts; ipts.reserve(gmm::sqr(gt_order+1));
      for (size_type ii=0; ii <= gt_order; ++ii) {
        for (size_type jj=0; jj <= gt_order; ++jj) {
          scalar_type r = R0 + scalar_type(j)*dR
	    + scalar_type(jj)*(dR/scalar_type(gt_order));
          scalar_type t = scalar_type(i)*dtheta
	    + scalar_type(ii)*dtheta/scalar_type(gt_order);
          ipts.push_back(mesh.add_point(base_node(r*cos(t),r*sin(t))));
        }
      }
      mesh.add_convex(pgt, ipts.begin());
    }
  }

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

  __wave_number = wave_number = complex_type
    (PARAM.real_value("WAVENUM_R", "Real part of the wave number"),
     PARAM.real_value("WAVENUM_I", "Imaginary part of the wave number"));

  with_mult = int(PARAM.int_value("DIRICHLET_VERSION",
				  "Dirichlet condition version"));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(ppi);
  mf_u.set_finite_element(pf_u);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(pf_u);
  } else {
    mf_rhs.set_finite_element(getfem::fem_descriptor(data_fem_name));
  }
  

  /* select boundaries */
  cout << "Selecting Robin and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    if (gmm::vect_norm2(mesh.points_of_face_of_convex(i.cv(),
							i.f())[0]) > 5.) {
      mesh.region(ROBIN_BOUNDARY_NUM).add(i.cv(),i.f());
    } else mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/


bool Helmholtz_problem::solve(plain_vector &U) {

  // Complex model.
  getfem::model Helmholtz_model(true);

  // Main unknown of the problem
  Helmholtz_model.add_fem_variable("u", mf_u);

  // Helmholtz term on u.
  plain_vector K(1); K[0] = wave_number;
  Helmholtz_model.add_initialized_fixed_size_data("k", K);
  add_Helmholtz_brick(Helmholtz_model, mim, "u", "k");

  // Fourier-Robin condition.
  plain_vector Q(1); Q[0] = wave_number * complex_type(0,1.);
  Helmholtz_model.add_initialized_fixed_size_data("Q", Q);
  add_Fourier_Robin_brick(Helmholtz_model, mim, "u", "Q", ROBIN_BOUNDARY_NUM);

  // Dirichlet condition
  plain_vector F(mf_rhs.nb_dof());
  getfem::interpolation_function(mf_rhs, F, incoming_field);
  Helmholtz_model.add_initialized_fem_data("DirichletData", mf_rhs, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (Helmholtz_model, mim, "u", mf_u,
     DIRICHLET_BOUNDARY_NUM, "DirichletData");

  // Helmholtz_model.listvar(cout);

  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(Helmholtz_model, iter);

  gmm::resize(U, mf_u.nb_dof());
  gmm::copy(Helmholtz_model.complex_variable("u"), U);

  cout << "U = " << U << endl;

  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    Helmholtz_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      getfem::stored_mesh_slice sl(p.mesh, p.mesh.nb_convex() < 2000 ? 8 : 6);
      exp.exporting(sl);
      exp.write_point_data(p.mf_u, gmm::real_part(U), "helmholtz_rfield");
      exp.write_point_data(p.mf_u, gmm::imag_part(U), "helmholtz_ifield");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi2 -d helmholtz.vtk -f WarpScalar -m Surface -m Outline"
	"\n";
    }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
