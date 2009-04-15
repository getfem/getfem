// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2009 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

/**
   @file helmholtz.cc
   @brief Helmholtz problem (Delta(u) + k^2 u = 0)

   Diffraction of a plane wave by a circular obstacle. 

   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"

/* some Getfem++ types that we will be using */
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

  scalar_type residual;        /* max residual for the iterative solvers         */
  int with_mult;

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  Helmholtz_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh) {}
};

complex_type incoming_field(const base_node &P, scalar_type k) {
  return complex_type(cos(k*P[1]+.2),sin(k*P[1]+.2));
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

  wave_number = complex_type
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

typedef getfem::standard_complex_model_state MODELSTATE;

bool Helmholtz_problem::solve(plain_vector &U) {

  // Helmholtz brick.
  getfem::mdbrick_Helmholtz<MODELSTATE> WAVE(mim, mf_u, wave_number);
  
  // (homogeneous) Robin condition
  getfem::mdbrick_QU_term<MODELSTATE> 
    ROBIN(WAVE, wave_number * complex_type(0,1.), ROBIN_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  plain_vector F(nb_dof_rhs);
  GMM_ASSERT1(!mf_rhs.is_reduced(),
	      "To be adapted, use interpolation_function");
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = incoming_field(mf_rhs.point_of_basic_dof(i), wave_number.real());

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<MODELSTATE> 
    final_model(ROBIN, DIRICHLET_BOUNDARY_NUM);
  final_model.set_constraints_type(getfem::constraints_type(with_mult));
  final_model.rhs().set(mf_rhs, F);
  

  // Generic solve.
  cout << "Number of variables : " << final_model.nb_dof() << endl;
  cout << "Number of constraints : " << final_model.nb_constraints() << endl;
  MODELSTATE MS;
  gmm::iteration iter(residual, 1, 400000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(WAVE.get_solution(MS), U);

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
	"mayavi -d helmholtz.vtk -f WarpScalar -m BandedSurfaceMap -m Outline"
	"\n";
    }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
