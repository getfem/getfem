/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Julien Pommier.                    */
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
 * Helmholtz problem (Delta(u) + k^2 u = 0)
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
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
using bgeot::complex_type; /* = std::complex<double> */
using bgeot::size_type;   /* = unsigned long */
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
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  complex_type wave_number;

  scalar_type residu;        /* max residu for the iterative solvers         */
  bool with_mult;

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  Helmholtz_problem(void) : mf_u(mesh), mf_rhs(mesh), mf_coef(mesh) {}
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
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  size_type Nt = PARAM.int_value("NTHETA", "Nomber of space steps "), 
    Nr=PARAM.int_value("NR", "Nomber of space steps ");
  size_type gt_order = PARAM.int_value("GTDEGREE",
		       "polynomial degree of geometric transformation");
  scalar_type dtheta=2*M_PI*1./Nt;
  scalar_type R0 = PARAM.real_value("R0","R0");
  scalar_type R1 = PARAM.real_value("R1","R1");
  scalar_type dR = (R1-R0)/(Nr-1);
  bgeot::pgeometric_trans pgt = bgeot::parallelepiped_geotrans(2, gt_order);
  for (size_type i=0; i < Nt; ++i) {
    for (size_type j=0; j < Nr-1; ++j) {
      std::vector<size_type> ipts; ipts.reserve(dal::sqr(gt_order+1));
      for (size_type ii=0; ii <= gt_order; ++ii) {
        for (size_type jj=0; jj <= gt_order; ++jj) {
          scalar_type r = R0 + j*dR + jj*(dR/gt_order);
          scalar_type t = i*dtheta + ii*dtheta/gt_order;
          ipts.push_back(mesh.add_point(base_node(r*cos(t),r*sin(t))));
        }
      }
      mesh.add_convex(pgt, ipts.begin());
    }
  }

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  wave_number = complex_type
    (PARAM.real_value("WAVENUM_R", "Real part of the wave number"),
     PARAM.real_value("WAVENUM_I", "Imaginary part of the wave number"));

  with_mult = PARAM.int_value("DIRICHLET_WITH_MULTIPLIERS",
			      "with multipliers");

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mf_u.set_finite_element(mesh.convex_index(), pf_u, ppi);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u, ppi);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name), ppi);
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0), ppi);

  /* select boundaries */
  cout << "Selecting Robin and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    if (bgeot::vect_norm2(mesh.points_of_face_of_convex(it->cv,
							it->f)[0]) > 5.) {
      mf_u.add_boundary_elt(ROBIN_BOUNDARY_NUM, it->cv, it->f);
    } else mf_u.add_boundary_elt(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

typedef getfem::standard_complex_model_state MODELSTATE;

bool Helmholtz_problem::solve(plain_vector &U) {

  // Helmholtz brick.
  getfem::mdbrick_Helmholtz<MODELSTATE> WAVE(mf_u, mf_coef, wave_number);
  
  // (homogeneous) Robin condition
  getfem::mdbrick_QU_term<MODELSTATE> 
    ROBIN(WAVE, mf_rhs, wave_number * complex_type(0,1.), ROBIN_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  plain_vector F(nb_dof_rhs);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = incoming_field(mf_rhs.point_of_dof(i), wave_number.real());

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<MODELSTATE> 
    final_model(ROBIN, mf_rhs, F, DIRICHLET_BOUNDARY_NUM, with_mult);
  
  // Generic solve.
  cout << "Number of variables : " << final_model.nb_dof() << endl;
  cout << "Number of constraints : " << final_model.nb_constraints() << endl;
  MODELSTATE MS;
  gmm::iteration iter(residu, 1, 400000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  WAVE.get_solution(MS, U);

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
    Helmholtz_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      getfem::stored_mesh_slice sl(p.mesh, 8);
      exp.exporting(sl);
      exp.write_point_data(p.mf_u, gmm::real_part(U), "helmholtz_rfield");
      exp.write_point_data(p.mf_u, gmm::imag_part(U), "helmholtz_ifield");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d helmholtz.vtk -f WarpScalar -m BandedSurfaceMap -m Outline"
	"\n";
    }
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
