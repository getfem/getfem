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

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  Helmholtz_problem(void) : mf_u(mesh), mf_rhs(mesh), mf_coef(mesh) {}
};

scalar_type incoming_field(const base_node &P) {
  scalar_type s = 0;
  for (size_type i=1; i < P.size(); ++i) s += P[i]*(1.-P[i]);
  return s;
}

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void Helmholtz_problem::init(void) {
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
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  wave_number = complex_type(PARAM.real_value("WAVENUM_R", "Real part of the wave number"),
			     PARAM.real_value("WAVENUM_I", "Imaginary part of the wave number"));

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
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM "
		<< data_fem_name << ". In that case you need to set "
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
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    base_node un = mesh.normal_of_face_of_convex(it->cv, it->f);
    un /= gmm::vect_norm2(un);
    if (dal::abs(un[0]) < 1.0E-7) { // new Neumann face
      mf_u.add_boundary_elt(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
    } else if (dal::abs(un[0] + 1.) < 1e-7) {
      mf_u.add_boundary_elt(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
    } else if (dal::abs(un[0] - 1.) < 1e-7) {
      mf_u.add_boundary_elt(ROBIN_BOUNDARY_NUM, it->cv, it->f);
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool Helmholtz_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();

  // Helmholtz brick.
  getfem::mdbrick_Helmholtz<> WAVE(mf_u, mf_coef, wave_number);
  
  // (homogeneous) Robin condition
  getfem::mdbrick_QU_term<getfem::standard_complex_model_state> 
    ROBIN(WAVE, mf_rhs, wave_number * complex_type(0,1.), 
	  ROBIN_BOUNDARY_NUM);
  
  // Defining the Dirichlet condition value.
  plain_vector F(nb_dof_rhs);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i] = incoming_field(mf_rhs.point_of_dof(i));

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<getfem::standard_complex_model_state> 
    final_model(ROBIN, mf_rhs,
		F, DIRICHLET_BOUNDARY_NUM);

  // Generic solve.
  cout << "Number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_complex_model_state MS;
  gmm::iteration iter(residu, 1, 40000);
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
    //p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");

    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.write_dataset(p.mf_u, gmm::real_part(U), "helmholtz_field");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d helmholtz.vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline\n";
    }
    // getfem::save_solution(p.datafilename + ".dataelt", p.mf_u, p.U, p.K);
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
