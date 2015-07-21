/*===========================================================================
 
 Copyright (C) 2002-2015 Yves Renard, Julien Pommier.
 
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
 * Goal : scalar Signorini problem with Xfem.
 *
 * Research program.
 */

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_contact_and_friction_nodal.h"
#include "getfem/getfem_import.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_vector;
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/* 
 * Exact solution 
 */
double Radius;

double u_exact(const base_node &p) {
  double sum = std::accumulate(p.begin(), p.end(), double(0));
  double norm_sqr = gmm::vect_norm2_sqr(p);
  return 5.0 * sin(sum) * (norm_sqr - Radius*Radius);
}

double rhs(const base_node &p) {
  double sum = std::accumulate(p.begin(), p.end(), double(0));
  double norm_sqr = gmm::vect_norm2_sqr(p);
  double N = double(gmm::vect_size(p));
  return 5.0 * (N * sin(sum) * (norm_sqr - Radius*Radius-2.0)
		- 4.0 * sum * cos(sum));
}

/*
 * Test procedure
 */

void test_mim(getfem::mesh_im_level_set &mim, getfem::mesh_fem &mf_rhs,
	      bool bound) {
  unsigned N =  mim.linked_mesh().dim();
  size_type nbdof = mf_rhs.nb_dof();
  plain_vector V(nbdof), W(1);
  std::fill(V.begin(), V.end(), 1.0);

  getfem::generic_assembly assem("u = data(#1); V()+=comp(Base(#1))(i).u(i);");
  assem.push_mi(mim); assem.push_mf(mf_rhs); assem.push_data(V);
  assem.push_vec(W);
  assem.assembly(getfem::mesh_region::all_convexes());
  double exact(0), R2 = Radius*Radius, R3 = R2*Radius;
  switch (N) {
  case 1: exact = bound ? 1.0 : Radius; break;
  case 2: exact = bound ? Radius*M_PI : R2*M_PI/2.0; break;
  case 3: exact = bound ? 2.0*M_PI*R2 : 2.0*M_PI*R3/3.0; break;
  default: assert(N <= 3);
  }
  if (bound) cout << "Boundary length: "; else cout << "Area: ";
  cout << W[0] << " should be " << exact << endl;
  assert(gmm::abs(exact-W[0])/exact < 0.01); 
}

/* 
 * Main program 
 */

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    
    // Read parameters.
    bgeot::md_param PARAM;
    PARAM.read_command_line(argc, argv);
    bool signorini = (PARAM.int_value("SIGNORINI", "Is signorini ?") != 0);

    // Load the mesh
    getfem::mesh mesh;
    std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
    getfem::import_mesh(MESH_FILE, mesh);
    unsigned N = mesh.dim();

    // center the mesh in (0, 0).
    base_node Pmin(N), Pmax(N);
    mesh.bounding_box(Pmin, Pmax);
    Pmin += Pmax; Pmin /= -2.0;
    Pmin[N-1] = -Pmax[N-1];
    mesh.translation(Pmin);

    // Level set definition
    unsigned lsdeg = unsigned(PARAM.int_value("LEVEL_SET_DEGREE", "level set degree"));
    Radius = PARAM.real_value("RADIUS", "Domain radius");
    getfem::level_set ls(mesh, bgeot::dim_type(lsdeg));
    const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
    for (unsigned i = 0; i < lsmf.nb_dof(); ++i)
      ls.values()[i] = gmm::vect_norm2_sqr(lsmf.point_of_basic_dof(i))-Radius*Radius;
    getfem::mesh_level_set mls(mesh);
    mls.add_level_set(ls);
    mls.adapt();
    
    getfem::mesh mcut;
    mls.global_cut_mesh(mcut);
    mcut.write_to_file("cut.mesh");

    // Integration method on the domain
    std::string IM = PARAM.string_value("IM", "Mesh file");
    std::string IMS = PARAM.string_value("IM_SIMPLEX", "Mesh file");
    int intins = getfem::mesh_im_level_set::INTEGRATE_INSIDE;
    getfem::mesh_im uncutmim(mesh);
    uncutmim.set_integration_method(mesh.convex_index(),
				    getfem::int_method_descriptor(IM));
    getfem::mesh_im_level_set mim(mls, intins,
				  getfem::int_method_descriptor(IMS));
    mim.set_integration_method(mesh.convex_index(),
			       getfem::int_method_descriptor(IM));
    mim.adapt();


    // Integration method on the boudary
    int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
    getfem::mesh_im_level_set mimbound(mls, intbound,
				       getfem::int_method_descriptor(IMS));
    mimbound.set_integration_method(mesh.convex_index(),
				    getfem::int_method_descriptor(IM));
    mimbound.adapt();
    

    // Finite element method for the unknown
    getfem::mesh_fem pre_mf(mesh);
    std::string FEM = PARAM.string_value("FEM", "finite element method");
    pre_mf.set_finite_element(mesh.convex_index(),
			      getfem::fem_descriptor(FEM));
    getfem::partial_mesh_fem mf(pre_mf);
    dal::bit_vector kept_dof = select_dofs_from_im(pre_mf, mim);
    dal::bit_vector rejected_elt;
    for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv)
      if (mim.int_method_of_element(cv) == getfem::im_none())
	rejected_elt.add(cv);
    mf.adapt(kept_dof, rejected_elt);
    size_type nb_dof = mf.nb_dof();
    
    // Finite element method for the rhs
    getfem::mesh_fem mf_rhs(mesh);
    std::string FEMR = PARAM.string_value("FEM_RHS", "finite element method");
    mf_rhs.set_finite_element(mesh.convex_index(),
			      getfem::fem_descriptor(FEMR));
    size_type nb_dof_rhs = mf_rhs.nb_dof();
    cout << "nb_dof_rhs = " << nb_dof_rhs << endl;
    
    // Finite element method for the multipliers
    getfem::mesh_fem pre_mf_mult(mesh);
    std::string FEMM = PARAM.string_value("FEM_MULT", "fem for multipliers");
    pre_mf_mult.set_finite_element(mesh.convex_index(),
				   getfem::fem_descriptor(FEMM));
    getfem::partial_mesh_fem mf_mult(pre_mf_mult);
    dal::bit_vector kept_dof_mult
      = select_dofs_from_im(pre_mf_mult, mimbound,N-1);
    mf_mult.adapt(kept_dof_mult, rejected_elt);
    size_type nb_dof_mult = mf_mult.nb_dof();
    cout << "nb_dof_mult = " << nb_dof_mult << endl;

    // Tests
    test_mim(mim, mf_rhs, false);
    test_mim(mimbound, mf_rhs, true);

    // Selecting Dirichlet boundary
    unsigned dirichlet_boundary = 1;
    getfem::mesh_region border_faces;
    getfem::outer_faces_of_mesh(mesh, border_faces);
    for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
      base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
      un /= gmm::vect_norm2(un);
      if (gmm::abs(un[N-1] - 1.0) < 0.01)
	mesh.region(dirichlet_boundary).add(i.cv(), i.f());
    }

 
    // Mass matrix on the boundary
    getfem::CONTACT_B_MATRIX B(nb_dof_mult, nb_dof);
    getfem::asm_mass_matrix(B, mimbound, mf_mult, mf);

    // Brick system
    getfem::model model;
    model.add_fem_variable("u", mf);
    model.add_initialized_scalar_data("a", 1.);
    getfem::add_generic_elliptic_brick(model, mim, "u", "a");

    plain_vector F(nb_dof_rhs);
    getfem::interpolation_function(mf_rhs, F, u_exact, dirichlet_boundary);
    model.add_initialized_fem_data("DData", mf_rhs, F);
    getfem::add_Dirichlet_condition_with_multipliers
      (model, mim, "u", mf_mult, dirichlet_boundary, "DData");
    
    getfem::interpolation_function(mf_rhs, F, rhs);
    model.add_initialized_fem_data("VolumicData", mf_rhs, F);
    getfem::add_source_term_brick(model, mim, "u", "VolumicData");
    
    model.add_fixed_size_variable("mult_n", nb_dof_mult);

    if (signorini) {
      model.add_initialized_scalar_data("r", 1);
      getfem::add_basic_contact_brick(model, "u", "mult_n", "r", B);
    } else {
      getfem::add_constraint_with_multipliers(model, "u", "mult_n", B,
                                              plain_vector(nb_dof_mult));
    }

    // Solving the problem
    cout << "Total number of unknown: " << model.nb_dof() << endl;
    gmm::iteration iter(1e-9, 1, 40000);
    getfem::standard_solve(model, iter);
    plain_vector U(nb_dof);
    gmm::copy(model.real_variable("u"), U);

    // interpolation of the solution on mf_rhs
    plain_vector Uint(nb_dof_rhs), Vint(nb_dof_rhs);
    getfem::interpolation(mf, mf_rhs, U, Uint);
    getfem::interpolation_function(mf_rhs, Vint, u_exact);

    // computation of max error.
    if (!signorini) {
      GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
      double errmax = 0.0;
      for (size_type i = 0; i < nb_dof_rhs; ++i)
	if (gmm::vect_norm2(mf_rhs.point_of_basic_dof(i)) < Radius)
	  errmax = std::max(errmax, gmm::abs(Uint[i]-Vint[i]));
      cout << "Linfty error: " << errmax << endl;
      cout << "L2 error: " << getfem::asm_L2_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	   << endl;
      cout << "H1 error: " << getfem::asm_H1_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	   << endl;
    }
      
    // export de la solution au format vtk.
    getfem::vtk_export exp("xfem_contact.vtk", (2==1));
    exp.exporting(mf); 
    exp.write_point_data(mf, U, "solution");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi2 -d xfem_contact.vtk -f WarpScalar -m Surface -m Outline\n";
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
