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
 * Goal : scalar Signorini problem with Xfem.
 *
 * Research program.
 */

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_derivatives.h>
#include <getfem_regular_meshes.h>
#include <getfem_model_solvers.h>
#include <getfem_mesh_im_level_set.h>
#include <getfem_partial_mesh_fem.h>
#include <getfem_import.h>
#include <gmm.h>

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
  return sin(sum) * (norm_sqr - Radius*Radius);
}

double rhs(const base_node &p) {
  double sum = std::accumulate(p.begin(), p.end(), double(0));
  double norm_sqr = gmm::vect_norm2_sqr(p);
  double N = double(gmm::vect_size(p));
  return N * sin(sum) * (norm_sqr - Radius*Radius-2.0) - 4.0 * sum * cos(sum);
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
  case 1: exact = bound ? 2.0 : 2*Radius; break;
  case 2: exact = bound ? 2.0*Radius*M_PI : R2*M_PI; break;
  case 3: exact = bound ? 4.0*M_PI*R2 : 4.0*M_PI*R3/3.0; break;
  default: assert(N <= 3);
  }
  if (bound) cout << "Boundary length: "; else cout << "Area: ";
  cout << W[0] << " should be " << exact << endl;
  assert(gmm::abs(exact-W[0])/exact < 0.1); 
}

/* 
 * Main program 
 */

int main(int argc, char *argv[]) {

  DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  // getfem::getfem_mesh_level_set_noisy();

  try {
    
    // Read parameters.
    ftool::md_param PARAM;
    PARAM.read_command_line(argc, argv);

    // Load the mesh
    getfem::mesh mesh;
    std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
    getfem::import_mesh(MESH_FILE, mesh);
    unsigned N = mesh.dim();

    // center the mesh in (0, 0).
    base_node Pmin(N), Pmax(N);
    mesh.bounding_box(Pmin, Pmax);
    Pmin += Pmax; Pmin /= -2.0;
    mesh.translation(Pmin);

    // Level set definition
    unsigned lsdeg = PARAM.int_value("LEVEL_SET_DEGREE", "level set degree");
    Radius = PARAM.real_value("RADIUS", "Domain radius");
    getfem::level_set ls(mesh, lsdeg);
    const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
    for (unsigned i = 0; i < lsmf.nb_dof(); ++i)
      ls.values()[i] = gmm::vect_norm2_sqr(lsmf.point_of_dof(i))-Radius*Radius;
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
    getfem::mesh_fem mf(mesh);
    std::string FEM = PARAM.string_value("FEM", "finite element method");
    mf.set_finite_element(mesh.convex_index(), getfem::fem_descriptor(FEM));
    size_type nb_dof = mf.nb_dof();
    
    // Finite element method for the rhs
    getfem::mesh_fem mf_rhs(mesh);
    std::string FEMR = PARAM.string_value("FEM_RHS", "finite element method");
    mf_rhs.set_finite_element(mesh.convex_index(),
			      getfem::fem_descriptor(FEMR));
    size_type nb_dof_rhs = mf_rhs.nb_dof();
    cout << "nb_dof_rhs = " << nb_dof_rhs << endl;
    
    // Finite element method for the multipliers
    getfem::mesh_fem mf_mult(mesh);
    std::string FEMM = PARAM.string_value("FEM_MULT", "fem for multipliers");
    mf_mult.set_finite_element(mesh.convex_index(),
			       getfem::fem_descriptor(FEMM));
    size_type nb_dof_mult = mf_mult.nb_dof();
    cout << "nb_dof_mult = " << nb_dof_mult << endl;
    

    // Partial Finite element methods
    getfem::partial_mesh_fem mf2(mf);
    dal::bit_vector kept_dof = select_dofs_from_im(mf, mim);
    mf2.adapt(kept_dof);
    nb_dof = mf2.nb_dof();

    getfem::partial_mesh_fem mf_mult2(mf_mult);
    dal::bit_vector kept_dof_mult = select_dofs_from_im(mf_mult, mimbound,N-1);
    mf_mult2.adapt(kept_dof_mult);
    nb_dof_mult = mf_mult2.nb_dof();


    // Stiffness matrix for the Poisson problem
    sparse_matrix K(nb_dof, nb_dof);
    getfem::asm_stiffness_matrix_for_homogeneous_laplacian(K, mim, mf2);

    

    // Mass matrix on the boundary
    sparse_matrix B(nb_dof_mult, nb_dof);
    getfem::asm_mass_matrix(B, mimbound, mf_mult2, mf2);

    // rhs
    plain_vector F(nb_dof), FD(nb_dof_rhs);
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      FD[i] = rhs(mf_rhs.point_of_dof(i));
    getfem::asm_source_term(F, mim, mf2, mf_rhs, FD);

    // Tests
    test_mim(mim, mf_rhs, false);
    test_mim(mimbound, mf_rhs, true);

    // global system
    sparse_matrix A(nb_dof + nb_dof_mult,
		    nb_dof + nb_dof_mult);
    gmm::sub_interval II(0, nb_dof), JJ(nb_dof, nb_dof_mult);
    gmm::copy(K, gmm::sub_matrix(A, II));
    gmm::copy(gmm::transposed(B), gmm::sub_matrix(A, II, JJ));
    gmm::copy(B, gmm::sub_matrix(A, JJ, II));
    plain_vector gF(nb_dof + nb_dof_mult);
    gmm::copy(F, gmm::sub_vector(gF, II));

    double rcond; 
    plain_vector XX(nb_dof + nb_dof_mult);
    SuperLU_solve(A, XX, gF, rcond);
    cout << "condition number : " << 1.0 / rcond << endl;

    plain_vector U(nb_dof);
    gmm::copy(gmm::sub_vector(XX, II), U);


    // interpolation of the solution on mf_rhs
    plain_vector Uint(nb_dof_rhs);
    getfem::interpolation(mf2, mf_rhs, U, Uint);
    plain_vector Vint(nb_dof_rhs);
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      Vint[i] = u_exact(mf_rhs.point_of_dof(i));

    // computation of max error.
    double errmax = 0.0;
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      if (gmm::vect_norm2(mf_rhs.point_of_dof(i)) < Radius)
	errmax = std::max(errmax, gmm::abs(Uint[i]-Vint[i]));
    cout << "Linfty error: " << errmax << endl;
    cout << "L2 error: " << getfem::asm_L2_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	 << endl;
    cout << "H1 error: " << getfem::asm_H1_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	 << endl;
      
    // export de la solution au format vtk.
    getfem::vtk_export exp("xfem_contact.vtk", (2==1));
    exp.exporting(mf2); 
    exp.write_point_data(mf2, U, "solution");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d xfem_contact.vtk  -f WarpScalar -m BandedSurfaceMap -m Outline\n";
    
    
    
    

//     // rhs for the problem on multipliers
//     plain_vector X(nb_real_dof), FF(nb_dof_mult);
//     SuperLU_solve(K2, X, F2, rcond);
//     gmm::mult(B2, X, FF);

//     // cg on the problem on multipliers
//     gmm::SMatrix SK(K2, B2);
//     plain_vector L(nb_dof_mult);
//     gmm::iteration iter(1e-10);
//     gmm::cg(SK, L, FF, gmm::identity_matrix(), iter);
   
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
