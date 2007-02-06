// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**
 * Goal : scalar Dirichlet problem with Xfem.
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
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_import.h"
#include "gmm/gmm.h"

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
int u_version;
double u_alpha = 1.5;
double u_B = 20;
double u_n = 7.0;
double dtheta = M_PI/36;
double u_exact(const base_node &p) {
  double norm_sqr = gmm::vect_norm2_sqr(p);
  switch (u_version) {
    case 0: {
      double sum = std::accumulate(p.begin(), p.end(), double(0));
      
      return 5.0 * sin(sum) * (norm_sqr - Radius*Radius);
    }
    case 1: {
      double r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
      return Radius*Radius - r*r *(1+A*(1.0 + sin(n*T)));
    }
    case 2: {
      double r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n,B=u_B;
      return (Radius*Radius - r*r *(1+A*(1.0 + sin(n*T)))) * cos(B*r);      
    }
  }
  GMM_ASSERT1(false, "Invalid exact solution");
}

double g_exact(const base_node &p) {
  double norm_sqr = gmm::vect_norm2_sqr(p);
  switch (u_version) {
    case 0: {
      double sum = std::accumulate(p.begin(), p.end(), double(0));
      if (norm_sqr < 1e-10) norm_sqr = 1e-10;
      return 5.0 * (sum * cos(sum) * (norm_sqr - Radius*Radius)
		    + 2.0 * norm_sqr * sin(sum)) / sqrt(norm_sqr);
    }
    case 1: {
      double r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
      return - sqrt(r*r*pow(2.0*sin(T)+2.0*sin(T)*A+2.0*sin(T)*A*sin(n*T)+cos(T)*A*cos(n*T)*n,2.0)+r*r*pow(-2.0*cos(T)-2.0*cos(T)*A-2.0*cos(T)*A*sin(n*T)+sin(T)*A*cos(n*T)*n,2.0));
    } 
    case 2: {
      double R=Radius, r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n,B=u_B;
      if (gmm::abs(r) < 1e-10) r = 1e-10;
      return -(4.0*r*cos(B*r)+8.0*r*cos(B*r)*A+8.0*r*cos(B*r)*A*A+2.0*sin(B*r)*B*R*R-2.0*sin(B*r)*B*r*r+r*A*A*pow(cos(n*T),2.0)*n*n*cos(B*r)+8.0*r*cos(B*r)*A*A*sin(n*T)-4.0*sin(B*r)*B*r*r*A*sin(n*T)+2.0*sin(B*r)*B*r*r*A*A*pow(cos(n*T),2.0)-4.0*r*cos(B*r)*A*A*pow(cos(n*T),2.0)+8.0*r*cos(B*r)*A*sin(n*T)-4.0*sin(B*r)*B*r*r*A*A*sin(n*T)-4.0*sin(B*r)*B*r*r*A*A-4.0*sin(B*r)*B*r*r*A+2.0*sin(B*r)*B*R*R*A+2.0*sin(B*r)*B*R*R*A*sin(n*T))/sqrt(A*A*pow(cos(n*T),2.0)*n*n+4.0+8.0*A+8.0*A*sin(n*T)+8.0*A*A+8.0*A*A*sin(n*T)-4.0*A*A*pow(cos(n*T),2.0));
    }
  }
  return 0;
}

double rhs(const base_node &p) {
  switch (u_version) {
    case 0: {
      double sum = std::accumulate(p.begin(), p.end(), double(0));
      double norm_sqr = gmm::vect_norm2_sqr(p);
      double N = double(gmm::vect_size(p));
      return 5.0 * (N * sin(sum) * (norm_sqr - Radius*Radius-2.0)
		    - 4.0 * sum * cos(sum));
    }
    case 1: {
      double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
      return -(-4.0-4.0*A-4.0*A*sin(n*T)+A*sin(n*T)*n*n);
    }
    case 2: {
      double R=Radius, r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n,B=u_B;
      if (gmm::abs(r) < 1e-10) r = 1e-10;
      return (4.0*r*cos(B*r)*A*sin(n*T)-5.0*sin(B*r)*B*r*r*A+4.0*r*cos(B*r)+4.0*
r*cos(B*r)*A+sin(B*r)*B*R*R-5.0*sin(B*r)*B*r*r-5.0*sin(B*r)*B*r*r*A*sin(n*T)-r*
r*r*cos(B*r)*B*B-r*A*sin(n*T)*n*n*cos(B*r)+r*cos(B*r)*B*B*R*R-r*r*r*cos(B*r)*B*
	      B*A-r*r*r*cos(B*r)*B*B*A*sin(n*T))/r;
    }
  }
  return 0;
}

double ls_value(const base_node &p) {
  switch (u_version) {
    case 0: return gmm::vect_norm2_sqr(p)-Radius*Radius;
    case 1: case 2: {
      double r=gmm::vect_norm2(p), A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
      return (r*r*(1+A*(1.0 + sin(n*T))) - Radius*Radius) / 15.0;
    }
  }
  return 0;
}

/*
 * Test procedure
 */

void test_mim(getfem::mesh_im_level_set &mim, getfem::mesh_fem &mf_rhs,
	      bool bound) {
  if (!u_version) {
    unsigned N =  mim.linked_mesh().dim();
    size_type nbdof = mf_rhs.nb_dof();
    plain_vector V(nbdof), W(1);
    std::fill(V.begin(), V.end(), 1.0);
    
    getfem::generic_assembly assem("u=data(#1); V()+=comp(Base(#1))(i).u(i);");
    assem.push_mi(mim); assem.push_mf(mf_rhs); assem.push_data(V);
    assem.push_vec(W);
    assem.assembly(getfem::mesh_region::all_convexes());
    double exact(0), R2 = Radius*Radius, R3 = R2*Radius;
    switch (N) {
      case 1: exact = bound ? 1.0 : 2.0*Radius; break;
      case 2: exact = bound ? Radius*M_PI : R2*M_PI; break;
      case 3: exact = bound ? 2.0*M_PI*R2 : 4.0*M_PI*R3/3.0; break;
      default: assert(N <= 3);
    }
    if (bound) cout << "Boundary length: "; else cout << "Area: ";
    cout << W[0] << " should be " << exact << endl;
    assert(gmm::abs(exact-W[0])/exact < 0.01); 
  }
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
    u_version = PARAM.int_value("EXACT_SOL", "Which exact solution");
    
    // Load the mesh
    getfem::mesh mesh;
    std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
    getfem::import_mesh(MESH_FILE, mesh);
    unsigned N = mesh.dim();

    // center the mesh in (0, 0).
    base_node Pmin(N), Pmax(N);
    mesh.bounding_box(Pmin, Pmax);
    Pmin += Pmax; Pmin /= -2.0;
    // Pmin[N-1] = -Pmax[N-1];
    mesh.translation(Pmin);
    scalar_type h = mesh.minimal_convex_radius_estimate();
    cout << "h = " << h << endl;

    // Level set definition
    unsigned lsdeg = PARAM.int_value("LEVEL_SET_DEGREE", "level set degree");
    Radius = PARAM.real_value("RADIUS", "Domain radius");
    getfem::level_set ls(mesh, lsdeg);
    getfem::level_set lsup(mesh, lsdeg, true), lsdown(mesh, lsdeg, true);
    const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
    for (unsigned i = 0; i < lsmf.nb_dof(); ++i) {
      lsup.values()[i] = lsdown.values()[i] = ls.values()[i]
	= ls_value(lsmf.point_of_dof(i));
      lsdown.values(1)[i] = lsmf.point_of_dof(i)[1];
      lsup.values(1)[i] = -lsmf.point_of_dof(i)[1];
    }

    scalar_type simplify_rate = std::min(0.03, 0.05 * sqrt(h));
    cout << "simplify rate = " << simplify_rate << endl;
    
    ls.simplify(simplify_rate); lsup.simplify(simplify_rate); lsdown.simplify(simplify_rate); 
    getfem::mesh_level_set mls(mesh), mlsup(mesh), mlsdown(mesh);
    mls.add_level_set(ls);
    mls.adapt();
    mlsup.add_level_set(lsup);
    mlsup.adapt();
    mlsdown.add_level_set(lsdown);
    mlsdown.adapt();
    
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


    // Integration methods on the boudary
    int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
    getfem::mesh_im_level_set mimbounddown(mlsdown, intbound,
					   getfem::int_method_descriptor(IMS));
    mimbounddown.set_integration_method(mesh.convex_index(),
					getfem::int_method_descriptor(IM));
    mimbounddown.adapt();
    getfem::mesh_im_level_set mimboundup(mlsup, intbound,
					 getfem::int_method_descriptor(IMS));
    mimboundup.set_integration_method(mesh.convex_index(),
				      getfem::int_method_descriptor(IM));
    mimboundup.adapt();
    
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
      = select_dofs_from_im(pre_mf_mult, mimbounddown, N-1);
    mf_mult.adapt(kept_dof_mult, rejected_elt);
    size_type nb_dof_mult = mf_mult.nb_dof();
    cout << "nb_dof_mult = " << nb_dof_mult << endl;

 
    // Mass matrix on the boundary
    sparse_matrix B2(mf_rhs.nb_dof(), nb_dof);
    getfem::asm_mass_matrix(B2, mimboundup, mf_rhs, mf);

    sparse_matrix B(nb_dof_mult, nb_dof);
    getfem::asm_mass_matrix(B, mimbounddown, mf_mult, mf);

    // Tests
    test_mim(mim, mf_rhs, false);
    test_mim(mimbounddown, mf_rhs, true);

    // Brick system
    getfem::mdbrick_generic_elliptic<> brick_laplacian(mim, mf);
    
    getfem::mdbrick_source_term<> brick_volumic_rhs(brick_laplacian);
    plain_vector F(nb_dof_rhs);
    getfem::interpolation_function(mf_rhs, F, rhs);
    brick_volumic_rhs.source_term().set(mf_rhs, F);

    // Neumann condition
    getfem::interpolation_function(mf_rhs, F, g_exact);
    plain_vector R(nb_dof);
    gmm::mult(gmm::transposed(B2), F, R);
    brick_volumic_rhs.set_auxF(R);

    // Dirichlet condition
    getfem::mdbrick_constraint<> brick_constraint(brick_volumic_rhs);
    brick_constraint.set_constraints(B, plain_vector(nb_dof_mult));
    brick_constraint.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);

    getfem::mdbrick_abstract<> *final_brick = &brick_constraint;
    
    // Solving the problem
    cout << "Total number of unknown: " << final_brick->nb_dof() << endl;
    getfem::standard_model_state MS(*final_brick);
    gmm::iteration iter(1e-9, 1, 40000);
    getfem::standard_solve(MS, *final_brick, iter);
    plain_vector U(nb_dof);
    gmm::copy(brick_laplacian.get_solution(MS), U);
    plain_vector LAMBDA(nb_dof_mult);
    gmm::copy(brick_constraint.get_mult(MS), LAMBDA);

    // interpolation of the solution on mf_rhs
    plain_vector Uint(nb_dof_rhs), Vint(nb_dof_rhs), Eint(nb_dof_rhs);
    getfem::interpolation(mf, mf_rhs, U, Uint);
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      { Vint[i] = u_exact(mf_rhs.point_of_dof(i)); Eint[i] = gmm::abs(Uint[i] - Vint[i]); }
    

    // computation of max error.
    double errmax = 0.0;
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      if (ls_value(mf_rhs.point_of_dof(i)) <= 0.0) 
	errmax = std::max(errmax, Eint[i]);
      else Eint[i] = 0.0;
    cout << "Linfty error: " << errmax << endl;
    cout << "L2 error: " << getfem::asm_L2_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	 << endl;
    cout << "H1 error: " << getfem::asm_H1_dist(mim,mf_rhs,Uint,mf_rhs,Vint)
	 << endl;
    
    // exporting solution in vtk format.
    {
      getfem::vtk_export exp("xfem_dirichlet.vtk", (2==1));
      exp.exporting(mf); 
      exp.write_point_data(mf, U, "solution");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d xfem_dirichlet.vtk -f WarpScalar -m BandedSurfaceMap "
	"-m Outline\n";
    }
    // exporting error in vtk format.
    {
      getfem::vtk_export exp("xfem_dirichlet_error.vtk", (2==1));
      exp.exporting(mf_rhs); 
      exp.write_point_data(mf_rhs, Eint, "error");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d xfem_dirichlet_error.vtk -f WarpScalar -m BandedSurfaceMap "
	"-m Outline\n";
    }
    // exporting multipliers in vtk format.
    {
      getfem::vtk_export exp("xfem_dirichlet_mult.vtk", (2==1));
      exp.exporting(mf_mult); 
      exp.write_point_data(mf_mult, LAMBDA, "multipliers");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d xfem_dirichlet_mult.vtk -f WarpScalar -m BandedSurfaceMap "
	"-m Outline\n";
    }

    lsmf.write_to_file("xfem_dirichlet_ls.mf", true);
    gmm::vecsave("xfem_dirichlet_ls.U", ls.values());

    unsigned nrefine = mf.linked_mesh().convex_index().card() < 200 ? 32 : 4;
    if (1)
    {
      cout << "saving the slice solution for matlab\n";
      getfem::stored_mesh_slice sl, sl0,sll;
    
      
      getfem::mesh_slicer slicer(mf.linked_mesh());
      getfem::slicer_build_stored_mesh_slice sbuild(sl);
      getfem::mesh_slice_cv_dof_data<plain_vector> mfU(mf,U);
      getfem::slicer_isovalues iso(mfU, 0.0, 0);
      getfem::slicer_build_stored_mesh_slice sbuild0(sl0);

      slicer.push_back_action(sbuild);  // full slice in sl
      slicer.push_back_action(iso);     // extract isosurface 0
      slicer.push_back_action(sbuild0); // store it into sl0
      slicer.exec(nrefine, mf.convex_index());

      getfem::mesh_slicer slicer2(mf.linked_mesh());
      getfem::mesh_slice_cv_dof_data<plain_vector> 
	mfL(ls.get_mesh_fem(), ls.values());
      getfem::slicer_isovalues iso2(mfL, 0.0, 0);
      getfem::slicer_build_stored_mesh_slice sbuildl(sll);
      slicer2.push_back_action(iso2);     // extract isosurface 0
      slicer2.push_back_action(sbuildl); // store it into sl0
      slicer2.exec(nrefine, mf.convex_index());
      
      
      sl.write_to_file("xfem_dirichlet.sl", true);
      sl0.write_to_file("xfem_dirichlet.sl0", true);
      sll.write_to_file("xfem_dirichlet.sll", true);
      plain_vector UU(sl.nb_points()), LL(sll.nb_points()); 
      sl.interpolate(mf, U, UU);
      gmm::vecsave("xfem_dirichlet.slU", UU);
      sll.interpolate(mf_mult, LAMBDA, LL);
      gmm::vecsave("xfem_dirichlet.slL", LL);
    }
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
