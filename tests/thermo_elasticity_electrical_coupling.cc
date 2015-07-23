/*===========================================================================

 Copyright (C) 2015-2015 Yves Renard.

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

/**@file thermo_elasticity_electrical_coupling.cc

   @brief Deformation of a plate under the coupling of thermal, elasticity, and
          electric effects.

     ______________________________________
   /|         __       __       __         |->
   /|        /  \     /  \     /  \        |->
   /|       |    |   |    |   |    |       |-> F
   /|        \__/     \__/     \__/        |->
   /|______________________________________|->
    

  Elastic problem: The plate is clamped at rhe left boundary and a
    traction density of force F is prescribed at the right boundary.
  Electric problem: The potential is prescribed to be 0V at the right
    boundary and 0.1V at the left boundary.
  Thermal problem: A thermal insulation condition is prescribed at the
    left and hole boudnaries. The remaining boundary and the plate itself
    is supposed to be submitted to an heat transfert with respect to the
    air at 20oC.
  Coupling terms:
    - Joule heating: source term  sigma|Grad_V|^2
    - Dependance of the thermal conductivity in temperature :
      sigma = 1/(rho_0(1+alpha(theta-T0)))
      with T0 = 20oC, rho_0 the resistance temperature coefficient at T0
      and alpha the second resistance temperature coefficient.
    - Thermal expansion:
      stress_tensor = clambdastar div(u) I + 2 cmu epsilon(u) - beta theta I
      with beta = alpha_th E/(1-2nu), alpha_th being the thermal
      expansion coefficient.
  The first two coupling terms are nonlinear ones.
*/

#include "getfem/getfem_model_solvers.h" // Include Getfem models and solvers
#include "getfem/getfem_export.h"  // Export in various format, including vtk
#include "gmm/gmm.h"                  // Include Gmm matrix interface library
#include "getfem/getfem_mesher.h"     // Experimental meshing facilities
#include "getfem/getfem_generic_assembly.h"


using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some Getfem++ types that we will be using */
using bgeot::size_type;
using bgeot::base_node;
using bgeot::base_small_vector;
typedef getfem::model_real_plain_vector plain_vector;


int main(int argc, char *argv[]) {
  bgeot::md_param PARAM; // Small tool which reads a parameter file
  PARAM.read_command_line(argc, argv);


  //
  // Physical parameters
  //
  double epsilon = PARAM.real_value("epsilon", "Thickness of the plate (cm)");
  double E = PARAM.real_value("E", "Young Modulus (N/cm^2)");
  double nu = PARAM.real_value("nu", "Poisson ratio");
  double clambda = E*nu/((1+nu)*(1-2*nu)); // First Lame coefficient (N/cm^2)
  double cmu = E/(2*(1+nu));               // Second Lame coefficient (N/cm^2)
  double clambdastar = 2*clambda*cmu/(clambda+2*cmu); // Lame coefficient
                                                 // for Plane stress (N/cm^2)
  double F = PARAM.real_value("F",
                              "Force density at the right boundary (N/cm^2)");
  double kappa = PARAM.real_value("kappa", "Thermal conductivity (W/(cm K))");
  double D = PARAM.real_value("D", "Heat transfert coefficient (W/(K cm^2))");
  double air_temp = PARAM.real_value("air_temp",
                                     "Temperature of the air in oC");
  double alpha_th = PARAM.real_value("alpha_th",
                                     "Thermal expansion coefficient (/K)");
  double T0 = PARAM.real_value("T0", "Reference temperature in oC");
  double rho_0 = PARAM.real_value("rho_0",
                                  "Resistance temperature coefficient at T0");
  double alpha = PARAM.real_value("alpha",
                                  "Second resistance temperature coefficient");

  //
  // Numerical parameters
  //
  double h = PARAM.real_value("h", "Approximate mesh size");
  bgeot::dim_type elements_degree = 
    bgeot::dim_type(PARAM.int_value("elements_degree",
                                    "Degree of the finite element methods"));
  bool export_mesh =
    (PARAM.int_value("export_mesh",
                     "Draw the mesh after mesh generation or not") != 0);
  bool solve_in_two_steps =
    (PARAM.int_value("solve_in_two_steps",
                     "Solve the elasticity pb separately or not") != 0);

  //
  // Mesh generation. Meshes can also been imported from several formats.
  //
  getfem::mesh mesh;
  getfem::mesher_rectangle mo1(base_node(0., 0.), base_node(100., 25.));
  getfem::mesher_ball mo2(base_node(25., 12.5), 8.);
  getfem::mesher_ball mo3(base_node(50., 12.5), 8.);
  getfem::mesher_ball mo4(base_node(75., 12.5), 8.);
  getfem::mesher_union mo5(mo2, mo3, mo4);
  getfem::mesher_setminus mo(mo1, mo5);

  cout << "Mesh generation" << endl;
  std::vector<getfem::base_node> fixed;
  getfem::build_mesh(mesh, mo, h, fixed, 2, -2);

  //
  // Boundary selection.
  //
  
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  getfem::mesh_region fb1
    = getfem::select_faces_in_box(mesh, border_faces, base_node(1., 1.),
                                 base_node(99., 24.));
  getfem::mesh_region fb2
    = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector( 1., 0.), 0.01);
  getfem::mesh_region fb3
    = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(-1., 0.), 0.01);
  getfem::mesh_region fb4
    = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(0.,  1.), 0.01);
  getfem::mesh_region fb5
    = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(0., -1.), 0.01);

  size_type RIGHT_BOUND = 1, LEFT_BOUND = 2, TOP_BOUND = 3, BOTTOM_BOUND = 4;
  mesh.region( RIGHT_BOUND) = getfem::mesh_region::subtract(fb2, fb1);
  mesh.region(  LEFT_BOUND) = getfem::mesh_region::subtract(fb3, fb1);
  mesh.region(   TOP_BOUND) = getfem::mesh_region::subtract(fb4, fb1);
  mesh.region(BOTTOM_BOUND) = getfem::mesh_region::subtract(fb5, fb1);
 
  if (export_mesh) {
    getfem::vtk_export exp("mesh.vtk", false);
    exp.exporting(mesh);
    exp.write_mesh();
    cout << "\nYou can view the mesh for instance with\n";
    cout << "mayavi2 -d mesh.vtk -f ExtractEdges -m Surface\n" << endl;
  }

  //
  // Definition of finite elements methods and integration method
  //

  getfem::mesh_fem mfu(mesh, 2); // Finite element for the elastic displacement
  mfu.set_classical_finite_element(elements_degree);
  getfem::mesh_fem mft(mesh, 1); // Finite element for temperature
                                 // and electrical field
  mft.set_classical_finite_element(elements_degree);
  getfem::mesh_fem mfvm(mesh, 1); // Finite element for Von Mises stress
                                  // interpolation
  mfvm.set_classical_discontinuous_finite_element(elements_degree);

  getfem::mesh_im  mim(mesh);     // Integration method
  mim.set_integration_method(bgeot::dim_type(gmm::sqr(elements_degree)));


  //
  // Model definition
  //

  getfem::model md;
  md.add_fem_variable("u", mfu);       // Displacement of the structure
  md.add_fem_variable("theta", mft);   // Temperature
  md.add_fem_variable("V", mft);       // Electric potential

  // Membrane elastic deformation
  md.add_initialized_scalar_data("cmu", cmu);
  md.add_initialized_scalar_data("clambdastar", clambdastar);
  getfem::add_isotropic_linearized_elasticity_brick
    (md, mim, "u", "clambdastar", "cmu");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", bgeot::dim_type(elements_degree-1), LEFT_BOUND);
  md.add_initialized_fixed_size_data("Fdata", base_small_vector(F*epsilon,0.));
  getfem::add_source_term_brick(md, mim, "u", "Fdata", RIGHT_BOUND);

  // Electrical field
  std::string sigmaeps = "(eps/(rho_0*(1+alpha*(theta-T0))))";
  md.add_initialized_scalar_data("eps", epsilon);
  md.add_initialized_scalar_data("rho_0", rho_0);
  md.add_initialized_scalar_data("alpha", alpha);
  md.add_initialized_scalar_data("T0", T0);
  getfem::add_nonlinear_generic_assembly_brick
    (md, mim, sigmaeps+"*(Grad_V.Grad_Test_V)");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "V", bgeot::dim_type(elements_degree-1), RIGHT_BOUND);
  md.add_initialized_scalar_data("DdataV", 0.1);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "V", bgeot::dim_type(elements_degree-1), LEFT_BOUND, "DdataV");
  
  // Thermal problem
  md.add_initialized_scalar_data("kaeps", kappa*epsilon);
  getfem::add_generic_elliptic_brick(md, mim, "theta", "kaeps");
  md.add_initialized_scalar_data("D2", D*2);
  md.add_initialized_scalar_data("D2airt", air_temp*D*2);
  getfem::add_mass_brick(md, mim, "theta", "D2");
  getfem::add_source_term_brick(md, mim, "theta", "D2airt");
  md.add_initialized_scalar_data("Deps", D/epsilon);
  md.add_initialized_scalar_data("Depsairt", air_temp*D/epsilon);
  getfem::add_Fourier_Robin_brick(md, mim, "theta", "Deps", TOP_BOUND);
  getfem::add_source_term_brick(md, mim, "theta", "Depsairt", TOP_BOUND);
  getfem::add_Fourier_Robin_brick(md, mim, "theta", "Deps", BOTTOM_BOUND);
  getfem::add_source_term_brick(md, mim, "theta", "Depsairt", BOTTOM_BOUND);



  // Joule heating term
  getfem::add_nonlinear_generic_assembly_brick
     (md, mim, "-"+sigmaeps+"*Norm_sqr(Grad_V)*Test_theta");

  // Thermal expansion term
  md.add_initialized_scalar_data("beta", alpha_th*E/(1-2*nu));
  getfem::add_linear_generic_assembly_brick
    (md, mim, "beta*(T0-theta)*Trace(Grad_Test_u)");


  //
  // Model solve
  //
  gmm::iteration iter(1E-9, 1, 100);

  if (solve_in_two_steps) {
    md.disable_variable("u");
    cout << "First problem with " << md.nb_dof() << " dofs" << endl;
    getfem::standard_solve(md, iter);
    md.enable_variable("u");
    md.disable_variable("theta");
    md.disable_variable("V");
    cout << "Second problem with " << md.nb_dof() << " dofs" << endl;
    iter.init();
    getfem::standard_solve(md, iter);
  } else {
    cout << "Global problem with " << md.nb_dof() << " dofs" << endl;
    getfem::standard_solve(md, iter);
  }


  //
  // Solution export
  //
  plain_vector U(mfu.nb_dof()); gmm::copy(md.real_variable("u"), U);
  plain_vector V(mft.nb_dof()); gmm::copy(md.real_variable("V"), V);
  plain_vector THETA(mft.nb_dof()); gmm::copy(md.real_variable("theta"),THETA);
  plain_vector VM(mfvm.nb_dof());
  getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
    (md, "u", "clambdastar", "cmu", mfvm, VM, false);
  plain_vector CO(mfvm.nb_dof() * 2);
  getfem::ga_interpolation_Lagrange_fem(md, "-"+sigmaeps+"*Grad_V",  mfvm, CO);
  
  getfem::vtk_export exp("displacement_with_von_mises.vtk", false);
  exp.exporting(mfu);
  exp.write_point_data(mfu, U, "elastostatic displacement");
  exp.write_point_data(mfvm, VM, "Von Mises stress");
  cout << "\nYou can view solutions with for instance:\n\nmayavi2 "
    "-d displacement_with_von_mises.vtk -f WarpVector -m Surface\n" << endl;
  
  getfem::vtk_export exp2("temperature.vtk", false);
  exp2.exporting(mft);
  exp2.write_point_data(mft, THETA, "Temperature");
  cout << "mayavi2 -d temperature.vtk -f WarpScalar -m Surface\n" << endl;

  getfem::vtk_export exp3("electric_potential.vtk", false);
  exp3.exporting(mft);
  exp3.write_point_data(mft, V, "Electric potential");
  cout << "mayavi2 -d electric_potential.vtk -f WarpScalar -m Surface\n"
       << endl;

  cout << "L2 norm of temperature = " << getfem::asm_L2_norm(mim, mft, THETA) << endl;

  return 0; 
}


