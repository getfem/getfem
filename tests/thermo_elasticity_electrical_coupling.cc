/*===========================================================================

 Copyright (C) 2015-2026 Yves Renard.

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
    Plane stress conditions are assumed.
  Electric problem: The potential is prescribed to be 0V at the right
    boundary and 0.1V at the left boundary.
  Thermal problem: A thermal insulation condition is prescribed at the
    left, right, and hole boudnaries. The remaining boundary and the
    plate front and back surfaces are supposed to transfer heat by
    convection with respect to the surrounding air at 20 deg C.
  Coupling terms:
    - Joule heating: source term  1/rho ||Grad_V||^2
    - Dependance of the thermal resistivity on temperature :
      rho = rho_0(1+alpha(T-T0))
      with T0 = 20 deg C, rho_0 the resistivity at T0
      and alpha the resistivity-temperature coefficient.
    - Thermal expansion:
      stress_tensor = E/(1+nu) ( nu/(1-nu) (div(u) - 2 alpha_th DT) I
                                + (epsilon(u) - alpha_th DT I) )
      with alpha_th being the thermal expansion coefficient.
  The first two coupling terms are nonlinear ones.
*/

#include "getfem/getfem_model_solvers.h" // Include Getfem models and solvers
#include "getfem/getfem_export.h"        // Export in various format, including vtu
#include "getfem/getfem_mesher.h"        // Experimental meshing facilities
#include "getfem/getfem_generic_assembly.h"


using std::cout;
using std::endl;

/* some GetFEM types that we will be using */
using bgeot::dim_type;
using bgeot::size_type;
using bgeot::base_node;
using bgeot::base_small_vector;


int main(int argc, char *argv[]) {
  GETFEM_MPI_INIT(argc, argv);


  //
  // Physical parameters
  //

  bgeot::md_param PARAM; // Small tool which reads a parameter file
  PARAM.read_command_line(argc, argv);
  double t = PARAM.real_value("t", "Thickness of the plate (cm)"),
         E = PARAM.real_value("E", "Young Modulus (N/cm^2)"),
         nu = PARAM.real_value("nu", "Poisson ratio"),
         F = PARAM.real_value("F",
                              "Force density at the right boundary (N/cm^2)"),
         kappa = PARAM.real_value("kappa", "Thermal conductivity (W/(cm K))"),
         D = PARAM.real_value("D", "Heat transfer coefficient (W/(K cm^2))"),
         air_temp = PARAM.real_value("air_temp",
                                     "Temperature of the air in deg C"),
         alpha_th = PARAM.real_value("alpha_th",
                                     "Thermal expansion coefficient (1/K)"),
         T0 = PARAM.real_value("T0", "Reference temperature in deg C"),
         rho_0 = PARAM.real_value("rho_0", "Resistivity at T0"),
         alpha = PARAM.real_value("alpha", "Resistivity-temperature coefficient");


  //
  // Numerical parameters
  //

  double h = PARAM.real_value("h", "Approximate mesh size");
  dim_type elements_degree
             = dim_type(PARAM.int_value("elements_degree",
                                        "Degree of the finite element methods"));
  bool export_mesh
         = (PARAM.int_value("export_mesh",
                            "Export the mesh after mesh generation or not") != 0),
       solve_in_two_steps
         = (PARAM.int_value("solve_in_two_steps",
                            "Solve the elasticity pb separately or not") != 0);


  //
  // Mesh generation. Meshes can also be imported in various formats.
  //

  getfem::mesh mesh;
  getfem::pmesher_signed_distance
    mo1 = getfem::new_mesher_rectangle(base_node(0., 0.), base_node(100., 25.)),
    mo2 = getfem::new_mesher_ball(base_node(25., 12.5), 8.),
    mo3 = getfem::new_mesher_ball(base_node(50., 12.5), 8.),
    mo4 = getfem::new_mesher_ball(base_node(75., 12.5), 8.),
    mo5 = getfem::new_mesher_union(mo2, mo3, mo4),
    mo = getfem::new_mesher_setminus(mo1, mo5);

  cout << "Mesh generation" << endl;
  std::vector<getfem::base_node> fixed;
  getfem::build_mesh(mesh, mo, h, fixed, 2, -2);


  //
  // Boundary selection
  //

  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  getfem::mesh_region
    fb1 = getfem::select_faces_in_box(mesh, border_faces,
                                      base_node(1., 1.), base_node(99., 24.)),
    fb2 = getfem::select_faces_of_normal(mesh, border_faces,
                                         base_small_vector( 1., 0.), 0.01),
    fb3 = getfem::select_faces_of_normal(mesh, border_faces,
                                         base_small_vector(-1., 0.), 0.01),
    fb4 = getfem::select_faces_of_normal(mesh, border_faces,
                                         base_small_vector(0.,  1.), 0.01),
    fb5 = getfem::select_faces_of_normal(mesh, border_faces,
                                         base_small_vector(0., -1.), 0.01),
    fb6 = getfem::select_faces_in_ball(mesh, border_faces,
                                       base_node(25., 12.5), 8.+0.01*h),
    fb7 = getfem::select_faces_in_ball(mesh, border_faces,
                                       base_node(50., 12.5), 8.+0.01*h),
    fb8 = getfem::select_faces_in_ball(mesh, border_faces,
                                       base_node(75., 12.5), 8.+0.01*h);

  size_type RIGHT_BOUND=1, LEFT_BOUND=2, TOP_BOUND=3, BOTTOM_BOUND=4,
            HOLE_BOUND=5, HOLE1_BOUND=6, HOLE2_BOUND=7, HOLE3_BOUND=8;
  mesh.region( RIGHT_BOUND) = getfem::mesh_region::subtract(fb2, fb1);
  mesh.region(  LEFT_BOUND) = getfem::mesh_region::subtract(fb3, fb1);
  mesh.region(   TOP_BOUND) = getfem::mesh_region::subtract(fb4, fb1);
  mesh.region(BOTTOM_BOUND) = getfem::mesh_region::subtract(fb5, fb1);
  mesh.region(  HOLE_BOUND) = fb1;
  mesh.region( HOLE1_BOUND) = fb6;
  mesh.region( HOLE2_BOUND) = fb7;
  mesh.region( HOLE3_BOUND) = fb8;

  dal::bit_vector nn = mesh.convex_index();
  bgeot::size_type i;
  for (i << nn; i != bgeot::size_type(-1); i << nn) {
    bgeot::pconvex_structure cvs = mesh.structure_of_convex(i);
    for (bgeot::short_type f = 0; f < cvs->nb_faces(); ++f) {
      if(mesh.region(HOLE_BOUND).is_in(i, f))
        GMM_ASSERT1(
          mesh.region(HOLE1_BOUND).is_in(i, f) ||
          mesh.region(HOLE2_BOUND).is_in(i, f) ||
          mesh.region(HOLE3_BOUND).is_in(i, f),
          "Error in select region"
        );
    }
  }

  if (export_mesh) {
    getfem::vtu_export exp("mesh.vtu", false);
    exp.exporting(mesh);
    exp.write_mesh();
    cout << "\nYou can view the mesh for instance with\n";
    cout << "mayavi2 -d mesh.vtu -f ExtractEdges -m Surface\n" << endl;
  }


  //
  // Definition of finite element methods and integration method
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
  mim.set_integration_method(2*elements_degree);


  //
  // Model definition
  //

  getfem::model md;
  md.add_fem_variable("u", mfu);   // Displacement of the structure
  md.add_fem_variable("T", mft);   // Temperature
  md.add_fem_variable("V", mft);   // Electric potential
  md.add_initialized_scalar_data("t", t);

  // Membrane elastic deformation and thermal expansion
  md.add_initialized_scalar_data("E", E);
  md.add_initialized_scalar_data("nu", nu);
  md.add_initialized_scalar_data("alpha_th", alpha_th);
  md.add_initialized_scalar_data("T0", T0);
  md.add_macro("sigma", "E/(1+nu)*( nu/(1-nu)*(Div(u)-2*alpha_th*T)*Id(2)"
                        "+(Sym(Grad(u))-alpha_th*T*Id(2)) )");
  getfem::add_linear_term(md, mim, "t*sigma:Grad(Test_u)");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", elements_degree-1, LEFT_BOUND);
  md.add_initialized_scalar_data("F", F);
  getfem::add_linear_term(md, mim, "-t*F*Test_u(1)", RIGHT_BOUND);

  // Electric field
  md.add_initialized_scalar_data("rho_0", rho_0);
  md.add_initialized_scalar_data("alpha", alpha);
  md.add_macro("rho", "rho_0*(1+alpha*(T-T0))");
  getfem::add_nonlinear_term(md, mim, "t/rho * Grad(V).Grad(Test_V)");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "V", elements_degree-1, RIGHT_BOUND);
  md.add_initialized_scalar_data("DdataV", 0.1);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "V", elements_degree-1, LEFT_BOUND, "DdataV");

  // Thermal problem
  md.add_initialized_scalar_data("kappaT", kappa);
  md.add_initialized_scalar_data("D", D);
  md.add_initialized_scalar_data("T_air", air_temp);
  getfem::add_linear_term(md, mim,
                          "t*kappaT*Grad(T).Grad(Test_T) + 2*D*(T-T_air)*Test_T");
  getfem::add_linear_term(md, mim, "t*D*(T-T_air).Test_T", TOP_BOUND);
  getfem::add_linear_term(md, mim, "t*D*(T-T_air).Test_T", BOTTOM_BOUND);
  // Joule heating term
  getfem::add_nonlinear_term(md, mim, "-t/rho * Norm_sqr(Grad(V))*Test_T");


  //
  // Model solve
  //

  gmm::iteration iter(1E-9, 1, 100);
  if (solve_in_two_steps) {
    md.disable_variable("u");
    cout << "First problem with " << md.nb_dof() << " dofs" << endl;
    getfem::standard_solve(md, iter);
    md.enable_variable("u");
    md.disable_variable("T");
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

  getfem::model_real_plain_vector VM(mfvm.nb_dof()), CO(mfvm.nb_dof() * 2);
  getfem::ga_local_projection // needs a discontinuous mesh_fem
    (md, mim, "sqrt(Norm_sqr(sigma)+sqr(sigma(1,2))-sigma(1,1)*sigma(2,2))",
     mfvm, VM);
  getfem::ga_interpolation_Lagrange_fem(md, "-t/rho * Grad(V)", mfvm, CO);
  
  getfem::vtu_export exp("displacement_with_von_mises.vtu", false);
  exp.exporting(mfu);
  exp.write_point_data(mfu, md.real_variable("u"), "elastostatic displacement");
  exp.write_point_data(mfvm, VM, "Von Mises stress");
  cout << "\nYou can view solutions with for instance:\n\nmayavi2 "
    "-d displacement_with_von_mises.vtu -f WarpVector -m Surface\n" << endl;

  getfem::vtu_export exp2("temperature.vtu", false);
  exp2.exporting(mft);
  exp2.write_point_data(mft, md.real_variable("T"), "Temperature");
  cout << "mayavi2 -d temperature.vtu -f WarpScalar -m Surface\n" << endl;

  getfem::vtu_export exp3("electric_potential.vtu", false);
  exp3.exporting(mft);
  exp3.write_point_data(mft, md.real_variable("V"), "Electric potential");
  cout << "mayavi2 -d electric_potential.vtu -f WarpScalar -m Surface\n"
       << endl;

  cout << "L2 norm of temperature = "
       << getfem::asm_L2_norm(mim, mft, md.real_variable("T")) << endl;

  GETFEM_MPI_FINALIZE;

  return 0; 
}
