// Copyright (C) 2015-2020 Yves Renard.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 3 of the License,  or
// (at your option) any later version along with the GCC Runtime Library
// Exception either version 3.1 or (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License and GCC Runtime Library Exception for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

// Deformation of a plate under the coupling of thermal, elasticity, and
// electric effects.
//
//
//     ______________________________________
//   /|         __       __       __         |->
//   /|        /  \     /  \     /  \        |->
//   /|       |    |   |    |   |    |       |-> F
//   /|        \__/     \__/     \__/        |->
//   /|______________________________________|->
//     
//
// Elastic problem: The plate is clamped at rhe left boundary and a
//   traction density of force F is prescribed at the right boundary.
// Electric problem: The potential is prescribed to be 0V at the right
//   boundary and 0.1V at the left boundary.
// Thermal problem: A thermal insulation condition is prescribed at the
//   left and hole boudnaries. The remaining boundary and the plate itself
//   is supposed to be submitted to heat transfer with respect to the
//   air at 20°C.
// Coupling terms:
//   - Joule heating: source term  sigma|Grad_V|^2
//   - Dependance of the thermal conductivity in temperature :
//     sigma = 1/(rho_0(1+alpha(theta-T0)))
//     with T0 = 20°C, rho_0 the resistance temperature coefficient at T0
//     and alpha the second resistance temperature coefficient.
//   - Thermal expansion:
//     stress_tensor = clambdastar div(u) I + 2 cmu epsilon(u) - beta theta I
//     with beta = alpha_th E/(1-2nu), alpha_th being the thermal
//     expansion coefficient.
// The first two coupling terms are nonlinear ones.


lines(0);
stacksize('max');

path = get_absolute_file_path('demo_thermo_elasticity_electrical_coupling.sce');

printf('demo thermo elasticity electrical coupling  started\n');

// trace on;

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end

gf_workspace('clear all');

//
// Physical parameters
//
epsilon = 1.;       // Thickness of the plate (cm)
E = 21E6;           // Young Modulus (N/cm^2)
nu = 0.3;           // Poisson ratio
clambda = E*nu/((1+nu)*(1-2*nu)); // First Lame coefficient (N/cm^2)
cmu = E/(2*(1+nu));               // Second Lame coefficient (N/cm^2)
clambdastar = 2*clambda*cmu/(clambda+2*cmu); // Lame coefficient for Plane stress (N/cm^2)
F = 100E2;          // Force density at the right boundary (N/cm^2)
kappa = 4.;         // Thermal conductivity (W/(cm K))
D = 10;             // Heat transfer coefficient (W/(K cm^2))
air_temp = 20;      // Temperature of the air in ??C.
alpha_th = 16.6E-6; // Thermal expansion coefficient (/K).
T0 = 20;            // Reference temperature in ??C.
rho_0 = 1.754E-8;   // Resistance temperature coefficient at T0 = 20??C
alpha = 0.0039;     // Second resistance temperature coefficient.

//
// Numerical parameters
//
h = 2;                      // Approximate mesh size
elements_degree = 2;        // Degree of the finite element methods
draw_mesh = %t;             // Draw the mesh after mesh generation or not
solve_in_two_steps = %t;    // Solve the elasticity problem separately or not

//
// Mesh generation. Meshes can also been imported from several formats.
//
mo1 = gf_mesher_object('rectangle', [0 0], [100 25]);
mo2 = gf_mesher_object('ball', [25 12.5], 8);
mo3 = gf_mesher_object('ball', [50 12.5], 8);
mo4 = gf_mesher_object('ball', [75 12.5], 8);
mo5 = gf_mesher_object('union', mo2, mo3, mo4);
mo  = gf_mesher_object('set minus', mo1, mo5);

disp('Mesh generation');
gf_util('trace level', 2);   // No trace for mesh generation
mesh = gf_mesh('generate', mo, h, 2);

//
// Boundary selection
//
fb1 = gf_mesh_get(mesh, 'outer faces in box', [1 1], [99 24]);  // Boundary of the holes
fb2 = gf_mesh_get(mesh, 'outer faces with direction', [ 1 0], 0.01); // Right boundary
fb3 = gf_mesh_get(mesh, 'outer faces with direction', [-1 0], 0.01); // Left boundary
fb4 = gf_mesh_get(mesh, 'outer faces with direction', [0  1], 0.01); // Top boundary
fb5 = gf_mesh_get(mesh, 'outer faces with direction', [0 -1], 0.01); // Bottom boundary

RIGHT_BOUND = 1; LEFT_BOUND = 2; TOP_BOUND = 3; BOTTOM_BOUND = 4; HOLE_BOUND = 5; 
gf_mesh_set(mesh, 'region',  RIGHT_BOUND, fb2);
gf_mesh_set(mesh, 'region',   LEFT_BOUND, fb3);
gf_mesh_set(mesh, 'region',    TOP_BOUND, fb4);
gf_mesh_set(mesh, 'region', BOTTOM_BOUND, fb5);
gf_mesh_set(mesh, 'region',   HOLE_BOUND, fb1);
gf_mesh_set(mesh, 'region subtract',  RIGHT_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region subtract',   LEFT_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region subtract',    TOP_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region subtract', BOTTOM_BOUND, HOLE_BOUND);

if (draw_mesh)
  scf(1);
  gf_plot_mesh(mesh, 'refine', 8, 'curved', 'on', 'regions', [RIGHT_BOUND LEFT_BOUND TOP_BOUND BOTTOM_BOUND]);
  title('Mesh');
  sleep(1000);
end

//
// Definition of finite element methods and integration method
//

mfu = gf_mesh_fem(mesh, 2); // Finite element for the elastic displacement
gf_mesh_fem_set(mfu, 'classical fem', elements_degree);
mft = gf_mesh_fem(mesh, 1); // Finite element for the temperature and the electrical field
gf_mesh_fem_set(mft, 'classical fem', elements_degree);
mfvm = gf_mesh_fem(mesh, 1); // Finite element for Von Mises stress interpolation
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', elements_degree-1);
mim = gf_mesh_im(mesh, elements_degree*2);   // Integration method


//
// Model definition
//

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);       // Displacement of the structure
gf_model_set(md, 'add fem variable', 'theta', mft);   // Temperature
gf_model_set(md, 'add fem variable', 'V', mft);       // Electric potential

// Membrane elastic deformation
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambdastar', [clambdastar]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'clambdastar', 'cmu');

gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', elements_degree-1, LEFT_BOUND);
gf_model_set(md, 'add initialized data', 'Fdata', [F*epsilon, 0]);
gf_model_set(md, 'add source term brick', mim, 'u', 'Fdata', RIGHT_BOUND);

// Electrical field
sigmaps = '(eps/(rho_0*(1+alpha*(theta-T0))))';
gf_model_set(md, 'add initialized data', 'eps', [epsilon]);
gf_model_set(md, 'add initialized data', 'rho_0', [rho_0]);
gf_model_set(md, 'add initialized data', 'alpha', [alpha]);
gf_model_set(md, 'add initialized data', 'T0', [T0]);
gf_model_set(md, 'add nonlinear generic assembly brick', mim, sigmaeps+'*(Grad_V.Grad_Test_V)');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'V', elements_degree-1, RIGHT_BOUND);
gf_model_set(md, 'add initialized data', 'DdataV', [0.1]);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'V', elements_degree-1, LEFT_BOUND, 'DdataV');

// Thermal problem
gf_model_set(md, 'add initialized data', 'kaeps', [kappa*epsilon]);
gf_model_set(md, 'add generic elliptic brick', mim, 'theta', 'kaeps');
gf_model_set(md, 'add initialized data', 'D2', [D*2]);
gf_model_set(md, 'add initialized data', 'D2airt', [air_temp*D*2]);
gf_model_set(md, 'add mass brick', mim, 'theta', 'D2');
gf_model_set(md, 'add source term brick', mim, 'theta', 'D2airt');
gf_model_set(md, 'add initialized data', 'Deps', [D/epsilon]);
gf_model_set(md, 'add initialized data', 'Depsairt', [air_temp*D/epsilon]);
gf_model_set(md, 'add Fourier Robin brick', mim, 'theta', 'Deps', TOP_BOUND);
gf_model_set(md, 'add source term brick', mim, 'theta', 'Depsairt', TOP_BOUND);
gf_model_set(md, 'add Fourier Robin brick', mim, 'theta', 'Deps', BOTTOM_BOUND);
gf_model_set(md, 'add source term brick', mim, 'theta', 'Depsairt', BOTTOM_BOUND);

// Joule heating term
gf_model_set(md, 'add nonlinear generic assembly brick', mim, '-'+sigmaeps+'*Norm_sqr(Grad_V)*Test_theta');

// Thermal expansion term
gf_model_set(md, 'add initialized data', 'beta', [alpha_th*E/(1-2*nu)]);
gf_model_set(md, 'add linear generic assembly brick', mim, 'beta*(T0-theta)*Trace(Grad_Test_u)');


//
// Model solve and solution plot
//
if (solve_in_two_steps) 
  gf_model_set(md, 'disable variable', 'u');
  disp(sprintf('First problem with %d dofs', gf_model_get(md, 'nbdof')));
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
  gf_model_set(md, 'enable variable', 'u');
  gf_model_set(md, 'disable variable', 'theta');
  gf_model_set(md, 'disable variable', 'V');
  disp(sprintf('Second problem with %d dofs', gf_model_get(md, 'nbdof')));
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
else
  disp(sprintf('Second problem with %d dofs', gf_model_get(md, 'nbdof')));
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
end
  
  
U = gf_model_get(md, 'variable', 'u');
V = gf_model_get(md, 'variable', 'V');
THETA = gf_model_get(md, 'variable', 'theta');
VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u', 'clambdastar', 'cmu', mfvm);
CO = matrix(gf_model_get(md, 'interpolation', '-'+sigmaeps+'*Grad_V', mfvm), [2 gf_mesh_fem_get(mfvm, 'nbdof')]);
    
hh = scf(2);
hh.color_map = jetcolormap(255);
subplot(3,1,1);
gf_plot(mfvm, VM, 'mesh', 'off', 'deformed_mesh','off', 'deformation', U, 'deformation_mf', mfu, 'deformation_scale', 100, 'refine', 8); colorbar(min(VM),max(VM));
title('Von Mises stress in N/cm^2 (on the deformed configuration, scale factor x100)');
subplot(3,1,2);
drawlater;
gf_plot(mft, V, 'mesh', 'off', 'deformed_mesh','off', 'deformation', U, 'deformation_mf', mfu, 'deformation_scale', 100, 'refine', 8); colorbar(min(V),max(V));
// gf_plot(mfvm, CO, 'quiver', 'on', 'quiver_density', 0.1, 'mesh', 'off', 'deformed_mesh','off', 'deformation_mf', mfu, 'deformation', U, 'deformation_scale', 100, 'refine', 8);
title('Electric potential in Volt (on the deformed configuration, scale factor x100)');
drawnow;
subplot(3,1,3);
gf_plot(mft, THETA, 'mesh', 'off', 'deformed_mesh','off', 'deformation', U, 'deformation_mf', mfu, 'deformation_scale', 100, 'refine', 8); colorbar(min(THETA),max(THETA));
title('Temperature in °C (on the deformed configuration, scale factor x100)');


printf('demo thermo elasticity electrical coupling terminated\n');
