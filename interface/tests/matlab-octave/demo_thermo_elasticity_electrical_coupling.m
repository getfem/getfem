% Copyright (C) 2015-2026 Yves Renard.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 3 of the License,  or
% (at your option) any later version along with the GCC Runtime Library
% Exception either version 3.1 or (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License and GCC Runtime Library Exception for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
%

gf_workspace('clear all');
clear all;

% Deformation of a plate under the coupling of thermal, elasticity, and
% electric effects.
%
%
%     ______________________________________
%   /|         __       __       __         |->
%   /|        /  \     /  \     /  \        |->
%   /|       |    |   |    |   |    |       |-> F
%   /|        \__/     \__/     \__/        |->
%   /|______________________________________|->
%     
%
% Elastic problem: The plate is clamped at rhe left boundary and a
%   traction density of force F is prescribed at the right boundary.
%   Plane stress conditions are assumed.
% Electric problem: The potential is prescribed to be 0V at the right
%   boundary and 0.1V at the left boundary.
% Thermal problem: A thermal insulation condition is prescribed at the
%   left, right, and hole boudnaries. The remaining boundary and the
%   plate front and back surfaces are supposed to transfer heat by
%   convection with respect to the surrounding air at 20 deg C.
% Coupling terms:
%   - Joule heating: source term  1/rho ||Grad_V||^2
%   - Dependance of the thermal resistivity on temperature :
%     rho = rho_0(1+alpha(T-T0))
%     with T0 = 20 deg C, rho_0 the resistivity at T0
%     and alpha the resistivity-temperature coefficient.
%   - Thermal expansion:
%     stress_tensor = E/(1+nu) ( nu/(1-nu) (div(u) - 2 alpha_th DT) I
%                               + (epsilon(u) - alpha_th DT I) )
%     with alpha_th being the thermal expansion coefficient.
% The first two coupling terms are nonlinear ones.


%
% Physical parameters
%

t = 1.;             % Thickness of the plate (cm)
E = 21E6;           % Young Modulus (N/cm^2)
nu = 0.3;           % Poisson ratio
F = 100E2;          % Force density at the right boundary (N/cm^2)
kappa = 4.;         % Thermal conductivity (W/(cm K))
D = 10;             % Heat transfer coefficient (W/(K cm^2))
air_temp = 20;      % Temperature of the air in deg C
alpha_th = 16.6E-6; % Thermal expansion coefficient (1/K)
T0 = 20;            % Reference temperature in deg C
rho_0 = 1.754E-8;   % Resistivity at T0
alpha = 0.0039;     % Resistivity-temperature coefficient


%
% Numerical parameters
%

h = 2;                      % Approximate mesh size
elements_degree = 2;        % Degree of the finite element methods
draw_mesh = false;          % Draw the mesh after mesh generation or not
solve_in_two_steps = true;  % Solve the elasticity problem separately or not


%
% Mesh generation. Meshes can also be imported in various formats.
%

mo1 = gf_mesher_object('rectangle', [0 0], [100 25]);
mo2 = gf_mesher_object('ball', [25 12.5], 8);
mo3 = gf_mesher_object('ball', [50 12.5], 8);
mo4 = gf_mesher_object('ball', [75 12.5], 8);
mo5 = gf_mesher_object('union', mo2, mo3, mo4);
mo  = gf_mesher_object('set_minus', mo1, mo5);

disp('Mesh generation');
gf_util('trace level', 2);   % No trace for mesh generation
mesh = gf_mesh('generate', mo, h, 2);


%
% Boundary selection
%

fb1 = gf_mesh_get(mesh, 'outer_faces_in_box', [1 1], [99 24]);       % Boundary of the holes
fb2 = gf_mesh_get(mesh, 'outer_faces_with_direction', [ 1 0], 0.01); % Right boundary
fb3 = gf_mesh_get(mesh, 'outer_faces_with_direction', [-1 0], 0.01); % Left boundary
fb4 = gf_mesh_get(mesh, 'outer_faces_with_direction', [0  1], 0.01); % Top boundary
fb5 = gf_mesh_get(mesh, 'outer_faces_with_direction', [0 -1], 0.01); % Bottom boundary

RIGHT_BOUND = 1; LEFT_BOUND = 2; TOP_BOUND = 3; BOTTOM_BOUND = 4; HOLE_BOUND = 5; 
gf_mesh_set(mesh, 'region',  RIGHT_BOUND, fb2);
gf_mesh_set(mesh, 'region',   LEFT_BOUND, fb3);
gf_mesh_set(mesh, 'region',    TOP_BOUND, fb4);
gf_mesh_set(mesh, 'region', BOTTOM_BOUND, fb5);
gf_mesh_set(mesh, 'region',   HOLE_BOUND, fb1);
gf_mesh_set(mesh, 'region_subtract',  RIGHT_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region_subtract',   LEFT_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region_subtract',    TOP_BOUND, HOLE_BOUND);
gf_mesh_set(mesh, 'region_subtract', BOTTOM_BOUND, HOLE_BOUND);

if (draw_mesh)
  figure(1);
  gf_plot_mesh(mesh, 'refine', 8, 'curved', 'on', 'regions', [RIGHT_BOUND LEFT_BOUND TOP_BOUND BOTTOM_BOUND]);
  title('Mesh');
  pause(1);
end


%
% Definition of finite element methods and integration method
%

mfu = gf_mesh_fem(mesh, 2); % Finite element for the elastic displacement
gf_mesh_fem_set(mfu, 'classical_fem', elements_degree);
mft = gf_mesh_fem(mesh, 1); % Finite element for temperature and electrical field
gf_mesh_fem_set(mft, 'classical_fem', elements_degree);
mfvm = gf_mesh_fem(mesh, 1); % Finite element for Von Mises stress interpolation
gf_mesh_fem_set(mfvm, 'classical_discontinuous_fem', elements_degree);
mim = gf_mesh_im(mesh, elements_degree*2);   % Integration method


%
% Model definition
%

md=gf_model('real');
gf_model_set(md, 'add_fem_variable', 'u', mfu);   % Displacement of the structure
gf_model_set(md, 'add_fem_variable', 'T', mft);   % Temperature
gf_model_set(md, 'add_fem_variable', 'V', mft);   % Electric potential
gf_model_set(md, 'add_initialized_data', 't', t);

% Membrane elastic deformation and thermal expansion
gf_model_set(md, 'add_initialized_data', 'E', E);
gf_model_set(md, 'add_initialized_data', 'nu', nu);
gf_model_set(md, 'add_initialized_data', 'alpha_th', alpha_th);
gf_model_set(md, 'add_initialized_data', 'T0', T0);
gf_model_set(md, 'add_macro', 'sigma',...
                              ['E/(1+nu)*( nu/(1-nu)*(Div(u)-2*alpha_th*T)*Id(2)'...
                               '+(Sym(Grad(u))-alpha_th*T*Id(2)) )']);
gf_model_set(md, 'add_linear_term', mim, 't*sigma:Grad(Test_u)');
gf_model_set(md, 'add_Dirichlet_condition_with_multipliers', mim, 'u', elements_degree-1, LEFT_BOUND);
gf_model_set(md, 'add_initialized_data', 'F', F);
gf_model_set(md, 'add_linear_term', mim, '-t*F*Test_u(1)', RIGHT_BOUND);

% Electric field
gf_model_set(md, 'add_initialized_data', 'rho_0', rho_0);
gf_model_set(md, 'add_initialized_data', 'alpha', alpha);
gf_model_set(md, 'add_macro', 'rho', 'rho_0*(1+alpha*(T-T0))');
gf_model_set(md, 'add_nonlinear_term', mim, 't/rho * Grad(V).Grad(Test_V)');
gf_model_set(md, 'add_Dirichlet_condition_with_multipliers', mim, 'V', elements_degree-1, RIGHT_BOUND);
gf_model_set(md, 'add_initialized_data', 'DdataV', 0.1);
gf_model_set(md, 'add_Dirichlet_condition_with_multipliers', mim, 'V', elements_degree-1, LEFT_BOUND, 'DdataV');

% Thermal problem
gf_model_set(md, 'add_initialized_data', 'kappaT', kappa);
gf_model_set(md, 'add_initialized_data', 'D', D);
gf_model_set(md, 'add_initialized_data', 'T_air', air_temp);
gf_model_set(md, 'add_linear_term', mim, 't*kappaT*Grad(T).Grad(Test_T) + 2*D*(T-T_air)*Test_T');
gf_model_set(md, 'add_linear_term', mim, 't*D*(T-T_air).Test_T', TOP_BOUND);
gf_model_set(md, 'add_linear_term', mim, 't*D*(T-T_air).Test_T', BOTTOM_BOUND);
% Joule heating term
gf_model_set(md, 'add_nonlinear_term', mim, '-t/rho * Norm_sqr(Grad(V))*Test_T');


%
% Model solve and solution plot
%

if (solve_in_two_steps)
  gf_model_set(md, 'disable_variable', 'u');
  disp(['First problem with ', num2str(gf_model_get(md, 'nbdof')), ' dofs']);
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
  gf_model_set(md, 'enable_variable', 'u');
  gf_model_set(md, 'disable_variable', 'T');
  gf_model_set(md, 'disable_variable', 'V');
  disp(['Second problem with ', num2str(gf_model_get(md, 'nbdof')), ' dofs']);
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
else
  disp(['Global problem with ', num2str(gf_model_get(md, 'nbdof')), ' dofs']);
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', 100, 'noisy');
end
  
U = gf_model_get(md, 'variable', 'u');
VM = gf_model_get(md, 'local_projection', mim, 'sqrt(Norm_sqr(sigma)+sqr(sigma(1,2))-sigma(1,1)*sigma(2,2))', mfvm);
CO = reshape(gf_model_get(md, 'interpolation', '-t/rho * Grad(V)', mfvm),...
             [2 gf_mesh_fem_get(mfvm, 'nbdof')]);

figure(2);
subplot(3,1,1);
gf_plot(mfvm, VM, 'mesh', 'off', 'disp_options', 'off',...
                  'deformed_mesh','off', 'deformation', U,... 
                  'deformation_mf', mfu, 'deformation_scale', 100, 'refine', 8);
colorbar;
title('Von Mises stress in N/cm^2 (on the deformed configuration, scale factor x100)');
subplot(3,1,2);
hold on;
gf_plot(mft, gf_model_get(md, 'variable', 'V'), 'mesh', 'off', 'disp_options', 'off',...
                                                'deformed_mesh','off',...
                                                'deformation', U, 'deformation_mf', mfu,...
                                                'deformation_scale', 100, 'refine', 8);
colorbar;
gf_plot(mfvm, CO, 'quiver', 'on', 'quiver_density', 0.1, 'mesh', 'off', 'disp_options', 'off',...
                  'deformed_mesh', 'off', 'deformation', U, 'deformation_mf', mfu,
                  'deformation_scale', 100, 'refine', 8);
colorbar;
title('Electric potential in Volt (on the deformed configuration, scale factor x100)');
hold off;
subplot(3,1,3);
gf_plot(mft, gf_model_get(md, 'variable', 'T'), 'mesh', 'off', 'disp_options', 'off',...
                                                'deformed_mesh', 'off',...
                                                'deformation', U, 'deformation_mf', mfu,...
                                                'deformation_scale', 100, 'refine', 8);
colorbar;
title('Temperature in ^oC (on the deformed configuration, scale factor x100)');
