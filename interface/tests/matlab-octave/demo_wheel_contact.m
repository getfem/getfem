% Copyright (C) 2015-2020 Yves Renard.
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
% Contact of a deformable 'wheel' onto a plane deformable obstacle.
%

view_mesh = true;
Dirichlet_version = false; % A Dirichlet condition instead of the applied load

%
% Physical parameters
%
E = 21E6;                         % Young Modulus (N/cm^2)
nu = 0.3;                         % Poisson ratio
clambda = E*nu/((1+nu)*(1-2*nu)); % First Lame coefficient (N/cm^2)
cmu = E/(2*(1+nu));               % Second Lame coefficient (N/cm^2)
clambdastar = 2*clambda*cmu/(clambda+2*cmu); % Lame coefficient for Plane stress (N/cm^2)
applied_force = 1E7;              % Force at the hole boundary (N)


%
% Numerical parameters
%
h = 4;                   % Approximate mesh size
elements_degree = 2;     % Degree of the finite element methods
gamma0 = 1./E;           % Augmentation parameter for the augmented Lagrangian 

%
% Mesh generation. Meshes can also been imported from several formats.
%
mo1 = gfMesherObject('ball', [0., 15.], 15.);
mo2 = gfMesherObject('ball', [0., 15.], 8.);
mo3 = gfMesherObject('set minus', mo1, mo2);

disp('Meshes generation');
gf_util('trace level', 1)   % No trace for mesh generation nor for assembly
mesh1 = gfMesh('generate', mo3, h, 2);
mesh2 = gfMesh('import','structured', sprintf('GT="GT_PK(2,1)";SIZES=[30,10];NOISED=0;NSUBDIV=[%d,%d];', floor(30/h)+1, floor(10/h)+1));
mesh2.translate([-15.;-10.]);


if (view_mesh)
  figure(2);
  gf_plot_mesh(mesh1);
  hold on;
  gf_plot_mesh(mesh2);
  hold off;
  pause(1);
end

%
% Boundary selection
%
fb1 = mesh1.outer_faces_in_box([-8.1, 6.9], [8.1, 23.1]);  % Boundary of the hole
fb2 = mesh1.outer_faces_with_direction([0., -1.], pi/4.5);  % Contact boundary of the wheel
fb3 = mesh2.outer_faces_with_direction([0., -1.], 0.01);   % Bottom boundary of the foundation

HOLE_BOUND=1; CONTACT_BOUND=2; BOTTOM_BOUND=3;

mesh1.set_region(HOLE_BOUND, fb1);
mesh1.set_region(CONTACT_BOUND, fb2);
mesh1.region_subtract(CONTACT_BOUND, HOLE_BOUND);
mesh2.set_region(BOTTOM_BOUND, fb3);


%
% Definition of finite elements methods and integration method
%

mfu1 = gfMeshFem(mesh1, 2);
mfu1.set_classical_fem(elements_degree);
mflambda = gfMeshFem(mesh1, 2);
mflambda.set_classical_fem(elements_degree-1);
mflambda_C = gfMeshFem(mesh1, 1);
mflambda_C.set_classical_fem(elements_degree-1);
mfu2 = gfMeshFem(mesh2, 2);
mfu2.set_classical_fem(elements_degree);
mfvm1 = gfMeshFem(mesh1, 1);
mfvm1.set_classical_discontinuous_fem(elements_degree);
mfvm2 = gfMeshFem(mesh2, 1);
mfvm2.set_classical_discontinuous_fem(elements_degree);
mim1  = gfMeshIm(mesh1, 6);
mim1c = gfMeshIm(mesh1, gfInteg('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'));
mim2  = gfMeshIm(mesh2, 4);


%
% Model definition
%

md=gfModel('real');
md.add_fem_variable('u1', mfu1);       % Displacement of the structure 1
md.add_fem_variable('u2', mfu2);       % Displacement of the structure 2

md.add_initialized_data('cmu', [cmu]);
md.add_initialized_data('clambdastar', [clambdastar]);
md.add_isotropic_linearized_elasticity_brick(mim1, 'u1', 'clambdastar', 'cmu');
md.add_isotropic_linearized_elasticity_brick(mim2, 'u2', 'clambdastar', 'cmu');
md.add_Dirichlet_condition_with_multipliers(mim2, 'u2', elements_degree-1, BOTTOM_BOUND);



% Contact condition (augmented Lagrangian)
md.add_initialized_data('gamma0', [gamma0]);
md.add_interpolate_transformation_from_expression('Proj1', mesh1, mesh2, '[X(1);0]');
md.add_filtered_fem_variable('lambda1', mflambda_C, CONTACT_BOUND);
md.add_nonlinear_generic_assembly_brick(mim1c, 'lambda1*(Test_u1.[0;1])-lambda1*(Interpolate(Test_u2,Proj1).[0;1])', CONTACT_BOUND);
md.add_nonlinear_generic_assembly_brick(mim1c, '-(gamma0*element_size)*(lambda1 + neg_part(lambda1+(1/(gamma0*element_size))*((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).[0;1])))*Test_lambda1', CONTACT_BOUND);

% Prescribed force in the hole (or Dirichlet condition)
if (Dirichlet_version)
  md.add_initialized_data('DData', [0., -1.0]);
  md.add_Dirichlet_condition_with_multipliers(mim1, 'u1', elements_degree-1, HOLE_BOUND, 'DData');
else
  md.add_filtered_fem_variable('lambda_D', mflambda, HOLE_BOUND);
  md.add_initialized_data('F', [applied_force/(8.*2.*pi)]);
  md.add_variable('alpha_D', 1);
  md.add_linear_generic_assembly_brick(mim1, '-lambda_D.Test_u1 + (alpha_D*[0;1] - u1).Test_lambda_D + (lambda_D.[0;1] + F)*Test_alpha_D + 1E-6*alpha_D*Test_alpha_D', HOLE_BOUND);
  % The small penalization 1E-6*alpha_D*Test_alpha_D seems necessary to
  % have a convergence in all cases. Why ?
end

%
% Model solve
%

disp(sprintf('Solve problem with %d dofs', md.nbdof()));
md.solve('max_res', 1E-9, 'max_iter', 40, 'noisy'); % , 'lsearch', 'simplest',  'alpha min', 0.8)
if (~Dirichlet_version)
  disp(sprintf('alpha_D = %g', md.variable('alpha_D')));
end
% md.test_tangent_matrix(1e-4, 10, 0.001);
md.variable('lambda1')

%
% Solution draw
%  
U1 = md.variable('u1');
U2 = md.variable('u2');
VM1 = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u1', 'clambdastar', 'cmu', mfvm1);
VM2 = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u2', 'clambdastar', 'cmu', mfvm2);

figure(1);
gf_plot(mfvm1, VM1, 'deformed_mesh', 'on', 'deformation', U1, 'deformation_mf', mfu1, 'deformation_scale', 1, 'refine', 8);
hold on;
gf_plot(mfvm2, VM2, 'deformed_mesh', 'on', 'deformation', U2, 'deformation_mf', mfu2, 'deformation_scale', 1, 'refine', 8);
xlabel('x'); ylabel('y');
title('Deformed configuration');
hold off;

