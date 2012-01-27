% Matlab GetFEM++ interface
%
% Copyright (C) 2009-2011 Yves Renard.
%
% This file is a part of GetFEM++
%
% GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 2.1 of the License,  or
% (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
%
% Elastodynamic problem with unilateral contact with a rigid obstacle.
% Newmark scheme.
%
% This program is used to check that matlab-getfem is working. This is also
% a good example of use of GetFEM++.
%

gf_workspace('clear all');
clear all;

plot_mesh = true;




m=gf_mesh('cartesian', [0:0.1:1]);

% Import the mesh : disc
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h4.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h2.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h1.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h0_5.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h0_3.mesh');

% Import the mesh : sphere
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_8_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_80_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_2000_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_16000_elts.mesh');


d = gf_mesh_get(m, 'dim'); % Mesh dimension


% Parameters of the model

clambda = 1;             % Lame coefficient
cmu = 1;                 % Lame coefficient
vertical_force = 0.05;   % Volumic load in the vertical direction
r = 10;                  % Augmentation parameter
dt = 0.1;                % time step
beta = 0.5;              % Newmark scheme coefficient
gamma = 0.5;             % Newmark scheme coefficient

niter = 100;             % Maximum number of iterations for Newton's algorithm.

% Signed distance representing the obstacle
if (d == 1) obstacle = 'x'; elseif (d == 2) obstacle = 'y'; else obstacle = 'z'; end;

% Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2;

border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(d, :) < -0.01));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
contact_boundary=border(:, find(normals(d, :) > 0.01));
gf_mesh_set(m, 'region', GAMMAD, contact_boundary);




% Finite element methods
u_degree = 1;
lambda_degree = 1;

mfu=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', u_degree);
mfd=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', u_degree);
mflambda=gf_mesh_fem(m, 1); % used only by versions 5 to 13
gf_mesh_fem_set(mflambda, 'classical fem', lambda_degree);
mfvm=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', u_degree-1);

% Integration method
mim=gf_mesh_im(m, 4);
mim_sing=gf_mesh_im(m);

first_elem = -1;

for i = gf_mesh_get(m, 'cvid')
   
    PT = gf_mesh_get(m, 'pts from cvid', i);
    if (min(PT) ~= 0.)
        first_elem = i;
        disp(sprintf('First element = %d', i));
        gf_mesh_im_set(mim_sing, 'integ', 4, i);
    end
end


M = gf_asm('mass matrix', mim, mfu);
M_sing = gf_asm('mass matrix', mim_sing, mfu);

% Plot the mesh
if (plot_mesh)
  figure(1);
  gf_plot_mesh(m, 'regions', [GAMMAC]);
  title('Mesh and contact boundary (in red)');
  pause(0.1);
end;

nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');

% Volumic density of force
F = zeros(nbdofd*d, 1);
F(d:d:nbdofd*d) = -vertical_force;

% Elasticity model
md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
ind_el = gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*beta));
ind_rhs = gf_model_set(md, 'add explicit rhs', 'u', zeros(nbdofu,1));

gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', u_degree, GAMMAD);

ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
gf_model_set(md, 'add fem variable', 'lambda_n', mflambda_partial);
gf_model_set(md, 'add initialized data', 'r', [r]);
OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
gf_model_set(md, 'add continuous contact with rigid obstacle brick', ...
      mim, 'u', 'lambda_n', 'obstacle', 'r', GAMMAC, 1);

nbdofl = size(ldof,1);
U0 = zeros(nbdofu, 1);
V0 = zeros(nbdofu, 1);
LAMBDA0 = zeros(nbdofl, 1);
FF = gf_asm('volumic source', mim, mfu, mfd, F');
K = gf_asm('linear elasticity', mim, mfu, mfd, ones(nbdofd,1)*clambda, ones(nbdofd,1)*cmu);


for t = 0:dt:1
  
  LL = gf_asm('boundary', GAMMAC, 'l=data(#2);V(#1)+=comp(vBase(#1).Base(#2).Normal())(:,i,k,i).l(k)', mim, mfu, mflambda_partial, LAMBDA0);

  % calcul de B
  B = M*U0 + dt*M*V0 + dt*dt*(1-beta)*(-K*U0+FF+LL);
  
  gf_model_set(md, 'set private rhs', ind_rhs, B)
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'very noisy', 'max_iter', niter);

  U1 = gf_model_get(md, 'variable', 'u');
  LAMBDA1 = gf_model_get(md, 'variable', 'u');
  lambda_n = gf_model_get(md, 'variable', 'lambda_n');

  % TODO : computation of the velocity V1 ...

  U0 = U1;
  LAMBDA0 = LAMBDA1;
  V0 = V1;
end
    
if (d > 1)
  VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
 	            'u', 'clambda', 'cmu', mfvm);
end


% set a custom colormap
% r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r);

figure(2);
if (d == 3)
  c=[0.1;0;20]; x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
  % Whole boundary
  % sl2=gf_slice({'boundary',{'none'}}, m, 5);
  % Slice, 3 planes
  % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,5);
  % Slice, 2 planes
  % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y}}},m,5);
  % Slice, 1 plane
  sl2=gf_slice({'boundary',{'planar',+1,c,x}}, m, 5);

  P=gf_slice_get(sl2,'pts'); dP=gf_compute(mfu,U,'interpolate on',sl2);
  gf_slice_set(sl2, 'pts', P+dP);
  VMsl=gf_compute(mfvm,VM,'interpolate on',sl2);
  set(gcf,'renderer','zbuffer');
  h=gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl);
  view(-80,-15); axis on; camlight; gf_colormap('chouette');
  % map=[1:-1/10:0]'*[1 1 1]; colormap(map); % for NB
    
  % gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
  %        gf_mesh_get(mfu,'outer faces'), 'deformation', U, ...
  %        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
  % view(-5,-10); camlight; colormap(map);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Sliced deformed configuration (not really a small deformation of course ...)');
elseif (d == 2)
  gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
          'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
  xlabel('x'); ylabel('y');
  title('Deformed configuration (not really a small deformation of course ...)');
end;

colorbar;
pause(0.1);
