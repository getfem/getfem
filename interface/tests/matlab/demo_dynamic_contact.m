% Copyright (C) 2009-2012 Yves Renard.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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




m=gf_mesh('cartesian', [0:0.05:1]);

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
vertical_force = 1.0;    % Volumic load in the vertical direction
r = 10;                  % Augmentation parameter
dt = 0.0005;              % time step
beta = 0.25;             % Newmark scheme coefficient
gamma = 0.5;             % Newmark scheme coefficient
theta = 0.5;             % theta-method scheme coefficient
dirichlet_val = 0.45;
singular_mass = 1;
with_theta_method = 1;

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

if (singular_mass)
  M = M_sing;
end

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
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
if (with_theta_method)
    gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*theta*theta));
else
    gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*beta));
end
ind_rhs = gf_model_set(md, 'add explicit rhs', 'u', zeros(nbdofu,1));

gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

gf_model_set(md, 'add initialized data', 'dirichletdata', [dirichlet_val]);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, GAMMAD, 'dirichletdata');

ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
gf_model_set(md, 'add fem variable', 'lambda_n', mflambda_partial);
gf_model_set(md, 'add initialized data', 'r', [r]);
OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
gf_model_set(md, 'add continuous contact with rigid obstacle brick', ...
      mim, 'u', 'lambda_n', 'obstacle', 'r', GAMMAC, 1);

nbdofl = size(ldof,1);
U0 = (gf_mesh_fem_get(mfu, 'eval', { sprintf('%g+0.5-0.5*x', dirichlet_val)}))'

gf_plot(mfu, U0');
pause;

MV0 = zeros(nbdofu, 1);
V1 = zeros(nbdofu, 1);
% LAMBDA0 = zeros(nbdofl, 1);
FF = gf_asm('volumic source', mim, mfu, mfd, F');
K = gf_asm('linear elasticity', mim, mfu, mfd, ones(nbdofd,1)*clambda, ones(nbdofd,1)*cmu);
% LL0 = gf_asm('boundary', GAMMAC, 'l=data(#2);V(#1)+=comp(vBase(#1).Base(#2).Normal())(:,i,k,i).l(k)', mim, mfu, mflambda_partial, LAMBDA0);
MA0 = FF-K*U0;

nit = 0;
tplot = 0;
for t = 0:dt:10
  
  % calcul de B
  if (with_theta_method)
      B = (M*U0 + dt*MV0)/(dt*dt*theta*theta) + (1-theta)*MA0/theta;
  else
      B = (M*U0 + dt*MV0 + dt*dt*(1/2-beta)*MA0)/(beta*dt*dt);
  end
  
  gf_model_set(md, 'set private rhs', ind_rhs, B)
  gf_model_get(md, 'solve', 'max_res', 1E-9, 'very noisy', 'max_iter', niter);

  U1 = (gf_model_get(md, 'variable', 'u'))';
  lambda_n = gf_model_get(md, 'variable', 'lambda_n');
  disp(sprintf('u1(1) = %g', U1(1)));
  disp(sprintf('lambda_n = %g', lambda_n(1)));
  
  if (with_theta_method)
     MV1 = ((M*U1 - M*U0)/dt -(1-theta)*MV0)/theta;
     MA1 = ((MV1-MV0)/dt - (1-theta)*MA0)/theta;
  else
     MA1 = (M*U1-M*U0-dt*MV0-dt*dt*(1/2-beta)*MA0)/(dt*dt*beta);
     MV1 = MV0 + dt*(gamma*MA1 + (1-gamma)*MA0);
  end   
      
  Msize = size(M,1);
  % MV1(Msize) = 0;
  if (singular_mass)
    V1(2:Msize) = M(2:Msize, 2:Msize) \ MV1(2:Msize);
    V1(1) = 0;
  else
    V1 = M \ MV1;
    disp('ok...');
  end
  % V1(Msize) = 0;
  
  E = (V1'*MV1 + U1'*K*U1)/2 - FF'*U1;
  disp(sprintf('energy = %g', E));
  
  nit = nit + 1;
  disp(sprintf('nit = %d', nit));
  if (t >= tplot)
    gf_plot(mfu, U1');
    pause;
    tplot = tplot + 0.1;
  end;
  

  U0 = U1;
  MV0 = MV1;
  MA0 = MA1;
end
    


