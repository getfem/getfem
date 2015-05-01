% Copyright (C) 2009-2012 Yves Renard.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
% Elastodynamic problem with unilateral contact with a rigid obstacle.
% Newmark and theta-method schemes.
%
% This program is used to check that matlab-getfem is working. This is also
% a good example of use of GetFEM++.
%

gf_workspace('clear all');
clear all;
gf_util('trace level', 0);

% NX = 20; m=gf_mesh('cartesian', [0:1/NX:1]); % Cas 1D

% Import the mesh : disc
m=gf_mesh('load', '../../../tests/meshes/disc_P2_h4.mesh');
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

if (d == 1)
  clambda = 1;             % Lame coefficient
  cmu = 1;                 % Lame coefficient
  friction = 0;            % Friction coefficient
  vertical_force = 0.0;    % Volumic load in the vertical direction
  r = 10;                  % Augmentation parameter
  dt = 0.001;              % Time step
  T = 4;                   % Simulation time
  dt_plot = 0.01;          % Drawing step;
  beta = 0.5;              % Newmark scheme coefficient
  gamma = 1.0;             % Newmark scheme coefficient
  theta = 0.5;             % Theta-method scheme coefficient
  dirichlet = 1;           % Dirichlet condition or not
  dirichlet_val = 0.0;
  scheme = 4;              % 1 = theta-method, 2 = Newmark, 3 = Newmark with beta = 0, 4 = midpoint modified
  u_degree = 1;
  v_degree = 1;
  lambda_degree = 1;
  Nitsche = 1;             % Use Nitsche's method or not
  gamma0_N = 0.001;        % Parameter gamma0 for Nitsche's method
  theta_N = 1;             % Parameter theta for Nitsche's method
else
  clambda = 20;            % Lame coefficient
  cmu = 20;                % Lame coefficient
  friction = 0;            % Friction coefficient
  vertical_force = 0.1;    % Volumic load in the vertical direction
  r = 10;                  % Augmentation parameter
  dt = 0.001;               % Time step
  T = 40;                  % Simulation time
  dt_plot = 0.5;           % Drawing step;
  beta = 0.25;             % Newmark scheme coefficient
  gamma = 0.5;             % Newmark scheme coefficient
  theta = 0.5;             % Theta-method scheme coefficient
  dirichlet = 0;           % Dirichlet condition or not
  dirichlet_val = 0.45;
  scheme = 3;              % 1 = theta-method, 2 = Newmark, 3 = Newmark with beta = 0, 4 = midpoint modified
  u_degree = 2;
  v_degree = 1;
  lambda_degree = 1;
  Nitsche = 1;             % Use Nitsche's method or not
  gamma0_N = 0.01;          % Parameter gamma0 for Nitsche's method
  theta_N =  1;          % Parameter theta for Nitsche's method
end
  
singular_mass = 0;         % 0 = standard method
                           % 1 = Mass elimination on boundary
                           % 2 = Mixed displacement/velocity
niter = 100;               % Maximum number of iterations for Newton's algorithm.
plot_mesh = false;
make_movie = 0;
residual = 1E-8;

if (scheme >= 3 && (Nitsche ~= 1 || singular_mass ~= 0))
    error('Incompatibility');
end

if (friction ~= 0 && d == 1)
    error('Not taken into account');
end

% Signed distance representing the obstacle
if (d == 1) obstacle = 'x'; elseif (d == 2) obstacle = 'y'; else obstacle = 'z'; end;

% Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2;

border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(d, :) < -0.01));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
dirichlet_boundary=border(:, find(normals(d, :) > 0.01));
gf_mesh_set(m, 'region', GAMMAD, dirichlet_boundary);

% Finite element methods

mfu=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', u_degree);
mfv=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfv, 'classical fem', v_degree);
mfd=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', u_degree);
if (friction == 0)
  mflambda=gf_mesh_fem(m, 1);
else
  mflambda=gf_mesh_fem(m, d); 
end
gf_mesh_fem_set(mflambda, 'classical fem', lambda_degree);
mfvm=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', u_degree-1);

% Integration method
mim=gf_mesh_im(m, 4);
mim_sing=gf_mesh_im(m);
if (d == 1)
  mim_friction = mim;
elseif (d == 2)
  mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'));
elseif (d == 3)
   mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),4)')); 
end;


first_elem = -1;

M = gf_asm('mass matrix', mim, mfu);

if (singular_mass == 1) % Rought singular mass matrix (no redistribution)
  for i = gf_mesh_get(m, 'cvid')
      if (size(find(contact_boundary(1,:) == i), 2) == 0) 
          gf_mesh_im_set(mim_sing, 'integ', 4, i);
      end
  end
  M_sing = gf_asm('mass matrix', mim_sing, mfu);
  M = M_sing;
end

if (singular_mass == 2)
  B = gf_asm('mass matrix', mim, mfv, mfu);
  C = gf_asm('mass matrix', mim, mfv, mfv);
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
F = zeros(d, nbdofd);
F(d,:) = -vertical_force;

% Elasticity model
md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
if (singular_mass == 2)
  gf_model_set(md, 'add fem variable', 'v', mfv);
  switch(scheme)
    case 1
      gf_model_set(md, 'add explicit matrix', 'u', 'v', (B')/(dt*dt*theta*theta));
    case 2
      gf_model_set(md, 'add explicit matrix', 'u', 'v', (B')/(dt*dt*beta));
  end
  gf_model_set(md, 'add explicit matrix', 'v', 'v', C);
  gf_model_set(md, 'add explicit matrix', 'v', 'u', -B);
else
  switch(scheme)
    case 1
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*theta*theta));
    case 2
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*beta));
    case 4
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*0.25));
  end
end
ind_rhs = gf_model_set(md, 'add explicit rhs', 'u', zeros(nbdofu,1));

gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);

gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

if (dirichlet)
  dirichletdata = zeros(1,d); dirichletdata(d) = dirichlet_val;
  gf_model_set(md, 'add initialized data', 'dirichletdata', dirichletdata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, GAMMAD, 'dirichletdata');
end

OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);

if (Nitsche)
  gf_model_set(md, 'add initialized data', 'gamma0', [gamma0_N]);
  if (scheme == 4)
      if (friction ~= 0)
         error('To be adapted for friction');
      end
      gf_model_set(md, 'add initialized data', 'friction_coeff', [0]);
      gf_model_set(md, 'add initialized data', 'alpha_f', [0]);
      gf_model_set(md, 'add fem data', 'wt', mfu);
      
      % N??cessaire ?
      % gf_model_set(md, 'add Nitsche contact with rigid obstacle brick', mim_friction, 'u', ...
      %  'obstacle', 'gamma0', GAMMAC, theta_N);
      
      gf_model_set(md, 'add Nitsche midpoint contact with rigid obstacle brick', mim_friction, 'u', ...
       'obstacle', 'gamma0', GAMMAC, theta_N, 'friction_coeff', 'alpha_f', 'wt', 2);
      gf_model_set(md, 'add Nitsche midpoint contact with rigid obstacle brick', mim_friction, 'u', ...
       'obstacle', 'gamma0', GAMMAC, theta_N, 'friction_coeff', 'alpha_f', 'wt', 1);
  else
    if (friction == 0)
      gf_model_set(md, 'add Nitsche contact with rigid obstacle brick', mim_friction, 'u', ...
                   'obstacle', 'gamma0', GAMMAC, theta_N);
    else
      gf_model_set(md, 'add initialized data', 'friction_coeff', [friction]);
      gf_model_set(md, 'add initialized data', 'alpha_f', [1./dt]);
      gf_model_set(md, 'add fem data', 'wt', mfu);
      gf_model_set(md, 'add Nitsche contact with rigid obstacle brick', mim_friction, 'u', ...
           'obstacle', 'gamma0', GAMMAC, theta_N, 'friction_coeff', 'alpha_f', 'wt');
    end
  end
else
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda', mflambda_partial);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (friction == 0)
    gf_model_set(md, 'add integral contact with rigid obstacle brick', ...
                 mim_friction, 'u', 'lambda', 'obstacle', 'r', GAMMAC, 1);
  else
    gf_model_set(md, 'add initialized data', 'friction_coeff', [friction]);
    gf_model_set(md, 'add initialized data', 'alpha_f', [1./dt]);
    gf_model_set(md, 'add fem data', 'wt', mfu);
    gf_model_set(md, 'add integral contact with rigid obstacle brick', mim_friction, 'u', ...
                 'lambda', 'obstacle', 'r', 'friction_coeff', GAMMAC, 1, 'alpha_f', 'wt');
  end
end

if (d == 1)
    U0 = (gf_mesh_fem_get(mfu, 'eval', { sprintf('%g+0.5-0.5*x', dirichlet_val)}))';
else
    U0 = zeros(nbdofu, 1);
    U0(d:d:nbdofu) = 5;
end;
if (singular_mass == 2)
    VV1 = B*U0; VV2 = C\VV1; MU0 = (B')*VV2;
else
    MU0 = M*U0;
end

MV0 = zeros(nbdofu, 1);
V1 = zeros(nbdofu, 1);
FF = gf_asm('volumic source', mim, mfu, mfd, F);
K = gf_asm('linear elasticity', mim, mfu, mfd, ones(nbdofd,1)*clambda, ones(nbdofd,1)*cmu);
MA0 = FF-K*U0;
if (singular_mass == 1)
  if (d == 1)
    MA0(1) = 0;
  else
    error('Take it into account !');
  end
elseif (singular_mass == 2) % to be verified
  VV1 = (B') \ MA0; VV2 = C*VV1; A0 = B\VV2; VV1 = B*A0; VV2 = C\VV1; MA0 = (B')*VV2;
end
nit = 0; tplot = 0;
if (make_movie)
  nim = 0;
  figure('position', [100 100 800 600]); % Necessary for the movie to be read by vlc
  mov = avifile('toto.avi');
end

for t = 0:dt:T
  disp(sprintf('t=%g', t));
  % calcul de LL
  
  switch(scheme)
    case 1
      LL = (MU0 + dt*MV0)/(dt*dt*theta*theta) + (1-theta)*MA0/theta;
    case 2
      LL = (MU0 + dt*MV0 + dt*dt*(1/2-beta)*MA0)/(beta*dt*dt);
    case 3
      LL = 0*MU0;
    case 4
      LL = MU0/(dt*dt*0.25) + MV0/(dt*0.5);
  end
  
  if (friction ~= 0 || scheme == 4)
    gf_model_set(md, 'variable', 'wt', U0);
    % disp(gf_model_get(md, 'variable', 'wt'));
  end
  
  if (scheme == 3)
    A0 = M \ MA0;
    V0 = M \ MV0;
    U1 = U0 + dt*dt*A0/2 + dt*V0;
    gf_model_set(md, 'variable', 'u', U1);
    gf_model_get(md, 'assembly', 'build_rhs');
  else
    gf_model_set(md, 'set private rhs', ind_rhs, LL);
    gf_model_get(md, 'solve', 'max_res', residual, 'max_iter', niter); % , 'noisy');
    U1 = (gf_model_get(md, 'variable', 'u'))';
  end
  
  if (singular_mass == 2)
    MU1 = (B')*(gf_model_get(md, 'variable', 'v'))';
  else
    MU1 = M*U1;
  end
  
  if (d == 1)
    disp(sprintf('u1(1) = %g', U1(1)));
    Msize = size(M,1);
    if (Nitsche == 0)
      lambda = gf_model_get(md, 'variable', 'lambda');
      disp(sprintf('lambda_n = %g', lambda(1)));
    end
  
    disp(sprintf('U0(N) = %g', U0(Msize)));
    disp(sprintf('U1(N) = %g', U1(Msize)));
    disp(sprintf('MV0(N) = %g', MV0(Msize)));
  end
  
  switch(scheme)
    case 1
      MV1 = ((MU1 - MU0)/dt -(1-theta)*MV0)/theta;
      MA1 = ((MV1-MV0)/dt - (1-theta)*MA0)/theta;
    case 2
      MA1 = (MU1-MU0-dt*MV0-dt*dt*(1/2-beta)*MA0)/(dt*dt*beta);
      MV1 = MV0 + dt*(gamma*MA1 + (1-gamma)*MA0);
    case 3
      MA1 = (gf_model_get(md, 'rhs'))';
      MV1 = MV0 + dt*(gamma*MA1+(1-gamma)*MA0);
    case 4
      U1_2 = U1;
      U1 = 2*U1_2 - U0;
      V1_2 = (U1 - U0)/dt;
      MV1 = 2*M*V1_2 - MV0;
      MA1 = 0*MV1;
      MU1 = M*U1;
  end
      
  if (singular_mass == 1)
    V1 = cgs(M, MV1); % Pseudo inverse ...
  elseif (singular_mass == 2)
    VV1 = (B') \ MV1; VV2 = C*VV1; V1 = B\VV2; 
  else
    V1 = M \ MV1;
  end

  
  E = (V1'*MV1 + U1'*K*U1)/2 - FF'*U1;
  disp(sprintf('energy = %g', E));
  
  nit = nit + 1;
  if (t >= tplot)
      if (d >= 2)
        VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		    'u', 'clambda', 'cmu', mfvm);
      end
      if (d == 1)
        X = [0:1/NX:1]';
        plot(zeros(1, Msize)-0.05, X+U1, '-b');
        hold on;
        plot(zeros(1, Msize)+0.05, U1+X, '-b');
        for i = 1:NX+1
           plot([-0.05 0.05], (U1(i)+X(i))*[1 1], 'b'); 
        end
        hold off;
        axis([-0.4 0.4 0.0 1.5]);
      elseif (d == 2)
        gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U1', ...
            'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8, 'disp_options', 'off');
        xlabel('x'); ylabel('y');
        % title('Deformed configuration (not really a small deformation of course ...)');
        % gf_colormap('chouette');
        gg = [ .7 .9 .4; .5 .9 .3;   .3 .8 .2;    .1 .7 .4;     .2 0.7 1.0000; .3 0.3 1.0000;
               1.0 .8 .1;  1.0 .6 .1;   1.0 .45 .1;   1.0 0.3 .1];
        r = reshape(repmat(gg',6,1),3,60)';
        colormap(r);
        colorbar;
        caxis([0 32]);
        axis([-25 25 -1 50]);
      else
        c=[0.1;0;20]; x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
        % Whole boundary
        % sl2=gf_slice({'boundary',{'none'}}, m, 5);
        % Slice, 3 planes
        % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,5);
        % Slice, 2 planes
        sl2=gf_slice({'boundary',{'union',{'planar',+1,c,y},{'planar',+1,c,x}}},m,5);
        % Slice, 1 plane
        % sl2=gf_slice({'boundary',{'planar',+1,c,x}}, m, 5);

        P=gf_slice_get(sl2,'pts'); dP=gf_compute(mfu,U1','interpolate on',sl2);
        gf_slice_set(sl2, 'pts', P+dP);
        VMsl=gf_compute(mfvm,VM,'interpolate on',sl2);
        set(gcf,'renderer','zbuffer');
        h=gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl);
        view(-80,-15); axis on; camlight; gf_colormap('chouette');
        % map=[1:-1/10:0]'*[1 1 1]; colormap(map); % for NB
    
        % gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
        %        gf_mesh_get(mfu,'outer faces'), 'deformation', U, ...
        %        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
        % view(-5,-10); camlight; colormap(map);
        xlabel('x'); ylabel('y'); zlabel('z');
        axis([-25 25 -25 25 -1 50]);
        caxis([0 20]);
        colorbar;
        % title('Sliced deformed configuration');
      end
    pause(0.1);
    if (make_movie)
      nim = nim + 1;
      F = getframe(gcf);
      mov = addframe(mov,F);
      Mov(:,nim) = getframe;
    end
    tplot = tplot + dt_plot;
  end;
  

  U0 = U1;
  MU0 = MU1;
  MV0 = MV1;
  MA0 = MA1;
end
   

if (make_movie)
  mov = close(mov);
  mov = aviread('toto.avi');
  movie(mov);
end


