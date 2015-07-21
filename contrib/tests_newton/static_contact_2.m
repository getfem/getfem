% Copyright (C) 2012-2015 Yves Renard, Julien Pommier.
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


% The newton should converge every time with about 10 to 15 iterations
% (should not go upper than 15 iterations ...)


clear all;
is_automatic = false;

if (is_automatic) 
    disp('automatic version');
    draw = false;
    plot_mesh = false;
    vertical_force = 20.0; % Volumic load in the vertical direction
    niter = 100;   % Maximum number of iterations for Newton's algorithm.
    friction_coeff = 1.0;  % coefficient of friction
    gf_util('trace level', 1);
else
  disp('non-automatic version');
  clear all;
  % main parameters 
  expe = 3;              % Experiment number
  r = 1e-6;                  % Augmentation parameter
  dirichlet_translation = -15;
  vertical_force = 20.0; % Volumic load in the vertical direction
  niter = 100;           % Maximum number of iterations for Newton's algorithm.
  friction_coeff = 1.0;  % coefficient of friction

  draw = false;
  plot_mesh = false;
  version = 1; % 1 : frictionless contact and the basic contact brick
              % 2 : contact with 'static' Coulomb friction and basic contact brick
              % 3 : frictionless contact and the contact with a
              %     rigid obstacle brick, symmetric version
              % 4 : contact with 'static' Coulomb friction and the contact with a
              %     rigid obstacle brick, symmetric version
              % 5 : frictionless contact and the integral brick
              %     Newton and Alart-Curnier augmented lagrangian,
              %     unsymmetric version
              % 6 : frictionless contact and the integral brick
              %     Newton and Alart-Curnier augmented lagrangian, symmetric
              %     version.
              % 7 : frictionless contact and the integral brick
              %     Newton and Alart-Curnier augmented lagrangian,
              %     unsymmetric version with an additional augmentation.
              % 8 : frictionless contact and the integral brick
              %     New unsymmetric method.
              % 9 : frictionless contact and the integral brick : Uzawa
              %     on the Lagrangian augmented by the penalization term.
              % 10 : contact with 'static' Coulomb friction and the integral
              %     brick. Newton and Alart-Curnier augmented lagrangian,
              %     unsymmetric version.
              % 11 : contact with 'static' Coulomb friction and the integral
              %     brick. Newton and Alart-Curnier augmented lagrangian,
              %     nearly symmetric version.
              % 12 : contact with 'static' Coulomb friction and the integral
              %     brick. Newton and Alart-Curnier augmented lagrangian,
              %     unsymmetric version with an additional augmentation.
              % 13 : contact with 'static' Coulomb friction and the integral
              %     brick. New unsymmetric method.
              % 14 : "unsymmetric" De Saxce version
              % 15 : New unsymmetric method with De Saxce projection
              % 16 : contact with 'static' Coulomb friction and the integral
              %     brick : Uzawa on the Lagrangian augmented by the penalization term.
              % 17 : contact with 'static' Coulomb friction and the integral
              %     brick : Uzawa on De Saxce augmented Lagrangian.
              % 18 : penalized contact with 'static' Coulomb friction
              %     (r is the penalization coefficient).
end

% Import the mesh : 2D punch
switch (expe) 
    case 1
        m=gf_mesh('load', 'punch2D_h4.mesh');
        with_dirichlet = 1;
    case 2
        m=gf_mesh('load', 'punch2D_h2.mesh');
        with_dirichlet = 1;
    case 3
        m=gf_mesh('load', 'punch2D_h1.mesh');
        with_dirichlet = 1;
    case 4
        m=gf_mesh('load', 'punch2D_h0_5.mesh');
        with_dirichlet = 1;
    case 5
        m=gf_mesh('load', 'punch2D_h0_25.mesh');
        with_dirichlet = 1;
    case 6
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/disc_P2_h8.mesh');
        with_dirichlet = 0;
    case 7
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/disc_P2_h4.mesh');
        with_dirichlet = 0;
    case 8
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/disc_P2_h2.mesh');
        with_dirichlet = 0;
    case 9
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/disc_P2_h1.mesh');
        with_dirichlet = 0;
    case 10
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/disc_P2_h0_5.mesh');
        with_dirichlet = 0;
    case 11 
        m=gf_mesh('load', 'punch3D_h5_12.mesh'); with_dirichlet = 1;
    case 12
        m=gf_mesh('load', 'punch3D_h3_1.mesh'); with_dirichlet = 1;
    case 13
        m=gf_mesh('load', 'punch3D_h1_8.mesh'); with_dirichlet = 1;
    case 14
        error('unattributed experiment');
    case 15
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/sphere_with_quadratic_tetra_8_elts.mesh');
        with_dirichlet = 0; % h = 20
    case 16
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/sphere_with_quadratic_tetra_80_elts.mesh');
        with_dirichlet = 0; % h = 8
    case 17
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
        with_dirichlet = 0; % h = 6
    case 18
        m=gf_mesh('load', '../../../../source++/getfem/tests/meshes/sphere_with_quadratic_tetra_2000_elts.mesh');
        with_dirichlet = 0; % h = 3.5
    case 19
        m=gf_mesh('load', 'sphere_with_quadratic_tetra_6000_elts.mesh');
        with_dirichlet = 0; % h = 2.3
end

d = gf_mesh_get(m, 'dim'); % Mesh dimension
h = mean(gf_mesh_get(m, 'convex radius'));
disp(sprintf('h = %g', h));


% condition_type = 3; % 0 = No kill rigid motions (for frictional problems)
                    % 1 = Explicitely kill horizontal rigid displacements
                    % 2 = Kill rigid displacements using a global penalization
                    % 3 = Add a Dirichlet condition on the top of the structure
if (with_dirichlet)
   disp('With a clamped boundary');
   condition_type = 3;
elseif (version == 2 || version == 4 || version >= 10)
   disp('No treatment for rigid displacements');
   condition_type = 0;
else
  condition_type = 1;
  disp('Kill horizontal rigid displacements');
end

% Parameters of the model
clambda = 1000;           % Lame coefficient
cmu = 1000;               % Lame coefficient
real_r = r * clambda;
% real_r = r;
% condition_type = 3; % 0 = No kill rigid motions (for friction problems)
                    % 1 = Explicitely kill horizontal rigid displacements
                    % 2 = Kill rigid displacements using a global penalization
                    % 3 = Add a Dirichlet condition on the top of the structure
penalty_parameter = 1E-6;    % Penalization coefficient for the global penalization
                             % and residual for Uzawa methods.
residual = 6e-11 * clambda;
diverged_residual = 1e14 * clambda; % Gives up when the residual is too large.
uzawa_residual = 1e-4;
if (d == 2)
    cpoints = [0, 0];   % constraigned points for 2d
    cunitv  = [1, 0];   % corresponding constraigned directions for 2d
else
    cpoints = [0, 0, 0,   0, 0, 0,   5, 0, 5];  % constraigned points for 3d
    cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0];  % corresponding constraigned directions for 3d
end;




 % Signed distance representing the obstacle
if (d == 2) obstacle = 'y'; else obstacle = 'z'; end;

% Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2;

border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(d, :) < -0.01));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
% contact_boundary=border(:, find(normals(d, :) > 0.9));
% gf_mesh_set(m, 'region', GAMMAD, contact_boundary);


P=gf_mesh_get(m,'pts'); % get list of mesh points coordinates
pidtop=find(P(d,:) > 39.999); % find those on top of the object
ftop=gf_mesh_get(m,'faces from pid',pidtop); 
gf_mesh_set(m, 'region', GAMMAD, ftop);




% Finite element methods
u_degree = 2;
lambda_degree = 2;

mfu=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', u_degree);
mfd=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', u_degree);
mflambda=gf_mesh_fem(m, 1); % used only by versions 5 to 13
gf_mesh_fem_set(mflambda, 'classical fem', lambda_degree);
mfvm=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', u_degree-1);

nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');
disp(sprintf('Nb dof on u : %d', nbdofu)); 

% Integration method
mim=gf_mesh_im(m, 4);
if (d == 2)
  mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),2)'));
else
   mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),2)')); 
end;

% Plot the mesh
if (plot_mesh)
  figure(1);
  if (d <= 3)
    gf_plot_mesh(m, 'regions', [GAMMAC]);
    title('Mesh and contact boundary (in red)');
    axis([-21 21 0 41]);
  elseif (d == 3)
    D = zeros(1, nbdofd);
    gf_plot(mfd, D, 'mesh', 'on',  'cvlst', gf_mesh_get(mfd, 'outer faces'), 'refine', 8);
  end
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times');
  set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 18);
  set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  % pause; print(gcf,'-dpng','-r300', 'mesh.png'); return;
  pause(1);
end;

% Volumic density of force
nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');
F = zeros(nbdofd*d, 1);
F(d:d:nbdofd*d) = -vertical_force;

% Elasticity model
md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'variable', 'u', 0.01*(rand(1, gf_mesh_fem_get(mfu, 'nbdof'))-0.5));
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

if (condition_type == 3)
  Ddata = zeros(1, d); Ddata(d) = dirichlet_translation;
  gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', u_degree, GAMMAD, 'Ddata');
elseif (condition_type == 1)
  gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
  gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
  gf_model_set(md, 'add pointwise constraints with multipliers', 'u', 'cpoints', 'cunitv');
elseif (condition_type == 2)  
  % Small penalty term to avoid rigid motion (should be replaced by an
  % explicit treatment of the rigid motion with a constraint matrix)
  gf_model_set(md, 'add initialized data', 'penalty_param', ...
              [penalty_parameter]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
end;

% The contact condition

cdof = gf_mesh_fem_get(mfu, 'dof on region', GAMMAC);
nbc = size(cdof, 2) / d;

if (nbc <= 0)
    disp('No contact zone');
    return;
end;

solved = false; nb_uzawa_iter = 0; converged = false;
if (version >= 1 && version <= 4) % defining the matrices BN and BT by hand
  contact_dof = cdof(d:d:nbc*d);
  contact_nodes = gf_mesh_fem_get(mfu, 'basic dof nodes', contact_dof);
  BN = sparse(nbc, nbdofu);
  ngap = zeros(nbc, 1);
  for i = 1:nbc
    BN(i, contact_dof(i)) = -1.0;
    ngap(i) = contact_nodes(d, i);
  end;
  if (version == 2 || version == 4)
    BT = sparse(nbc*(d-1), nbdofu);
    for i = 1:nbc
      for j = 1:d-1
        BT(j+(i-1)*(d-1), contact_dof(i)-d+j) = 1.0;
      end;
    end;
  end;

  gf_model_set(md, 'add variable', 'lambda_n', nbc);
  gf_model_set(md, 'variable', 'lambda_n', 0.01*(rand(1, nbc)-0.5));
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  if (version == 2 || version == 4)
    gf_model_set(md, 'add variable', 'lambda_t', nbc*(d-1));
    gf_model_set(md, 'variable', 'lambda_t', 0.01*(rand(1, nbc*(d-1))-0.5));
    gf_model_set(md, 'add initialized data', 'friction_coeff', ...
                 [friction_coeff]);
  end;
  gf_model_set(md, 'add initialized data', 'ngap', ngap);
  gf_model_set(md, 'add initialized data', 'alpha', ones(nbc, 1));
  if (version == 1 || version == 3)
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', 'r', ...
        BN, 'ngap', 'alpha', 1+(version - 1)/2);
  else
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', ...
		 'lambda_t', 'r', BN, BT, 'friction_coeff', 'ngap', 'alpha', 1+(version - 2)/2);
  end;
% elseif (version == 3 || version == 4) % BN and BT defined by contact brick
% 
%   gf_model_set(md, 'add variable', 'lambda_n', nbc);
%   gf_model_set(md, 'add initialized data', 'r', [r]);
%   if (version == 3)
%     gf_model_set(md, 'add nodal contact with rigid obstacle brick', mim, 'u', ...
% 	         'lambda_n', 'r', GAMMAC, obstacle, 0);
%   else
%     gf_model_set(md, 'add variable', 'lambda_t', nbc * (d-1));
%     gf_model_set(md, 'add initialized data', 'friction_coeff', ...
% 		 [friction_coeff]);
%     gf_model_set(md, 'add nodal contact with rigid obstacle brick', mim, 'u', ...
% 	         'lambda_n', 'lambda_t', 'r', 'friction_coeff', GAMMAC, ...
% 		 obstacle, 0);
%   end;

elseif (version >= 5 && version <= 8) % The integral version, Newton
 
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda_n', mflambda_partial);
  gf_model_set(md, 'variable', 'lambda_n', 0.01*(rand(1, gf_mesh_fem_get(mflambda_partial, 'nbdof'))-0.5));
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add integral contact with rigid obstacle brick', ...
      mim_friction, 'u', 'lambda_n', 'obstacle', 'r', GAMMAC, version-4);
          
elseif (version == 9) % The integral version, Uzawa on the augmented Lagrangian
    
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  nbc = gf_mesh_fem_get(mflambda_partial, 'nbdof');
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  M = gf_asm('mass matrix', mim, mflambda_partial, mflambda_partial, GAMMAC);
  lambda_n = zeros(1, nbc);
  % lambda_n = (rand(1, nbc)-0.5) * 0.01;
  gf_model_set(md, 'add initialized fem data', 'lambda_n', mflambda_partial, lambda_n);
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', GAMMAC, 2, 'lambda_n');
  nb_newton_iter = 0; nb_uzawa_iter = 0;
        
  for ii=1:100
      disp(sprintf('Uzawa iteration %d', ii));
      nb_uzawa_iter = nb_uzawa_iter + 1;
[nbit, converged] = gf_model_get(md, 'solve', 'max_res', residual, 'diverged_res', diverged_residual, 'max_iter', niter, 'noisy'); % , 'very noisy');
      nb_newton_iter = nb_newton_iter  + nbit;
      if (nb_newton_iter >= niter || ~converged) 
          nb_newton_iter = niter;
          break;
      end
      U = gf_model_get(md, 'variable', 'u');
      lambda_n_old = lambda_n;
      lambda_n = (M\ gf_asm('integral contact Uzawa projection', GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda_n, mfd, OBS, real_r))';
      gf_model_set(md, 'variable', 'lambda_n', lambda_n);
      difff = max(abs(lambda_n-lambda_n_old)) / max(abs(lambda_n));
      disp(sprintf('diff: %g   threshold: %g', difff, uzawa_residual));
      % pause;
      if (difff < uzawa_residual) break; end;
  end;
  
  solved = true;
  
elseif (version >= 10 && version <= 15) % The integral version with friction, Newton
 
  gf_mesh_fem_set(mflambda, 'qdim', d);
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda', mflambda_partial);
  gf_model_set(md, 'variable', 'lambda', 0.01*(rand(1, gf_mesh_fem_get(mflambda_partial, 'nbdof'))-0.5));
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add integral contact with rigid obstacle brick', mim_friction, 'u', ...
	         'lambda', 'obstacle', 'r', 'friction_coeff', GAMMAC, version-9);

elseif (version == 16 || version == 17) % The integral version, Uzawa on the augmented Lagrangian with friction
  
  gf_mesh_fem_set(mflambda, 'qdim', d);
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  nbc = gf_mesh_fem_get(mflambda_partial, 'nbdof');
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  M = gf_asm('mass matrix', mim, mflambda_partial, mflambda_partial, GAMMAC);
  % lambda = (rand(1, nbc)-0.5) * 0.01;
  lambda = zeros(1, nbc);
  gf_model_set(md, 'add initialized fem data', 'lambda', mflambda_partial, lambda);
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', 'friction_coeff', GAMMAC, version - 14, 'lambda');
  nb_newton_iter = 0; nb_uzawa_iter = 0;
  for ii=1:100
      disp(sprintf('Uzawa iteration %d', ii));
      nb_uzawa_iter = nb_uzawa_iter + 1;
      [nbit, converged] = gf_model_get(md, 'solve', 'max_res', residual, 'diverged_res', diverged_residual, 'max_iter', niter, 'noisy'); % , 'very noisy');
      nb_newton_iter = nb_newton_iter  + nbit;
      if (nb_newton_iter >= niter || ~converged) 
          nb_newton_iter = niter;
          break;
      end
      U = gf_model_get(md, 'variable', 'u');
      lambda_old = lambda;
      lambda = (M\ gf_asm('integral contact Uzawa projection', GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda, mfd, OBS, real_r, friction_coeff, version-15))';
      gf_model_set(md, 'variable', 'lambda', lambda);
      difff = max(abs(lambda-lambda_old))/max(abs(lambda));
      disp(sprintf('diff: %g   threshold: %g', difff, uzawa_residual));
      
      % pause;
      if (difff < uzawa_residual) break; end;
  end;
  
  solved = true;

elseif (version == 18)
 
  gf_model_set(md, 'add initialized data', 'r', [real_r]);
  gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
  gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);
  gf_model_set(md, 'add penalized contact with rigid obstacle brick', mim_friction, 'u', ...
	         'obstacle', 'r', 'friction_coeff', GAMMAC);
    
else
  error('Inexistent version');
end

% Solve the problem
if (~solved)
  [nb_newton_iter, converged] = gf_model_get(md, 'solve', 'max_res', residual, 'diverged_res', diverged_residual, 'noisy', 'max_iter', niter); % , 'lsearch', 'simplest'); % , 'with pseudo potential');
end;

if (~converged)
    nb_newton_iter = niter;
end

U = gf_model_get(md, 'variable', 'u');
% lambda_n = gf_model_get(md, 'variable', 'lambda_n');
VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u', 'clambda', 'cmu', mfvm);
    

% set a custom colormap
% r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r);


if (draw)

  figure(2);
  if (d == 3)
    c=[0.1;0;20]; x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
    % Whole boundary
    % sl2=gf_slice({'boundary',{'none'}}, m, 5);
    % Slice, 3 planes
    % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,5);
    % Slice, 2 planes
    sl2=gf_slice({'boundary',{'union',{'planar',+1,c,y},{'planar',+1,c,x}}},m,5);
    % Slice, 1 plane
    % sl2=gf_slice({'boundary',{'planar',+1,c,x}}, m, 5);

    P=gf_slice_get(sl2,'pts'); dP=gf_compute(mfu,U,'interpolate on',sl2);
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
    % title('Sliced deformed configuration (not really a small deformation of course ...)');
  else
    gf_plot(mfvm, VM, 'deformed_mesh', 'off', 'deformation', U, ...
            'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
    xlabel('x'); ylabel('y');
    % title('Deformed configuration (not really a small deformation of course ...)');
    % gf_colormap('chouette');
    gg = [ .7 .9 .4; .5 .9 .3;   .3 .8 .2;    .1 .7 .4;     .2 0.7 1.0000; .3 0.3 1.0000;
	       1.0 .8 .1;  1.0 .6 .1;   1.0 .45 .1;   1.0 0.3 .1];
    r = reshape(repmat(gg',6,1),3,60)';
    colormap(r);
    % caxis([0 3]);
    % axis([-11 11 -1 36]); 
  end;

  % colorbar;
  axesobj = findobj('type', 'axes'); set(axesobj, 'fontname', 'times');
  set(axesobj, 'fontunits', 'points'); set(axesobj, 'fontsize', 18);
  set(axesobj, 'fontweight', 'bold'); set(axesobj, 'linewidth', 2);
  pause(1); print(gcf,'-dpng','-r300', 'deformation.png');
  pause(0.1);
end;
