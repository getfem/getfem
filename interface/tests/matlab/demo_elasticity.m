% Copyright (C) 2005-2020 Julien Pommier.
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

clear all;
% parameters
option = 2;            % 1 = use Lame coefficients brick
                       % 2 = use plane strain approximation brick
                       % 3 = use plane stress approximation brick
d = 2;                 % dimension (cannot be changed for the moment)
E = 1;                 % Young's modulus
nu = 0.3;              % Poisson ratio
cmu = E/(2*(1+nu));    % Lame coefficient
if (option == 3)
  clambda = 2*cmu*nu/(1-2*nu); % Lame coefficient, 3D and plane strain
else
  clambda = 2*cmu*nu/(1-nu);   % Lame coefficient, plane stress
end
dirichlet_version = 2; % 1 = With multipliers, 2 = Nitsche's method
theta = 1;             % Nitsche's method parameter theta
gamma0 = 0.001;        % Nitsche's method parameter gamma0 (gamma = gamma0*h)
incompressible = false;% Test with incompressibility or not
NX = 30;

% trace on;
gf_workspace('clear all');
m = gf_mesh('cartesian',[0:1/NX:1],[0:1/NX:1]);
%m=gf_mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')


% create a mesh_fem of for a field of dimension d (i.e. a vector field)
mf = gf_mesh_fem(m,d);
% assign the Q2 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'fem', gf_fem('FEM_QK(2,2)'));

if (incompressible) 
  mfp = gf_mesh_fem(m,1);
  gf_mesh_fem_set(mfp, 'fem', gf_fem('FEM_QK(2,1)'));
end

mf_H = gf_mesh_fem(m,1);
gf_mesh_fem_set(mf_H, 'fem', gf_fem('FEM_QK(2,1)'));

mfdu=gf_mesh_fem(m,1); gf_mesh_fem_set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));

% Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
%mim = gf_mesh_im(m, gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,5),4)'));
% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
GAMMAD = 1;
gf_mesh_set(m, 'boundary', GAMMAD, border);
% gf_plot_mesh(m, 'regions', [1]); % the boundary edges appears in red
% pause(1);


% Polynomial exact solution
% Uexact = gf_mesh_fem_get(mf, 'eval', { '(x.^4).*(y.^2)', 'x.*y'});
% F = gf_mesh_fem_get(mf, 'eval', {sprintf('-(%g)*(12*(x.^2).*(y.^2)+1) - (%g)*(24*(x.^2).*(y.^2)+2*(x.^4)+1)', clambda, cmu), sprintf('-8*((%g)+(%g))*((x.^3).*y)', clambda, cmu)});

% Exact incompressible solution in terms of trigonometric functions
a = 8;
Uexact = gf_mesh_fem_get(mf, 'eval', { sprintf('sin((%g)*x)', a), sprintf('-(%g)*y.*cos((%g)*x)',a,a)});
F = gf_mesh_fem_get(mf, 'eval', {sprintf('(%g)*((%g)^2)*sin((%g)*x)', cmu, a, a), sprintf('-(%g)*((%g)^3)*y.*cos((%g)*x)', cmu, a, a)});
if (incompressible)
  Pexact = 100*gf_mesh_fem_get(mfp, 'eval', {sprintf('sin((%g)*(x-y))', a)});
  F = F + 100*gf_mesh_fem_get(mf, 'eval', {sprintf('(%g)*cos((%g)*(x-y))', a, a), sprintf('-(%g)*cos((%g)*(x-y))', a, a)});
end

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
if (option == 1)
  gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
  gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
  % gf_model_set(md, 'add linear term', mim, ...
  %       '(clambda*Div_Test_u*Id(qdim(u)) +
  %       cmu*(Grad_Test_u''+Grad_Test_u)):Grad_Test2_u');
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'clambda', 'cmu');
elseif (option == 2)
  gf_model_set(md, 'add initialized data', 'E', [E]);
  gf_model_set(md, 'add initialized data', 'nu', [nu]);
  gf_model_set(md, 'add isotropic linearized elasticity brick pstrain', mim, 'u', 'E', 'nu');
elseif (option == 3)
  gf_model_set(md, 'add initialized data', 'E', [E]);
  gf_model_set(md, 'add initialized data', 'nu', [nu]);
  gf_model_set(md, 'add isotropic linearized elasticity brick pstress', mim, 'u', 'E', 'nu');
end

if (incompressible)
  gf_model_set(md, 'add fem variable', 'p', mfp);
  gf_model_set(md, 'add linear incompressibility brick', mim, 'u', 'p');
  % Not necessary to fix the pressure at a point ?
  % gf_model_set(md, 'add initialized data', 'cpoints', [0.5, 0.5]);
  % gf_model_set(md, 'add pointwise constraints with multipliers', 'p', 'cpoints');
end
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
if (dirichlet_version == 1)
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, GAMMAD, 'DirichletData');
else
  gf_model_set(md, 'add initialized data', 'gamma0', [gamma0]);
  expr = gf_model_get(md, 'Neumann term', 'u', GAMMAD);
  gf_model_set(md, 'add Dirichlet condition with Nitsche method', mim, 'u', expr, 'gamma0', GAMMAD, theta, 'DirichletData');
end

tic;    
gf_model_get(md, 'solve', 'noisy', 'max iter', 100);
U = gf_model_get(md, 'variable', 'u');
toc;

figure(1);
subplot(1+incompressible, 2, 1);
if (option == 1)
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'clambda', 'cmu', mfdu);
elseif (option == 2)
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises pstrain', 'u', 'E', 'nu', mfdu);
elseif (option == 3)
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises pstress', 'u', 'E', 'nu', mfdu);
end
    
gf_plot(mfdu, VM, 'deformed_mesh', 'on', 'deformation', U, 'deformation_mf', mf, 'refine', 4); 
colorbar;title('approximated solution');

subplot(1+incompressible, 2, 2);
gf_plot(mf,U-Uexact,'mesh','on', 'norm', 'on'); 
colorbar; title('difference with exact solution');

if (incompressible)
  P = gf_model_get(md, 'variable', 'p');
  P = P - (P(1) - Pexact(1));
  subplot(2, 2, 3);
  gf_plot(mfp, P);
  colorbar;title('approximated pressure');
  subplot(2, 2, 4);
  gf_plot(mfp, P-Pexact);
  colorbar;title('difference with exact pressure');
end

disp(sprintf('H1 norm of error: %g', gf_compute(mf,U-Uexact,'H1 norm',mim)));

