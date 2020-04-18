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

% The main interest of this file is to show how interpolate transformation
% can be used to prescribe a periodic boundary condition


% Options for prescribing the Dirichlet condition
dirichlet_version = 0; % 0 = simplification, 1 = with multipliers, 2 = penalization,  3 = Nitsche's method
theta = 1;       % Nitsche's method parameter theta
gamma0 = 0.001;  % Nitsche's method parameter gamma0 (gamma = gamma0*h)
r = 1e8;         % Penalization parameter
draw = true;

asize =  size(who('automatic_var654'));
if (asize(1)) draw = false; end;

% trace on;
gf_workspace('clear all');
NX = 20;
m = gf_mesh('cartesian',[0:1/NX:1],[0:1/NX:1]);
%m=gf_mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% assign the Q2 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

% create a mesh_fem of for a field of dimension 1 for the multiplier
mf_lambda = gf_mesh_fem(m,1);
% assign the Q1 fem to all elements of the mesh_fem (only the boundary will be used),
gf_mesh_fem_set(mf_lambda,'fem',gf_fem('FEM_QK(2,1)'));

% Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
%mim = gf_mesh_im(m, gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,5),4)'));
% detect the border of the mesh
GAMMAD = 1;
GAMMAP = 2;
border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
dirichlet_boundary=border(:, find(abs(normals(1, :)) < 0.01));
gf_mesh_set(m, 'region', GAMMAD, dirichlet_boundary);
periodic_boundary=border(:, find(normals(1, :) > 0.01));
gf_mesh_set(m, 'region', GAMMAP, periodic_boundary);

if (draw)
  gf_plot_mesh(m, 'regions', [1]); % , 'convexes', 'on'); % the boundary edges appears in red
  pause(1);
end

% interpolate the exact solution
Uexact = gf_mesh_fem_get(mf, 'eval', { 'cos(2*pi*(x+0.1)).*sin(2*pi*y)' });

% its second derivative
F      = gf_mesh_fem_get(mf, 'eval', { '8*pi*pi*cos(2*pi*(x+0.1)).*sin(2*pi*y)' });

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
% gf_model_set(md, 'add linear term', mim, 'Grad_Test2_u.Grad_Test_u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
switch (dirichlet_version)
  case 0,
    gf_model_set(md, 'add Dirichlet condition with simplification', 'u', GAMMAD, 'DirichletData');   
  case 1, 
    gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, GAMMAD, 'DirichletData');
  case 2,
    gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', r, GAMMAD, 'DirichletData');
  case 3,
    gf_model_set(md, 'add initialized data', 'gamma0', [gamma0]);
    gf_model_set(md, 'add Dirichlet condition with Nitsche method', mim, 'u', 'gamma0', GAMMAD, theta, 'DirichletData');
end
% periodic condition
gf_model_set(md, 'add filtered fem variable', 'lambda', mf_lambda, GAMMAP); % multiplier for the periodic condition
gf_model_set(md, 'add interpolate transformation from expression', 'transform', m, m, 'X-[1;0]');
gf_model_set(md, 'add linear term', mim, '(Interpolate(u,transform)-u)*lambda', GAMMAP);
gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

if (draw)
  subplot(2,1,1); gf_plot(mf,U,'mesh','on','contour',.01:.01:.1); 
  colorbar; title('computed solution');

  subplot(2,1,2); gf_plot(mf,U-Uexact,'mesh','on'); 
  colorbar;title('difference with exact solution');
end

err = gf_compute(mf, U-Uexact, 'H1 norm', mim);

disp(sprintf('H1 norm of error: %g', err));

if (err > 1E-3)
   error('Laplacian test: error to big');
end


