% Copyright (C) 2015-2015 Yves Renard.
%
% This file is a part of GetFEM++
%
% GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

% A Poisson problem solved with a Discontinuous Galerkin method (or
% interior penalty method).

% Options for prescribing the Dirichlet condition
dirichlet_version = 0; % 0 = simplification, 1 = with multipliers,
                       % 2 = penalization,  3 = Nitsche's method
theta = 1;       % Nitsche's method parameter theta
gamma0 = 0.001;  % Nitsche's method parameter gamma0 (gamma = gamma0*h)
r = 1e8;         % Penalization parameter
draw = true;
quadrangles = true;
K = 2;           % Degree of the discontinuous finite element method
interior_penalty_factor = 1E7; % Parameter of the interior penalty term

asize =  size(who('automatic_var654'));
if (asize(1)) draw = false; end;

% trace on;
gf_workspace('clear all');
NX = 20;
if (quadrangles)
  m = gf_mesh('cartesian',[0:1/NX:1],[0:1/NX:1]);
else
  m=gf_mesh('import','structured',sprintf('GT="GT_PK(2,1)";SIZES=[1,1];NOISED=0;NSUBDIV=[%d,%d];', NX, NX));
end

% Create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% Assign the QK or PK fem to all convexes of the mesh_fem, and define an
% integration method
if (quadrangles)
  gf_mesh_fem_set(mf,'fem',gf_fem(sprintf('FEM_QK_DISCONTINUOUS(2,%d)', K)));
  mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,6)'));
else
  gf_mesh_fem_set(mf,'fem',gf_fem(sprintf('FEM_PK_DISCONTINUOUS(2,%d)', K)));
  mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(6)'));
end

% Detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
GAMMAD=1;
gf_mesh_set(m, 'boundary', GAMMAD, border);
% Inner edges for the interior penalty terms
in_faces = gf_mesh_get(m,'inner faces');
INNER_FACES=2;
gf_mesh_set(m, 'boundary', INNER_FACES, in_faces);

% Mesh plot
if (draw)
  subplot(2,2,1); gf_plot_mesh(m, 'regions', [GAMMAD]);  % the boundary edges appear in red
  subplot(2,2,2); gf_plot_mesh(m, 'regions', [INNER_FACES]); % inner edges appear in red
  pause(1);
end

% Interpolate the exact solution
% Uexact = gf_mesh_fem_get(mf, 'eval', { '10*y.*(y-1).*x.*(x-1)+10*x.^5' });
Uexact = gf_mesh_fem_get(mf, 'eval', { 'x.*sin(2*pi*x).*sin(2*pi*y)' });
% Opposite of its laplacian
% F      = gf_mesh_fem_get(mf, 'eval', { '-(20*(x.^2+y.^2)-20*x-20*y+200*x.^3)' });
F      = gf_mesh_fem_get(mf, 'eval', { '4*pi*(2*pi*x.*sin(2*pi*x) - cos(2*pi*x)).*sin(2*pi*y)' });

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add linear generic assembly brick', mim, 'Grad_u.Grad_Test_u');
% gf_model_set(md, 'add Laplacian brick', mim, 'u');
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
    expr = gf_model_get(md, 'Neumann term', 'u', GAMMAD);
    gf_model_set(md, 'add Dirichlet condition with Nitsche method', mim, 'u', expr, 'gamma0', GAMMAD, theta, 'DirichletData');
end

% Interior penalty terms
gf_model_set(md, 'add initialized data', 'alpha', [interior_penalty_factor]);
gf_model_set(md, 'add linear generic assembly brick', mim, '(u-Interpolate(u,neighbour_elt))*(Grad_Test_u.Normal)/2 + (u-Interpolate(u,neighbour_elt))*(Interpolate(Grad_Test_u,neighbour_elt).Normal)/2', INNER_FACES);
gf_model_set(md, 'add linear generic assembly brick', mim, '-Test_u*(Grad_u.Normal + Interpolate(Grad_u,neighbour_elt).Normal)/2+Interpolate(Test_u,neighbour_elt)*(Grad_u.Normal + Interpolate(Grad_u,neighbour_elt).Normal)/2', INNER_FACES);
gf_model_set(md, 'add linear generic assembly brick', mim, 'alpha*(u-Interpolate(u,neighbour_elt))*Test_u - alpha*(u-Interpolate(u,neighbour_elt))*Interpolate(Test_u,neighbour_elt)', INNER_FACES);
% gf_model_set(md, 'add linear generic assembly brick', mim, 'alpha*(u-Interpolate(u,neighbour_elt))*(Test_u-Interpolate(Test_u,neighbour_elt))', INNER_FACES);


gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

if (draw)
  subplot(2,2,3); gf_plot(mf,U,'mesh','off'); 
  colorbar; title('computed solution');

  subplot(2,2,4); gf_plot(mf,Uexact-U,'mesh','on'); 
  colorbar;title('difference with exact solution');
end

err = gf_compute(mf, Uexact-U, 'H1 norm', mim);

disp(sprintf('H1 norm of error: %g', err));

if (err > 2E-4)
   error('Laplacian test: error to big');
end


