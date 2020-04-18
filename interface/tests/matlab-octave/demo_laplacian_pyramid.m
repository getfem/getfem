% Copyright (C) 2005-2020 Yves Renard, Julien Pommier.
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

% Options for prescribing the Dirichlet condition
clear all;
dirichlet_version = 1; % 0 = simplification, 1 = with multipliers,
                       % 2 = penalization,  3 = Nitsche's method
theta = 1;       % Nitsche's method parameter theta
gamma0 = 0.001;  % Nitsche's method parameter gamma0 (gamma = gamma0*h)
r = 1e8;         % Penalization parameter
draw = true;
draw_mesh = false;
with_pyramids = true;

K = 2;           % Degree of the finite element method

asize =  size(who('automatic_var654'));
if (asize(1)) draw = false; end;

% trace on;
NX = 10;
if (with_pyramids)
  m = gf_mesh('pyramidal', [0:1/NX:1], [0:1/NX:1], [0:1/NX:1]);
else
  m = gf_mesh('regular simplices', [0:1/NX:1], [0:1/NX:1], [0:1/NX:1]);
end

% Create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% Assign the pyramidal fem to all convexes of the mesh_fem, and define an
% integration method
if (with_pyramids)
  gf_mesh_fem_set(mf,'fem',gf_fem(sprintf('FEM_PYRAMID_LAGRANGE(%d)', K)));
  % mim = gf_mesh_im(m, gf_integ('IM_PYRAMID_COMPOSITE(IM_TETRAHEDRON(5))'));
  mim = gf_mesh_im(m, gf_integ('IM_PYRAMID(IM_GAUSS_PARALLELEPIPED(3,5))'));
  % mim = gf_mesh_im(m, gf_integ('IM_PYRAMID_COMPOSITE(IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),3))'));
else
  gf_mesh_fem_set(mf,'fem',gf_fem(sprintf('FEM_PK(3,%d)', K)));
  mim = gf_mesh_im(m, gf_integ('IM_TETRAHEDRON(5)'));
end


% Detect the border of the mesh
border = gf_mesh_get(m, 'outer faces');
% Mark it as boundary GAMMAD=1
GAMMAD=1;
gf_mesh_set(m, 'region', GAMMAD, border);
if (draw_mesh)
  figure(1);
  gf_plot_mesh(m, 'regions', [GAMMAD], 'convexes', 'on'); % the boundary edges appear in red
  pause(1);
end

% Interpolate the exact solution
% Uexact = gf_mesh_fem_get(mf, 'eval', { '10*y.*(y-1).*x.*(x-1)+10*x.^5' });
% Uexact = gf_mesh_fem_get(mf, 'eval', { 'x.*sin(2*pi*x).*sin(2*pi*y)' });
Uexact = gf_mesh_fem_get(mf, 'eval', { 'sin(pi*x/2).*sin(pi*y/2).*sin(pi*z/2)' });
% Opposite of its laplacian
% F      = gf_mesh_fem_get(mf, 'eval', { '-(20*(x.^2+y.^2)-20*x-20*y+200*x.^3)' });
% F      = gf_mesh_fem_get(mf, 'eval', { '4*pi*(2*pi*x.*sin(2*pi*x) - cos(2*pi*x)).*sin(2*pi*y)' });
F = gf_mesh_fem_get(mf, 'eval', { '3*pi*pi*sin(pi*x/2).*sin(pi*y/2).*sin(pi*z/2)/4' });

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add linear term', mim, 'Grad_u.Grad_Test_u');
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
gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');
if (draw)
  figure(2);
  subplot(2,1,1); gf_plot(mf,U,'mesh','off', 'cvlst', gf_mesh_get(m, 'outer faces'), 'refine', 4); 
  colorbar; title('computed solution');

  subplot(2,1,2); gf_plot(mf,Uexact,'mesh','on', 'cvlst', gf_mesh_get(m, 'outer faces'), 'refine', 4); 
  colorbar;title('exact solution');
end

err = gf_compute(mf, Uexact-U, 'H1 norm', mim);

disp(sprintf('H1 norm of error: %g', err));

% M = gf_asm('mass matrix', mim, mf);
% K = gf_asm('laplacian', mim, mf, mf, ones(1, gf_mesh_fem_get(mf, 'nbdof')));

if (0) % Drawing the shape functions on the reference element
  m2 = gf_mesh('empty', 3);
  gf_mesh_set(m2, 'add convex', gf_geotrans('GT_PYRAMID(1)'), [-1 -1 0;  1, -1, 0; -1,  1, 0;  1,  1, 0;  0,  0, 1]');
  % gf_mesh_set(m2, 'add convex', gf_geotrans('GT_PYRAMID(1)'), [-1 -1 2;  1, -1, 2; -1,  1, 2;  1,  1, 2;  0,  0, 1]');
  % gf_mesh_set(m2, 'add convex', gf_geotrans('GT_PYRAMID(1)'), [ 1 -1 0; 1,  1, 0;   1, -1, 2;  1,  1, 2;  0,  0, 1]');
  mf2 = gf_mesh_fem(m2,1);
  gf_mesh_fem_set(mf2,'fem',gf_fem('FEM_PYRAMID_LAGRANGE(2)'));
  Utest = zeros(1,gf_mesh_fem_get(mf2, 'nbdof'));
  % gf_mesh_fem_get(mf2, 'basic dof nodes')
  % mim2 = gf_mesh_im(m2, gf_integ('IM_PYRAMID_COMPOSITE(IM_TETRAHEDRON(2))'));
  mim2 = gf_mesh_im(m2, gf_integ('IM_PYRAMID(IM_GAUSS_PARALLELEPIPED(3,7))'));
  % mim2 = gf_mesh_im(m2, gf_integ('IM_PYRAMID_COMPOSITE(IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),5))'));
  format long
  M2 = gf_asm('mass matrix', mim2, mf2)
  K2 = gf_asm('laplacian', mim2, mf2, mf2, ones(1, gf_mesh_fem_get(mf2, 'nbdof')))
  for i = 1:size(Utest,2)
    Utest(i) = 1;
    figure(3);
    gf_plot(mf2,Utest,'mesh','on', 'cvlst', gf_mesh_get(m2, 'outer faces'), 'refine', 8); 
    colorbar; title('shape function');
    pause;
    Utest(i) = 0;
  end
end

if (0) % Drawing the shape functions on the whole mesh
  Utest = U*0;
  for i = 1:size(U,2)
    Utest(i) = 1;
    figure(3);
    gf_plot(mf,Utest,'mesh','on', 'cvlst', gf_mesh_get(m, 'outer faces'), 'refine', 8); 
    colorbar; title('shape function');
    Utest(i) = 0;
    pause;
  end
end

if (err > 0.033)
   error('Laplacian test: error to big');
end


