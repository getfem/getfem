% Copyright (C) 2005-2012 Julien Pommier.
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

% Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
%mim = gf_mesh_im(m, gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,5),4)'));
% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);
if (draw)
  gf_plot_mesh(m, 'regions', [1]); % the boundary edges appears in red
  pause(1);
end

% interpolate the exact solution
Uexact = gf_mesh_fem_get(mf, 'eval', { 'y.*(y-1).*x.*(x-1)+x.^5' });
% its second derivative
F      = gf_mesh_fem_get(mf, 'eval', { '-(2*(x.^2+y.^2)-2*x-2*y+20*x.^3)' });

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
switch (dirichlet_version)
  case 0,
    gf_model_set(md, 'add Dirichlet condition with simplification', 'u', 1, 'DirichletData');   
  case 1, 
    gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 1, 'DirichletData');
  case 2,
    gf_model_set(md, 'add Dirichlet condition with penalization', mim, 'u', r, 1, 'DirichletData');
  case 3,
    gf_model_set(md, 'add initialized data', 'gamma0', [gamma0]);
    gf_model_set(md, 'add Dirichlet condition with Nitsche method', mim, 'u', 'gamma0', 1, theta, 'DirichletData');
end
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

if (err > 2E-5)
   error('Laplacian test: error to big');
end


