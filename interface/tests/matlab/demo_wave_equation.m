% Copyright (C) 2008-2020 Yves Renard.
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


% Simple demo of a wave equation solved with a Newmark scheme with the
% Getfem tool for time integration schemes

gf_workspace('clear all');
m = gf_mesh('cartesian',[0:.2:1],[0:.2:1]);
% m=gf_mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% assign the Q2 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

% Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));

% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);

% interpolate the initial data
U0 = gf_mesh_fem_get(mf, 'eval', { 'y.*(y-1).*x.*(x-1).*x.*x' });
V0 = 0*U0;

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 1);

% transient part.
T = 15.0;
dt = 0.025;
beta = 0.25;
gamma = 0.5;

gf_model_set(md, 'add Newmark scheme', 'u', beta, gamma);
gf_model_set(md, 'add mass brick', mim, 'Dot2_u');
gf_model_set(md, 'set time step', dt);

% Initial data.
gf_model_set(md, 'variable', 'Previous_u',  U0);
gf_model_set(md, 'variable', 'Previous_Dot_u',  V0);

% Initialisation of the acceleration 'Previous_Dot2_u'
gf_model_set(md, 'perform init time derivative', dt/20.);
gf_model_get(md, 'solve');

% Iterations
for t = 0:dt:T

  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');
  V = gf_model_get(md, 'variable', 'Dot_u');

  subplot(2,1,1); gf_plot(mf, U, 'mesh', 'on', 'contour', .01:.01:.1); 
  colorbar; title(sprintf('computed solution u for t=%g', t));

  subplot(2,1,2); gf_plot(mf, V, 'mesh', 'on', 'contour', .01:.01:.1); 
  colorbar; title(sprintf('computed solution du/dt for t=%g', t));
  pause(0.1);

  gf_model_set(md, 'shift variables for time integration');
end;


