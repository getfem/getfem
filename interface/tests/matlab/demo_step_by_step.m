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

m = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);

% we enable vertices and convexes labels
gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

% assign the same integration method on all elements
mim=gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));

% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #42
gf_mesh_set(m, 'region', 42, border);
gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red

% empty real model
md = gf_model('real');

% declare that "u" is an unknown of the system
% on the finite element method `mf`
gf_model_set(md, 'add fem variable', 'u', mf);

% add generic elliptic brick on "u"
gf_model_set(md, 'add Laplacian brick', mim, 'u');

% add Dirichlet condition
Uexact = gf_mesh_fem_get(mf, 'eval', {'(x-.5).^2 + (y-.5).^2 + x/5 - y/3'});
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');

% add source term
f = gf_mesh_fem_get(mf, 'eval', { '2*(x.^2+y.^2)-2*(x+y)+20*(x.^3)' });
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, f);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

% solve the linear system
gf_model_get(md, 'solve');

% extracted solution
u = gf_model_get(md, 'variable', 'u');
% display
gf_plot(mf, u, 'mesh','on');
