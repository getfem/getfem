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


m = gf_mesh('cartesian', [0:.1:1], [0:.1:1]);
mf = gf_mesh_fem(m,1); % create a meshfem of for a field of dimension 1 (i.e. a scalar field)
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

disp(gf_fem_get(gf_fem('FEM_QK(2,2)'), 'poly_str'));

mim=gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2, 4)'));
% mim=gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)')); % not allowed with the high level generic assembly

border = gf_mesh_get(m,'outer faces');
gf_mesh_set(m, 'region', 42, border); % create the region (:#42

% the boundary edges appears in red
gf_plot_mesh(m, 'regions', [42], 'vertices','on','convexes','on'); 


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
% gf_model_set(md, 'add linear term', mim, 'Grad_u.Grad_Test_u');
gf_model_set(md, 'add Laplacian brick', mim, 'u');
R=gf_mesh_fem_get(mf, 'eval', {'(x-.5).^2 + (y-.5).^2 + x/5 - y/3'});
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, R);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');

gf_model_get(md, 'variable list');

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');
gf_plot(mf, U, 'mesh','on');
