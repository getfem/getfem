// Copyright (C) 2010-2020 Yann COLLETTE.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 3 of the License,  or
// (at your option) any later version along with the GCC Runtime Library
// Exception either version 3.1 or (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License and GCC Runtime Library Exception for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

m  = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);
mf = gf_mesh_fem(m,1); // create a meshfem of for a field of dimension 1 (i.e. a scalar field)
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

disp(gf_fem_get(gf_fem('FEM_QK(2,2)'), 'poly_str'));

mim = gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));

border = gf_mesh_get(m,'outer faces');
gf_mesh_set(m, 'region', 42, border); // create the region

// the boundary edges appears in red
h = scf();
drawlater;
gf_plot_mesh(m, 'regions', [42], 'vertices','on','convexes','on'); 
drawnow;

md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');

R = gf_mesh_fem_get_eval(mf, list(list('(x-.5).^2 + (y-.5).^2 + x/5 - y/3')));
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, R);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');
gf_model_get(md, 'variable list');
gf_model_get(md, 'solve');

U = gf_model_get(md, 'variable', 'u');

h = scf();
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mf, U, 'mesh','on');
h.color_map = jetcolormap(255);
drawnow;
