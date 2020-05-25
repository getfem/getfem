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

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_laplacian.sce');

printf('demo laplacian started\n');

// trace on;

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

m = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);
//m = gf_mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')

// create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);

// assign the Q2 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

// Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
//mim = gf_mesh_im(m, gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,5),4)'));

// detect the border of the mesh
border = gf_mesh_get(m,'outer faces');

// mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);

h = scf();
drawlater;
gf_plot_mesh(m, 'regions', [1]); // the boundary edges appears in red
drawnow;

// interpolate the exact solution
Uexact = gf_mesh_fem_get_eval(mf, list(list('y.*(y-1).*x.*(x-1)+x.^5')));

// its second derivative
F      = gf_mesh_fem_get_eval(mf, list(list('-(2*(x.^2+y.^2)-2*x-2*y+20*x.^3)')));

md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 1, 'DirichletData');

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

disp(sprintf('H1 norm of error: %g', gf_compute(mf,U-Uexact,'H1 norm',mim)));

h = scf();
drawlater;
h.color_map = jetcolormap(255);
subplot(2,1,1); 
gf_plot(mf,U,'mesh','on','contour',.01:.01:.1); 
colorbar(min(U),max(U)); 
title('computed solution');

subplot(2,1,2); 
gf_plot(mf,U-Uexact,'mesh','on'); 
colorbar(min(U-Uexact),max(U-Uexact));
title('difference with exact solution');
h.color_map = jetcolormap(255);
drawnow;

printf('demo laplacian terminated\n');
