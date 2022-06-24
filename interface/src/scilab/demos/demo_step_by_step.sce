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

path = get_absolute_file_path('demo_step_by_step.sce');

printf('demo step_by_step started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

// creation of a simple cartesian mesh
m = gf_mesh('cartesian', 0:0.1:1.1, 0:0.1:1.1);

// create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m, 1);

// assign the Q2 fem to all convexes of the MeshFem
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

// view the expression of its basis functions on the reference convex
printf('The expression of its basis functions on the reference convex\n');
disp(gf_fem_get(gf_fem('FEM_QK(2,2)'),'poly_str')); 

// an exact integration will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));

// detect the border of the mesh
border = gf_mesh_get(m,'outer_faces');

// mark it as boundary #42
gf_mesh_set(m,'region',42, border);


// empty real model
md = gf_model('real');

// declare that "u" is an unknown of the system
// on the finite element method `mf`
gf_model_set(md, 'add fem variable', 'u', mf);

// add generic elliptic brick on "u"
gf_model_set(md, 'add Laplacian brick', mim, 'u');

// add Dirichlet condition
Uexact = gf_mesh_fem_get_eval(mf, list('(x-.5).^2 + (y-.5).^2 + x/5 - y/3'));
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');

// add source term
f = gf_mesh_fem_get_eval(mf, list('2*(x.^2+y.^2)-2*(x+y)+20*x.^3'));
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, f);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

// solve the linear system
gf_model_get(md, 'solve');

// extracted solution
u = gf_model_get(md, 'variable', 'u');


// export computed solution
gf_mesh_fem_get(mf,'export_to_pos', path + '/sol.pos',u,'Computed solution');

// display
hh = scf();
hh.color_map = jetcolormap(255);
gf_plot(mf, u, 'mesh','on');

printf('demo step_by_step terminated\n');
