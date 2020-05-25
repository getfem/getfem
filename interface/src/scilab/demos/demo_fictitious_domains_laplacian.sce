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

printf('This demo use levelset to impose (weakly) a Dirichlet condition on a part of an \n');
printf('implicit boundary defined by the zero of the levelset and a Neumann condition on \n');
printf('the remaining part of that boundary. A Poisson problem\n');

gf_workspace('clear all');
NX=10;
N = 2;
ls_degree = 1;
R = 0.4;

if (N == 3) then
  m = gfMesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
  //m = gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
  mfu0 = gfMeshFem(m,1);
  mf_mult = gfMeshFem(m,1);
  set(mfu0, 'fem', gf_fem('FEM_QK(3,2)'));
  set(mf_mult, 'fem', gf_fem('FEM_QK(3,1)'));
  adapt_im = 'IM_TETRAHEDRON(6)'
elseif (N == 2) then
  m = gfMesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
  //m = gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
  mfu0 = gfMeshFem(m,1);
  mf_mult = gfMeshFem(m,1);
  set(mfu0, 'fem', gf_fem('FEM_QK(2,2)'));
  set(mf_mult, 'fem', gf_fem('FEM_QK(2,1)'));
  adapt_im = 'IM_TRIANGLE(6)'
else 
  error('Wrong dimension');
end

ls  = gf_levelset(m, ls_degree);
ls2 = gf_levelset(m, ls_degree, 'with_secondary');

mf_ls = gf_levelset_get(ls, 'mf');
P = gf_mesh_fem_get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);
if (N == 3) then
  z = P(3,:);
else
  z = 0 * x;
end
ULS = ((x.^2 + y.^2 + z.^2).^1.5 - R^3);
ULS2 = x;
gf_levelset_set(ls, 'values', ULS);
gf_levelset_set(ls2, 'values', ULS, ULS2);


mls = gfMeshLevelSet(m);
set(mls, 'add', ls);
set(mls, 'add',ls2);
set(mls, 'adapt');
mim_bound2 = gfMeshIm('levelset', mls, 'boundary(a)', gf_integ(adapt_im));
mim_bound  = gfMeshIm('levelset', mls, 'boundary(b)', gf_integ(adapt_im));
mim_int    = gfMeshIm('levelset', mls, 'inside(a)',   gf_integ(adapt_im));
set(mim_int, 'integ', 4);

// Some verifications
A1 = gf_asm('volumic','V()+=comp()',mim_bound);
A2 = gf_asm('volumic','V()+=comp()',mim_bound2);
V  = gf_asm('volumic','V()+=comp()',mim_int);
if (N == 2) then
  disp(sprintf('length : %g should be %g', A1, %pi*R));
  disp(sprintf('length : %g should be %g', A2, 2*%pi*R));
  disp(sprintf('area   : %g should be %g', V, %pi*R^2));
else
  disp(sprintf('area   : %g should be %g', A1, 4*%pi*R^2/2));
  disp(sprintf('area   : %g should be %g', A2, 4*%pi*R^2));
  disp(sprintf('volume : %g should be %g', V,  4*%pi*R^3/3.));   
end

// partial mesh fem
dof_out = get(mfu0, 'dof from im', mim_int);
cv_out  = get(mim_int, 'convex_index');
cv_in   = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu     = gfMeshFem('partial', mfu0, dof_out, cv_in);

// data
if (N == 2) then
  Volumic_data = gf_mesh_fem_get_eval(mfu0, list('45*sqrt(x.^2+y.^2)'));
  surface_data = gf_mesh_fem_get_eval(mfu0, list('-15*(x.^2+y.^2)'));
  Sol_U = gf_mesh_fem_get_eval(mfu0, list(sprintf('5*((%g)^3-(x.^2+y.^2).^1.5)', R)));
else
  Volumic_data = gf_mesh_fem_get_eval(mfu0, list('60*sqrt(x.^2+y.^2+z.^2)'));
  surface_data = gf_mesh_fem_get_eval(mfu0, list('-15*(x.^2+y.^2+z.^2)'));
  Sol_U = gf_mesh_fem_get_eval(mfu0, list(sprintf('5*((%g)^3-(x.^2+y.^2+z.^2).^1.5)', R)));
end

// getfem model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add Laplacian brick', mim_int, 'u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mfu0, Volumic_data);
gf_model_set(md, 'add source term brick', mim_int, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'SurfaceData', mfu0, surface_data);
gf_model_set(md, 'add source term brick', mim_bound2, 'u', 'SurfaceData');
gf_model_set(md, 'add initialized fem data', 'SurfaceData2', mfu0, -surface_data);
gf_model_set(md, 'add source term brick', mim_bound, 'u', 'SurfaceData2');
gf_model_set(md, 'add multiplier', 'mult_dir', mf_mult, 'u');
gf_model_set(md, 'add Dirichlet condition with multipliers', ...
	     mim_bound, 'u', 'mult_dir', -1);

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

// Comparison with the exaxt solution
ERRL2 = gf_compute(mfu, U, 'L2 dist',      mim_int, mfu0, Sol_U);
ERRH1 = gf_compute(mfu, U, 'H1 semi dist', mim_int, mfu0, Sol_U);
disp(sprintf('L2 error= %g', ERRL2));
disp(sprintf('semi H1 error = %g', ERRH1));

hh = scf();
hh.color_map = gf_colormap('chouette');
if (N == 2) then
  gf_plot(mfu, U, 'mesh','on', 'refine', 2);
else
  gf_plot(mfu, U, 'mesh','on', 'cvlst', gf_mesh_get(m, 'outer faces'), 'refine', 2);
end
