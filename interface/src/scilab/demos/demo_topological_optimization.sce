// Copyright (C) 2009 Alassane SY, Yves Renard.
// Copyright (C) 2010 Yann Collette.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//  Shape optimization with topological gradient
//   (with a fictitious domain approach).
//
//  This program is used to check that scilab-getfem is working. This is
//  also a good example of use of GetFEM.
//

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_topological_optimization.sce');

printf('demo topological_optimization started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

// parameters
NX          = 30;
ls_degree   = 1;
alpha       = 1;
beta        = 1;
rayon_trous = 0.2;
Index       = 1;

// Finite element and integration methods definition
m = gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
// m=gf_mesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
// mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
mf_basic = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mf_basic,'fem',gf_fem('FEM_QK(2,2)'));
mls = gf_mesh_levelset(m);
ls  = gf_levelset(m, ls_degree);
gf_mesh_levelset_set(mls, 'add', ls);
mf_ls = gf_levelset_get(ls, 'mf');
P     = gf_mesh_fem_get(mf_ls, 'basic dof nodes');

x   = P(1,:);
y   = P(2,:);
ULS = 1000*ones(1,length(x));

h = scf();
h.color_map = jetcolormap(255);
  
// Loop on the topological optimization
while(1) 
  gf_workspace('push');
  gf_levelset_set(ls, 'values', ULS);
  gf_mesh_levelset_set(mls, 'adapt');
  mim_bound = gf_mesh_im('levelset', mls, 'boundary', gf_integ('IM_TRIANGLE(6)'));
  mim       = gf_mesh_im('levelset', mls, 'outside',  gf_integ('IM_TRIANGLE(6)'));
  gf_mesh_im_set(mim, 'integ', 4);
  mf_mult = gf_mesh_fem(m);
  gf_mesh_fem_set(mf_mult, 'fem', gf_fem('FEM_QK(2,1)'));
  M   = gf_asm('mass matrix', mim, mf_basic);
  D   = abs(full(diag(M)));
  ind = find(D > 1E-8);
  mf  = gf_mesh_fem('partial', mf_basic, ind);
  S   = gf_asm('volumic','V()+=comp()',mim);
  printf('remaining surface:'); disp(S);

  // Problem definition (Laplace(u) + u = f)
  md = gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mf);
  gf_model_set(md, 'add Laplacian brick', mim, 'u');
  gf_model_set(md, 'add fem data', 'VolumicData', mf_basic);
  gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
  gf_model_set(md, 'add initialized data', 'rho', [1.]);
  gf_model_set(md, 'add mass brick', mim, 'u', 'rho');
  gf_model_set(md, 'add multiplier', 'mult_dir', mf_mult, 'u');

  // To be completely robust, a stabilization should be used on the Dirichlet
  // boundary to ensure the inf-sup condition (Nitsche or Barbosa-Hughes)
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim_bound, 'u', 'mult_dir', -1);
  // Solving the direct problem.
  U0 = gf_mesh_fem_get_eval(mf_basic, ...
                            list(list('0.4*(3.*sin(%pi*(x+y)) + ((x-0.5).^10 + (y-0.5).^10 + (x+0.5).^10 + (y+0.5).^10))')));
  gf_model_set(md, 'variable', 'VolumicData', U0);
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');

  drawlater;
  clf();
  subplot(2,1,1);
  gf_plot(mf, U);
  U_tmp = gf_levelset_get(ls,'values');
  [h1,h2]=gf_plot(mf_ls, U_tmp, 'contour', 0,'pcolor','off');
  // set(h2{1},'LineWidth',2);
  // set(h2{1},'Color','green');
  colorbar(min(U_tmp),max(U_tmp));
  title('u');

  // Solving the adjoint problem.
  UBASIC = gf_compute(mf, U, 'interpolate on', mf_basic);
  F = 2*(UBASIC-U0);
  gf_model_set(md, 'variable', 'VolumicData', F);
  gf_model_get(md, 'solve');
  W = gf_model_get(md, 'variable', 'u');

  // Computation of the topological gradient
  mf_g=gf_mesh_fem(m, 1);
  gf_mesh_fem_set(mf_g,'fem', ...
                  gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'));
  DU = gf_compute(mf, U, 'gradient', mf_g);
  DW = gf_compute(mf, W, 'gradient', mf_g);
  nbdof = gf_mesh_fem_get(mf_g, 'nbdof');
  DU = matrix(DU, 2, nbdof);
  DW = matrix(DW, 2, nbdof);
  UU = gf_compute(mf, U, 'interpolate on', mf_g);
  UU0 = gf_compute(mf_basic, U0, 'interpolate on', mf_g);
  LS = gf_compute(mf_ls, ULS, 'interpolate on', mf_g);
  G = (-4*%pi*( alpha*(DU(1,:).^2 + DU(2,:).^2 + DU(1,:).*DW(1,:) + ...
      DU(2,:).*DW(2,:)) + beta*(UU-UU0).^2)) .* (sign(LS)+1.)/2;

  subplot(2,1,2);
  gf_plot(mf_g, G);
  title('Topological gradient');
  colorbar(min(G),max(G));
  drawnow;
  xs2png(h.figure_id, path + sprintf('/topological_opt%03d.png',Index));
  Index = Index + 1;
  sleep(10);

  // Find the point where the topological gradient is minimum
  [val, i] = min(G);
  if (val >= -12) then
    disp('Topological optimization finished.');
    return;
  end

  point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);
  gf_workspace('pop');

  // Updating the level set to add the hole
  R = -(val+7) / 200;
  xc = point(1);
  yc = point(2);
  ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);
end
 
printf('demo topological_optimization terminated\n');
