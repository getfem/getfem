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

path = get_absolute_file_path('demo_fictitious_domains.sce');

printf('demo fictitious domains started\n');

printf('This demo uses levelset to impose (weakly) a Dirichlet condition on an\n');
printf('implicit boundary defined by the zero of the levelset\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

NX = 40;
ls_degree = 2;

m = gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
//m = gf_mesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
_ls = gf_levelset(m, ls_degree);
ls2 = gf_levelset(m, ls_degree, 'with_secondary');

mf_ls  = gf_levelset_get(_ls, 'mf');
mf_ls2 = gf_levelset_get(ls2, 'mf');

P = gf_mesh_fem_get(mf_ls, 'basic dof nodes');
x = P(1,:);
y = P(2,:);
//ULS = ((x + 0.25).^2 + (y - 0.4).^2) - 0.05^2;
//ULS = min(ULS, ((x - 0.25).^2 + (y - 0.4).^2) - 0.05^2);

ULS = 1000*ones(1,length(x));
rand('seed',1);

if 0 then
  for ix=1:5
    for iy=1:5
      xc = ((ix-1)/4) * 0.8 - 0.4;
      yc = ((iy-1)/4) * 0.8 - 0.4;
      if (modulo(iy,2)==0) then
	xc = xc + 0.05;
      else
	xc = xc - 0.05;
      end;
      R = 0.03 + 0.005*(iy-1);
      ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);
    end
  end
else
  for i=1:8
    xc  = rand() - 0.5;
    yc  = rand() - 0.5;
    R   = rand() * 0.09 + 0.02;
    ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);
  end
end

gf_levelset_set(_ls, 'values', ULS);

ULS2  = 1000*ones(1,length(x));
ULS2s = 1000*ones(1,length(x));

for i=1:1
  xc = 0; //rand() - 0.5;
  yc = 0.0; //rand() - 0.5;
  theta = %pi/3; //%pi*rand();
  n = [-sin(theta) cos(theta)];
  
  R = 0.19; //rand() * 0.09 + 0.02;
  ULS2 = min(ULS2, ((x-xc)*n(1) + (y-yc)*n(2)));
  //ULS2s = min(ULS2s, ((x - xc).^2 + (y - yc).^2) - R^2);
  ULS2s = min(ULS2s, (abs(y - yc)+abs(x-xc) - R));
end

gf_levelset_set(ls2, 'values', ULS2, ULS2s); //'-y-x+.2'); //, 'sqr(y-.2) - 0.04');

mls = gf_mesh_levelset(m);
gf_mesh_levelset_set(mls, 'add', _ls);
gf_mesh_levelset_set(mls, 'add', ls2);
gf_mesh_levelset_set(mls, 'adapt');

mim_bound = gf_mesh_im('levelset',mls,'boundary(a+b)', gf_integ('IM_TRIANGLE(6)')); //, gf_integ('IM_QUAD(5)'));
mim       = gf_mesh_im('levelset',mls,'all(a+b)', gf_integ('IM_TRIANGLE(6)'));
gf_mesh_im_set(mim, 'integ', 4);

mfu0 = gf_mesh_fem(m,2); gf_mesh_fem_set(mfu0, 'fem', gf_fem('FEM_QK(2,3)'));
mfdu = gf_mesh_fem(m,1); gf_mesh_fem_set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));
mf_mult = gf_mesh_fem(m,2); gf_mesh_fem_set(mf_mult, 'fem', gf_fem('FEM_QK(2,1)'));

A = gf_asm('volumic','V()+=comp()',mim_bound)

//printf('generating first plot\n');
//h = scf(); 
//h.color_map = jetcolormap(255);
//drawlater;
//gf_plot_mesh(gf_mesh_levelset_get(mls, 'cut mesh'));
//gf_plot_mesh(gf_mesh_levelset_get(mls, 'cut_mesh'), 'curved', 'on');
//gf_plot(mf_ls, ULS);
//h.color_map = jetcolormap(255);
//drawnow;

dof_out = gf_mesh_fem_get(mfu0, 'dof from im', mim);
cv_out  = gf_mesh_im_get(mim, 'convex_index');
cv_in   = _setdiff(gf_mesh_get(m, 'cvid'), cv_out);

// mfu = gf_mesh_fem('partial', mfu0, dof_out, cv_in);

md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu0);
gf_model_set(md, 'add initialized data', 'lambda', [1]);
gf_model_set(md, 'add initialized data', 'mu', [1]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
gf_model_set(md, 'add initialized data', 'VolumicData', [0; 10]);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add multiplier', 'mult_dir', mf_mult, 'u');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim_bound, 'u', 'mult_dir', -1);

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'lambda', 'mu', mfdu);

printf('generating second plot\n');
  
h = scf();
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mfdu, VM, 'deformed_mesh', 'on', 'deformation', U, 'deformation_mf', mfu0, 'refine', 8, 'cvlst', cv_out); 
//gf_plot(mfu0, U, 'norm', 'on', 'deformed_mesh', 'on', 'deformation', U,	'deformation_mf', mfu0, 'refine', 8, 'cvlst', cv_out); 
//// gf_mesh_fem_set(mfu0,'qdim',1); Unorm=sqrt(U(1:2:$).^2 + U(2:2:$).^2);
// [h1,h2] = gf_plot(mfu0, Unorm,'contour',0.00001,'pcolor','off');
[h1,h2]=gf_plot(mf_ls, gf_levelset_get(_ls,'values'), 'contour', 0,'pcolor','off');
//h.color_map = jetcolormap(255);
//[h1,h2]=gf_plot(mf_ls2, get(ls2,'values'), 'contour',0,'pcolor','off');
plot([xc + R*n(2); xc - R*n(2)],[yc - R*n(1), yc + R*n(1)],'b-');
colorbar(min(U),max(U));
drawnow;

printf('demo fictitious domains terminated\n');
