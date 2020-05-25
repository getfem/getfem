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

// Example of automatic refinement of the mesh
// In this example, the refinement will focus on the
// transition between the Dirichlet and the Neumann boundary.

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_refine.sce');

printf('demo refine started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

//clear all; clf;

L = 100;
H = 22;
N = 2;

if (N == 2) then // 2D beam
  m    = gf_mesh('regular simplices',0:10:L, 0:11:H);
  mim  = gf_mesh_im(m);    gf_mesh_im_set(mim,  'integ', gf_integ('IM_TRIANGLE(6)'));
  mfu  = gf_mesh_fem(m,N); gf_mesh_fem_set(mfu, 'fem',   gf_fem('FEM_PK(2,2)'));
  mfd  = gf_mesh_fem(m);   gf_mesh_fem_set(mfd, 'fem',   gf_fem('FEM_PK(2,1)'));
  mf0  = gf_mesh_fem(m);   gf_mesh_fem_set(mf0, 'fem',   gf_fem('FEM_PK(2,0)'));
  mfdu = gf_mesh_fem(m);   gf_mesh_fem_set(mfdu,'fem',   gf_fem('FEM_PK_DISCONTINUOUS(2,2)'));
else         // 3D beam
  m    = gf_mesh('regular simplices',0:10:L, 0:11:H, 0:11:H);
  mim  = gf_mesh_im(m);    gf_mesh_im_set(mim,  'integ', gf_integ('IM_TETRAHEDRON(5)'));
  mfu  = gf_mesh_fem(m,N); gf_mesh_fem_set(mfu, 'fem',   gf_fem('FEM_PK(3,2)'));
  mfd  = gf_mesh_fem(m);   gf_mesh_fem_set(mfd, 'fem',   gf_fem('FEM_PK(3,1)'));
  mf0  = gf_mesh_fem(m);   gf_mesh_fem_set(mf0, 'fem',   gf_fem('FEM_PK(3,0)'));
  mfdu = gf_mesh_fem(m);   gf_mesh_fem_set(mfdu,'fem',   gf_fem('FEM_PK_DISCONTINUOUS(3,1)'));
end

lambda = 121150; 
mu     = 80769;

P = gf_mesh_get(m,'pts');
fleft  = gf_mesh_get(m,'faces from pid',find(abs(P(1,:))<1e-6));
fright = gf_mesh_get(m,'faces from pid',find(abs(P(1,:) - L)<1e-6));

// assign boundary numbers
gf_mesh_set(m,'boundary',1,fleft);
gf_mesh_set(m,'boundary',2,fright);

F = zeros(N,1); F(2) = -20; // the external force

md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add initialized data', 'mu', [mu]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
gf_model_set(md, 'add initialized data', 'VolumicData', F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 1);

h = scf();
h.color_map = jetcolormap(255);

for step=1:8
  dd = gf_mesh_fem_get(mf0, 'basic dof from cvid');
  
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');

  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'lambda', 'mu', mfdu);

  if (N==3) then 
    opt = list('cvlst', get(m,'outer_faces')); 
  else 
    opt = list(); 
  end

  drawlater;
  clf();
  subplot(2,1,1);
  gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U, 'deformation_mf',mfu,'refine', 4, 'deformation_scale',1, opt(:)); 
  colorbar(min(U),max(U)); 
  title('Von Mises stress');
  
  ERR   = gf_compute(mfu,U,'error estimate', mim);
  E     = ERR;
  E(dd) = ERR;

  subplot(2,1,2);
  gf_plot(mf0, E, 'mesh','on', 'refine', 1, opt(:)); 
  colorbar(min(E),max(E));
  title('Error estimate')
  h.color_map = jetcolormap(255);
  drawnow;

  sleep(1000);
  
  Index = find(ERR > 1e-3);
  gf_mesh_set(m, 'refine', Index);
  gf_mesh_set(m, 'optimize structure');
end

printf('demo refine terminated\n');
