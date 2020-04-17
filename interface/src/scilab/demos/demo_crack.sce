// Copyright (C) 2009 Luis Saavedra, Yves Renard.
// Copyright (C) 2009-2020, Yann Collette
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
// Linear Elastostatic problem with a crack.
// A good example of use of GetFEM.

lines(0);

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

path = get_absolute_file_path('demo_crack.sce');

printf('demo crack started\n');

// Parameters:
nx         = 20;
DIRICHLET  = 101;
Lambda     = 1.25e10;   // Lame coefficient
Mu         = 1.875e10;  // Lame coefficient

// Global Functions
ck0    = gf_global_function('crack',0);
ck1    = gf_global_function('crack',1);
ck2    = gf_global_function('crack',2);
ck3    = gf_global_function('crack',3);
ckoff  = gf_global_function('cutoff',2,0.4,0.01,0.4);
ckoff0 = gf_global_function('product', ck0, ckoff);
ckoff1 = gf_global_function('product', ck1, ckoff);
ckoff2 = gf_global_function('product', ck2, ckoff);
ckoff3 = gf_global_function('product', ck3, ckoff);

// Mesh in action:
m = gf_mesh('regular_simplices', -0.5:1.0/nx:0.5+1.0/nx, -0.5:1.0/nx:0.5+1.0/nx);
// m = gf_mesh('import','gmsh',path + 'data/quad.msh')

// boundary set:
gf_mesh_set(m,'region',DIRICHLET, gf_mesh_get(m,'outer_faces'));

// MeshFem in action:
mf_pre_u = gf_mesh_fem(m);
gf_mesh_fem_set(mf_pre_u,'fem',gf_fem('FEM_PK(2,1)'));

// Levelset in action:
_ls = gf_levelset(m,1,'y','x');
mls = gf_mesh_levelset(m);
gf_mesh_levelset_set(mls,'add',_ls);
gf_mesh_levelset_set(mls,'adapt');
mfls_u    = gf_mesh_fem('levelset',mls,mf_pre_u);
mf_sing_u = gf_mesh_fem('global function', m, _ls, list(ckoff0,ckoff1,ckoff2,ckoff3),1);
mf_u      = gf_mesh_fem('sum',mf_sing_u,mfls_u);
gf_mesh_fem_set(mf_u,'qdim',2);

// exact solution:
mf_ue = gf_mesh_fem('global function', m, _ls, list(ck0,ck1,ck2,ck3));
A = 2+2*Mu/(Lambda+2*Mu);
B =-2*(Lambda+Mu)/(Lambda+2*Mu);
Ue = zeros(2,4);
Ue(1,1) =   0; Ue(2,1) = A-B; // sin(theta/2)
Ue(1,2) = A+B; Ue(2,2) = 0;   // cos(theta/2)
Ue(1,3) =  -B; Ue(2,3) = 0;   // sin(theta/2)*sin(theta)
Ue(1,4) =   0; Ue(2,4) = B;   // cos(theta/2)*cos(theta)
Ue = Ue / 2*%pi;
Ue = matrix(Ue,1,8);

// MeshIm in action:
mim = gf_mesh_im('levelset', mls, 'all', ...
                 gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'), ...
                 gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'), ...
                 gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'));

// Model in action:
md = gf_model('real');
gf_model_set(md,'add_fem_variable', 'u', mf_u);

// data
gf_model_set(md,'add_initialized_data','lambda', [Lambda]);
gf_model_set(md,'add_initialized_data','mu', [Mu]);
gf_model_set(md,'add_isotropic_linearized_elasticity_brick',mim,'u','lambda','mu');
// gf_model_set(md,'add_variable','mult_spec',6);
// BB = gf_spmat('empty',6,gf_mesh_fem_get(mf_u,'nbdof'));
// gf_model_set(md,'add_constraint_with_multipliers','u','mult_spec',BB,zeros(6,1));
gf_model_set(md,'add_initialized_fem_data','DirichletData', mf_ue, Ue);
gf_model_set(md,'add_Dirichlet_condition_with_penalization',mim,'u', 1e12, DIRICHLET, 'DirichletData');

// assembly of the linear system and solve:
gf_model_get(md,'solve');
U = gf_model_get(md,'variable','u');

// Interpolation of the solution on a cut mesh for the drawing purpose
cut_mesh = gf_mesh_levelset_get(mls,'cut_mesh');
mfv = gf_mesh_fem(cut_mesh,2);
gf_mesh_fem_set(mfv,'classical_discontinuous_fem',2,0.001);
gf_mesh_fem_set(mf_ue,'qdim',2);
V  = gf_compute(mf_u,U,'interpolate_on',mfv);
Ve = gf_compute(mf_ue,Ue,'interpolate_on',mfv);

// computation of the Von Mises stress
mfvm = gf_mesh_fem(cut_mesh);
gf_mesh_fem_set(mfvm,'classical_discontinuous_fem',2,0.001);
gf_model_set(md,'add initialized fem data', 'u_cut', mfv, V);
VM = gf_model_get(md,'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u_cut', 'lambda', 'mu', mfvm);

// export to pos
gf_mesh_fem_get(mfv,'export_to_pos', path + '/crack.pos',V,'V',Ve,'Ve', mfvm, VM,'Von Mises');
printf('You can view the solution with (for example): gmsh %scrack.pos\n', path);

// drawing the solution
h = scf();
h.color_map = jetcolormap(255);
drawlater;
subplot(2,1,1);
title('the mesh before cracking');
gf_plot_mesh(mf_pre_u);
subplot(2,1,2);
title('the mesh after cracking');
gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation_mf', mfv, 'deformation', V, 'deformation_scale', 0.20);
colorbar(min(VM),max(VM));
drawnow;

printf('demo crack terminated\n');
