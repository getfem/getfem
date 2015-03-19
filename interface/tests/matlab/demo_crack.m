% Copyright (C) 2009-2012 Luis Saavedra, Yves Renard.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
%
%%%  Linear Elastostatic problem with a crack. Use of XFem.
%
%  This program is used to check that matlab-getfem is working. This is
%  also a good example of use of GetFEM++.
%
%%%

gf_workspace('clear all');

% Parameters
nx         = 20;
DIRICHLET  = 101;
Lambda     = 1.25e10;   % Lame coefficient
Mu         = 1.875e10;  % Lame coefficient

% Global Functions
ck0 = gf_global_function('crack',0);
ck1 = gf_global_function('crack',1);
ck2 = gf_global_function('crack',2);
ck3 = gf_global_function('crack',3);
coff = gf_global_function('cutoff',2,0.4,0.01,0.4);
ckoff0 = gf_global_function('product', ck0, coff);
ckoff1 = gf_global_function('product', ck1, coff);
ckoff2 = gf_global_function('product', ck2, coff);
ckoff3 = gf_global_function('product', ck3, coff);

% Mesh in action:
m = gf_mesh('regular_simplices', -0.5:1.0/nx:0.5+1.0/nx, -0.5:1.0/nx:0.5+1.0/nx);
% m = gf_mesh('import','gmsh','quad.msh')

% boundary set:
gf_mesh_set(m,'region',DIRICHLET, gf_mesh_get(m,'outer_faces'));
% MeshFem in action:
mf_pre_u = gf_mesh_fem(m);
gf_mesh_fem_set(mf_pre_u,'fem',gf_fem('FEM_PK(2,1)'));
% Levelset in action:
ls  = gf_levelset(m,1,'y','x');
mls = gf_mesh_levelset(m);
gf_mesh_levelset_set(mls,'add',ls);
gf_mesh_levelset_set(mls,'adapt');
mfls_u    = gf_mesh_fem('levelset',mls,mf_pre_u);
mf_sing_u = gf_mesh_fem('global function',m,ls, {ckoff0,ckoff1,ckoff2,ckoff3},1);
mf_u      = gf_mesh_fem('sum',mf_sing_u,mfls_u);
gf_mesh_fem_set(mf_u,'qdim',2);
% exact solution:
mf_ue = gf_mesh_fem('global function',m,ls,{ck0,ck1,ck2,ck3});
A = 2+2*Mu/(Lambda+2*Mu);
B=-2*(Lambda+Mu)/(Lambda+2*Mu)
Ue = zeros(2,4);
Ue(1,1) =   0; Ue(2,1) = A-B; % sin(theta/2)
Ue(1,2) = A+B; Ue(2,2) = 0;   % cos(theta/2)
Ue(1,3) =  -B; Ue(2,3) = 0;   % sin(theta/2)*sin(theta)
Ue(1,4) =   0; Ue(2,4) = B;   % cos(theta/2)*cos(theta)
Ue = Ue / 2*pi;
Ue=reshape(Ue,1,8);
% MeshIm in action:
mim = gf_mesh_im('levelset', mls, 'all', ...
		gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'), ...
		gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'), ...
		gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'));
% Model in action:
md = gf_model('real');
gf_model_set(md,'add_fem_variable', 'u', mf_u);
% data
gf_model_set(md,'add_initialized_data','lambda', [Lambda]);
gf_model_set(md,'add_initialized_data','mu', [Mu]);
gf_model_set(md,'add_isotropic_linearized_elasticity_brick',mim,'u','lambda','mu');
gf_model_set(md,'add_initialized_fem_data','DirichletData', mf_ue, Ue);
gf_model_set(md,'add_Dirichlet_condition_with_penalization',mim,'u', 1e12, DIRICHLET, 'DirichletData');

% assembly of the linear system and solve:
gf_model_get(md,'solve');
U = gf_model_get(md,'variable','u');

% Interpolation of the solution on a cut mesh for the drawing purpose
cut_mesh = gf_mesh_levelset_get(mls,'cut_mesh');
mfv = gf_mesh_fem(cut_mesh,2);
gf_mesh_fem_set(mfv,'classical_discontinuous_fem',2,0.001);
gf_mesh_fem_set(mf_ue,'qdim',2);
V  = gf_compute(mf_u,U,'interpolate_on',mfv);
Ve = gf_compute(mf_ue,Ue,'interpolate_on',mfv);

% computation of the Von Mises stress
mfvm = gf_mesh_fem(cut_mesh);
gf_mesh_fem_set(mfvm,'classical_discontinuous_fem',2,0.001);
gf_model_set(md,'add initialized fem data', 'u_cut', mfv, V);
VM = gf_model_get(md,'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u_cut', 'lambda', 'mu', mfvm);

% export to pos
gf_mesh_fem_get(mfv,'export_to_pos','crack.pos',V,'V',Ve,'Ve', mfvm, VM,'Von Mises');
disp('You can view the solution with (for example): gmsh crack.pos\n');

% drawing the solution
gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation_mf', mfv, 'deformation', V, 'deformation_scale', 0.20);
colorbar;
