// Matlab GetFEM++ interface
//
// Copyright (C) 2009 Yves Renard.
//
// This file is a part of GetFEM++
//
// GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
// Static equilibrium of an elastic solid in contact with a rigid foundation
//
// This program is used to check that matlab-getfem is working. This is also
// a good example of use of GetFEM++.
//

gf_workspace('clear all');

// parameters of the model
lambda = 1;  // Lame coefficient
mu = 1;      // Lame coefficient

// set a custom colormap
h = scf();
h.color_map = jetcolormap(255);

//r = [0.7 .7 .7];
//l = r($,:); s=63; s1=20; s2=25; s3=48; s4=55; 
//for i=1:s
//  c1 = max(min((i-s1)/(s2-s1),1),0);
//  c2 = max(min((i-s3)/(s4-s3),1),0); 
//  r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; 
//end
//h.color_map = r;

// Import the mesh
m = gf_mesh('load', 'data/disc_P2_h1.mesh');

d = gf_mesh_get(m, 'dim');

// Selection of the contact boundary
border  = gf_mesh_get(m, 'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary = border(:,find(normals(d, :) < 0));
GAMMAC = 1;
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);

// Plot the mesh
drawlater;
gf_plot_mesh(m, 'regions', [GAMMAC]);
title('Mesh and contact boundary (in blue)');
drawnow;
sleep(100);

// Finite element methods
mfu = gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', 2);
mfd = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', 1);
mfvm = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', 1);
mim = gf_mesh_im(m, 4);

// Volumic density of force
nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');
F = zeros(nbdofd*d, 1);
F(d:d:nbdofd*d) = -0.01;

// Elasticity model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'mu', [mu]);
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

// Small penalty term to avoid rigid motion (should be replaced by an explicit
// treatment of the rigid motion with a constraint matrix)
gf_model_set(md, 'add initialized data', 'penalty_param', [1E-8]);          
gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');

// The contact condition

cdof = gf_mesh_fem_get(mfu, 'dof on region', GAMMAC);
nbc = size(cdof, 2) / d;
contact_dof = cdof(d:d:nbc*d);
contact_nodes = gf_mesh_fem_get(mfu, 'basic dof nodes', contact_dof);
BN = spzeros(nbc, nbdofu);
for i = 1:nbc
  BN(i, contact_dof(i)) = -1.0;
  gap(i) = contact_nodes(d, i);
end

printf('nbc = %d\n', nbc);

gf_model_set(md, 'add variable', 'lambda_n', nbc);
gf_model_set(md, 'add initialized data', 'r', [1.0]);
gf_model_set(md, 'add initialized data', 'gap', gap);
gf_model_set(md, 'add initialized data', 'alpha', ones(nbc, 1));
gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', 'r', BN, 'gap', 'alpha', 0);
// gf_model_set(md, 'add contact with rigid obstacle brick', mim, 'u', 'lambda_n', 'r', GAMMAC, 'y', 0);

// Solve the problem

gf_model_get(md, 'solve', 'max_res', 1E-7, 'very noisy', 'max_iter', 40);

U = gf_model_get(md, 'variable', 'u');
lambda_n = gf_model_get(md, 'variable', 'lambda_n');
VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u', 'lambda', 'mu', mfvm);

h = scf();
drawlater;
h.color_map = jetcolormap(255);
if (d == 3) then
  gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
          gf_mesh_get(mfdu,'outer faces'), 'deformation', U, ...
          'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
else
  gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
          'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
end

title('Deformed configuration');
colorbar(min(VM),max(VM));
drawnow;
sleep(100);
