// Copyright (C) 2009 Yves Renard.
// Copyright (C) 2009-2010 Yann Collette.
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

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_static_contact.sce');

printf('demo static_contact started\n');

gf_workspace('clear all');

// Import the mesh
//m = gf_mesh('load', path + '/data/disc_P2_h1.mesh');
m = gf_mesh('load', path + '/data/disc_P2_h2.mesh');
//m = gf_mesh('load', path + '/data/disc_P2_h0.5.mesh');

d = gf_mesh_get(m, 'dim');

// parameters of the model
lambda         = 1;   // Lame coefficient
mu             = 1;   // Lame coefficient
friction_coeff = 0.2; // coefficient of friction
r              = 5;   // Augmentation parameter
version        = 4;   // 1 : frictionless contact and the basic contact brick
                      // 2 : contact with 'static' friction and the basic contact brick
                      // 3 : frictionless contact and the contact with a
                      //     rigid obstacle brick
                      // 4 : contact with 'static' friction and the contact with a
                      //     rigid obstacle brick
if (d == 2) then
  obstacle = 'y';  // Signed distance representing the obstacle
else
  obstacle = 'z';  // Signed distance representing the obstacle
end

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
F(d:d:nbdofd*d) = -0.02;

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

if (version == 1 | version == 2) then // defining the matrices BN and BT by hand
  contact_dof = cdof(d:d:nbc*d);
  contact_nodes = gf_mesh_fem_get(mfu, 'basic dof nodes', contact_dof);
  BN = spzeros(nbc, nbdofu);
  for i = 1:nbc
    BN(i, contact_dof(i)) = -1.0;
    gap(i) = contact_nodes(d, i);
  end
  if (version == 2) then
    BT = spzeros(nbc*(d-1), nbdofu);
    for i = 1:nbc
      BT(i*(d-1), contact_dof(i)-d+1) = 1.0;
      if (d > 2) then
        BT(i*(d-1)+1, contact_dof(i)-d+2) = 1.0;
      end
    end
  end

  gf_model_set(md, 'add variable', 'lambda_n', nbc);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (version == 2) then
    gf_model_set(md, 'add variable', 'lambda_t', nbc * (d-1));
    gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
  end
  gf_model_set(md, 'add initialized data', 'gap', gap);
  gf_model_set(md, 'add initialized data', 'alpha', ones(nbc, 1));
  if (version == 1) then
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', 'r', BN, 'gap', 'alpha', 0);
  else
    gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', 'lambda_t', 'r', BN, BT, 'friction_coeff', 'gap', 'alpha', 0);
  end
elseif (version == 3 | version == 4) then // BN and BT defined by the contact brick
  gf_model_set(md, 'add variable', 'lambda_n', nbc);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (version == 3) then
    gf_model_set(md, 'add contact with rigid obstacle brick', mim, 'u', 'lambda_n', 'r', GAMMAC, obstacle, 0);
  else
    gf_model_set(md, 'add variable', 'lambda_t', nbc * (d-1));
    gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
    gf_model_set(md, 'add contact with rigid obstacle brick', mim, 'u', ...
	         'lambda_n', 'lambda_t', 'r', 'friction_coeff', GAMMAC, obstacle, 0);
  end
else
  error('Inexistent version');
end

// Solve the problem

gf_model_get(md, 'solve', 'max_res', 1E-7, 'very noisy', 'max_iter', 100);

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

title('Deformed configuration (not really a small deformation of course ...)');
colorbar(min(VM),max(VM));
drawnow;
sleep(100);

printf('demo static_contact terminated\n');
