% Copyright (C) 2012-2012 Yves Renard.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 2.1 of the License,  or
% (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

gf_workspace('clear all');
clear all;

lambda = 1.; mu = 1.;   % Elasticity parameters
r = 1.0;                % Augmentation parameter
f_coeff = 1.;           % Friction coefficient
vf = 0.01;                % Vertical force
penalty_parameter = 0.1;


mesh1 = gf_mesh('load', '../../../tests/meshes/disc_with_a_hole.mesh');
mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]');

N = gf_mesh_get(mesh1, 'dim');

mfu1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(mfu1, 'classical fem', 2);
pre_mflambda1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(pre_mflambda1, 'classical fem', 1);
mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', 1);
fb1 = gf_mesh_get(mesh1, 'outer faces');
CONTACT_BOUNDARY1 = 1;
gf_mesh_set(mesh1,'boundary', CONTACT_BOUNDARY1, fb1);
dol1 = gf_mesh_fem_get(pre_mflambda1, 'basic dof on region', CONTACT_BOUNDARY1);
mflambda1 = gf_mesh_fem('partial',  pre_mflambda1, dol1);
mim1 = gf_mesh_im(mesh1, 4);


mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', 1);
pre_mflambda2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(pre_mflambda2, 'classical fem', 1);
mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', 1);
fb2 = gf_mesh_get(mesh2, 'outer faces');
CONTACT_BOUNDARY2 = 2;
gf_mesh_set(mesh2,'boundary', CONTACT_BOUNDARY2, fb2);
dol2 = gf_mesh_fem_get(pre_mflambda2, 'basic dof on region', CONTACT_BOUNDARY2);
mflambda2 = gf_mesh_fem('partial',  pre_mflambda2, dol2);
mim2 = gf_mesh_im(mesh2, 8);

two_bodies = 1;

md=gf_model('real');
gf_model_set(md, 'add initialized data', 'lambda', lambda);
gf_model_set(md, 'add initialized data', 'mu', mu);

if (two_bodies) 
  gf_model_set(md, 'add fem variable', 'u1', mfu1);
  gf_model_set(md, 'add fem variable', 'lambda1', mflambda1);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1', 'lambda', 'mu');
%   gf_model_set(md, 'add initialized data', 'cpoints1', [0 0.5 0 1.5 0 0.5 0 1.5]);
%   gf_model_set(md, 'add initialized data', 'cunitv1', [1 0 1 0 0 1 0 1]);
%   gf_model_set(md, 'add initialized data', 'cdata', [0 0 -0.01 -0.01]);
%   gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints1', 'cunitv1', 'cdata');
  gf_model_set(md, 'add initialized data', 'cpoints1', [0 0.5 0 1.5]);
  gf_model_set(md, 'add initialized data', 'cunitv1', [1 0 1 0]);
  gf_model_set(md, 'add initialized data', 'cdata', [0 0]);
  gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints1', 'cunitv1', 'cdata');
  gf_model_set(md, 'add initialized data', 'data1', [0 -vf]);
  gf_model_set(md, 'add source term brick', mim1, 'u1', 'data1');
  gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);          
  gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');
end;

gf_model_set(md, 'add fem variable', 'u2', mfu2);
gf_model_set(md, 'add fem variable', 'lambda2', mflambda2);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2', 'lambda', 'mu');
gf_model_set(md, 'add initialized data', 'cpoints2', [0 0]);
gf_model_set(md, 'add initialized data', 'cunitv2', [1 0]);
gf_model_set(md, 'add pointwise constraints with multipliers', 'u2', 'cpoints2', 'cunitv2');
gf_model_set(md, 'add initialized data', 'data2', [0 -vf]);
gf_model_set(md, 'add source term brick', mim2, 'u2', 'data2');
gf_model_set(md, 'add initialized data', 'penalty_param2', [penalty_parameter]);          
gf_model_set(md, 'add mass brick', mim2, 'u2', 'penalty_param2');

gf_model_set(md, 'add initialized data', 'r', r);
gf_model_set(md, 'add initialized data', 'f', f_coeff);

indb = gf_model_set(md, 'add integral large sliding contact brick', mim2, 'u2', 'lambda2', 'r', 'f', CONTACT_BOUNDARY2);

if (two_bodies) 
  gf_model_set(md, 'add boundary to large sliding contact brick', indb, mim1, 'u1', 'lambda1', CONTACT_BOUNDARY1);
end;

gf_model_set(md, 'add rigid obstacle to large sliding contact brick', indb, 'y');

% gf_model_get(md, 'test tangent matrix', 1E-6, 10, 0.00001);



for i=1:100

   
    
gf_model_get(md, 'solve', 'noisy', 'max_iter', 50, 'max_res', 1e-8); % , 'lsearch', 'simplest');

U2 = gf_model_get(md, 'variable', 'u2');
VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u2', 'lambda', 'mu', mfvm2);

gf_plot(mfvm2,VM2,'mesh','off', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;

if (two_bodies)
   hold on
   U1 = gf_model_get(md, 'variable', 'u1');
   VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u1', 'lambda', 'mu', mfvm1);
   gf_plot(mfvm1,VM1,'mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;
   hold off
end;

axis([-2, 2, -0.2, 3]);
pause(1);

 vf = vf + 0.001;
 gf_model_set(md, 'variable', 'data1', [0 -vf]);
 gf_model_set(md, 'variable', 'data2', [0 -vf]);

end;






