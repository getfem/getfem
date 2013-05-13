% Copyright (C) 2012-2013 Yves Renard.
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


clear all;
% gf_workspace('clear all');

test_case = 0; % 0 = 2D punch on a rigid obstacle
               % 1 = 2D with two differente meshes
               % 2 = 2D with multi-body and only one mesh
               % 3 = 3D case (sphere / parallelepiped) (two meshes)

lambda = 1.; mu = 1.;   % Elasticity parameters
r = 100;                 % Augmentation parameter
f_coeff = 0;            % Friction coefficient
test_tangent_matrix = false;
nonlinear_elasticity = true;
max_iter = 50;

if (test_case == 2)
  vf = 0.01;              % Vertical force
  vf_mult = 1.05;
  penalty_parameter = 0.1;
  release_dist = 0.05;
  max_res = 1E-9;
elseif (test_case == 0)
  vf = 0.001;
  vf_mult = 1.1;
  penalty_parameter = 0;
  dirichlet_translation = 0.5;
  max_res = 1E-8;
  release_dist = 5;
else
  vf = 0.01;              % Vertical force
  vf_mult = 1.5;
  penalty_parameter = 0.01;
  max_res = 1E-8;
  if (test_case == 1)
    release_dist = 0.1;
  else
    release_dist = 5;
  end
end;    

if (test_case == 0)
  % mesh1=gf_mesh('load', '../../../tests/meshes/punch2D_1.mesh');
  mesh1=gf_mesh('load', '../../../tests/meshes/punch2D_2.mesh');
elseif (test_case == 1)
  mesh1 = gf_mesh('load', '../../../tests/meshes/disc_with_a_hole.mesh');
  % mesh1 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0.1];SIZES=[1,0.1];NSUBDIV=[20,2]');
  mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]');
elseif (test_case == 2)
  mesh1 = gf_mesh('load', '../../../tests/meshes/multi_body.mesh');
else
  mesh1 = gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
  mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(3,1)";ORG=[-15,-15,-4];SIZES=[30,30,4];NSUBDIV=[10,10,2]');
end

N = gf_mesh_get(mesh1, 'dim');

mfu1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(mfu1, 'classical fem', 2);
pre_mflambda1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(pre_mflambda1, 'classical fem', 1);
mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', 1);
CONTACT_BOUNDARY1 = 1;
DIRICHLET_BOUNDARY1 = 3;
if (test_case ~= 0)
  fb1 = gf_mesh_get(mesh1, 'outer faces');
  gf_mesh_set(mesh1,'boundary', CONTACT_BOUNDARY1, fb1);
else
  border = gf_mesh_get(mesh1,'outer faces');
  normals = gf_mesh_get(mesh1, 'normal of faces', border);
  contact_boundary=border(:, find(normals(N, :) < -0.01));
  gf_mesh_set(mesh1, 'region', CONTACT_BOUNDARY1, contact_boundary);
  P=gf_mesh_get(mesh1,'pts'); % get list of mesh points coordinates
  pidtop=find(P(N,:) > 39.999); % find those on top of the object
  ftop=gf_mesh_get(mesh1,'faces from pid',pidtop); 
  gf_mesh_set(mesh1, 'region', DIRICHLET_BOUNDARY1, ftop);
end




% dol1 = gf_mesh_fem_get(pre_mflambda1, 'basic dof on region', CONTACT_BOUNDARY1);
% mflambda1 = gf_mesh_fem('partial',  pre_mflambda1, dol1);
mim1 = gf_mesh_im(mesh1, 4);

if (test_case ~= 2 && test_case ~= 0) 
  mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', 2);
  pre_mflambda2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(pre_mflambda2, 'classical fem', 1);
  mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', 1);
  fb2 = gf_mesh_get(mesh2, 'outer faces');
  CONTACT_BOUNDARY2 = 2;
  gf_mesh_set(mesh2,'boundary', CONTACT_BOUNDARY2, fb2);
  % dol2 = gf_mesh_fem_get(pre_mflambda2, 'basic dof on region', CONTACT_BOUNDARY2);
  % mflambda2 = gf_mesh_fem('partial',  pre_mflambda2, dol2);
  mim2 = gf_mesh_im(mesh2, 4);
end

md=gf_model('real');
gf_model_set(md, 'add initialized data', 'lambda', lambda);
gf_model_set(md, 'add initialized data', 'mu', mu);

F = zeros(1, N); F(N) = -vf;

gf_model_set(md, 'add fem variable', 'u1', mfu1);
gf_model_set(md, 'add filtered fem variable', 'lambda1', pre_mflambda1, CONTACT_BOUNDARY1);

if (nonlinear_elasticity)
  lawname = 'Ciarlet Geymonat';
  params = [lambda;mu;mu/2-lambda/8];
  gf_model_set(md,'add initialized data','params', params);
  gf_model_set(md, 'add nonlinear elasticity brick', mim1, 'u1', lawname, 'params');
else
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1', 'lambda', 'mu');
end
if (test_case == 1)
  %   gf_model_set(md, 'add initialized data', 'cpoints1', [0 0.5 0 1.5 0 0.5 0 1.5]);
  %   gf_model_set(md, 'add initialized data', 'cunitv1', [1 0 1 0 0 1 0 1]);
  %   gf_model_set(md, 'add initialized data', 'cdata', [0 0 -0.01 -0.01]);
  %   gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints1', 'cunitv1', 'cdata');
  gf_model_set(md, 'add initialized data', 'cpoints1', [0 0.5 0 1.5]);
  gf_model_set(md, 'add initialized data', 'cunitv1', [1 0 1 0]);
  gf_model_set(md, 'add initialized data', 'cdata', [0 0]);
  gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints1', 'cunitv1', 'cdata');
end
gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);
indmass = gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');
gf_model_set(md, 'add initialized data', 'data1', F);
gf_model_set(md, 'add source term brick', mim1, 'u1', 'data1');

if (test_case ~= 2 && test_case ~= 0)
  gf_model_set(md, 'add fem variable', 'u2', mfu2);
  gf_model_set(md, 'add filtered fem variable', 'lambda2', pre_mflambda2, CONTACT_BOUNDARY2);
  if (nonlinear_elasticity)
    gf_model_set(md, 'add nonlinear elasticity brick', mim2, 'u2', lawname, 'params');
  else
    gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2', 'lambda', 'mu');
  end
  if (test_case == 1)
    gf_model_set(md, 'add initialized data', 'cpoints2', [0 0]);
    gf_model_set(md, 'add initialized data', 'cunitv2', [1 0]);
    gf_model_set(md, 'add pointwise constraints with multipliers', 'u2', 'cpoints2', 'cunitv2');
  end
  gf_model_set(md, 'add initialized data', 'data2', F);
  gf_model_set(md, 'add source term brick', mim2, 'u2', 'data2');
  gf_model_set(md, 'add initialized data', 'penalty_param2', [penalty_parameter]);          
  gf_model_set(md, 'add mass brick', mim2, 'u2', 'penalty_param2');
end

if (test_case == 0)
  Ddata = zeros(1, N); Ddata(N) = dirichlet_translation;
  gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim1, 'u1', 1, DIRICHLET_BOUNDARY1, 'Ddata');
end
  
gf_model_set(md, 'add initialized data', 'r', r);
gf_model_set(md, 'add initialized data', 'f', f_coeff);



% delaunay: after dist 
mcff=gf_multi_contact_frame(md, N, release_dist, false, true, 0.2, true, 0, false);
gf_multi_contact_frame_set(mcff, 'add master boundary', mim1, CONTACT_BOUNDARY1, 'u1', 'lambda1');
if (test_case == 1) 
  gf_multi_contact_frame_set(mcff, 'add master boundary', mim2, CONTACT_BOUNDARY2, 'u2', 'lambda2');
  gf_multi_contact_frame_set(mcff, 'add obstacle', 'y');
elseif (test_case == 0) 
  gf_multi_contact_frame_set(mcff, 'add obstacle', '80-sqrt(x^2+(y-80)^2)');
elseif (test_case == 2) 
  gf_multi_contact_frame_set(mcff, 'add obstacle', '2-sqrt(x^2+(y-1)^2)');      
else
  gf_multi_contact_frame_set(mcff, 'add master boundary', mim2, CONTACT_BOUNDARY2, 'u2', 'lambda2');
  gf_multi_contact_frame_set(mcff, 'add obstacle', 'z+5');
end

gf_model_set(md, 'add integral large sliding contact brick', mcff, 'r', 'f');


for nit=1:100

  if (test_tangent_matrix) 
    errmax = gf_model_get(md, 'test tangent matrix', 1E-6, 20, 0.0001);
    disp(sprintf('errmax = %g', errmax));
    if (errmax > 1E-3) error('bad tangent matrix'); end;
    pause;
  end
    
  gf_model_get(md, 'solve', 'noisy', 'max_iter', max_iter, 'max_res', max_res); % , 'lsearch', 'simplest');

  U1 = gf_model_get(md, 'variable', 'u1');
  VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u1', 'lambda', 'mu', mfvm1);

  gf_plot(mfvm1,VM1,'mesh', 'off', 'deformed_mesh','on', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;

  hold on % quiver plot of the multiplier
  lambda1 = gf_model_get(md, 'variable', 'lambda1');
  mf_lambda1 = gf_model_get(md, 'mesh fem of variable', 'lambda1');
  sl=gf_slice({'boundary'}, mf_lambda1, CONTACT_BOUNDARY1);
  bound_lambda1=gf_compute(mf_lambda1, lambda1,'interpolate on', sl);
  bound_u1=gf_compute(mfu1, U1,'interpolate on', sl);
  pts = gf_slice_get(sl, 'pts');
  quiver(bound_u1(1,:)+pts(1,:), bound_u1(2,:)+pts(2,:), bound_lambda1(1,:), bound_lambda1(2,:))
  hold off
  
  % hold on
  % gf_plot(mf_lambda1, lambda1,'mesh', 'off', 'deformed_mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8);
  % hold off
  
  if (test_case ~= 2 && test_case ~= 0)
     hold on
     U2 = gf_model_get(md, 'variable', 'u2');
     VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		  'u2', 'lambda', 'mu', mfvm2);
     gf_plot(mfvm2,VM2,'mesh', 'off', 'deformed_mesh','on', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;
     hold off
  end;

  hold on
  % tic;
  % gf_multi_contact_frame_get(mcff, 'compute pairs');
  % toc
  slpt = gf_multi_contact_frame_get(mcff, 'slave points');
  mapt = gf_multi_contact_frame_get(mcff, 'master points');
  if (N == 2)
    line([slpt(1,:); mapt(1,:)], [slpt(2,:); mapt(2,:)], 'Color', 'blue');
    scatter(slpt(1,:), slpt(2, :), 20, 'red');
    scatter(mapt(1,:), mapt(2, :), 20, 'cyan');
  elseif (N == 3)
    line([slpt(1,:); mapt(1,:)], [slpt(2,:); mapt(2,:)],  [slpt(3,:); mapt(3,:)], 'Color', 'blue');
    scatter3(slpt(1,:), slpt(2, :), slpt(3, :), 20, 'red');
    scatter3(mapt(1,:), mapt(2, :), mapt(3, :), 20, 'cyan');
  end
  if (test_case == 0)
   rectangle('position', [-80, 0, 160, 160], 'Curvature', [1 1]);  % draw the obstacle
   axis([-15 15 -3 44]);
  end
  if (test_case == 2)
   rectangle('position', [-2, -1, 4, 4], 'Curvature', [1 1]);  % draw the obstacle
   axis([-1.3 1.3 -1.1 0.8]);
  end
  hold off

  pause;

  vf = vf * vf_mult; F(N) = -vf;
  gf_model_set(md, 'variable', 'data1', F);
  if (test_case ~= 2 && test_case ~= 0)
    gf_model_set(md, 'variable', 'data2', F);
  end
  
  if (test_case == 0)
    Ddata(N) = Ddata(N) - 0.25;
    gf_model_set(md, 'variable', 'Ddata', Ddata);
  end
  

end;
