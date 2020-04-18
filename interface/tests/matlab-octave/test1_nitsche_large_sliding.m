% Copyright (C) 2012-2020 Yves Renard.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
gf_workspace('clear all');

test_case =1; % 0 = 2D punch on a rigid obstacle
               % 1 = 2D punch on a deformable obstacle (one slave, one master)
               % 2 = 2D with two different meshes
               % 3 = 2D with multi-body and only one mesh
               % 4 = 3D case (sphere / parallelepiped) (two meshes)

clambda1 = 1.; cmu1 = 1.;   % Elasticity parameters
clambda2 = 1.; cmu2 = 1.;   % Elasticity parameters
r = 1;                    % Augmentation parameter
alpha = 1.0;                % Alpha coefficient for "sliding velocity"
f_coeff = 0;              % Friction coefficient

test_tangent_matrix = false;
nonlinear_elasticity = true;
max_iter = 30;
draw_mesh = true;
do_plot = true;
generic_assembly_contact_brick = true;

switch(test_case)
  case {0,1}
    vf = 0.0;
    vf_mult = 1.0;
    penalty_parameter = 0;
    dirichlet_translation = -1;
    max_res = 1E-8;
    release_dist = 1.5;
    self_contact = false;
    if (test_case)
        two_meshes = true;
    else
        two_meshes = false;
    end
  case 3
    vf = 0.01;              % Vertical force
    vf_mult = 1.05;
    penalty_parameter = 0.1;
    release_dist = 0.05;
    max_res = 1E-8;
    self_contact = true;
    two_meshes = false;
  case {2,4}
    vf = 0.01;              % Vertical force
    vf_mult = 1.5;
    penalty_parameter = 0.01;
    max_res = 1E-8;
    if (test_case == 2)
      release_dist = 0.2;
    else
      release_dist = 5;
    end
    self_contact = true;
    two_meshes = true;
end;    

switch (test_case) 
  case 0
    % mesh1 = gf_mesh('load', '../../../tests/meshes/punch2D_1.mesh');
    mesh1 = gf_mesh('load', '../../../tests/meshes/punch2D_2.mesh');
  case 1
    % mesh1 = gf_mesh('load', '../../../tests/meshes/punch2D_1.mesh');
    mesh1 = gf_mesh('load', '../../../tests/meshes/punch2D_2.mesh');
    mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-14,-5];SIZES=[28,5];NSUBDIV=[28,5]');
  case 2
    mesh1 = gf_mesh('load', '../../../tests/meshes/disc_with_a_hole.mesh');
    % mesh1 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0.1];SIZES=[1,0.1];NSUBDIV=[20,2]');
    mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]');
  case 3
    mesh1 = gf_mesh('load', '../../../tests/meshes/multi_body.mesh');
  case 4
    mesh1 = gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
    mesh2 = gf_mesh('import', 'structured', 'GT="GT_PK(3,1)";ORG=[-15,-15,-4];SIZES=[30,30,4];NSUBDIV=[10,10,2]');
end

N = gf_mesh_get(mesh1, 'dim');

mfu1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(mfu1, 'classical fem', 2);
pre_mflambda1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(pre_mflambda1, 'classical fem', 1);
mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', 1);
CONTACT_BOUNDARY1 = 1;
DIRICHLET_BOUNDARY1 = 3;
if (test_case >= 2)
  fb1 = gf_mesh_get(mesh1, 'outer faces');
  gf_mesh_set(mesh1,'region', CONTACT_BOUNDARY1, fb1);
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
mim1_contact = gf_mesh_im(mesh1, 4);
im1_nodes = gf_mesh_im_get(mim1_contact, 'im nodes', gf_mesh_get(mesh1, 'region', CONTACT_BOUNDARY1));
im1_nodes = im1_nodes(1:N,:);

if (two_meshes) 
  mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', 2);
  pre_mflambda2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(pre_mflambda2, 'classical fem', 1);
  mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', 1);
  
  CONTACT_BOUNDARY2 = 2;
  if (test_case ~= 1)
    fb2 = gf_mesh_get(mesh2, 'outer faces');
    gf_mesh_set(mesh2,'region', CONTACT_BOUNDARY2, fb2);
  else
    border = gf_mesh_get(mesh2,'outer faces');
    normals = gf_mesh_get(mesh2, 'normal of faces', border);
    contact_boundary=border(:, find(normals(N, :) > 0.01));
    gf_mesh_set(mesh2, 'region', CONTACT_BOUNDARY2, contact_boundary);
    dirichlet_boundary=border(:, find(normals(N, :) < -0.01));
    DIRICHLET_BOUNDARY2 = 5;
    gf_mesh_set(mesh2, 'region', DIRICHLET_BOUNDARY2, dirichlet_boundary);
  end
  mim2 = gf_mesh_im(mesh2, 4);
  mim2_contact = gf_mesh_im(mesh2, 4);
  im2_nodes = gf_mesh_im_get(mim2_contact, 'im nodes', gf_mesh_get(mesh2, 'region', CONTACT_BOUNDARY2));
  im2_nodes = im2_nodes(1:N,:);
end

if (draw_mesh)
  gf_plot_mesh(mesh1, 'regions', CONTACT_BOUNDARY1);
  if (test_case ~= 3 && test_case ~= 0) 
    hold on
   % gf_plot_mesh(mesh2, 'regions', CONTACT_BOUNDARY2);
    hold off
  end
  pause(2);
end


md=gf_model('real');

F = zeros(1, N); F(N) = -vf;

gf_model_set(md, 'add fem variable', 'u1', mfu1);
% gf_model_set(md, 'add filtered fem variable', 'lambda1', pre_mflambda1, CONTACT_BOUNDARY1);

if (nonlinear_elasticity)
  lawname = 'Ciarlet Geymonat';
  params1 = [clambda1;cmu1;cmu1/2-clambda1/8];
  gf_model_set(md,'add initialized data','params1', params1);
  gf_model_set(md, 'add finite strain elasticity brick', mim1, lawname, 'u1', 'params1');
else
  gf_model_set(md, 'add initialized data', 'clambda1', clambda1);
  gf_model_set(md, 'add initialized data', 'cmu1', cmu1);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1', 'clambda1', 'cmu1');
end
if (test_case == 2)
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

if (two_meshes)
  gf_model_set(md, 'add fem variable', 'u2', mfu2);
  % if (self_contact)
  %   gf_model_set(md, 'add filtered fem variable', 'lambda2', pre_mflambda2, CONTACT_BOUNDARY2);
  % end
  
  if (nonlinear_elasticity)
    lawname = 'Ciarlet Geymonat';
    params2 = [clambda2;cmu2;cmu2/2-clambda2/8];
    gf_model_set(md,'add initialized data','params2', params2);
    gf_model_set(md, 'add finite strain elasticity brick', mim2, lawname, 'u2', 'params2');
  else
    gf_model_set(md, 'add initialized data', 'clambda2', clambda2);
    gf_model_set(md, 'add initialized data', 'cmu2', cmu2);
    gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2', 'clambda2', 'cmu2');
  end
  if (test_case == 2)
    gf_model_set(md, 'add initialized data', 'cpoints2', [0 0]);
    gf_model_set(md, 'add initialized data', 'cunitv2', [1 0]);
    gf_model_set(md, 'add pointwise constraints with multipliers', 'u2', 'cpoints2', 'cunitv2');
  end
  gf_model_set(md, 'add initialized data', 'penalty_param2', [penalty_parameter]);          
  gf_model_set(md, 'add mass brick', mim2, 'u2', 'penalty_param2');
  dat='data2';
  gf_model_set(md, 'add initialized data', dat, F);
  gf_model_set(md, 'add source term brick', mim2, 'u2', 'data2');
  if (test_case == 1)
    Ddata = zeros(1, N);
    gf_model_set(md, 'add initialized data', 'Ddata2', Ddata);
    gf_model_set(md, 'add Dirichlet condition with multipliers', mim2, 'u2', 1, DIRICHLET_BOUNDARY2, 'Ddata2');
  end
end

if (test_case <= 1)
  Ddata = zeros(1, N); Ddata(N) = dirichlet_translation;
  gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim1, 'u1', 1, DIRICHLET_BOUNDARY1, 'Ddata');
end
  


if (generic_assembly_contact_brick)
  direct_generic_assembly = true;
  if (direct_generic_assembly)  % Direct use of high-level generic assembly
    gf_model_set(md, 'add raytracing transformation', 'contact_trans', release_dist);
    if (two_meshes) % The definition of a variable group is not mandatory. Just for test.
      gf_model_set(md, 'define variable group', 'u', 'u1', 'u2');
    else
      gf_model_set(md, 'define variable group', 'u', 'u1');
    end
  
    if (self_contact)
      gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans', mesh1, 'u', CONTACT_BOUNDARY1);
    else
      gf_model_set(md, 'add slave contact boundary to raytracing transformation', 'contact_trans', mesh1, 'u', CONTACT_BOUNDARY1);
    end

    switch(test_case)
      case 0
        gf_model_set(md, 'add rigid obstacle to raytracing transformation', 'contact_trans', '80-sqrt(sqr(x)+sqr(y-80))', N);
      case 1
       gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans', mesh2, 'u', CONTACT_BOUNDARY2);
      case 2
        gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans', mesh2, 'u', CONTACT_BOUNDARY2);
        gf_model_set(md, 'add rigid obstacle to raytracing transformation', 'contact_trans', 'y+1', N);
      case 3
        gf_model_set(md, 'add rigid obstacle to raytracing transformation', 'contact_trans', '2-sqrt(sqr(x)+sqr(y-1))', N);
      case 4
        gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans', mesh2, 'u', CONTACT_BOUNDARY2);
        gf_model_set(md, 'add rigid obstacle to raytracing transformation', 'contact_trans', 'z+5', N);
    end

   gamma='(r*element_size)';
   sigma_chap_u1 =gf_model_get(md, 'Neumann term', 'u1', CONTACT_BOUNDARY1)
   g1_n = '((Interpolate(u2,contact_trans)-u1+Interpolate(X,contact_trans)-X).Transformed_unit_vector(Grad_u1, Normal))';
   g1_t='u1';
   gf_model_set(md, 'add initialized data', 'r', r);
   gf_model_set(md, 'add initialized data', 'f', f_coeff);
   %gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,(Test_u1).Transformed_unit_vector(Grad_u1, Normal)*pos_part(',g1_n,'-',gamma,'*(',sigma_chap_u1,'.Transformed_unit_vector(Grad_u1, Normal))),1)'),CONTACT_BOUNDARY1);
   %gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,-(Interpolate(Test_u2,contact_trans)).Transformed_unit_vector(Grad_u1, Normal)*pos_part(',g1_n,'-',gamma,'*',sigma_chap_u1,'.Transformed_unit_vector(Grad_u1, Normal)),1)'),CONTACT_BOUNDARY1);
   
   % assembly_string1 = strcat('Interpolate_filter(contact_trans,-(Test_u1).(Coulomb_friction_coupled_projection((',sigma_chap_u1,'),Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',f,1/',gamma,')),1)')
   gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,-(Test_u1).(Coulomb_friction_coupled_projection((',sigma_chap_u1,'),Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',f,1/',gamma,')),1)'),CONTACT_BOUNDARY1);
   % assembly_string2 = strcat('Interpolate_filter(contact_trans,(Interpolate(Test_u2,contact_trans)).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',f,','1/',gamma,')),1)')
   gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,(Interpolate(Test_u2,contact_trans)).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',f,','1/',gamma,')),1)'),CONTACT_BOUNDARY1);
  end
end
    % errmax = gf_model_get(md, 'test tangent matrix', 1E-9, 10, 0.0001);
    % errmax = gf_model_get(md, 'test tangent matrix term', 'lambda1', 'u1', 1E-8, 20, 0.00001);
   % disp(sprintf('errmax = %g', errmax));
  
  %if (errmax > 1E-3) error('bad tangent matrix'); end;
    % return;
   % pause(0.1);
 

  gf_model_get(md, 'solve', 'noisy', 'max_iter', max_iter, 'max_res', max_res); % , 'lsearch', 'simplest');

  
  if (do_plot)
    U1 = gf_model_get(md, 'variable', 'u1');
    if (nonlinear_elasticity)
      VM1 = gf_model_get(md, 'compute Von Mises or Tresca', 'u1', lawname, 'params1', mfvm1);
    else
      VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
	    	  'u1', 'clambda1', 'cmu1', mfvm1);
    end
    gf_plot(mfvm1,VM1,'mesh', 'off', 'deformed_mesh','on', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;

    % hold on % quiver plot of the multiplier
    % lambda1 = gf_model_get(md, 'variable', 'lambda1');
    % mf_lambda1 = gf_model_get(md, 'mesh fem of variable', 'lambda1');
    % sl=gf_slice({'boundary'}, mf_lambda1, CONTACT_BOUNDARY1);
    % bound_lambda1=gf_compute(mf_lambda1, lambda1,'interpolate on', sl);
    % bound_u1=gf_compute(mfu1, U1,'interpolate on', sl);
    % pts = gf_slice_get(sl, 'pts');
    % quiver(bound_u1(1,:)+pts(1,:), bound_u1(2,:)+pts(2,:), bound_lambda1(1,:), bound_lambda1(2,:))
    % hold off
  
    % hold on
    % gf_plot(mf_lambda1, lambda1,'mesh', 'off', 'deformed_mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8);
    % hold off
  
    if (two_meshes)
       hold on
       U2 = gf_model_get(md, 'variable', 'u2');
       if (nonlinear_elasticity)
         VM2 = gf_model_get(md, 'compute Von Mises or Tresca', 'u2', lawname, 'params2', mfvm2);
       else
         VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		      'u2', 'clambda2', 'cmu2', mfvm2);
       end
       gf_plot(mfvm2,VM2,'mesh', 'off', 'deformed_mesh','on', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;
       hold off
    end;
   end;
