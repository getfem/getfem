

asize =  size(who('automatic_var654777'));
if (asize(1)) is_automatic = true; else is_automatic = false; end

gf_workspace('clear all');

approximation_type = 2; % 0 : Augmentend Lagrangian
                        % 1 : Nitsche (biased)
                        % 2 : Unbiased Nitsche method

draw_mesh = true;
ref_sol = 0
max_res = 1E-8;
 max_iter = 50;   %N = 2;% 2 ou 3 dimensions
R=0.25;
dirichlet_val = 0;
  f_coeff=0;  
  clambda = 1;           % Lame coefficient
cmu = 1;               % Lame coefficient
vertical_force = -0.3;
    penalty_parameter = 0.00001;
if (ref_sol == 0)
    method = [-1]; % theta
    gamma00 = [0.1]; %1/200
    nxy =[30]; % 2D ->400 and 3D -> 30
    ls_degree = 2;

    
end
if (ref_sol == 1)
    method = [0];
    gamma00 = [1/200];
    if (N==2)
        nxy=[2 5 10 15 27 37 50 60 75 80];
    else
        nxy = [5 10 12 15 20 25];
    end
    ls_degree = 2;


end
if (ref_sol == 2)
    method = [0 1 -1];
    gamma00 = [400 200 100 50 25 10 1 1/10 1/25 1/50 1/100 1/200 1/400];
    nxy = [21];
    ls_degree = 2;


end
for xx = 1:1:size(method,2)
for yy = 1:1:size(gamma00,2)
for zz = 1:1:size(nxy,2)
theta = method(xx);
gamma0 = gamma00(yy);
NX = 4*floor(nxy(zz)/4);

mo =  gf_mesher_object('ball',[0 0],0.25)
mesh1 = gf_mesh('generate', mo, 1/NX ,2) ;
%mesh1 = gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
mesh2 = gf_mesh('regular simplices', -.5:(1/NX):.5, -0.5:(1/NX):-0.25);


N = gf_mesh_get(mesh1, 'dim');

mfu1 = gf_mesh_fem(mesh1, N); gf_mesh_fem_set(mfu1, 'classical fem', 2);
mflambda1 = gf_mesh_fem(mesh1, 1); gf_mesh_fem_set(mflambda1, 'classical fem', 1);

mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', 1);

mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', 2);

mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', 1);

CONTACT_BOUNDARY1 = 1;

fb1 = gf_mesh_get(mesh1, 'outer faces');
normals1 = gf_mesh_get(mesh1, 'normal of faces', fb1);
contact_boundary1=fb1(:,find(normals1(N,:) <= 2/NX)); 
gf_mesh_set(mesh1,'region', CONTACT_BOUNDARY1,contact_boundary1);

mim1 = gf_mesh_im(mesh1, 4);
mim1_contact = gf_mesh_im(mesh1, 6);
  
CONTACT_BOUNDARY2 = 2;
    
border = gf_mesh_get(mesh2,'outer faces');
    normals2 = gf_mesh_get(mesh2, 'normal of faces', border);
       
   P = gf_mesh_get(mesh2, 'pts');
   P=P(:,find(P(1,:)>=-0.25));
   P=P(:,find(P(1,:)<=0.25)) ;
   P=P(:,find(P(2,:)>=-0.25)) ;
   Id=gf_mesh_get(mesh2, 'pid from coords', P);
   contact_boundary2=gf_mesh_get(mesh2,'faces from pid',Id);

   gf_mesh_set(mesh2, 'region', CONTACT_BOUNDARY2,  contact_boundary2);
    
    
dirichlet_boundary=border(:, find(normals2(N, :) < -0.01));
DIRICHLET_BOUNDARY2 = 3;
gf_mesh_set(mesh2, 'region', DIRICHLET_BOUNDARY2, dirichlet_boundary);

mim2 = gf_mesh_im(mesh2, 4);
mim2_contact = gf_mesh_im(mesh2, 6);

    if (draw_mesh)
    gf_plot_mesh(mesh2);
    hold on
    gf_plot_mesh(mesh1,  'regions', CONTACT_BOUNDARY1);
    hold off
    % pause;
    end
    
    %Elastic model 
md=gf_model('real');
gf_model_set(md,'add fem variable', 'u1', mfu1);
gf_model_set(md,'add fem variable', 'u2', mfu2);

gf_model_set(md,'add initialized data', 'gamma0', gamma0);
gf_model_set(md, 'add initialized data', 'theta', theta);
gf_model_set(md, 'add initialized data', 'friction_coeff', [f_coeff]);

gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1','clambda', 'cmu');
gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2','clambda', 'cmu');

gf_model_set(md, 'add initialized data', 'Fdata', [0 vertical_force]);
gf_model_set(md, 'add source term brick', mim1, 'u1', 'Fdata');

Ddata = zeros(1, N);
gf_model_set(md, 'add initialized data', 'Ddata2', Ddata);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim2, 'u2', 1, DIRICHLET_BOUNDARY2, 'Ddata2');
   
if (N==2)
  cpoints = [0,0    0,0.1]; % constrained points for 2d
  cunitv  = [1,0   1,0];   % corresponding constrained directions for 2d, mieux avec [0, 0.1]
else
  cpoints = [0, 0, 0,    0, 0, 0,   0, 0, 0.1]; % constrained points for 3d
  cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0];   % corresponding constrained directions for 3d, mieux avec [0, 0.1]
end
    
gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints', 'cunitv');

gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);
indmass = gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');



if (approximation_type == 0)
  gf_model_set(md, 'add initialized data', 'N1', [0;-1]);
  gf_model_set(md, 'add_interpolate_transformation_from_expression', 'Proj1', mesh1, mesh2, '[X(1);-0.25]');
  
  gf_model_set(md, 'add filtered fem variable', 'lambda1', mflambda1, CONTACT_BOUNDARY1);
  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, '-lambda1*(Test_u1.N1)+lambda1*(Interpolate(Test_u2,Proj1).N1)', CONTACT_BOUNDARY1);
  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, '-(gamma0*element_size)*(lambda1 + neg_part(lambda1-(1/(gamma0*element_size))*((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).N1)))*Test_lambda1', CONTACT_BOUNDARY1);
    
else
  
  
  gf_model_set(md, 'add initialized data', 'N1', [0;-1]);
  gf_model_set(md, 'add_interpolate_transformation_from_expression', 'Proj1', mesh1, mesh2, '[X(1);-0.25]');
  gamma='(gamma0*element_size)';
  sigma_vh1_n ='(((clambda*Trace(Grad_Test_u1)*Id(qdim(u1)) + cmu*(Grad_Test_u1 + Grad_Test_u1''))*Normal).N1)';
  sigma_uh1_n ='(((clambda*Trace(Grad_u1)*Id(qdim(u1)) + cmu*(Grad_u1 + Grad_u1''))*Normal).N1)';
  Pn_gamma_u1=strcat('((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).N1) - ',gamma,'*(',sigma_uh1_n,')');

  if (approximation_type == 2) coeff = '0.5*'; else coeff = ''; end

  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, strcat('-',coeff,'theta*',gamma,'*(',sigma_uh1_n,')*(',sigma_vh1_n,')'), CONTACT_BOUNDARY1);
  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, strcat('-',coeff,'theta*',sigma_vh1_n,'*pos_part(',Pn_gamma_u1,')'), CONTACT_BOUNDARY1);
  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, strcat(coeff,'(1/',gamma,')*(Test_u1.N1)*pos_part(',Pn_gamma_u1,')'), CONTACT_BOUNDARY1);
  gf_model_set(md, 'add nonlinear generic assembly brick', mim1_contact, strcat('-',coeff,'(1/',gamma,')*(Interpolate(Test_u2,Proj1).N1)*pos_part(',Pn_gamma_u1,')'),CONTACT_BOUNDARY1);

  if (approximation_type == 2)
    gf_model_set(md, 'add initialized data', 'N2', [0;1]);
    gf_model_set(md,'add_interpolate_transformation_from_expression', 'Proj2', mesh2,mesh1, '[ X(1) ; 0.0001-sqrt(0.0625-sqr(X(1)))]');
    sigma_vh2_n ='(((clambda*Trace(Grad_Test_u2)*Id(qdim(u2)) + cmu*(Grad_Test_u2 + Grad_Test_u2'')).Normal).N2)';
    sigma_uh2_n ='(((clambda*Trace(Grad_u2)*Id(qdim(u2)) + cmu*(Grad_u2 + Grad_u2'')).Normal).N2)';
    Pn_gamma_u2=strcat('((u2-Interpolate(u1,Proj2)).N2) - ',gamma,'*(',sigma_uh2_n,')-(Interpolate(X,Proj2)-X).N2');

    gf_model_set(md, 'add nonlinear generic assembly brick', mim2_contact, strcat('-0.5*theta*',gamma,'*(',sigma_uh2_n,')*(',sigma_vh2_n,')'), CONTACT_BOUNDARY2);

    gf_model_set(md, 'add nonlinear generic assembly brick', mim2_contact, strcat('-0.5*theta*',sigma_vh2_n,'*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
    gf_model_set(md, 'add nonlinear generic assembly brick', mim2_contact, strcat('0.5*(1/',gamma,')*(Test_u2.N2)*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
    gf_model_set(md, 'add nonlinear generic assembly brick', mim2_contact, strcat('-0.5*(1/',gamma,')*(Interpolate(Test_u1,Proj2).N2)*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
  end
end














solve=true;
% gf_model_get(md, 'test tangent matrix term', 'u1', 'u2', 1e-6, niter, 10.0);
% gf_model_get(md, 'test tangent matrix', 1e-6, 10, 10);

niter= 100;
gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8);

U1 = gf_model_get(md, 'variable', 'u1');
UU1 = gf_model_get(md, 'variable', 'u1');
VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                          'u1', 'clambda', 'cmu', mfvm1);
U2 = gf_model_get(md, 'variable', 'u2');  
UU2 = gf_model_get(md, 'variable', 'u2');
VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                          'u2', 'clambda', 'cmu', mfvm2);

% plot figure
    if (ref_sol==1)
        sl1=gf_slice({'none'}, mesh1, 5);
        P1=gf_slice_get(sl1,'pts'); dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
        gf_slice_set(sl1, 'pts', P1 + dP1);
        VMsl1=gf_compute(mfvm1,VM1,'interpolate on',sl1);
        sl2=gf_slice({'none'}, mesh2, 5);
        P2=gf_slice_get(sl2,'pts'); dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
        gf_slice_set(sl2, 'pts', P2+dP2);
        VMsl2=gf_compute(mfvm2,VM2,'interpolate on',sl2);
        gf_plot_slice(sl1,'mesh','off','mesh_slice_edges','off','data',VMsl1);
        hold on;
        gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl2);
        hold off;   
    %else
        sl1=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R},{'planar',1,[0;0;0],[1;0;0]}}}, m, 5);
        sl2=gf_slice({'boundary',{'intersection',{'planar',-1,[0;0;zc],[0;0;1]},{'planar',1,[0;0;0],[1;0;0]}}}, m, 5); 
        P1=gf_slice_get(sl1,'pts');P2=gf_slice_get(sl2,'pts');
        dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
        dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
        gf_slice_set(sl1, 'pts', P1+dP1);
        gf_slice_set(sl2, 'pts', P2+dP2);           
        VMsl1=gf_compute(mfvm,VM1,'interpolate on',sl1);
        VMsl2=gf_compute(mfvm,VM2,'interpolate on',sl2);
        set(gcf,'renderer','zbuffer');
        h=gf_plot_slice(sl1,'mesh','off','mesh_faces','on','mesh_slice_edges','on','data',VMsl1);
        hold on;
        h=gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl2);
        hold off;
        view(-55,10); axis on; camlight('headlight'); 
        xlabel('x'); ylabel('y'); zlabel('z');
        %title('3D');       
        
    end
  
gf_plot(mfvm1,VM1,'mesh', 'off', 'deformed_mesh','on', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;
hold on
gf_plot(mfvm2,VM2,'mesh', 'off', 'deformed_mesh','on', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;

return

    if (ref_sol == 0)
        map(1,:)=[0 200 255]; % bleu fonc??
        map(2,:)=[0 255 255];
        map(3,:)=[128 255 128];
        map(4,:)=[255 255 0];
        map(5,:)=[255 128 0];
        map(6,:)=[255 0 0];
        colormap(map./255)
    gf_mesh_fem_get(mfu1, 'save', 'sol_ref_mesh_fem1','with_mesh');
    gf_mesh_fem_get(mfu2, 'save', 'sol_ref_mesh_fem2','with_mesh');
    save sol_de_reference1 UU1;
    save sol_de_reference2 UU2;
  else
    meshref = gf_mesh('load', 'sol_ref_mesh_fem');
   end
end
end
end


























