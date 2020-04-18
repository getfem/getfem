% Copyright (C) 2015-2020 Rabii Mlika, Yves Renard.
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

gf_workspace('clear all');
clear all;
draw_mesh = false;       
                  
% The test case: The numerical tests in two dimensions (resp. three dimensions) are performed on a domain Ω =]−0.5, 0.5[^2 (resp. Ω =]−0.5, 0.5[^3 
% containing the first body: Ω1 , a disk of radius R and center (0,0) (resp. a sphere of radius R and center (0,0,0)), and the second: Ω2 =]−0.5, 0.5[×]−0.5, −R[ 
% (resp. Ω2 =]−0.5, 0.5[2 ×]−0.5, R[). The contact surface Γ_c1 is the lower semicircle and Γ_c2 is the top surface of Ω2 (i.e.Γ1 = {x ∈ ∂Ω1 ; x2 <=0} and 
% Γ_c2 = {x ∈ ∂Ω2 ; x2 = −R}. A Dirichlet condition is prescribed on the bottom of the rectangle (resp. cuboid).Since no Dirichlet condition is applied on Ω1 the problem is only
% semi-coercive,so we apply a penalisation on it and to overcome the non-definiteness coming from the free rigid motions, the horizontal displacement is prescribed to be zero on the two points of coordinates (0,0) and
% (0,0.1) which blocks the horizontal translation and the rigid rotation.The projector Π1 is defined from Γ1 to Γ2 in the vertical direction. All remaining parts of the boundaries are
% considered traction free. The Lame coefficients are λ and μ and we apply a vertical volume density of force on Ω1.
                        
N =2             % 2 or 3 dimensions

R=0.25;                 % Radiaus of Ω1.
dirichlet_val = 0;      % Dirchelet condition.
f_coeff=0;              % friction coefficient.
clambda = 1;            % Lame coefficient λ.
cmu = 1;                % Lame coefficient μ.
vertical_force = -0.1;  % Verticvertical volume density of force on Ω1.
penalty_parameter = 1E-7;    % penalisation parmeter on Ω1.
elelments_degre =2        %  degre of elments (1 or 2).
release_dist = 0.5;
nonlinear_elasticity=true;

approximation_type =1% 0 : Augmentend Lagrangian
                        % 1 : Nitsche (biased)
                        % 2 : Unbiased Nitsche method


 theta = 0;       %   theta
 gamma0 =1/100;   %   Nitsche's parmeter gamma0
 Nxy =50;         %   mesh size (=1/nxy) 2D ->250 and 3D -> 32  
 

 
 

R= 1/(floor(1/R))    
NX = floor(Nxy*R)/R

%mesh constuction
    if     (N==2) 
        mo1 =  gf_mesher_object('ball',[0 0],R); % Ω1
        mesh1 = gf_mesh('generate', mo1, 1/NX ,2) ; 
        mesh2=gf_mesh('regular simplices', -.5:(1/NX):.5, -.5:(1/NX):-R);
    elseif (N==3)
        mo1 =  gf_mesher_object('ball',[0 0 0],R); % Ω1
        mesh1 = gf_mesh('generate', mo1, 1/NX ,2) ; 
        mesh2=gf_mesh('regular simplices', -.5:(1/NX):.5,-.5:(1/NX):.5, -.5:(1/NX):-R);
    end
    

    mfu1 = gf_mesh_fem(mesh1, N) ;gf_mesh_fem_set(mfu1, 'classical fem', elelments_degre);
    mflambda1 = gf_mesh_fem(mesh1, 1); gf_mesh_fem_set(mflambda1, 'classical fem', elelments_degre);

    mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', elelments_degre);

    mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', elelments_degre);

    mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', elelments_degre);

    mim1 = gf_mesh_im(mesh1, 4);
    mim1_contact = gf_mesh_im(mesh1, 6);
    mim2 = gf_mesh_im(mesh2, 4);
    mim2_contact = gf_mesh_im(mesh2, 6);
   
%Contact and dirchelet boundaries.
    CONTACT_BOUNDARY1 = 1;

    fb1 = gf_mesh_get(mesh1, 'outer faces');
    normals1 = gf_mesh_get(mesh1, 'normal of faces', fb1);
    contact_boundary1=fb1(:,find(normals1(N,:) < 0.2)); 
    gf_mesh_set(mesh1,'region', CONTACT_BOUNDARY1,contact_boundary1);
    
    CONTACT_BOUNDARY2 = 2;
    
    border = gf_mesh_get(mesh2,'outer faces');
    normals2 = gf_mesh_get(mesh2, 'normal of faces', border);
       
   P = gf_mesh_get(mesh2, 'pts');
   if (N==2) 
    P=P(:,find(P(1,:)>=-R));
    P=P(:,find(P(1,:)<=R)) ;
    P=P(:,find(P(2,:)>=-R)) ;
   elseif (N==3)
    P=P(:,find((P(1,:)).^2+(P(2,:)).^2<=R^2));
    P=P(:,find(P(3,:)>=-R)) ;
   end
   Id=gf_mesh_get(mesh2, 'pid from coords', P);
   contact_boundary2=gf_mesh_get(mesh2,'faces from pid',Id);

   gf_mesh_set(mesh2, 'region', CONTACT_BOUNDARY2,  contact_boundary2);
    
    
    dirichlet_boundary=border(:, find(normals2(N, :) < -0.01));
    DIRICHLET_BOUNDARY2 = 3;
    gf_mesh_set(mesh2, 'region', DIRICHLET_BOUNDARY2, dirichlet_boundary);

   
%Drawing the Mesh
    if (draw_mesh)
        figure(2);
    gf_plot_mesh(mesh1);
    hold on
    gf_plot_mesh(mesh2, 'regions',CONTACT_BOUNDARY2); % ,  'regions', CONTACT_BOUNDARY1);
    hold off
    pause;
    end
    
%Elastic model 
    md=gf_model('real');
    gf_model_set(md,'add fem variable', 'u1', mfu1);
    gf_model_set(md,'add fem variable', 'u2', mfu2);

    gf_model_set(md,'add initialized data', 'gamma0', gamma0);
    gf_model_set(md, 'add initialized data', 'theta', theta);
    gf_model_set(md, 'add initialized data', 'friction_coeff', [f_coeff]);
     gf_model_set(md, 'add initialized data', 'R', R);

    gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
    gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
    if (nonlinear_elasticity)
      lawname = 'Ciarlet Geymonat';
      params = [clambda;cmu;cmu/2-clambda/8];
      gf_model_set(md,'add initialized data','params', params);
      gf_model_set(md, 'add finite strain elasticity brick', mim1, lawname, 'u1', 'params');
      gf_model_set(md, 'add finite strain elasticity brick', mim2, lawname, 'u2', 'params');
    else
      gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1','clambda', 'cmu');
      gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2','clambda', 'cmu');
    end
    
    gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);
    indmass = gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');
% Boundary conditions: 
    if N==2
       gf_model_set(md, 'add initialized data', 'Fdata', [0 vertical_force]);
       Ddata = zeros(1, 2);
       cpoints = [0,0    0,0.1]; % constrained points for 2d
       cunitv  = [1,0   1,0];   % corresponding constrained directions for 2d.
    elseif N==3
       gf_model_set(md, 'add initialized data', 'Fdata', [0 0 vertical_force]);
       Ddata = zeros(1, 3);
       cpoints = [0, 0, 0,    0, 0, 0,   0, 0, 0.1]; % constrained points for 3d
       cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0];   % corresponding constrained directions for 3d.
    end
    gf_model_set(md, 'add source term brick', mim1, 'u1', 'Fdata');
    gf_model_set(md, 'add initialized data', 'Ddata2', Ddata);
    gf_model_set(md, 'add Dirichlet condition with multipliers', mim2, 'u2', 1, DIRICHLET_BOUNDARY2, 'Ddata2');  
    gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
    gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
    gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints', 'cunitv');

%Contact method:
gf_model_set(md, 'add raytracing transformation', 'contact_trans', release_dist);
gf_model_set(md, 'add slave contact boundary to raytracing transformation', 'contact_trans', mesh1, 'u1', CONTACT_BOUNDARY1);
 gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans', mesh2, 'u2', CONTACT_BOUNDARY2);
    if (approximation_type == 0)
      gf_model_set(md, 'add filtered fem variable', 'lambda1', mflambda1, CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, '-lambda1*(Test_u1.Transformed_unit_vector(Grad_u1, Normal))', CONTACT_BOUNDARY1);
       gf_model_set(md, 'add nonlinear term', mim1_contact, 'Interpolate_filter(contact_trans,lambda1*(Interpolate(Test_u2,contact_trans).Transformed_unit_vector(Grad_u1, Normal)),1)', CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, 'Interpolate_filter(contact_trans,-(gamma0*element_size)*(lambda1 + neg_part(lambda1-(1/(gamma0*element_size))*( (u1-Interpolate(u2,contact_trans))+(X-Interpolate(X,contact_trans)) ).Transformed_unit_vector(Grad_u1, Normal)))*Test_lambda1,1)', CONTACT_BOUNDARY1);

    else
      gamma='(gamma0*element_size)';
      %sigma_vh1 ='((clambda*Trace(Grad_Test_u1)*Id(qdim(u1)) + cmu*(Grad_Test_u1 + Grad_Test_u1''))*Normal)';
      expr = gf_model_get(md, 'Neumann term', 'u1', CONTACT_BOUNDARY1)
      sigma_chap_u1 =expr;
      g1_n = '-((u1-Interpolate(u2,contact_trans)+X-Interpolate(X,contact_trans)).Transformed_unit_vector(Grad_u1, Normal))';
      g1_t='(u1) ';
     
      if (approximation_type == 1) coeff = ''; elseif (approximation_type == 2) coeff = '0.5*'; end
     % gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('-',coeff,'Interpolate_filter(contact_trans,(Test_u1).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Normal,',g1_t,',',g1_n,',[0;friction_coeff;friction_coeff],1/',gamma,')),1)'),CONTACT_BOUNDARY1);
     % gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('',coeff,'Interpolate_filter(contact_trans,Interpolate(Test_u2,contact_trans).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Normal,',g1_t,',',g1_n,',[0;friction_coeff;friction_coeff],1/',gamma,')),1)'),CONTACT_BOUNDARY1);
     %gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,',coeff,'(Test_u1).Transformed_unit_vector(Grad_u1, Normal)*(1/',gamma,')*pos_part(',g1_n,'-',gamma,'*(',sigma_chap_u1,'.Transformed_unit_vector(Grad_u1, Normal))),1)'),CONTACT_BOUNDARY1);
      %gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,-',coeff,'(Interpolate(Test_u2,contact_trans)).Transformed_unit_vector(Grad_u1, Normal)*(1/',gamma,')*pos_part(',g1_n,'-',gamma,'*',sigma_chap_u1,'.Transformed_unit_vector(Grad_u1, Normal)),1)'),CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans, -',coeff,'(Test_u1).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',friction_coeff,1/',gamma,')),1)'),CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('Interpolate_filter(contact_trans,',coeff,'(Interpolate(Test_u2,contact_trans)).(Coulomb_friction_coupled_projection(',sigma_chap_u1,',Transformed_unit_vector(Grad_u1, Normal),',g1_t,',',g1_n,',friction_coeff,1/',gamma,')),1)'),CONTACT_BOUNDARY1);
      if (approximation_type == 2)
       
        sigma_chap_u2=gf_model_get(md, 'Neumann term', 'u2', CONTACT_BOUNDARY2);
        %sigma_vh2 ='((clambda*Trace(Grad_Test_u2)*Id(qdim(u2)) + cmu*(Grad_Test_u2 + Grad_Test_u2'')).Normal)';
        %sigma_uh2 ='((clambda*Trace(Grad_u2)*Id(qdim(u2)) + cmu*(Grad_u2 + Grad_u2'')).Normal)';
        g2_n = '-(( u2-Interpolate(u1,contact_trans2)+X-Interpolate(X,contact_trans2) ).Transformed_unit_vector(Grad_u2, Normal))';
        g2_t='(u2) ';

      %gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('-',coeff,'theta*',gamma,'*((',sigma_uh2,').(',sigma_vh2,'))'), CONTACT_BOUNDARY2);
      %gf_model_set(md, 'add nonlinear term',mim2_contact, strcat('',coeff,'theta*',gamma,'*',sigma_vh2,'.(Coulomb_friction_coupled_projection(',sigma_uh2,',N2,',g2_t,',',g2_n,',[0;friction_coeff*',gamma,';friction_coeff*',gamma,'],1/',gamma,'))'),CONTACT_BOUNDARY2);
      gf_model_set(md, 'add raytracing transformation', 'contact_trans2', release_dist);
      gf_model_set(md, 'add slave contact boundary to raytracing transformation', 'contact_trans2', mesh2, 'u2', CONTACT_BOUNDARY2);
       gf_model_set(md, 'add master contact boundary to raytracing transformation', 'contact_trans2', mesh1, 'u1', CONTACT_BOUNDARY1); 
       
     %gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('Interpolate_filter(contact_trans2,',coeff,'(Test_u2).Transformed_unit_vector(Grad_u2, Normal)*(1/',gamma,')*pos_part(',g2_n,'-',gamma,'*',sigma_chap_u2,'.Transformed_unit_vector(Grad_u2, Normal)),1)'),CONTACT_BOUNDARY2);
     %gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('Interpolate_filter(contact_trans2,-',coeff,'(Interpolate(Test_u1,contact_trans2)).Transformed_unit_vector(Grad_u2, Normal)*(1/',gamma,')*pos_part(',g2_n,'-',gamma,'*',sigma_chap_u2,'.Transformed_unit_vector(Grad_u2, Normal)),1)'),CONTACT_BOUNDARY2);
        gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('Interpolate_filter(contact_trans2, -',coeff,'(Test_u2).(Coulomb_friction_coupled_projection(',sigma_chap_u2,',Transformed_unit_vector(Grad_u2, Normal),',g2_t,',',g2_n,',friction_coeff,1/',gamma,')),1)'),CONTACT_BOUNDARY2);
      gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('Interpolate_filter(contact_trans2,',coeff,'(Interpolate(Test_u1,contact_trans2)).(Coulomb_friction_coupled_projection(', sigma_chap_u2,',Transformed_unit_vector(Grad_u2, Normal),',g2_t,',',g2_n,',friction_coeff,1/',gamma,')),1)'),CONTACT_BOUNDARY2);
      end
    end


%resolution:

    solve=true;
    niter= 100;
    max_res = 1E-8;
    max_iter = 100;

    gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'lsolver', 'mumps');

    U1 = gf_model_get(md, 'variable', 'u1');
    UU1 = gf_model_get(md, 'variable', 'u1');
    VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u1', 'clambda', 'cmu', mfvm1);
    U2 = gf_model_get(md, 'variable', 'u2');  
    UU2 = gf_model_get(md, 'variable', 'u2');
    VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u2', 'clambda', 'cmu', mfvm2);                         
    

% Von Mises stress contour plot (for ref_sol=0):   

         if N==2
            
            gf_plot(mfvm1,VM1,'mesh', 'off', 'deformed_mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;
            hold on
            gf_plot(mfvm2,VM2,'mesh', 'off', 'deformed_mesh','off', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;
         elseif N==3
            sl1=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R},{'planar',1,[0;0;0],[1;0;0]}}}, mesh1, 5);
            sl2=gf_slice({'boundary',{'intersection',{'planar',-1,[0;0;-R],[0;0;1]},{'planar',1,[0;0;0],[1;0;0]}}}, mesh2, 5); 
            P1=gf_slice_get(sl1,'pts');P2=gf_slice_get(sl2,'pts');
            dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
            dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
            gf_slice_set(sl1, 'pts', P1+dP1);
            gf_slice_set(sl2, 'pts', P2+dP2);           
            VMsl1=gf_compute(mfvm1,VM1,'interpolate on',sl1);
            VMsl2=gf_compute(mfvm2,VM2,'interpolate on',sl2);
            set(gcf,'renderer','zbuffer');
            h=gf_plot_slice(sl1,'mesh','off','mesh_faces','on','mesh_slice_edges','on','data',VMsl1);
            hold on;
            h=gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl2);
            hold off;
            view(-55,10); axis on; camlight('headlight'); %gf_colormap('tank');
            xlabel('x'); ylabel('y'); zlabel('z');
            %title('3D');      
         end
