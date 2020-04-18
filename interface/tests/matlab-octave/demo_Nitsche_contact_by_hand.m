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



asize =  size(who('automatic_var654777'));
if (asize(1)) is_automatic = true; else is_automatic = false; end



gf_workspace('clear all');
clear all;
approximation_type = 0  % 0 : Augmentend Lagrangian
                        % 1 : Nitsche (biased)
                        % 2 : Unbiased Nitsche method

draw_mesh = true;       
ref_sol = 0             % 0 : Reference solution (Von Mises)
                        % 1 : Convergence curves in L2 and H1 norms on Omega_1and Omega_2.
                        % 2 : Error as fonction of gamma0 for different values of theta 
                        
% The test case: The numerical tests in two dimensions (resp. three dimensions) are performed on a domain Omega =]-0.5, 0.5[^2 (resp. Omega =]-0.5, 0.5[^3 
% containing the first body: Omega_1 , a disk of radius R and center (0,0) (resp. a sphere of radius 0.25 and center (0,0,0)), and the second: Omega_2 =]-0.5, 0.5[ times ]-0.5, -0.25[
% (resp. Omega_2 =]-0.5, 0.5[2 times ]-0.5, 0.25[). The contact surface Gamma_c1 is the lower semicircle and Gamma_c2 is the top surface of Omega_2 (i.e.Gamma_1 = {x in partial Omega_1 ; x2 <=0} and 
% Gamma_c2 = {x in partial Omega_2 ; x2 = -0.25}. A Dirichlet condition is prescribed on the bottom of the rectangle (resp. cuboid).Since no Dirichlet condition is applied on Omega_1 the problem is only
% semi-coercive,so we apply a penalisation on it and to overcome the non-definiteness coming from the free rigid motions, the horizontal displacement is prescribed to be zero on the two points of coordinates (0,0) and
% (0,0.1) which blocks the horizontal translation and the rigid rotation.The projector PI_1 is defined from Gamma_1 to Gamma_2 in the vertical direction. All remaining parts of the boundaries are
% considered traction free. The Lame coefficients are lambda and mu and we apply a vertical volume density of force on Omega_1.
                        
N = 2                   % 2 or 3 dimensions

R=0.25;                 % Radius of Omega_1.
dirichlet_val = 0;      % Dirchelet condition.
f_coeff=0;              % friction coefficient.
clambda = 1;            % Lame coefficient lambda.
cmu = 1;                % Lame coefficient mu.
vertical_force = -0.1;  % Vertical volume density of force on Omega_1.
penalty_parameter = 1E-7;    % penalisation parmeter on Omega_1.
elements_degree = 2          %  degree of elements (1 or 2).
    
 if (ref_sol == 0)
    Theta = [-1];       %   theta
    Gamma0 = [1/100];   %   Nitsche's parmeter gamma0
    Nxy =[50];          %   mesh size (=1/nxy) 2D ->250 and 3D -> 25  
 

 elseif (ref_sol == 1)  
    Theta= [1];
    Gamma0 = [1/100];
     if (N==2)
        Nxy=[10 15 27 37 40 50 60 65 70 75];
     else
        nxy = [5 10 12 15 20 25];
     end
 elseif (ref_sol == 2)
    Theta = [0 1 -1];
    Gamma0 = [400 200 100 50 25 10 1 1/10 1/25 1/50 1/100 1/200 1/400];
    Nxy = [21];
 end

for xx = 1:1:size(Theta,2)
for yy = 1:1:size(Gamma0,2)
for zz = 1:1:size(Nxy,2)
    
theta = Theta(xx);
gamma0 = Gamma0(yy);
NX = Nxy(zz)

%mesh constuction
    if     (N==2) 
        mo1 =  gf_mesher_object('ball',[0 0],R); % Omega_1
        mesh1 = gf_mesh('generate', mo1, 1/NX ,4) ; 
        mo2=gf_mesher_object('rectangle', [-0.5 -0.5], [0.5 -0.25]); %  Omega_2
        mesh2 = gf_mesh('generate', mo2, 1/NX ,2) ; 
    elseif (N==3)
        mo1 =  gf_mesher_object('ball',[0 0 0],R); % Omega_1
        mesh1 = gf_mesh('generate', mo1, 1/NX ,2) ; 
        mo2=gf_mesher_object('rectangle', [-0.5 -0.5 -0.5], [0.5  0.5 -0.25]); %  Omega_2
        mesh2 = gf_mesh('generate', mo2, 1/NX ,2) ; 
    end
    

    mfu1 = gf_mesh_fem(mesh1, N) ;gf_mesh_fem_set(mfu1, 'classical fem', elements_degree);
    mflambda1 = gf_mesh_fem(mesh1, 1); gf_mesh_fem_set(mflambda1, 'classical fem', elements_degree);

    mfvm1 = gf_mesh_fem(mesh1); gf_mesh_fem_set(mfvm1, 'classical discontinuous fem', elements_degree);

    mfu2 = gf_mesh_fem(mesh2, N); gf_mesh_fem_set(mfu2, 'classical fem', elements_degree);

    mfvm2 = gf_mesh_fem(mesh2); gf_mesh_fem_set(mfvm2, 'classical discontinuous fem', elements_degree);

    mim1 = gf_mesh_im(mesh1, 4);
    mim1_contact = gf_mesh_im(mesh1, 6);
    mim2 = gf_mesh_im(mesh2, 4);
    mim2_contact = gf_mesh_im(mesh2, 6);
   
%Contact and dirchelet boundaries.
    CONTACT_BOUNDARY1 = 1;

    fb1 = gf_mesh_get(mesh1, 'outer faces');
    normals1 = gf_mesh_get(mesh1, 'normal of faces', fb1);
    contact_boundary1=fb1(:,find(normals1(N,:) < 0.1)); 
    gf_mesh_set(mesh1,'region', CONTACT_BOUNDARY1,contact_boundary1);
    
    CONTACT_BOUNDARY2 = 2;
    
    border = gf_mesh_get(mesh2,'outer faces');
    normals2 = gf_mesh_get(mesh2, 'normal of faces', border);
       
   P = gf_mesh_get(mesh2, 'pts');
   if (N==2) 
    P=P(:,find(P(1,:)>=-0.25));
    P=P(:,find(P(1,:)<=0.25)) ;
    P=P(:,find(P(2,:)>=-0.25)) ;
   elseif (N==3)
    P=P(:,find(P(1,:)>=-0.25));
    P=P(:,find(P(1,:)<=0.25)) ;
    P=P(:,find(P(2,:)>=-0.25));
    P=P(:,find(P(2,:)<=0.25)) ;
    P=P(:,find(P(3,:)>=-0.25)) ;
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
    gf_plot_mesh(mesh2);
    hold on
    gf_plot_mesh(mesh1, 'refine', 8, 'curved', 'on'); % ,  'regions', CONTACT_BOUNDARY1);
    hold off
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
    gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1','clambda', 'cmu');
    gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2','clambda', 'cmu');
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
    if N==2 
         gf_model_set(md, 'add initialized data', 'N1', [0;-1]);
         gf_model_set(md, 'add_interpolate_transformation_from_expression', 'Proj1', mesh1, mesh2, '[X(1);-0.25]');
    elseif N==3
        gf_model_set(md, 'add initialized data', 'N1', [0;0;-1]);
        gf_model_set(md, 'add_interpolate_transformation_from_expression', 'Proj1', mesh1, mesh2, '[X(1);X(2);-0.25]');
    end
  
    if (approximation_type == 0)
      gf_model_set(md, 'add filtered fem variable', 'lambda1', mflambda1, CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, '-lambda1*(Test_u1.N1)+lambda1*(Interpolate(Test_u2,Proj1).N1)', CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, '-(gamma0*element_size)*(lambda1 + neg_part(lambda1-(1/(gamma0*element_size))*((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).N1)))*Test_lambda1', CONTACT_BOUNDARY1);

    else
      gamma='(gamma0*element_size)';
      sigma_vh1_n ='(((clambda*Trace(Grad_Test_u1)*Id(qdim(u1)) + cmu*(Grad_Test_u1 + Grad_Test_u1''))*Normal).N1)';
      sigma_uh1_n ='(((clambda*Trace(Grad_u1)*Id(qdim(u1)) + cmu*(Grad_u1 + Grad_u1''))*Normal).N1)';
      Pn_gamma_u1=strcat('((u1-Interpolate(u2,Proj1)+X-Interpolate(X,Proj1)).N1) - ',gamma,'*(',sigma_uh1_n,')');
 
      if (approximation_type == 1) coeff = ''; elseif (approximation_type == 2) coeff = '0.5*'; end
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('-',coeff,'theta*',gamma,'*(',sigma_uh1_n,')*(',sigma_vh1_n,')'), CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('-',coeff,'theta*',sigma_vh1_n,'*pos_part(',Pn_gamma_u1,')'), CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat(coeff,'(1/',gamma,')*(Test_u1.N1)*pos_part(',Pn_gamma_u1,')'), CONTACT_BOUNDARY1);
      gf_model_set(md, 'add nonlinear term', mim1_contact, strcat('-',coeff,'(1/',gamma,')*(Interpolate(Test_u2,Proj1).N1)*pos_part(',Pn_gamma_u1,')'),CONTACT_BOUNDARY1);
     
      if (approximation_type == 2)
         if N==2 
         gf_model_set(md, 'add initialized data', 'N2', [0;1]);
         gf_model_set(md,'add_interpolate_transformation_from_expression', 'Proj2', mesh2,mesh1, '[ X(1) ; 0.0001-sqrt(sqr(R)-sqr(X(1)))]');
        elseif N==3
        gf_model_set(md, 'add initialized data', 'N2', [0;0;1]);
         gf_model_set(md,'add_interpolate_transformation_from_expression', 'Proj2', mesh2,mesh1, '[ X(1) ; X(2) ; 0.01-sqrt(sqr(R)-(sqr(X(1))+sqr(X(2))) )]');
         end
        sigma_vh2_n ='(((clambda*Trace(Grad_Test_u2)*Id(qdim(u2)) + cmu*(Grad_Test_u2 + Grad_Test_u2'')).Normal).N2)';
        sigma_uh2_n ='(((clambda*Trace(Grad_u2)*Id(qdim(u2)) + cmu*(Grad_u2 + Grad_u2'')).Normal).N2)';
        Pn_gamma_u2=strcat('((u2-Interpolate(u1,Proj2)).N2) - ',gamma,'*(',sigma_uh2_n,')-(Interpolate(X,Proj2)-X).N2');
        
        gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('-0.5*theta*',gamma,'*(',sigma_uh2_n,')*(',sigma_vh2_n,')'), CONTACT_BOUNDARY2);
        gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('-0.5*theta*',sigma_vh2_n,'*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
        gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('0.5*(1/',gamma,')*(Test_u2.N2)*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
        gf_model_set(md, 'add nonlinear term', mim2_contact, strcat('-0.5*(1/',gamma,')*(Interpolate(Test_u1,Proj2).N2)*pos_part(',Pn_gamma_u2,')'), CONTACT_BOUNDARY2);
      end
    end


%resolution:
    solve=true;
    niter= 100;
    max_res = 1E-8;
    max_iter = 100;

    % gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'lsolver', 'mumps');
    gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy', 'lsolver', 'mumps');

    U1 = gf_model_get(md, 'variable', 'u1');
    UU1 = gf_model_get(md, 'variable', 'u1');
    VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                              'u1', 'clambda', 'cmu', mfvm1);
    U2 = gf_model_get(md, 'variable', 'u2');  
    UU2 = gf_model_get(md, 'variable', 'u2');
    VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
                              'u2', 'clambda', 'cmu', mfvm2);
%Van Mises stress contour plot (for ref_sol=0):                          
    if ref_sol==0                      
         if N==2
            gf_plot(mfvm1,VM1,'mesh', 'off', 'deformed_mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;
            hold on
            gf_plot(mfvm2,VM2,'mesh', 'off', 'deformed_mesh','off', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;
         elseif N==3
            sl1=gf_slice({'boundary',{'intersection',{'ball',-1,[0;0;0],R},{'planar',1,[0;0;0],[1;0;0]}}}, mesh1, 5);
            sl2=gf_slice({'boundary',{'intersection',{'planar',-1,[0;0;-0.25],[0;0;1]},{'planar',1,[0;0;0],[1;0;0]}}}, mesh2, 5); 
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

%Saving the reference solution:
        gf_mesh_get(mesh1, 'save', 'sol_ref_mesh1');
        gf_mesh_get(mesh2, 'save', 'sol_ref_mesh2');

        gf_mesh_fem_get(mfu1, 'save', 'sol_ref_mesh_fem1');
        gf_mesh_fem_get(mfu2, 'save', 'sol_ref_mesh_fem2');
        %gf_mesh_im_get(mim1,'save','sol_ref_mim1');
        %gf_mesh_im_get(mim2,'save','sol_ref_mim2');
        save sol_de_reference1 UU1;
        save sol_de_reference2 UU2;
       
    else
% Reconstruction of reference mesh: mesh_ref1, mesh_ref2, mfu_ref1, mfu_ref2, min_ref1, min_ref2
        
        mesh_ref1 = gf_mesh('load', 'sol_ref_mesh1');
        mesh_ref2 = gf_mesh('load', 'sol_ref_mesh2');
       mfu_ref1 = gf_mesh_fem('load', 'sol_ref_mesh_fem1',mesh_ref1);
       mfu_ref2 = gf_mesh_fem('load', 'sol_ref_mesh_fem2',mesh_ref2);

        N =gf_mesh_get(mesh_ref2,'dim');
        %mfu_ref1 = gf_mesh_fem(mesh_ref1, N); gf_mesh_fem_set(mfu_ref1, 'classical fem', elements_degree);
        %mfu_ref2 = gf_mesh_fem(mesh_ref2, N);gf_mesh_fem_set(mfu_ref2, 'classical fem', elements_degree);
        mim_ref1 = gf_mesh_im(mesh_ref1, 4);
        mim_ref2 = gf_mesh_im(mesh_ref2, 4);


        U1ref = load('sol_de_reference1');
        U2ref = load('sol_de_reference2');

 %U1ee = gf_compute(mfu_ref1, U1ref.UU1, 'interpolate on', mfu_ref1);

        U1e = gf_compute(mfu1, U1, 'interpolate on', mfu_ref1);
        U2e = gf_compute(mfu2, U2, 'interpolate on', mfu_ref2);
        
%Calculation of the displacement relatif error (in percentage) for H1 and L2 norms:
        n_tot1 = gf_compute(mfu_ref1, U1e-U1ref.UU1, 'L2 norm', mim_ref1);
        n_tot2 = gf_compute(mfu_ref2, U2e-U2ref.UU2, 'L2 norm', mim_ref2);
        n_ref1 = gf_compute(mfu_ref1, U1ref.UU1, 'L2 norm',  mim_ref1);
        n_ref2 = gf_compute(mfu_ref2, U2ref.UU2, 'L2 norm',  mim_ref2);
        m_tot1 = gf_compute(mfu_ref1, U1e-U1ref.UU1, 'H1 norm', mim_ref1);
        m_tot2 = gf_compute(mfu_ref2, U2e-U2ref.UU2, 'H1 norm', mim_ref2);
        m_ref1 = gf_compute(mfu_ref1, U1ref.UU1, 'H1 norm',  mim_ref1);
        m_ref2 = gf_compute(mfu_ref2, U2ref.UU2, 'H1 norm',  mim_ref2);
        n1 = 100*n_tot1/n_ref1; 
        n2 = 100*n_tot2/n_ref2;
        m1 = 100*m_tot1/m_ref1;
        m2 = 100*m_tot2/m_ref2;
        Y11(yy,zz,xx)=n1;
        Y12(yy,zz,xx)=n2;
        Y21(yy,zz,xx)=m1;
        Y22(yy,zz,xx)=m2;
    end
end
end
end

if (ref_sol == 1 ) % Curves of error depending on h
    Y11(1,:,1)
    Y12(1,:,1)
    Y21(1,:,1)
    Y22(1,:,1)
    X=1./Nxy 

   
    figure(2);
    msize= size(X,2);
    loglog(X,Y11(1,:,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y12(1,:,1),'+--k', 'linewidth', 2, 'MarkerSize', 15);
    hold off;
    P1 = polyfit(log(X),log(Y11(1,:,1)),1);
    P2 = polyfit(log(X),log(Y12(1,:,1)),1);
    legend(strcat('norm on Omega 1  (slope=',num2str(P1(1)), ')'), ...
           strcat('norm on Omega 2 (slope=',num2str(P2(1)), ')'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('h');
    ylabel('L^2 relative error (in %)');
    set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
    

    figure(3);
 
    loglog(X,Y21(1,:,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y22(1,:,1),'+--k', 'linewidth', 2, 'MarkerSize', 15);
    hold off;
    P3 = polyfit(log(X),log(Y21(1,:,1)),1);
    P4 = polyfit(log(X),log(Y22(1,:,1)),1);
    legend(strcat('norm on Omega 1  (slope=',num2str(P3(1)), ')'), ...
           strcat('norm on Omega 2 (slope=',num2str(P4(1)), ')'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('h');
    ylabel('H^1 relative error (in %)');
    set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
       
end

if (ref_sol == 2 )  % Curve of error depending on gamma0
    
    % Theta = [0 1 -1];
    
    Y11(:,1,:)
    Y12(:,1,:)
    Y21(:,1,:)
    Y22(:,1,:)
    X = gamma 

    
    figure(2);
    

    loglog(X,Y11(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y11(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y11(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y11(:,1,1)'),1);
    P2 = polyfit(log(X),log(Y11(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y11(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('L^2(Omega 1) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'})  
    
    
    figure(3);
    

    loglog(X,Y12(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y12(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y12(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y12(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y12(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y12(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('L^2(Omega 2) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 
   
    
    
     figure(4);
    

    loglog(X,Y21(:,1,1)','o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y21(:,1,2)','+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y21(:,1,3)','x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y21(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y21(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y21(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('H^1(Omega 1) relative error (in %)');
    % set(gca,'XTickLabel',{'0.01';'0.1';'1';'...'}) 

    
    figure(5);
    

    loglog(X,Y22(:,1,1),'o-k', 'linewidth', 2, 'MarkerSize', 15 )
    hold on;
    loglog(X,Y22(:,1,2),'+-k', 'linewidth', 2, 'MarkerSize', 15 );
    loglog(X,Y22(:,1,3),'x-k', 'linewidth', 2, 'MarkerSize', 15 );
    hold off;
        
    P1 = polyfit(log(X),log(Y22(:,1,1)'),1); % the first and second are too bad;
    P2 = polyfit(log(X),log(Y22(:,1,2)'),1);
    P3 = polyfit(log(X),log(Y22(:,1,3)'),1);
    legend(strcat('norm for theta = 0 '), ...
           strcat('norm for theta = 1 '), ...
           strcat('norm for theta = -1'), ...
           'Location', 'NorthWest');
    grid on;
    axesobj = findobj('type', 'axes');
    set(axesobj, 'fontname', 'times'); set(axesobj, 'fontunits', 'points');
    set(axesobj, 'fontsize', 18); set(axesobj, 'fontweight', 'bold');
    set(axesobj, 'linewidth', 4);
    xlabel('gamma0');
    ylabel('H^1(Omega 2) relative error (in %)');
   
end

