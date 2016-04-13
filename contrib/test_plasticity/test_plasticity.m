% Copyright (C) 2010-2015 Yves Renard, Farshid Dabaghi.
%
% This file is a part of GetFEM++
%
% GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
clc;


% This program perform a convergence test for several implementations of
% small strain isotropic plasticity in GetFEM++, in particular different
% complementarity functions 


option = 3; % 1 : without hardening, without im_data, without plastic multiplier
            % 2 : without hardening, with plastic multiplier
            % 3 : with hardening, with plastic multiplier
            % 4 : with kinematic hardening, with im_data, without plastic multiplier
            
load_type = 2; % 1 : vertical
               % 2 : horizontal
               
bi_material = false;
test_tangent_matrix = false;
do_plot = true;


if (load_type == 2 && option < 3)
    error('Will not work with this load : will break the body');
end
               





% Initialize used data
T = 10;
NT = 40;
LX = 100;
LY = 20;
NX = 40;

NY = ceil(NX * LY / (2 * LX))*2;
DT = T/NT;

% theta is the parameter of the generalized mid point scheme.
% theta = 1/2 for mid point scheme and theta = 1 for backward Euler scheme.
theta = 1;


if (load_type == 1)
    f = [0 -600]';
    t = sin(2*pi*(0:DT:T)/20);
else
    f = [15000 0]';
    t = sin(2*pi*(0:DT:T)/10) + 0.1;
end;


% Create the mesh
% m = gfMesh('triangles grid', [0:(LX/NX):LX], [0:(LY/NY):LY]);
m = gfMesh('import','structured',sprintf('GT="GT_PK(2,1)";SIZES=[%d,%d];NOISED=0;NSUBDIV=[%d,%d];', LX, LY, NX, NY));
N = gf_mesh_get(m, 'dim');
  
% Plotting
% gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% Define used MeshIm
mim=gfMeshIm(m);  set(mim, 'integ', gfInteg('IM_TRIANGLE(6)')); % Gauss methods on triangles

% Define used MeshFem
if (option >= 2)
  mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,2)'));
else
  mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,1)'));
  mf_sigma=gfMeshFem(m,4); set(mf_sigma, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
end
% mf_xi = gfMeshFem(m); set(mf_xi, 'fem', gfFem('FEM_PK(2,2)'));
mf_xi = gfMeshFem(m); set(mf_xi, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_data=gfMeshFem(m); set(mf_data, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_vm = gfMeshFem(m); set(mf_vm, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));

% Find the border of the domain
P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6); % Retrieve index of points which x near to 0
pidright=find(abs(P(1,:) - LX)<1e-6); % Retrieve index of points which x near to L
fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);
set(m,'boundary',1,fleft); % for Dirichlet condition
set(m,'boundary',2,fright); % for Neumann condition

% Decomposed the mesh into 2 regions with different values of Lamé coeff
if (bi_material) separation = LY/2; else separation = 0; end
pidtop    = find(P(2,:)>=separation-1E-6); % Retrieve index of points of the top part
pidbottom = find(P(2,:)<=separation+1E-6); % Retrieve index of points of the bottom part
cvidtop   = get(m, 'cvid from pid', pidtop);
cvidbottom= get(m, 'cvid from pid', pidbottom);
CVtop     = sort(get(mf_data, 'basic dof from cvid', cvidtop));
CVbottom  = sort(get(mf_data, 'basic dof from cvid', cvidbottom));

% Definition of Lame coeff
lambda(CVbottom,1) = 121150; % Steel
lambda(CVtop,1) = 84605; % Iron
mu(CVbottom,1) = 80769; %Steel
mu(CVtop,1) = 77839; % Iron
% Definition of plastic threshold
von_mises_threshold(CVbottom) = 7000;
von_mises_threshold(CVtop) = 8000;
% Definition of Kinematic hardening parameter
if (option >= 3)
  Hk = mu(1)/5; Hi = Hk;
else
  Hk = 0; Hi = Hk;
end

% Create the model
md = gfModel('real');
set(md, 'add fem variable', 'u', mf_u);
set(md, 'add initialized fem data', 'lambda', mf_data, lambda);
set(md, 'add initialized fem data', 'mu', mf_data, mu);
set(md, 'add initialized fem data', 'von_mises_threshold', mf_data, von_mises_threshold);


switch (option)
  case 1
    gf_model_set(md, 'add fem data', 'Previous_u', mf_u);
    % Declare that sigma is a data of the system on mf_sigma
    set(md, 'add fem data', 'sigma', mf_sigma);
    % Add plasticity brick on u
    set(md, 'add elastoplasticity brick', mim, 'VM', 'u', 'Previous_u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
  case 2
    set(md, 'add fem variable', 'xi', mf_xi);
    set(md, 'add fem data', 'Previous_u', mf_u);
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add initialized data', 'r', [1e-8]);
    mim_data = gf_mesh_im_data(mim, -1, [N, N]);
    mim_data_scal = gf_mesh_im_data(mim, -1, 1);
    set(md, 'add im data', 'Epn', mim_data);
    
    if (theta == 1)
      Etheta = '(Sym(Grad_u))';
      Eptheta = strcat('((Epn+2*mu*xi*Deviator(',Etheta,'))/(1+2*mu*xi))');
      Epnp1 = Eptheta;
      sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
      sigma_theta = sigma_np1;
    else
      Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
      Eptheta = strcat('((Epn+2*mu*theta*xi*Deviator(',Etheta,'))/(1+2*mu*theta*xi))');
      Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
      sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
      sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
    end
    
    fbound = strcat('(Norm(Deviator(',sigma_theta,')) - von_mises_threshold)');
    % fbound = strcat('(Norm(',Eptheta,'-Epn)-theta*xi*von_mises_threshold)');
    expr = strcat(sigma_np1, ':Grad_Test_u + (1/r)*(xi - pos_part(xi+r*',fbound,'))*Test_xi');
    % expr = strcat(sigma_np1, ':Grad_Test_u + (',fbound,' + pos_part(-xi/r-',fbound,'))*Test_xi');
    gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
        
  case 3
    set(md, 'add fem variable', 'xi', mf_xi);
    set(md, 'add fem data', 'Previous_u', mf_u);
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add initialized data', 'r', [1e-8]);
    mim_data = gf_mesh_im_data(mim, -1, [N, N]);
    mim_data_scal = gf_mesh_im_data(mim, -1, 1);
    set(md, 'add im data', 'Epn', mim_data);
    set(md, 'add im data', 'alphan', mim_data_scal);
    set(md, 'add initialized data', 'Hk', [Hk]);
    set(md, 'add initialized data', 'Hi', [Hi]);
    
    if (theta == 1)
      Etheta = '(Sym(Grad_u))';
      Eptheta = strcat('((Epn+2*mu*xi*Deviator(',Etheta,'))/(1+(2*mu+Hk)*xi))');
      Epnp1 = Eptheta;
      sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
      sigma_theta = sigma_np1;
      alpha_theta = strcat('(alphan+xi*(Norm(2*mu*Deviator(',Etheta,')-(2*mu+Hk)*',Eptheta,')))');
      alpha_np1 = alpha_theta;
      % alpha_theta = strcat('(alphan+Norm(',Eptheta,'-Epn))'); % do not work
      % alpha_np1 = alpha_theta; % do not work
    else
      Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
      Eptheta = strcat('((Epn+2*mu*theta*xi*Deviator(',Etheta,'))/(1+(2*mu+Hk)*theta*xi))');
      Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
      sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
      sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
      % alpha_theta = strcat('alphan+theta*xi*(Norm(Deviator(',sigma_theta,')-Hk*',Eptheta,'))');
      % alpha_np1 = strcat('alphan+xi*(Norm(Deviator(',sigma_theta,')-Hk*',Eptheta,'))');
      % alpha_theta = strcat('(alphan+Norm(',Eptheta,'-Epn))');  % do not work
      % alpha_np1 = strcat('(alphan+Norm(',Eptheta,'-Epn)/theta)'); % do not work
    end
    
    % fbound = strcat('(Norm(Deviator(',sigma_theta,')-Hk*',Eptheta,') - von_mises_threshold - Hi*',alpha_theta,')');
    fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+Hk)*',Eptheta,') - von_mises_threshold - Hi*',alpha_theta,')');
    expr = strcat(sigma_np1, ':Grad_Test_u + (1/r)*(xi - pos_part(xi+r*',fbound,'))*Test_xi');
    % expr = strcat(sigma_np1, ':Grad_Test_u + (',fbound,' + pos_part(-xi/r-',fbound,'))*Test_xi');
    gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    
      
  case 4
    gf_model_set(md, 'add fem data', 'Previous_u', mf_u);
    mim_data = gf_mesh_im_data(mim, -1, [N, N]);
    gf_model_set(md, 'add im data', 'sigma', mim_data);
  
    % Declare that theta is a data of the system 
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add initialized data', 'Hk', [Hk]);

    Is = 'Reshape(Id(meshdim*meshdim),meshdim,meshdim,meshdim,meshdim)';
    IxI = 'Id(meshdim)@Id(meshdim)';
    coeff_long = '((lambda)*(Hk))/((2*(mu)+(Hk))*(meshdim*(lambda)+2*(mu)+(Hk)))';
    B_inv = sprintf('((2*(mu)/(2*(mu)+(Hk)))*(%s) + (%s)*(%s))', Is, coeff_long, IxI);
    B = sprintf('((1+(Hk)/(2*(mu)))*(%s) - (((lambda)*(Hk))/(2*(mu)*(meshdim*(lambda)+2*(mu))))*(%s))', Is, IxI);
    ApH = sprintf('((2*(mu)+(Hk))*(%s) + (lambda)*(%s))', Is, IxI);
    Enp1 = '((Grad_u+Grad_u'')/2)';
    En = '((Grad_Previous_u+Grad_Previous_u'')/2)';
  
    % Expression de sigma for Implicit Euler method
    % expr_sigma = strcat('(', B_inv, '*(Von_Mises_projection((-(Hk)*', Enp1, ')+(', ApH, '*(',Enp1,'-',En,')) + (', B, '*sigma), von_mises_threshold) + Hk*', Enp1, '))');
  
    % Expression de sigma for generalized mid-point scheme
    expr_sigma = strcat('(', B_inv, '*(Von_Mises_projection((',B,'*((1-theta)*sigma))+(-(Hk)*(((1-theta)*',En ...
        ,')+(theta*', Enp1, ')))+(theta*', ApH, '*(',Enp1,'-',En,')) + (theta*', B, '*sigma), von_mises_threshold) + (Hk)*(((1-theta)*',En,')+(theta*', Enp1, '))))');
  
    gf_model_set(md, 'add nonlinear generic assembly brick', mim, strcat(expr_sigma, ':Grad_Test_u'));
end

% Add homogeneous Dirichlet condition to u on the left hand side of the domain
set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);

% Add a source term to the system
set(md,'add initialized fem data', 'VolumicData', mf_data, get(mf_data, 'eval',{f(1,1)*t(1);f(2,1)*t(1)}));
set(md, 'add source term brick', mim, 'u', 'VolumicData', 2);

VM=zeros(1,get(mf_vm, 'nbdof'));



for step=1:size(t,2),
    disp(sprintf('step %d / %d, coeff = %g', step, size(t,2), t(step)));
    set(md, 'variable', 'VolumicData', get(mf_data, 'eval',{f(1,1)*t(step);f(2,1)*t(step)}));

    if (test_tangent_matrix)
      gf_model_get(md, 'test tangent matrix', 1E-8, 10, 0.000001);
    end;
   
    % Solve the system
    get(md, 'solve', 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'max_iter', 50, 'max_res', 1e-6);
    % get(md, 'solve', 'noisy', 'max_iter', 80);

    % Retrieve the solution U
    U = get(md, 'variable', 'u');
    
    % Compute new plasticity constraints used to compute 
    % the Von Mises or Tresca stress
    switch (option)
      case 1
        get(md, 'elastoplasticity next iter', mim, 'u', 'Previous_u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
        plast = get(md, 'compute plastic part', mim, mf_vm, 'u', 'Previous_u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
        % Compute Von Mises or Tresca stress
        VM = get(md, 'compute elastoplasticity Von Mises or Tresca', 'sigma', mf_vm, 'Von Mises');
      case 2
        von_mises = strcat('sqrt(3/2)*Norm(Deviator(', sigma_np1, '))');
        VMGP = get(md, 'interpolation', von_mises, mim_data_scal);
        M = gf_asm('mass matrix', mim, mf_vm);
        L = gf_asm('generic', mim, 1, 'vmgp*Test_vm', -1, 'vmgp', 0, mim_data_scal, VMGP, 'vm', 1, mf_vm, zeros(get(mf_vm, 'nbdof'),1));
        VM = (M \ L)';
        PLGP = get(md, 'interpolation', strcat('Norm(',Epnp1,')'), mim_data_scal);
        L = gf_asm('generic', mim, 1, 'plgp*Test_vm', -1, 'plgp', 0, mim_data_scal, PLGP, 'vm', 1, mf_vm, zeros(get(mf_vm, 'nbdof'),1));
        plast = (M \ L)';
        NewEpn = get(md, 'interpolation', Epnp1, mim_data);
        set(md, 'variable', 'Epn', NewEpn);
        gf_model_set(md, 'variable', 'Previous_u', U);
        % gf_model_set(md, 'variable', 'xi', zeros(get(mf_xi, 'nbdof'), 1));
        % xi = get(md, 'variable', 'xi')
        
      case 3
        von_mises = strcat('sqrt(3/2)*Norm(Deviator(', sigma_np1, '))');
        VMGP = get(md, 'interpolation', von_mises, mim_data_scal);
        M = gf_asm('mass matrix', mim, mf_vm);
        L = gf_asm('generic', mim, 1, 'vmgp*Test_vm', -1, 'vmgp', 0, mim_data_scal, VMGP, 'vm', 1, mf_vm, zeros(get(mf_vm, 'nbdof'),1));
        VM = (M \ L)';
        PLGP = get(md, 'interpolation', strcat('Norm(',Epnp1,')'), mim_data_scal);
        L = gf_asm('generic', mim, 1, 'plgp*Test_vm', -1, 'plgp', 0, mim_data_scal, PLGP, 'vm', 1, mf_vm, zeros(get(mf_vm, 'nbdof'),1));
        plast = (M \ L)';
        
        
        sigma = get(md, 'interpolation', sigma_np1, mim_data);
        Epsilon_u = gf_model_get(md, 'interpolation', '((Grad_u+Grad_u'')/2)', mim_data);
        ind_gauss_pt = 22500;
        if (size(sigma, 2) <= N*N*(ind_gauss_pt + 1))
          ind_gauss_pt = floor(3*size(sigma, 2) / (4*N*N));
        end
        sigma_fig(1,step)=sigma(N*N*ind_gauss_pt + 1);
        Epsilon_u_fig(1,step)=Epsilon_u(N*N*ind_gauss_pt + 1);
        
        
        NewEpn = get(md, 'interpolation', Epnp1, mim_data);
        set(md, 'variable', 'Epn', NewEpn);
        Newalphan = get(md, 'interpolation', alpha_np1, mim_data_scal);
        set(md, 'variable', 'alphan', Newalphan);
        gf_model_set(md, 'variable', 'Previous_u', U);
        
      case 4
        sigma_0 = gf_model_get(md, 'variable', 'sigma');
        sigma = gf_model_get(md, 'interpolation', expr_sigma, mim_data);
        U_0 = gf_model_get(md, 'variable', 'Previous_u');
        U_ntheta = theta*U + (1-theta)*U_0;
      
        M = gf_asm('mass matrix', mim, mf_vm);
        L = gf_asm('generic', mim, 1, 'sqrt(3/2)*Norm(Deviator(sigma))*Test_vm', -1, 'sigma', 0, mim_data, sigma, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1));
        VM = (M \ L)';
        coeff1='-lambda/(2*mu*(meshdim*lambda+2*mu))';
        coeff2='1/(2*mu)';
        Ainv=sprintf('(%s)*(%s) + (%s)*(%s)', coeff1, IxI, coeff2, Is);
        Ep = sprintf('(Grad_u+Grad_u'')/2 - (%s)*sigma', Ainv);
        L = gf_asm('generic', mim, 1, sprintf('Norm(%s)*Test_vm', Ep), -1, 'sigma', 0, mim_data, sigma, 'u', 0, mf_u, U, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1), 'mu', 0, mf_data, mu, 'lambda', 0, mf_data, lambda);
        plast = (M \ L)';
      
        gf_model_set(md, 'variable', 'u', U_ntheta);
        Epsilon_u = gf_model_get(md, 'interpolation', '((Grad_u+Grad_u'')/2)', mim_data);
        gf_model_set(md, 'variable', 'u', U);
        ind_gauss_pt = 22500;
        if (size(sigma, 2) <= N*N*(ind_gauss_pt + 1))
          ind_gauss_pt = floor(3*size(sigma, 2) / (4*N*N));
        end
        sigma_fig(1,step)=sigma(N*N*ind_gauss_pt + 1);
        Epsilon_u_fig(1,step)=Epsilon_u(N*N*ind_gauss_pt + 1);
      
        sigma = (sigma - (1-theta)*sigma_0)/theta;
        gf_model_set(md, 'variable', 'sigma', sigma);
        gf_model_set(md, 'variable', 'Previous_u', U);
    end
      
       
    if (do_plot)
      figure(2)
      subplot(3,1,1);
      gf_plot(mf_vm,VM, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0); % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Von Mises criterion for t = ', num2str(step)]);

      subplot(3,1,2);
      gf_plot(mf_vm,plast, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0);  % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Plastification for t = ', num2str(step)]);
    
      if (option >= 3)
        subplot(3,1,3);
        plot(Epsilon_u_fig, sigma_fig,'r','LineWidth',2)
        xlabel('Strain');
        ylabel('Stress')
        axis([-0.15 0.35 -16000 25000 ]);
        % hold on;
      end;
      
      pause(0.1);
    end
 
end;








