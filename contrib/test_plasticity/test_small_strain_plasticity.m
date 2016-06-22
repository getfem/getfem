% Copyright (C) 2010-2016 Yves Renard, Farshid Dabaghi.
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

gf_util('trace level', 1);


% This program performs a test for several implementations of
% small strain isotropic plasticity in GetFEM++


option = 2; % 1 : without hardening, without im_data, without plastic multiplier
            % 2 : without hardening, with plastic multiplier
            % 3 : with kinematic and isotropic hardening, with plastic multiplier
            % 4 : with kinematic and isotropic hardening, with im_data, without plastic multiplier
            % 5 : Souza-Auricchio model with plastic multiplier
            
load_type = 1; % 1 : vertical
               % 2 : horizontal
               
bi_material = false;
test_tangent_matrix = false;
do_plot = true;
use_small_strain_pl_brick = 1; % Use the (new) small strain plasticity brick when possible


if (load_type == 2 && option < 3)
    error('Will not work with this load : will break the body');
end
               

% Physical parameters
% E = 70000;
% nu = 0.33;
% lambda_top = E*nu/((1+nu)*(1-2*nu));
% mu_top = E / (2*(1+nu));

lambda_top = 84605;      % Iron
mu_top = 77839;          % Iron
lambda_bottom = 121150;  % Steel
mu_bottom = 80769;       % Steel
sigma_y_top = 8000;
sigma_y_bottom = 7000;

Hk = mu_top/5; Hi = 0; % Kinematic and isotropic hardening parameters


% Numerica parameters
T = 10;
NT = 40;
LX = 100;
LY = 20;
NX = 40;

NY = ceil(NX * LY / (2 * LX))*2;
DT = T/NT;

theta = 1.0; % Parameter for the generalized mid point scheme.


if (load_type == 1)
    f = [0 -500]';
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
mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,2)'));
mf_sigma=gfMeshFem(m,2,2); set(mf_sigma, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
% mf_xi = gfMeshFem(m); set(mf_xi, 'fem', gfFem('FEM_PK(2,2)'));
mf_xi = gfMeshFem(m); set(mf_xi, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_delta = gfMeshFem(m); set(mf_delta, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_data=gfMeshFem(m); set(mf_data, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_vm = gfMeshFem(m); set(mf_vm, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));

% Find the boundary of the domain
P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6); % Retrieve index of points which x near to 0
pidright=find(abs(P(1,:) - LX)<1e-6); % Retrieve index of points which x near to L
fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);
set(m,'boundary',1,fleft); % for Dirichlet condition
set(m,'boundary',2,fright); % for Neumann condition

% Decompose the mesh into two regions with different values of Lamé coefficients
if (bi_material) separation = LY/2; else separation = 0; end
pidtop    = find(P(2,:)>=separation-1E-6); % Retrieve index of points of the top part
pidbottom = find(P(2,:)<=separation+1E-6); % Retrieve index of points of the bottom part
cvidtop   = get(m, 'cvid from pid', pidtop);
cvidbottom= get(m, 'cvid from pid', pidbottom);
CVtop     = sort(get(mf_data, 'basic dof from cvid', cvidtop));
CVbottom  = sort(get(mf_data, 'basic dof from cvid', cvidbottom));

% Definition of Lame coeff
lambda(CVbottom) = lambda_bottom;
lambda(CVtop) = lambda_top;
mu(CVbottom) = mu_bottom;
mu(CVtop) = mu_top;
sigma_y(CVbottom) = sigma_y_bottom;
sigma_y(CVtop) = sigma_y_top;

% Create the model
md = gfModel('real');
set(md, 'add fem variable', 'u', mf_u);
set(md, 'add fem data', 'Previous_u', mf_u);
set(md, 'add initialized fem data', 'lambda', mf_data, lambda);
set(md, 'add initialized fem data', 'mu', mf_data, mu);
set(md, 'add initialized fem data', 'sigma_y', mf_data, sigma_y);
set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);
set(md, 'add initialized fem data', 'VolumicData', mf_data, zeros(get(mf_data, 'nbdof')*N,1))
set(md, 'add source term brick', mim, 'u', 'VolumicData', 2);
set(md, 'set time step', DT);
mim_data = gf_mesh_im_data(mim, -1, [N, N]);
mim_data_scal = gf_mesh_im_data(mim, -1, 1);

switch (option)
  case 1
    if (use_small_strain_pl_brick)
      set(md, 'add fem data', 'xi', mf_xi);
      set(md, 'add fem data', 'Previous_xi', mf_xi);
      set(md, 'add initialized data', 'theta', [theta]);
      set(md, 'add im data', 'Epn', mim_data);
      set(md, 'add small strain elastoplasticity brick', mim, 'Prandtl Reuss', false, 'u', 'xi', 'Epn', 'lambda', 'mu', 'sigma_y', 'theta', 'timestep');
    else      
      % Declare that sigma is a data of the system on mf_sigma
      set(md, 'add fem data', 'sigma', mf_sigma);
      % Add old plasticity brick on u
      set(md, 'add elastoplasticity brick', mim, 'VM', 'u', 'Previous_u', 'lambda', 'mu', 'sigma_y', 'sigma');
    end
    
  case 2
    if (use_small_strain_pl_brick)
      set(md, 'add fem variable', 'xi', mf_xi);
      set(md, 'add fem data', 'Previous_xi', mf_xi);
      set(md, 'add initialized data', 'theta', [theta]);
      set(md, 'add im data', 'Epn', mim_data);
      % set(md, 'add fem data', 'Epn', mf_sigma);
      
      set(md, 'add small strain elastoplasticity brick', mim, 'Prandtl Reuss', true, 'u', 'xi', 'Epn', 'lambda', 'mu', 'sigma_y', 'theta', 'timestep');
    else
      set(md, 'add fem variable', 'xi', mf_xi);
      set(md, 'add initialized data', 'theta', [theta]);
      set(md, 'add initialized data', 'r', [1e-8]);
      set(md, 'add im data', 'Epn', mim_data);
    
      if (theta == 1 || true)
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
    
      fbound = strcat('(Norm(Deviator(',sigma_theta,')) - sqrt(2/3)*sigma_y)');
      % fbound = strcat('(Norm(',Eptheta,'-Epn)-theta*xi*sigma_y)');
      expr = strcat(sigma_np1, ':Grad_Test_u+(1/r)*(xi-pos_part(xi+r*',fbound,'))*Test_xi');
      % expr = strcat(sigma_np1, ':Grad_Test_u + (',fbound,' + pos_part(-xi/r-',fbound,'))*Test_xi');
      gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    end  
  case 3
    set(md, 'add fem variable', 'xi', mf_xi);
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add initialized data', 'r', [1e-8]);
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
    else
      Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
      Eptheta = strcat('((Epn+2*mu*theta*xi*Deviator(',Etheta,'))/(1+(2*mu+Hk)*theta*xi))');
      Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
      sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
      sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
      alpha_theta = strcat('(alphan+theta*xi*(Norm(2*mu*Deviator(',Etheta,')-(2*mu+Hk)*',Eptheta,')))');
      % alpha_np1 = strcat('((', alpha_theta, ' - (1-theta)*alphan)/theta)');
      alpha_np1 = strcat('(alphan+xi*(Norm(2*mu*Deviator(',Etheta,')-(2*mu+Hk)*',Eptheta,')))');
      
      % alpha_theta = strcat('(alphan+Norm(',Eptheta,'-Epn))');  % do not work
      % alpha_np1 = strcat('(alphan+Norm(',Eptheta,'-Epn)/theta)'); % do not work
    end
    
    % fbound = strcat('(Norm(Deviator(',sigma_theta,')-Hk*',Eptheta,') - sigma_y - Hi*',alpha_theta,')');
    fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+Hk)*',Eptheta,') - sqrt(2/3)*(sigma_y + Hi*',alpha_theta,'))');
    expr = strcat(sigma_np1, ':Grad_Test_u + (1/r)*(xi - pos_part(xi+r*',fbound,'))*Test_xi');
    % expr = strcat(sigma_np1, ':Grad_Test_u + (',fbound,' + pos_part(-xi/r-',fbound,'))*Test_xi');
    gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    
  case 4
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add im data', 'Epn', mim_data);
    set(md, 'add im data', 'alphan', mim_data_scal);
    set(md, 'add initialized data', 'Hk', [Hk]);
    set(md, 'add initialized data', 'Hi', [Hi]);
    
    Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
    Btheta = strcat('((2*mu)*Deviator(',Etheta,')-(2*mu+Hk)*Epn)');
    alpha_theta = strcat('(max(alphan, ((2*mu+Hk)*alphan+Norm(',Btheta,') - sqrt(2/3)*sigma_y)/(2*mu+Hk+sqrt(2/3)*Hi)))');
    alpha_np1 = strcat('((', alpha_theta, ' - (1-theta)*alphan)/theta)');
    Eptheta = strcat('(Epn+(1/(2*mu+Hk))*pos_part(1-sqrt(2/3)*(sigma_y+Hi*',alpha_theta,')/(Norm(',Btheta,')+1e-25))*',Btheta,')');
    Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
    sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
    sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
    
    expr = strcat(sigma_np1, ':Grad_Test_u');
    gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    
  case 5
    set(md, 'add fem variable', 'xi', mf_xi);
    set(md, 'add initialized data', 'theta', [theta]);
    set(md, 'add initialized data', 'r1', [1e-8]);
    set(md, 'add initialized data', 'r2', [1]);
    set(md, 'add im data', 'Epn', mim_data);
    set(md, 'add initialized data', 'c1', [0]); % [7.5*3]);
    set(md, 'add initialized data', 'c2', [Hk]);
    set(md, 'add initialized data', 'c3', [0.15]); % [0.03]);
    
    
    
    if (1) % Version with two multipliers
        set(md, 'add fem variable', 'delta', mf_delta);
        Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
        Btheta = strcat('(Epn+theta*xi*2*mu*Deviator(',Etheta,'))');
        Eptheta = strcat('((',Btheta,')/(1+(2*mu+c2+delta)*theta*xi))'); % version without c1
        % Eptheta = strcat('(',Btheta,'*pos_part(1-theta*xi*c1/(Norm(',Btheta,')+1E-6))/(1+(2*mu+c2)*theta*xi + theta*delta))');
        Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
        sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
        sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
    
        
        fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2+delta)*',Eptheta,') - sqrt(2/3)*sigma_y)'); % version without c1
        
        
        % fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2+delta)*',Eptheta,'-c1*Normalized_reg(',Eptheta,',1E-6)) - sqrt(2/3)*sigma_y)');
        fbound_delta = strcat('(Norm(',Eptheta,')-c3)');
        expr = strcat(sigma_np1, ':Grad_Test_u + (10/r1)*(xi - pos_part(xi+r1*',fbound,'))*Test_xi - (100/r2)*(delta - pos_part(delta+r2*',fbound_delta,'))*Test_delta');
        gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    
    else
    
        Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))';
        Btheta = strcat('(Epn+theta*xi*2*mu*Deviator(',Etheta,'))');
        Eptheta = strcat('((',Btheta,')*min(c3/(max(Norm(',Btheta,'), c3/2)), pos_part(1-theta*xi*c1/(Norm(',Btheta,')+0.001))/(1+(2*mu+c2)*theta*xi)))');
        Eptheta = strcat('(',Btheta,'*min(c3/(max(Norm(',Btheta,'),c3/2)), 1/(1+(2*mu+c2)*theta*xi)))'); % version without c1
        % Eptheta = strcat('(',Btheta,'*min(c3/(Norm(',Btheta,')+1e-10), 1/(1+(2*mu+c2)*theta*xi)))');
        Epnp1 = strcat('((', Eptheta, ' - (1-theta)*Epn)/theta)');
        sigma_np1 = strcat('(lambda*Trace(Sym(Grad_u)-',Epnp1, ')*Id(meshdim) + 2*mu*(Sym(Grad_u)-', Epnp1,'))');
        sigma_theta = strcat('(lambda*Trace(',Etheta,'-',Eptheta, ')*Id(meshdim) + 2*mu*(',Etheta,'-', Eptheta,'))');
      
        % fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2)*',Eptheta,'-max(c1, (Norm(',Btheta,')-c3)/(theta*xi+1e-25)-(2*mu+c2)*c3)*Normalized(',Eptheta,')) - sqrt(2/3)*sigma_y)');
        fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2)*',Eptheta,'-pos_part( pos_part(Norm(',Btheta,')/(theta*xi+1e-10) - c1)/c3 - (1/(theta*xi+1e-10)+2*mu+c2))*',Eptheta,') - sqrt(2/3)*sigma_y)');
        % fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2)*',Eptheta,'-pos_part( Norm(',Btheta,')/(c3*(theta*xi+1e-10)) - (1/(theta*xi+1e-10)+2*mu+c2))*',Eptheta,') - sqrt(2/3)*sigma_y)');
        % fbound = strcat('(Norm(2*mu*Deviator(',Etheta,')-(2*mu+c2)*',Eptheta,') - sqrt(2/3)*sigma_y)');
        expr = strcat(sigma_np1, ':Grad_Test_u + (1/r1)*(xi - pos_part(xi+r1*',fbound,'))*Test_xi');
        % expr = strcat(sigma_np1, ':Grad_Test_u + (',fbound,' + pos_part(-xi/r1-',fbound,'))*Test_xi');
        gf_model_set(md, 'add nonlinear generic assembly brick', mim, expr);
    end
    
      
end

VM=zeros(1,get(mf_vm, 'nbdof'));



for step=1:size(t,2),
    disp(sprintf('step %d / %d, coeff = %g', step, size(t,2), t(step)));
    source = get(mf_data, 'eval',{f(1,1)*t(step);f(2,1)*t(step)});
    set(md, 'variable', 'VolumicData', source);

    if (test_tangent_matrix)
      gf_model_get(md, 'test tangent matrix', 1E-8, 10, 0.000001);
    end;
    
    if (option == 5)
       set(md, 'variable', 'delta', zeros(1, get(mf_delta, 'nbdof')));
       set(md, 'variable', 'xi', zeros(1, get(mf_xi, 'nbdof')));
    end
   
    % Solve the system
    get(md, 'solve', 'noisy', 'max_iter', 50, 'lsearch', 'simplest',  'alpha min', 0.5, 'max_res', 1e-6);
    
    if (option == 5)
       delta = get(md, 'variable', 'delta');
       norm(delta)
    end

    U = get(md, 'variable', 'u');
    
    % Compute new plastic internal variables
    switch (option)
      case 1
        if (use_small_strain_pl_brick)
            get(md, 'small strain elastoplasticity next iter', mim, 'Prandtl Reuss', false, 'u', 'xi', 'Epn', 'lambda', 'mu', 'sigma_y', 'theta', 'timestep');
        else
          get(md, 'elastoplasticity next iter', mim, 'u', 'Previous_u', 'VM', 'lambda', 'mu', 'sigma_y', 'sigma');
          plast = get(md, 'compute plastic part', mim, mf_vm, 'u', 'Previous_u', 'VM', 'lambda', 'mu', 'sigma_y', 'sigma');
          % Compute Von Mises or Tresca stress
          VM = get(md, 'compute elastoplasticity Von Mises or Tresca', 'sigma', mf_vm, 'Von Mises');
        end
      case 2
        if (use_small_strain_pl_brick)
          get(md, 'small strain elastoplasticity next iter', mim, 'Prandtl Reuss', true, 'u', 'xi', 'Epn', 'lambda', 'mu', 'sigma_y', 'theta', 'timestep');         
        else
          NewEpn = get(md, 'interpolation', Epnp1, mim_data);
          set(md, 'variable', 'Epn', NewEpn);
          set(md, 'variable', 'Previous_u', U);
        end
      case 3       
        NewEpn = get(md, 'interpolation', Epnp1, mim_data);
        Newalphan = get(md, 'interpolation', alpha_np1, mim_data_scal);
        norm(Newalphan)
        set(md, 'variable', 'Epn', NewEpn);
        set(md, 'variable', 'alphan', Newalphan);
        set(md, 'variable', 'Previous_u', U);
        
      case 4
        NewEpn = get(md, 'interpolation', Epnp1, mim_data);
        Newalphan = get(md, 'interpolation', alpha_np1, mim_data_scal);
        set(md, 'variable', 'Epn', NewEpn);
        set(md, 'variable', 'alphan', Newalphan);
        gf_model_set(md, 'variable', 'Previous_u', U);
        
      case 5
        NewEpn = get(md, 'interpolation', Epnp1, mim_data);
        set(md, 'variable', 'Epn', NewEpn);
        gf_model_set(md, 'variable', 'Previous_u', U);
    end
    
    % Compute Von Mises and plastic part for graphical post-treatment
    if (do_plot)
      if (option == 1 && ~use_small_strain_pl_brick) 
        sigma1 = 'sigma';
      else
        Ep = strcat('Norm(Epn)');
        sigma1 = '(lambda*Trace(Sym(Grad_u))*Id(meshdim) + 2*mu*(Sym(Grad_u)-Epn))';
        von_mises = strcat('sqrt(3/2)*Norm(Deviator(', sigma1, '))');
        VM = get(md, 'local projection', mim, von_mises, mf_vm);
        plast = get(md, 'local projection', mim, Ep, mf_vm);
      end
      
      sigma = get(md, 'interpolation', sigma1, mim_data);
      Epsilon_u = gf_model_get(md, 'interpolation', 'Sym(Grad_u)', mim_data);
      ind_gauss_pt = 22500;
      if (size(sigma, 2) <= N*N*(ind_gauss_pt + 1))
        ind_gauss_pt = floor(3*size(sigma, 2) / (4*N*N));
      end
      sigma_fig(1,step)=sigma(N*N*ind_gauss_pt + 1);
      Epsilon_u_fig(1,step)=Epsilon_u(N*N*ind_gauss_pt + 1);
    
      figure(2)
      subplot(3,1,1);
      gf_plot(mf_vm,VM, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0); % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Von Mises criterion for t = ', num2str(step)]);

      subplot(3,1,2);
      % gf_plot(mf_vm,plast, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0);  % 'deformed_mesh', 'on')
      gf_plot(mf_vm,plast,'refine', 4, 'disp_options', 0);  % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Plastification for t = ', num2str(step)]);
    
      subplot(3,1,3);
      plot(Epsilon_u_fig, sigma_fig,'r','LineWidth',2)
      xlabel('Strain');
      ylabel('Stress')
      axis([-0.15 0.35 -16000 25000]);
            
      pause(0.1);
    end
 
end;








