% Copyright (C) 2009-2020 Yves Renard, Franz Chouly.
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
%
% Elastodynamic problem with unilateral contact with a rigid obstacle.
% Newmark, theta-method schemes with Nitsche's method or singular
% mass method.
%
% This program is used to check that matlab-getfem is working. This is also
% a good example of use of GetFEM.
%

gf_workspace('clear all');
clear all;
gf_util('trace level', 0);
figure(1);

NX = 10; m=gf_mesh('cartesian', [0:1/NX:1]); % Cas 1D

% Import the mesh : disc
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h4.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h2.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h1.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h0_5.mesh');
% m=gf_mesh('load', '../../../tests/meshes/disc_P2_h0_3.mesh');

% Import the mesh : sphere
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_8_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_80_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_2000_elts.mesh');
% m=gf_mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_16000_elts.mesh');


d = gf_mesh_get(m, 'dim'); % Mesh dimension


% Parameters of the model
if (d == 1)
  clambda = 1;             % Lame coefficient
  cmu = 0;                 % Lame coefficient
  vertical_force = 0.0;    % Volumic load in the vertical direction
  r = 10;                  % Augmentation parameter
  gamma0_N = 200;          % Parameter gamma0 for Nitsche's method
  % gamma0_N = 1e+4;       % Parameter gamma0 for Nitsche's method
  % penalty_coeff = gamma0_N*NX; % Penalty coefficient for penalty method
  penalty_coeff = clambda*NX;      % Penalty coefficient for penalty method
  int_penalty_coeff = 1;   % Penalty coefficient for experimental method
  dt = 0.001;              % Time step
  T = 12;                   % Simulation time
  dt_plot = 0.04;          % Drawing step;
  dirichlet = 1;           % Dirichlet condition or not
  dirichlet_val = 0.0;
  u_degree = 1;
  v_degree = 1;
  lambda_degree = 1;
else
  clambda = 20;            % Lame coefficient
  cmu = 20;                % Lame coefficient
  vertical_force = 0.1;    % Volumic load in the vertical direction
  r = 10;                  % Augmentation parameter
  gamma0_N = 1000;         % Parameter gamma0 for Nitsche's method
  penalty_coeff = 100;     % Penalty coefficient for penalty method
  dt = 0.01;               % Time step
  T = 40;                  % Simulation time
  dt_plot = 0.5;           % Drawing step;
  dirichlet = 0;           % Dirichlet condition or not
  dirichlet_val = 0.45;
  u_degree = 2;
  v_degree = 1;
  lambda_degree = 1;
end
  
friction = 0;              % Friction coefficient

beta = 0.25;               % Newmark scheme coefficient
gamma = 0.5;               % Newmark scheme coefficient
theta = 0.5;               % Theta-method scheme coefficient

Contact_option = 3;        % 0 : Lagrange multiplier, 1 : Nitsche, 2 : Penalty method, 3 : experimental method
theta_N = 1;               % Parameter theta for Nitsche's method

singular_mass = 0;         % 0 = standard method
                           % 1 = Mass elimination on boundary
                           % 2 = Mixed displacement/velocity

scheme = 3;                % 1 = theta-method, 2 = Newmark, 3 = Newmark with beta = 0, 4 = midpoint modified (Experimental: compile GetFEM with --enable-experimental)
niter = 100;               % Maximum number of iterations for Newton's algorithm.
plot_mesh = false;
make_movie = 0;
residual = 1E-8;

if ((scheme >= 3 && (Contact_option < 1 || singular_mass ~= 0)) || (scheme == 3 && Contact_option < 1))
    error('Incompatibility');
end

if (friction ~= 0 && d == 1)
    error('Not taken into account');
end

% To store 
% - the energy / augmented energy / augmentation term
% - the displacement / velocity / velocity at midpoint
% - the normal and contact stress / contact stress at midpointMass elimination on boundary
Etime = zeros(1,round(T/dt)+1);
Eatime = zeros(1,round(T/dt)+1);
Rtime = zeros(1,round(T/dt)+1);
Utime = zeros(1,round(T/dt)+1);
Vtime = zeros(1,round(T/dt)+1);
Vmtime = zeros(1,round(T/dt)+1);
Stime = zeros(1,round(T/dt)+1);  % Direct computation of the normal stress
SRtime = zeros(1,round(T/dt)+1); % Computation with VR of the normal stress
CStime = zeros(1,round(T/dt)+1); % Contact stress
CSmtime = zeros(1,round(T/dt)+1); % Contact stress at midpoint
errL2 = zeros(1,round(T/dt)+1);
errH1 = zeros(1,round(T/dt)+1);

% Compute and display the wave speed and the Courant number
c0 = sqrt((clambda+2*cmu)); % rho = 1 ...
courantN = c0 * dt * NX; % h = 1/NX; 
fprintf('Wave speed c0 = %g\n', c0);
fprintf('Courant number nuC = %g\n', courantN);

% Signed distance representing the obstacle
if (d == 1) obstacle = 'x'; elseif (d == 2) obstacle = 'y'; else obstacle = 'z'; end;

% Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2;

border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(d, :) < -0.01));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
dirichlet_boundary=border(:, find(normals(d, :) > 0.01));
gf_mesh_set(m, 'region', GAMMAD, dirichlet_boundary);

% Finite element methods
mfu=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfu, 'classical fem', u_degree);
mfv=gf_mesh_fem(m, d);
gf_mesh_fem_set(mfv, 'classical fem', v_degree);
mfd=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfd, 'classical fem', u_degree);
if (friction == 0)
  mflambda=gf_mesh_fem(m, 1);
else
  mflambda=gf_mesh_fem(m, d); 
end
gf_mesh_fem_set(mflambda, 'classical fem', lambda_degree);
mfvm=gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm, 'classical discontinuous fem', u_degree-1);

% Integration method
mim=gf_mesh_im(m, 4);
mim_sing=gf_mesh_im(m);
if (d == 1)
  mim_friction = mim;
elseif (d == 2)
  mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'));
elseif (d == 3)
   mim_friction=gf_mesh_im(m, ...
      gf_integ('IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),4)')); 
end;


first_elem = -1;

M = gf_asm('mass matrix', mim, mfu);

if (singular_mass == 1) % Rought singular mass matrix (no redistribution)
  for i = gf_mesh_get(m, 'cvid')
      if (size(find(contact_boundary(1,:) == i), 2) == 0) 
          gf_mesh_im_set(mim_sing, 'integ', 4, i);
      end
  end
  M_sing = gf_asm('mass matrix', mim_sing, mfu);
  M = M_sing;
end

if (singular_mass == 2)
  B = gf_asm('mass matrix', mim, mfv, mfu);
  C = gf_asm('mass matrix', mim, mfv, mfv);
end

if (dirichlet == 1 && scheme == 3) % penalisation of homogeneous Dirichlet condition for explicit scheme
    GD = gf_asm('mass matrix', mim, mfu, mfu, GAMMAD);
    M = M + 1e13*GD'*GD;
end

% Plot the mesh
if (plot_mesh) hold off;
  figure(1);
  gf_plot_mesh(m, 'regions', [GAMMAC]);
  title('Mesh and contact boundary (in red)');
  pause(0.1);
end;

nbdofd = gf_mesh_fem_get(mfd, 'nbdof');
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');

% Volumic density of force
F = zeros(d, nbdofd);
F(d,:) = -vertical_force;

% Elasticity model
md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add fem data', 'wt', mfu);
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');
            
if (singular_mass == 2)
  gf_model_set(md, 'add fem variable', 'v', mfv);
  switch(scheme)
    case 1
      gf_model_set(md, 'add explicit matrix', 'u', 'v', (B')/(dt*dt*theta*theta));
    case 2
      gf_model_set(md, 'add explicit matrix', 'u', 'v', (B')/(dt*dt*beta));
  end
  gf_model_set(md, 'add explicit matrix', 'v', 'v', C);
  gf_model_set(md, 'add explicit matrix', 'v', 'u', -B);
else
  switch(scheme)
    case 1
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*theta*theta));
    case 2
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*beta));
    case 4
      gf_model_set(md, 'add explicit matrix', 'u', 'u', M/(dt*dt*0.25));
  end
end
ind_rhs = gf_model_set(md, 'add explicit rhs', 'u', zeros(nbdofu,1));

gf_model_set(md, 'add initialized fem data', 'volumicload', mfd, F);

gf_model_set(md, 'add source term brick', mim, 'u', 'volumicload');

if (dirichlet)
  dirichletdata = zeros(1,d); dirichletdata(d) = dirichlet_val;
  gf_model_set(md, 'add initialized data', 'dirichletdata', dirichletdata);
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, GAMMAD, 'dirichletdata');
end

OBS = gf_mesh_fem_get(mfd, 'eval', { obstacle });
gf_model_set(md, 'add initialized fem data', 'obstacle', mfd, OBS);

if (Contact_option == 1)
  gf_model_set(md, 'add initialized data', 'gamma0', [gamma0_N]);
  gf_model_set(md, 'add initialized data', 'theta_N', [theta_N]);
  expr_Neumann = gf_model_get(md, 'Neumann term', 'u', GAMMAC);
  if (scheme == 4)
      if (friction ~= 0)
         error('To be adapted for friction');
      end
      gf_model_set(md, 'add initialized data', 'friction_coeff', [0]);
      gf_model_set(md, 'add initialized data', 'alpha_f', [0]);
      
      % Very experimental brick : compile GetFEM with the option --enable-experimental
      expr_Neumann_wt = '((clambda*Div_wt)*Normal + (cmu*(Grad_wt+(Grad_wt'')))*Normal)';
      gf_model_set(md, 'add Nitsche midpoint contact with rigid obstacle brick', mim_friction, 'u', ...
      expr_Neumann, expr_Neumann_wt, 'obstacle', 'gamma0', GAMMAC, theta_N, 'friction_coeff', 'alpha_f', 'wt');
      
  else
    if (friction == 0)
      gf_model_set(md, 'add Nitsche contact with rigid obstacle brick', mim_friction, 'u', ...
                   expr_Neumann, 'obstacle', 'gamma0', GAMMAC, theta_N);
    else
      gf_model_set(md, 'add initialized data', 'friction_coeff', [friction]);
      gf_model_set(md, 'add initialized data', 'alpha_f', [1./dt]);
      gf_model_set(md, 'add Nitsche contact with rigid obstacle brick', mim_friction, 'u', ...
           expr_Neumann, 'obstacle', 'gamma0', GAMMAC, theta_N, 'friction_coeff', 'alpha_f', 'wt');
    end
  end
elseif (Contact_option == 0)
  ldof = gf_mesh_fem_get(mflambda, 'dof on region', GAMMAC);
  mflambda_partial = gf_mesh_fem('partial', mflambda, ldof);
  gf_model_set(md, 'add fem variable', 'lambda', mflambda_partial);
  gf_model_set(md, 'add initialized data', 'r', [r]);
  if (friction == 0)
    gf_model_set(md, 'add integral contact with rigid obstacle brick', ...
                 mim_friction, 'u', 'lambda', 'obstacle', 'r', GAMMAC, 1);
  else
    gf_model_set(md, 'add initialized data', 'friction_coeff', [friction]);
    gf_model_set(md, 'add initialized data', 'alpha_f', [1./dt]);
    gf_model_set(md, 'add integral contact with rigid obstacle brick', mim_friction, 'u', ...
                 'lambda', 'obstacle', 'r', 'friction_coeff', GAMMAC, 1, 'alpha_f', 'wt');
  end
elseif (Contact_option == 2)
    gf_model_set(md, 'add initialized data', 'penalty_coeff', [penalty_coeff]);
    gf_model_set(md, 'add nonlinear term', mim_friction, 'penalty_coeff*sqr(neg_part(obstacle + u.(Normalized(Grad_obstacle))))', GAMMAC);
    if (friction ~= 0)
        error('Sorry friction to be taken into account for penalty method')
    end
elseif (Contact_option == 3)
    gf_model_set(md, 'add nonlinear term', mim_friction, '0', GAMMAC); % In order to have at list a nonlinear term ...
    gf_model_set(md, 'add fem data', 'uel', mfu);
    if (scheme ~= 3)
        error('experimental method only implemented for explicit scheme, sorry')
    end
    if (friction ~= 0)
        error('Sorry friction to be taken into account for decoupled dynamic method')
    end
else
    error('Invalid contact option')
end

if (d == 1)
    U0 = (gf_mesh_fem_get(mfu, 'eval', { sprintf('%g+0.5-0.5*x', dirichlet_val)}))';
else
    U0 = zeros(nbdofu, 1);
    U0(d:d:nbdofu) = 5;
end;
if (singular_mass == 2)
    VV1 = B*U0; VV2 = C\VV1; MU0 = (B')*VV2;
else
    MU0 = M*U0;
end

MV0 = zeros(nbdofu, 1);
V1 = zeros(nbdofu, 1);
FF = gf_asm('volumic source', mim, mfu, mfd, F);
% K = gf_asm('linear elasticity', mim, mfu, mfd, ones(nbdofd,1)*clambda, ones(nbdofd,1)*cmu);
K = gf_asm('generic', mim, 2, '(clambda*Div_u*Id(meshdim)+2*cmu*Sym(Grad_u)):Grad_Test_u', -1, md);
Iu = gf_model_get(md, 'interval of variable', 'u'); Iu = Iu(1):(Iu(1)+Iu(2)-1);
K = K(Iu, Iu); K2 = K;
if (Contact_option == 1)
  KK = gf_asm('generic', mim, 2, '(-theta_N*element_size/gamma0)*((clambda*Div_u*Id(meshdim)+2*cmu*Sym(Grad_u))*Normal).((clambda*Div_Test_u*Id(meshdim)+2*cmu*Sym(Grad_Test_u))*Normal)', GAMMAC, md);
  K2 = K + KK(Iu,Iu);
end
MA0 = FF-K2*U0;

if (singular_mass == 1)
  if (d == 1)
    MA0(1) = 0;
  else
    error('Take it into account !');
  end
elseif (singular_mass == 2) % to be verified
  VV1 = (B') \ MA0; VV2 = C*VV1; A0 = B\VV2; VV1 = B*A0; VV2 = C\VV1; MA0 = (B')*VV2;
end
nit = 0; tplot = 0;
if (make_movie)
  nim = 0;
  figure('position', [100 100 800 600]); % Necessary for the movie to be read by vlc
  mov = avifile('toto.avi');
end

% Init the vector of physical quantities
E = (U0'*K*U0)/2 - FF'*U0;
Etime(1) = E;
Utime(1) = U0(1);
Vtime(1) = 0;
Stime(1) = -(clambda + 2*cmu)*(NX/1)*(U0(1)-U0(2));
SRtime(1) = - ( U0'*K(:,1) + MA0(1) - FF(1) );
CStime(1) = 0; % Contact stress (initial) -> 0 in this case
               % (in fact, no contact stress)

Rtime(1) = (0.5/gamma0_N)*(1/NX)*( Stime(1)^2 - CStime(1)^2 );
Eatime(1) = Etime(1) - theta_N * Rtime(1);
Vmtime(1) = NaN; % No mid time-step already
CSmtime(1) = NaN; % No mid time-step already

errL2(1) = 0;
errH1(1) = 0;

for t = dt:dt:T
  disp(sprintf('t=%g', t));
  % calcul de LL
  
  switch(scheme)
    case 1
      LL = (MU0 + dt*MV0)/(dt*dt*theta*theta) + (1-theta)*MA0/theta;
    case 2
      LL = (MU0 + dt*MV0 + dt*dt*(1/2-beta)*MA0)/(beta*dt*dt);
    case 3
      LL = 0*MU0;
    case 4
      LL = MU0/(dt*dt*0.25) + MV0/(dt*0.5);
  end
  
  gf_model_set(md, 'variable', 'wt', U0);
  
  if (scheme == 3)
    A0 = M \ MA0;
    V0 = M \ MV0;
    U1 = U0 + dt*dt*A0/2 + dt*V0;
    gf_model_set(md, 'variable', 'u', U1);
    gf_model_get(md, 'assembly', 'build_rhs');
  else
    gf_model_set(md, 'set private rhs', ind_rhs, LL);
    gf_model_get(md, 'solve', 'max_res', residual, 'max_iter', niter); % , 'noisy');
    U1 = (gf_model_get(md, 'variable', 'u'))';
  end
  
  if (singular_mass == 2)
    MU1 = (B')*(gf_model_get(md, 'variable', 'v'))';
  else
    MU1 = M*U1;
  end
  
  if (d == 1)
    fprintf('u1(1) = %g\n', U1(1));
    Msize = size(M,1);
    if (Contact_option == 0)
      lambda = gf_model_get(md, 'variable', 'lambda');
      fprintf('lambda_n = %g\n', lambda(1));
    end
  
    fprintf('U0(N) = %g\n', U0(Msize));
    fprintf('U1(N) = %g\n', U1(Msize));
    fprintf('MV0(N) = %g\n', MV0(Msize));
  end
  
  switch(scheme)
    case 1
      MV1 = ((MU1 - MU0)/dt -(1-theta)*MV0)/theta;
      MA1 = ((MV1-MV0)/dt - (1-theta)*MA0)/theta;
    case 2
      MA1 = (MU1-MU0-dt*MV0-dt*dt*(1/2-beta)*MA0)/(dt*dt*beta);
      MV1 = MV0 + dt*(gamma*MA1 + (1-gamma)*MA0);
    case 3 
      MA1 = (gf_model_get(md, 'rhs'))';
      MA1 = MA1(Iu);
      
      if (Contact_option == 3)
        if (0) % New version  
        % Computation of U_el 
        U_el = U1;
        U_el2 = (gf_model_get(md, 'interpolation', 'u + pos_part(obstacle-u.(Normalized(Grad_obstacle)))*Normalized(Grad_obstacle)', mfu, GAMMAC))';
        I = gf_mesh_fem_get(mfu, 'basic dof on region', GAMMAC);
        U_el(I) = U_el2(I);
        gf_model_set(md, 'variable', 'u', U_el);
        gf_model_get(md, 'assembly', 'build_rhs');
        MA2 = (gf_model_get(md, 'rhs'))';
        MA2 = MA2(Iu);
        gf_model_set(md, 'variable', 'u', U1);
        gf_model_get(md, 'assembly', 'build_rhs');
        MA3=(MA1+MA2)/2;
        MA3(I) = MA1(I);
        MA1 = MA3;
        % MA2(I) = MA1(I);
        % MA1 = MA2;
        
        elseif(0) % Old approximately working version
        % Computation of U_el 
        U_el = U1;
        U_el2 = (gf_model_get(md, 'interpolation', 'u + pos_part(obstacle-u.(Normalized(Grad_obstacle)))*Normalized(Grad_obstacle)', mfu, GAMMAC))';
        I = gf_mesh_fem_get(mfu, 'basic dof on region', GAMMAC);
        U_el(I) = U_el2(I);
        % Additional term
        gf_model_set(md, 'variable', 'uel', U_el);
        % F = gf_asm('generic', mim_friction, 1, 'Test_u*Normalized(Grad_obstacle)', GAMMAC, md)
        % pause
        F = gf_asm('generic', mim_friction, 1, '(((clambda*(Div_uel-Div_u)*Id(meshdim)+2*cmu*Sym(Grad_uel-Grad_u))*Normal).(Normalized(Grad_obstacle)))*(Test_u.(Normalized(Grad_obstacle)))', GAMMAC, md);
        MA1 = MA1 + F(Iu);
        else % penalty on the second node, very experimental
            if (U1(1) <= 0 && U1(2) <= 0)
              MA1(2) = MA1(2) - NX*int_penalty_coeff*(U1(2));
            end
        end
      end
      MV1 = MV0 + dt*(gamma*MA1+(1-gamma)*MA0);
    case 4
      U1_2 = U1;
      U1 = 2*U1_2 - U0;
      V1_2 = (U1 - U0)/dt;
      MV1 = 2*M*V1_2 - MV0;
      MA1 = 0*MV1;
      MU1 = M*U1;
  end
      
  if (singular_mass == 1)
    V1 = cgs(M, MV1); % Pseudo inverse ...
  elseif (singular_mass == 2)
    VV1 = (B') \ MV1; VV2 = C*VV1; V1 = B\VV2; 
  else
    V1 = M \ MV1;
  end

  if (d == 1)
    % Compute the analytical solution
    X = [0:1/(NX*u_degree):1]';
    UA = []; GUA = [];
    for i=1:length(X) [ UA(i) GUA(i) ] = uAnalytic(X(i),t+dt); end % Should be improved
  
    % Compute the norm here with U1
    % disp ('**** Compute the L2/H1 errors ****');
    errL2(nit+2) = gf_compute(mfu,U1'-UA,'L2 norm',mim);
    errH1(nit+2) = gf_compute(mfu,U1'-UA,'H1 norm',mim);
  end
  
  E = (V1'*MV1 + U1'*K*U1)/2 - FF'*U1;
  fprintf('energy = %g\n', E);
  
  % Save relevant physical quantities
  Etime(nit+2) = E;
  Utime(nit+2) = U1(1);
  if (d == 1)
    if (singular_mass == 1)
      Vtime(nit+2) = V1(2);
    else
      Vtime(nit+2) = V1(1);
    end
    Stime(nit+2) = -(clambda + 2*cmu)*(NX/1)*(U1(1)-U1(2));
    SRtime(nit+2) = - ( U1'*K(:,1) + MA1(1) - FF(1) );
    % Contact stress: depend on Nitsche or not
    if (Contact_option == 0)
      CStime(nit+2) = lambda(1);
    else
      CStime(nit+2) = -gamma0_N*max(0,-U1(1)-(1/gamma0_N)*Stime(nit+2)); 
    end
  
    Rtime(nit+2) = (0.5/gamma0_N)*(1/NX)*( Stime(nit+2)^2 - CStime(nit+2)^2 );
    Eatime(nit+2) = Etime(nit+2) - theta_N * Rtime(nit+2);
    Vmtime(nit+2) = 0.5*(Vtime(nit+1)+Vtime(nit+2));
    CSmtime(nit+2) = 0.5*(CStime(nit+1)+CStime(nit+2)); 
  end
    
  nit = nit + 1;
  if (t >= tplot)
      if (d >= 2)
        VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
		    'u', 'clambda', 'cmu', mfvm);
      end
      if (d == 1)
        figure(1);
        X = [0:1/(NX*u_degree):1]';
        plot(zeros(1, Msize)-0.05, X+U1, '-b');
        hold on;
        plot(zeros(1, Msize)+0.05, U1+X, '-b');
        for i = 1:(NX*u_degree)+1
           plot([-0.05 0.05], (U1(i)+X(i))*[1 1], 'b'); 
        end
        hold off;
        axis([-0.4 0.4 0.0 1.5]);
      elseif (d == 2)
        figure(1);
        gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U1', ...
            'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8, 'disp_options', 'off');
        xlabel('x'); ylabel('y');
        % title('Deformed configuration (not really a small deformation of course ...)');
        % gf_colormap('chouette');
        gg = [ .7 .9 .4; .5 .9 .3;   .3 .8 .2;    .1 .7 .4;     .2 0.7 1.0000; .3 0.3 1.0000;
               1.0 .8 .1;  1.0 .6 .1;   1.0 .45 .1;   1.0 0.3 .1];
        r = reshape(repmat(gg',6,1),3,60)';
        colormap(r);
        colorbar;
        caxis([0 32]);
        axis([-25 25 -1 50]);
      else
        figure(1);
        c=[0.1;0;20]; x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
        % Whole boundary
        % sl2=gf_slice({'boundary',{'none'}}, m, 5);
        % Slice, 3 planes
        % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,5);
        % Slice, 2 planes
        sl2=gf_slice({'boundary',{'union',{'planar',+1,c,y},{'planar',+1,c,x}}},m,5);
        % Slice, 1 plane
        % sl2=gf_slice({'boundary',{'planar',+1,c,x}}, m, 5);

        P=gf_slice_get(sl2,'pts'); dP=gf_compute(mfu,U1','interpolate on',sl2);
        gf_slice_set(sl2, 'pts', P+dP);
        VMsl=gf_compute(mfvm,VM,'interpolate on',sl2);
        set(gcf,'renderer','zbuffer');
        h=gf_plot_slice(sl2,'mesh','off','mesh_slice_edges','off','data',VMsl);
        view(-80,-15); axis on; camlight; gf_colormap('chouette');
        % map=[1:-1/10:0]'*[1 1 1]; colormap(map); % for NB
    
        % gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
        %        gf_mesh_get(mfu,'outer faces'), 'deformation', U, ...
        %        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
        % view(-5,-10); camlight; colormap(map);
        xlabel('x'); ylabel('y'); zlabel('z');
        axis([-25 25 -25 25 -1 50]);
        caxis([0 20]);
        colorbar;
        % title('Sliced deformed configuration');
      end
    pause(0.1);
    if (make_movie)
      nim = nim + 1;
      F = getframe(gcf);
      mov = addframe(mov,F);
      Mov(:,nim) = getframe;
    end
    tplot = tplot + dt_plot;
  end
  

  U0 = U1;
  MU0 = MU1;
  MV0 = MV1;
  MA0 = MA1;
end
   

if (make_movie)
  mov = close(mov);
  mov = aviread('toto.avi');
  movie(mov);
end

% Final display

if (d == 1)
  figure(2), grid on;
  t=0:dt:T; tp = rem(t,3);
  u = (tp<=1) .* (0.5-0.5*tp) + (tp>1).*(tp<=2) .* 0 + (tp>2) .* (0.5*tp-1);
  plot(0:dt:T,Utime,'linewidth',1); hold on;
  plot(0:dt:T,u,'linewidth',1,'color',[1 0 1]); hold off;
  legend('numerical','analytical');
  xlabel('time'); ylabel('displacement u');

  figure(3), grid on;
  tp = rem(t,3);
  u = (tp<=1) .* (-0.5) + (tp>1).*(tp<=2) .* 0 + (tp>2) .* (0.5);
  plot(0:dt:T,Vtime,'linewidth',1); hold on;
  plot(0:dt:T,u,'linewidth',1,'color',[1 0 1]); hold off;
  legend('numerical','analytical');
  xlabel('time'); ylabel('velocity v');

  figure(4), grid on;
  %plot(0:dt:T,Stime,'linewidth',1);
  %plot(0:dt:T,SRtime,'linewidth',1,'color',[0 1 0]);
  tp = rem(t,3);
  u = (tp<=1) .* 0 + (tp>1).*(tp<=2) .* (-0.5) + (tp>2) .* 0;
  plot(0:dt:T,CStime,'linewidth',1,'color',[1 0 0]); hold on;
  plot(0:dt:T,u,'linewidth',1,'color',[1 0 1]); hold off;
  legend(...%'sigma_n num (D)','sigma_n num (VR)',
      'contact stress','analytical');
  xlabel('time'); ylabel('normal / contact stress \sigma');
  
  figure(5), grid on;
  plot(0:dt:T,Etime,'linewidth',1);
  xlabel('time'); ylabel('energy E');

  if (0) 
  figure(6), grid on;
  tp = rem(t,3);
  u = (tp<=1) .* (-0.5) + (tp>1).*(tp<=2) .* 0 + (tp>2) .* (0.5);
  plot(0:dt:T,Vmtime,'linewidth',1); hold on;
  plot(0:dt:T,u,'linewidth',1,'color',[1 0 1]); hold off;
  legend('numerical','analytical');
  xlabel('time'); ylabel('velocity v (midpoint)');

  figure(7), grid on;
  %plot(0:dt:T,Stime,'linewidth',1);
  %plot(0:dt:T,SRtime,'linewidth',1,'color',[0 1 0]);
  tp = rem(t,3);
  u = (tp<=1) .* 0 + (tp>1).*(tp<=2) .* (-0.5) + (tp>2) .* 0;
  plot(0:dt:T,CSmtime,'linewidth',1,'color',[1 0 0]); hold on;
  plot(0:dt:T,u,'linewidth',1,'color',[1 0 1]); hold off;
  legend(...%'sigma_n num (D)','sigma_n num (VR)',
      'contact stress','analytical');
  xlabel('time'); ylabel('normal / contact stress \sigma (midpoint)');

  figure(8), grid on;
  plot(0:dt:T,Rtime,'linewidth',1);
  xlabel('time'); ylabel('R(t)');

  figure(9), grid on;
  plot(0:dt:T,Eatime,'linewidth',1);
  xlabel('time'); ylabel('augmented energy E_\Theta');
  end
end


