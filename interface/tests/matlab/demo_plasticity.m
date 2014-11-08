% Copyright (C) 2010-2012 Amandine Cottaz, Yves Renard.
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

clear all
gf_workspace('clear all');
clc

%%

% We try to compute a plasticity problem with a Von Mises criterion
% For convenience we consider an homogenous Dirichlet condition on the left
% of the domain and an easy computed Neumann Condition on the right

%%

with_hardening = 0;
bi_material = true;
test_tangent_matrix = 0;
do_plot = true;

% Initialize used data
L = 100;
H = 20;
NX = 50;
NY = 20;

f = [0 -330]';
t = [0 0.9032 1 1.1 1.3 1.5 1.7 1.74 1.7 1.5 1.3 1.1 1 0.9032 0.7 0.5 0.3 0.1 0];
if (with_hardening == 1)
  f = [0 -330]';
  t = [0 0.9032 1 1.1 1.3 1.5 1.7 1.74 1.7 1.5 1.3 1.1 1 0.9032 0.7 0.5 0.3 0.1 0];
  % t = [0 0.5 0.9 1.2 1.5 1.8 1.5 1.2 0.9 0.5 0.3 0];
end

% Create the mesh
% m = gfMesh('triangles grid', [0:(L/NX):L], [0:(H/NY):H]);
m = gfMesh('import','structured',sprintf('GT="GT_PK(2,1)";SIZES=[%d,%d];NOISED=0;NSUBDIV=[%d,%d];', L, H, NX, NY));

% Plotting
% gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% Define used MeshIm
mim=gfMeshIm(m);  set(mim, 'integ', gfInteg('IM_TRIANGLE(6)')); % Gauss methods on triangles

% Define used MeshFem
if (with_hardening == 1)
  mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,2)'));
else
  mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,1)'));
end
mf_data=gfMeshFem(m); set(mf_data, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
% mf_sigma=gfMeshFem(m,4); set(mf_sigma, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
mf_sigma=gfMeshFem(m,4); set(mf_sigma, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_vm = gfMeshFem(m); set(mf_vm, 'fem', gfFem('FEM_PK_DISCONTINUOUS(2,1)'));

% Find the border of the domain
P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6); % Retrieve index of points which x near to 0
pidright=find(abs(P(1,:) - L)<1e-6); % Retrieve index of points which x near to L
fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);
set(m,'boundary',1,fleft); % for Dirichlet condition
set(m,'boundary',2,fright); % for Neumann condition

% Decomposed the mesh into 2 regions with different values of Lamé coeff
if (bi_material) separation = H/2; else separation = 0; end
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
% Definition of hardening parameter
H = mu(1)/5;

% Create the model
md = gfModel('real');

% Declare that u is the unknown of the system on mf_u
% 2 is the number of version of the data stored, for the time integration scheme 
set(md, 'add fem variable', 'u', mf_u, 2);

% Declare that lambda is a data of the system on mf_data
set(md, 'add initialized fem data', 'lambda', mf_data, lambda);

% Declare that mu is a data of the system on mf_data
set(md, 'add initialized fem data', 'mu', mf_data, mu);

% Declare that von_mises_threshold is a data of the system on mf_data
set(md, 'add initialized fem data', 'von_mises_threshold', mf_data, von_mises_threshold);


if (with_hardening)
  N = gf_mesh_get(m, 'dim');
  gf_model_set(md, 'add fem data', 'Previous_u', mf_u);
  mim_data = gf_mesh_im_data(mim, -1, [N, N]);
  gf_model_set(md, 'add im data', 'sigma', mim_data);
   
  set(md, 'add initialized data', 'H', [H]);

  Is = 'Reshape(Id(meshdim*meshdim),meshdim,meshdim,meshdim,meshdim)';
  IxI = 'Id(meshdim)@Id(meshdim)';
  coeff_long = '((lambda)*(H))/((2*(mu)+(H))*(meshdim*(lambda)+2*(mu)+(H)))';
  B_inv = sprintf('((2*(mu)/(2*(mu)+(H)))*(%s) + (%s)*(%s))', Is, coeff_long, IxI);
  B = sprintf('((1+(H)/(2*(mu)))*(%s) + (-(lambda)*(H)/(2*(mu)*(meshdim*(lambda)+2*(mu))))*(%s))', Is, IxI);
  ApH = sprintf('((2*(mu)+(H))*(%s) + (lambda)*(%s))', Is, IxI);
  Enp1 = '((Grad_u+Grad_u'')/2)';
  En = '((Grad_Previous_u+Grad_Previous_u'')/2)';
  expr_sigma = strcat('(', B_inv, '*(Von_Mises_projection(((H)*', Enp1, ')+(', ApH, '*(',Enp1,'-',En,')) + (', B, '*sigma), von_mises_threshold) + H*', Enp1, '))');
  
  gf_model_set(md, 'add nonlinear generic assembly brick', mim, strcat(expr_sigma, ':Grad_Test_u'));
  % gf_model_set(md, 'add finite strain elasticity brick', mim, 'u', 'SaintVenant Kirchhoff', '[lambda; mu]');
else
    
  % Declare that sigma is a data of the system on mf_sigma
  set(md, 'add fem data', 'sigma', mf_sigma);
  % Add plasticity brick on u
  set(md, 'add elastoplasticity brick', mim, 'VM', 'u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
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
    get(md, 'solve','lsolver', 'superlu', 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'max_iter', 100, 'max_res', 1e-6);
    % get(md, 'solve', 'noisy', 'max_iter', 80);

    % Retrieve the solution U
    U = get(md, 'variable', 'u', 0);
    
    % Compute new plasticity constraints used to compute 
    % the Von Mises or Tresca stress
    if (with_hardening)
      sigma = gf_model_get(md, 'interpolation', expr_sigma, mim_data);
      gf_model_set(md, 'variable', 'sigma', sigma);
      gf_model_set(md, 'variable', 'Previous_u', U);
      
      M = gf_asm('mass matrix', mim, mf_vm);
      L = gf_asm('generic', mim, 1, 'sqrt(3/2)*Norm(Deviator(sigma))*Test_vm', -1, 'sigma', 0, mim_data, sigma, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1));
      VM = (M \ L)';
      coeff1='-lambda/(2*mu*(meshdim*lambda+2*mu))';
      coeff2='1/(2*mu)';
      Ainv=sprintf('(%s)*(%s) + (%s)*(%s)', coeff1, IxI, coeff2, Is);
      Ep = sprintf('(Grad_u+Grad_u'')/2 - (%s)*sigma', Ainv);
      L = gf_asm('generic', mim, 1, sprintf('Norm(%s)*Test_vm', Ep), -1, 'sigma', 0, mim_data, sigma, 'u', 0, mf_u, U, 'vm', 1, mf_vm, zeros(gf_mesh_fem_get(mf_vm, 'nbdof'),1), 'mu', 0, mf_data, mu, 'lambda', 0, mf_data, lambda);
      plast = (M \ L)';
    else
      get(md, 'elastoplasticity next iter', mim, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
      plast = get(md, 'compute plastic part', mim, mf_vm, 'u', 'VM', 'lambda', 'mu', 'von_mises_threshold', 'sigma');
      % Compute Von Mises or Tresca stress
      VM = get(md, 'compute elastoplasticity Von Mises or Tresca', 'sigma', mf_vm, 'Von Mises');
    end
  
    if (do_plot)
      figure(2)
      subplot(2,1,1);
      gf_plot(mf_vm,VM, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0); % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Von Mises criterion for t = ', num2str(step)]);

      subplot(2,1,2);
      gf_plot(mf_vm,plast, 'deformation',U,'deformation_mf',mf_u,'refine', 4, 'deformation_scale',1, 'disp_options', 0);  % 'deformed_mesh', 'on')
      colorbar;
      axis([-20 120 -20 40]);
      % caxis([0 10000]);
      n = t(step);
      title(['Plastification for t = ', num2str(step)]);
    
      pause(0.1);
    end

end;








