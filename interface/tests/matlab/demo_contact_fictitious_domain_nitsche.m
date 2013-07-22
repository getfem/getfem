disp('Resolution of a contact problem in 2D with two elastics bodies');
disp('with a fictitious domain method and Nitsche s method');


clear all;
% gf_workspace('clear all');
NX=5;
ls_degree = 1; % pour 2 tous les matrices ne sont pas nulles
R=0.25;
dirichlet_val = 0;
gamma0 = 1;
theta = 0; %Pb theta = 1;
%N = 2 %la dimension


%definition of fictitious domain's mesh with quadrangles and order 1 of level-set


m=gf_mesh('regular simplices', -.5:(1/NX):.5, -.5:(1/NX):.5);
%m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
ls1=gf_levelset(m, ls_degree);
ls2=gf_levelset(m, ls_degree);
mf_ls1=gfObject(gf_levelset_get(ls1, 'mf'));
mf_ls2=gfObject(gf_levelset_get(ls2, 'mf'));
mfu=gfMeshFem(m,2);
set(mfu, 'fem', gf_fem('FEM_PK(2,1)'));
mfvm=gfMeshFem(m,1);
set(mfvm, 'fem', gf_fem('FEM_PK(2,0)'));

% set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));
mls1=gfMeshLevelSet(m);
mls2=gfMeshLevelSet(m);

%definition of Omega 1 (circle)
 
P=get(mf_ls1, 'basic dof nodes');
x = P(1,:); y = P(2,:);
ULS1=1000*ones(1,numel(x));
ULS1 = min(ULS1, sqrt(x.^2 + y.^2) - R);
gf_levelset_set(ls1, 'values', ULS1);

%definition of Omega 2 (rectangle)

P=get(mf_ls2, 'basic dof nodes');
x = P(1,:); y = P(2,:);
ULS2=1000*ones(1,numel(x));
yc = -0.25; xc=0; 
% R2=0.125;R1=0.5;
ULS2=min(ULS2,y-yc);
gf_levelset_set(ls2, 'values', ULS2); 

%figure

set(mls1, 'add', ls1);
set(mls1, 'adapt');

set(mls2, 'add', ls2);
set(mls2, 'adapt');

%Dirichlet's boundary
 
GAMMAC = 1; GAMMAD = 2;


border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(2, :) < -0.01));%normal dans la direction -e2
gf_mesh_set(m, 'region', GAMMAD, contact_boundary);
%gf_model_set(md,'add inialized data', 'dirichlet data',[dirichlet_val])


% figure 1 : plot figure

clf; gf_plot_mesh(get(mls1,'cut mesh'));
hold on; gf_plot_mesh(get(mls2,'cut mesh')); hold off;

%gf_plot_mesh(get(mls, 'cut_mesh'), 'curved', 'on');
%hold on; gf_plot(mf_ls,ULS);

hold on; gf_plot_mesh(m, 'regions', GAMMAD, 'convexes', 'on'); %plot de bord avec condition de type Dirichlet
title('boundary with Dirichlet condition in red');hold off;


%Finites elements' method on mls1 and mls2

mim_bound = gfMeshIm('levelset',mls1,'boundary', gf_integ('IM_TRIANGLE(5)'));
mim = gfMeshIm('levelset',mls1,'all', gf_integ('IM_TRIANGLE(5)')); 
set(mim, 'integ', 4);


%mfu=gfMeshFem(m,2); set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));
%mfdu=gfMeshFem(m,1); set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));

%Elastic model 

md=gf_model('real');
gf_model_set(md,'add fem variable', 'u1',mfu);
gf_model_set(md,'add fem variable', 'u2',mfu);
gf_model_set(md,'add initialized fem data', 'd1', mf_ls1, ULS1);
gf_model_set(md,'add initialized fem data', 'd2', mf_ls2, ULS2);
gf_model_set(md,'add initialized data', 'gamma0', gamma0);





 clambda = 1;           % Lame coefficient
 cmu = 1;               % Lame coefficient
 gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
 gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
 gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u1','clambda', 'cmu');
 gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u2','clambda', 'cmu');
  
 
 
  cpoints = [0, 0];   % constrained points for 2d
 cunitv  = [0.1, 0];   % corresponding constrained directions for 2d, mieux avec [0, 0.1]
 gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
 gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
 gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints', 'cunitv');

 
 
 
 
 
 gf_model_set(md, 'add initialized data', 'Fdata', [0 -0.0001]); % initiale [0 -1]
 gf_model_set(md, 'add source term brick', mim, 'u1', 'Fdata');
 Ddata = zeros(1, 2); u1_degree=2; u2_degree=2;%Dimension 2
 gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
 % gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u1', u1_degree, GAMMAD, 'Ddata'); %neccessaire?
 % gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u2', u2_degree, GAMMAD, 'Ddata'); %neccessaire?
 gf_model_set(md, 'add Dirichlet condition with simplification', 'u2', GAMMAD, 'Ddata'); %neccessaire?



 
 % marche pas

% cpoints2 = [0, 0.4];   % constrained points for 2d
%  cunitv2  = [0.1, 0.4];   % corresponding constrained directions for 2d, mieux avec [0, 0.1]
%  gf_model_set(md, 'add initialized data', 'cpoints2', cpoints2);
%  gf_model_set(md, 'add initialized data', 'cunitv2', cunitv2);
%  gf_model_set(md, 'add pointwise constraints with multipliers', 'u2', 'cpoints2', 'cunitv2');
 
 
 

gf_model_set(md,'add Nitsche fictitious domain contact brick', mim_bound, 'u1', 'u2', 'd1', 'd2', 'gamma0', theta); 


disp('solve');
niter= 10; solve=true;



gf_model_get(md, 'test tangent matrix term','u1','u1', 1e-6, niter, 2);

% gf_model_get(md, 'test tangent matrix', 1e-6, niter, 02);

pause;

niter= 1;

gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'very noisy');






figure(2);

U1 = gf_model_get(md, 'variable', 'u1');

sl1=gf_slice({'isovalues', -1, mf_ls1, ULS1, 0}, m, 5);
P1=gf_slice_get(sl1,'pts'); dP1=gf_compute(mfu,U1,'interpolate on',sl1);
gf_slice_set(sl1, 'pts', P1 + dP1);
VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
    		      'u1', 'clambda', 'cmu', mfvm);
VMsl1=gf_compute(mfvm,VM1,'interpolate on',sl1);



U2 = gf_model_get(md, 'variable', 'u2');

sl2=gf_slice({'isovalues', -1, mf_ls2, ULS2, 0}, m, 5);
P2=gf_slice_get(sl2,'pts'); dP2=gf_compute(mfu,U2,'interpolate on',sl2);
gf_slice_set(sl2, 'pts', P2+dP2);
VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
    		      'u2', 'clambda', 'cmu', mfvm);
VMsl2=gf_compute(mfvm,VM2,'interpolate on',sl2);


hold on;
gf_plot_slice(sl1,'mesh','on','mesh_slice_edges','off','data',VMsl1);
gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl2);
hold off;






