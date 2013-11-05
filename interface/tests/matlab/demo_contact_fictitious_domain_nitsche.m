disp('Resolution of a contact problem in 2D with two elastics bodies');
disp('with a fictitious domain method and Nitsche s method');


clear all;
% gf_workspace('clear all');
NX=100;
ls_degree = 2; % pour 2 tous les matrices ne sont pas nulles
R=0.25;
dirichlet_val = 0;
gamma0 = 1/10;
theta = -1;
%N = 2 %la dimension
penalty_parameter = 10E-4;
vertical_force = -0.1;


%nxy = [100];%  18 26 34 42 50 58 66 74 82 90 98
 %zz = 1;
%for zz=1:1:size(nxy,2)
%NX=nxy(zz);

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

clf; %gf_plot_mesh(get(mls1,'cut mesh'));
%hold on; gf_plot_mesh(get(mls2,'cut mesh')); hold off;

%gf_plot_mesh(get(mls, 'cut_mesh'), 'curved', 'on');
%hold on; gf_plot(mf_ls,ULS);

%hold on; gf_plot_mesh(m, 'regions', GAMMAD, 'convexes', 'on'); %plot de bord avec condition de type Dirichlet
%title('boundary with Dirichlet condition in red');hold off;


%Finites elements' method on mls1 and mls2

mim_bound = gfMeshIm('levelset',mls1,'boundary', gf_integ('IM_TRIANGLE(5)'));
mim = gfMeshIm('levelset',mls1,'all', gf_integ('IM_TRIANGLE(5)')); 
mim1 = gfMeshIm('levelset', mls1, 'inside', gf_integ('IM_TRIANGLE(5)')); 
mim2 = gfMeshIm('levelset', mls2, 'inside', gf_integ('IM_TRIANGLE(5)')); 
set(mim, 'integ', 4);
set(mim1, 'integ', 4);
set(mim2, 'integ', 4);


dof_out = get(mfu, 'dof from im', mim1);
cv_out = get(mim1, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu1 = gfMeshFem('partial', mfu, dof_out, cv_in);

dof_out = get(mfu, 'dof from im', mim2);
cv_out = get(mim2, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu2 = gfMeshFem('partial', mfu, dof_out, cv_in);

%mfu=gfMeshFem(m,2); set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));
%mfdu=gfMeshFem(m,1); set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));

%Elastic model 

md=gf_model('real');
gf_model_set(md,'add fem variable', 'u1', mfu1);
gf_model_set(md,'add fem variable', 'u2', mfu2);
gf_model_set(md,'add initialized fem data', 'd1', mf_ls1, ULS1);
gf_model_set(md,'add initialized fem data', 'd2', mf_ls2, ULS2);
gf_model_set(md,'add initialized data', 'gamma0', gamma0);





clambda = 1;           % Lame coefficient
cmu = 1;               % Lame coefficient
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim1, 'u1','clambda', 'cmu');
gf_model_set(md, 'add isotropic linearized elasticity brick', mim2, 'u2','clambda', 'cmu');
  
 
gf_model_set(md, 'add initialized data', 'Fdata', [0 vertical_force]); % initiale [0 -1]
gf_model_set(md, 'add source term brick', mim1, 'u1', 'Fdata');
Ddata = zeros(1, 2); u1_degree=2; u2_degree=2; %Dimension 2
gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
% gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u1', u1_degree, GAMMAD, 'Ddata'); %neccessaire?
% gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u2', u2_degree, GAMMAD, 'Ddata'); %neccessaire?
gf_model_set(md, 'add Dirichlet condition with simplification', 'u2', GAMMAD, 'Ddata'); %neccessaire?
 
  
cpoints = [0, 0,   0, 0.1]; % constrained points for 2d
cunitv  = [1, 0,   1, 0];   % corresponding constrained directions for 2d, mieux avec [0, 0.1]
gf_model_set(md, 'add initialized data', 'cpoints', cpoints);
gf_model_set(md, 'add initialized data', 'cunitv', cunitv);
% gf_model_set(md, 'add pointwise constraints with multipliers', 'u1', 'cpoints', 'cunitv');
% gf_model_set(md, 'add pointwise constraints with penalization', 'u1', 100, 'cpoints', 'cunitv');
gf_model_set(md, 'add initialized data', 'penalty_param1', [penalty_parameter]);
indmass = gf_model_set(md, 'add mass brick', mim1, 'u1', 'penalty_param1');
gf_model_set(md, 'add initialized data', 'penalty_param2', [penalty_parameter]);
indmass = gf_model_set(md, 'add mass brick', mim2, 'u2', 'penalty_param2');

gf_model_set(md,'add Nitsche fictitious domain contact brick', mim_bound, 'u1', 'u2', 'd1', 'd2', 'gamma0', theta); 


disp('solve');
%niter= 50; solve=true;

% gf_model_get(md, 'test tangent matrix term', 'u1', 'u2', 1e-6, niter, 10.0);

%gf_model_get(md, 'test tangent matrix', 1e-6, niter, 10);

% pause;

niter= 100;

gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter, 'noisy');


%figure(2);

U1 = gf_model_get(md, 'variable', 'u1');

sl1=gf_slice({'isovalues', -1, mf_ls1, ULS1, 0}, m, 5);
P1=gf_slice_get(sl1,'pts'); dP1=gf_compute(mfu1,U1,'interpolate on',sl1);
gf_slice_set(sl1, 'pts', P1 + dP1);
VM1 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
    		      'u1', 'clambda', 'cmu', mfvm);
VMsl1=gf_compute(mfvm,VM1,'interpolate on',sl1);



U2 = gf_model_get(md, 'variable', 'u2');

sl2=gf_slice({'isovalues', -1, mf_ls2, ULS2, 0}, m, 5);
P2=gf_slice_get(sl2,'pts'); dP2=gf_compute(mfu2,U2,'interpolate on',sl2);
gf_slice_set(sl2, 'pts', P2+dP2);
VM2 = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
    		      'u2', 'clambda', 'cmu', mfvm);
VMsl2=gf_compute(mfvm,VM2,'interpolate on',sl2);


hold on;
gf_plot_slice(sl1,'mesh','on','mesh_slice_edges','off','data',VMsl1);
gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl2);
hold off;




% save the reference solution
% 
% V1e = VM1;
% V2e = VM2;
% gf_mesh_fem_get(mfvm, 'save', 'sol_ref_mesh_fem','with_mesh');
% save sol_de_reference V1e;
% save sol_de_reference V2e;
%  
% % errors in L2 and H1  
% 
% Discrétisation de réf N > n/4
% % 
meshref = gf_mesh('load', 'sol_ref_mesh_fem');
mfref = gf_mesh_fem('load', 'sol_ref_mesh_fem',meshref)
Vref1 = load('sol_de_reference', 'V1e');
Vref2 = load('sol_de_reference', 'V2e');
% mim2 = gf_mesh_im(meshref, 2);
% Ve1 = gf_compute(mfu1, U1, 'interpolate_on', mfref);
% Ve2 = gf_compute(mfu2, U2, 'interpolate_on', mfref);
% n_tot = gf_compute(mfref, Ve1-Vref1.VMsl1, 'L2
% norm',mim2)+gf_compute(mfref, Ve2-Vref2.VMsl2, 'L2 norm',  mim2);% plutot sur chaque corps
% n_ref = gf_compute(mfref, Vref1.VMsl1, 'L2 norm',  mim2)+gf_compute(mfref, Vref2.VMsl2, 'L2 norm',  mim2);
% m_tot = gf_compute(mfref, Ve1-Vref1.VMsl1, 'H1 norm',  mim2)+gf_compute(mfref, Ve2-Vref2.VMsl2, 'H1 norm',  mim2);
% m_ref = gf_compute(mfref, Vref1.VMsl1, 'H1 norm',  mim2)+gf_compute(mfref, Vref2.VMsl2, 'H1 norm',  mim2);
% n = 100*n_tot/n_ref
% m = 100*m_tot/m_ref
% nddl(zz)= gf_model_get(md,'nbdof')
% Y1(zz)=n % tab contenant les err en norme L2
% Y2(zz)=m % tab contenant les err en norme H1
% end
% X12=1./nxy; % [1/10 1/16 1/22 1/28 1/34 1/40];
% 
% % figure(2);
% loglog(X12,Y1);
% polyfit(log(X12),log(Y1),1)
% 
% % figure(3);
% loglog(X12,Y2);
% polyfit(log(X12),log(Y2),1)




