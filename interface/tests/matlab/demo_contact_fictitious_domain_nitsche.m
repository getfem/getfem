disp('Resolution of a contact problem in 2D with two elastics bodies');
disp('with a fictitious domain method and Nitsche s method');


clear all;
% gf_workspace('clear all');
NX=5;
ls_degree = 1;
R=0.25;
dirichlet_val = 0;
gamma0 = 1;
theta = 0;
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


clf; gf_plot_mesh(get(mls1,'cut mesh'));
hold on; gf_plot_mesh(get(mls2,'cut mesh')); hold off;

%gf_plot_mesh(get(mls, 'cut_mesh'), 'curved', 'on');
%hold on; gf_plot(mf_ls,ULS);

hold on; gf_plot_mesh(m, 'regions', GAMMAD, 'convexes', 'on'); %plot de bord avec condition de type Dirichlet
title('boundary with Dirichlet condition in red');
hold off;

%Finites elements' method on mls1 and mls2

mim_bound = gfMeshIm('levelset',mls1,'boundary', gf_integ('IM_TRIANGLE(5)'));
mim = gfMeshIm('levelset',mls1,'all', gf_integ('IM_TRIANGLE(5)')); 
set(mim, 'integ', 4);

%c'est suffisant?





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
  
 Ddata = zeros(1, 2); Ddata(2) = -5; u1_degree=2; u2_degree=2;%Dimension 2
 gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
 gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u1', u1_degree, GAMMAD, 'Ddata'); %neccessaire?
 gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u2', u2_degree, GAMMAD, 'Ddata'); %neccessaire?



 



gf_model_set(md,'add Nitsche fictitious domain contact brick', mim_bound, 'u1', 'u2', 'd1', 'd2', 'gamma0', theta); 








niter= 10; solve=true; 
disp('solve');

%pause;

gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter); % , 'very noisy');

% 
% % Solve the problem
% if (~solved)
%   gf_model_get(md, 'test tangent matrix', 1e-6, 10, 0.1);
%   gf_model_get(md, 'solve', 'max_res', 1E-9, 'very noisy', 'max_iter', niter); % ,  'lsearch', 'simplest'); % , 'with pseudo potential');
% end;
% 
% U = gf_model_get(md, 'variable', 'u');
% % lambda_n = gf_model_get(md, 'variable', 'lambda_n');
% VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca', ...
% 		  'u', 'clambda', 'cmu', mfvm);
%     
% 
% % set a custom colormap
% % r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r);
% 
% figure(2);
% if (d == 3)
%   c=[0.1;0;20]; x=[1;0;0]; y=[0;1;0]; z=[0;0;1];
%   % Whole boundary
%   % sl2=gf_slice({'boundary',{'none'}}, m, 5);
%   % Slice, 3 planes
%   % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,5);
%   % Slice, 2 planes
%   % sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y}}},m,5);
%   % Slice, 1 plane
%   sl2=gf_slice({'boundary',{'planar',+1,c,x}}, m, 5);
% 
%   P=gf_slice_get(sl2,'pts'); dP=gf_compute(mfu,U,'interpolate on',sl2);
%   gf_slice_set(sl2, 'pts', P+dP);
%   VMsl=gf_compute(mfvm,VM,'interpolate on',sl2);
%   set(gcf,'renderer','zbuffer');
%   h=gf_plot_slice(sl2,'mesh','on','mesh_slice_edges','off','data',VMsl);
%   view(-80,-15); axis on; camlight; gf_colormap('chouette');
%   % map=[1:-1/10:0]'*[1 1 1]; colormap(map); % for NB
%     
%   % gf_plot(mfvm, VM, 'mesh', 'off', 'cvlst', ...
%   %        gf_mesh_get(mfu,'outer faces'), 'deformation', U, ...
%   %        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
%   % view(-5,-10); camlight; colormap(map);
%   xlabel('x'); ylabel('y'); zlabel('z');
%   title('Sliced deformed configuration (not really a small deformation of course ...)');
% else
%   gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
%           'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 8);
%   xlabel('x'); ylabel('y');
%   title('Deformed configuration (not really a small deformation of course ...)');
% end;
% 
% colorbar;
% pause(0.1);





