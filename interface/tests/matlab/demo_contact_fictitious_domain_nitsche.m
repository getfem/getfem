disp('Résolution d un problème de contact en 2D de deux corps élastiques');
disp('par une méthode de domain fictif et Méthode de Nitsche');


clear all;
% gf_workspace('clear all');
NX=25;
ls_degree = 1;
R=0.25;
dirichlet_val = 0;
gamma0 = 1;
theta = 0;
%N = 2 %la dimension

%définition du maillage du Domaine fictif avec des quadrangles avec une ls d'ordre 1

m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
ls1=gf_levelset(m, ls_degree);
ls2=gf_levelset(m, ls_degree);
mf_ls1=gfObject(gf_levelset_get(ls1, 'mf'));
mf_ls2=gfObject(gf_levelset_get(ls2, 'mf'));
mfu=gfMeshFem(m,2);
set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));
mls1=gfMeshLevelSet(m);
mls2=gfMeshLevelSet(m);

%définition de Omega 1 (cercle)
 
P=get(mf_ls1, 'basic dof nodes');
x = P(1,:); y = P(2,:);
ULS1=1000*ones(1,numel(x));
ULS1 = min(ULS1,((x.^2 + y.^2 ) - R^2));
gf_levelset_set(ls1, 'values', ULS1);

%définition de Omega 2 (rectangle)

ULS2=1000*ones(1,numel(x));
yc = -0.375;xc=0; 
R2=0.125;R1=0.5;
ULS2=min(ULS2,(y-yc) - R2);
gf_levelset_set(ls2, 'values', ULS2); 

%figure

set(mls1, 'add', ls1);
set(mls1, 'adapt');

set(mls2, 'add', ls2);
set(mls2, 'adapt');

%bord de Dirichlet
 
GAMMAC = 1; GAMMAD = 2;


border = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);
contact_boundary=border(:, find(normals(2, :) < -0.01));%normal dans la direction -e2
gf_mesh_set(m, 'region', GAMMAD, contact_boundary);
%gf_model_set(md,'add inialized data', 'dirichlet data',[dirichlet_val])


clf; gf_plot_mesh(get(mls1,'cut mesh'));
hold on; gf_plot_mesh(get(mls2,'cut mesh'));

%gf_plot_mesh(get(mls, 'cut_mesh'), 'curved', 'on');
%hold on; gf_plot(mf_ls,ULS);

gf_plot_mesh(m, 'regions', GAMMAD); %plot de bord avec condition de type Dirichlet
title('bord avec condition de type Dirichlet en rouge');

%méthode d'élements finis sur mls1 et mls2

mim_bound = gfMeshIm('levelset',mls1,'boundary', gf_integ('IM_QUAD(5)'));
mim = gfMeshIm('levelset',mls1,'all', gf_integ('IM_QUAD(5)'));
set(mim, 'integ', 4);







%mfu=gfMeshFem(m,2); set(mfu, 'fem', gf_fem('FEM_QK(2,1)'));
%mfdu=gfMeshFem(m,1); set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));

%modèle élastique

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








niter= 1000; solve=true;
disp('solve');
pause;
gf_model_get(md, 'solve', 'max_res', 1E-9, 'max_iter', niter); % , 'very noisy');


















%Appelle à la fonction




