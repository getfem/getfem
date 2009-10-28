
clear all;
gf_workspace('clear all');

% parameters

ls_degree = 1;
lambda = 1;
mu = 1;
rayon_trous = 0.025;
% threshold = 0.25; % NX = 10
threshold = 1.25;
NX = 40;
DX = 1./NX;
if (0.2 * NX ~= round(0.2 * NX))
  disp('Bad value for NX');
  return;
end;

% Mesh definition
% m=gfMesh('cartesian', -1:DX:1, -.5:DX:.5);
m=gfMesh('regular simplices', -1:(1/NX):1, -.5:(1/NX):.5);
N = gf_mesh_get(m, 'dim');
pts =  gf_mesh_get(m, 'pts');

% Find the boundary GammaD and GammaN (to Omega)
pidleft = find((abs(pts(1, :)+1.0) < 1E-7) .* (abs(pts(2, :)) < 0.50000001));
fidleft = gf_mesh_get(m, 'faces from pid', pidleft);
normals = gf_mesh_get(m, 'normal of faces', fidleft);
fidleft=fidleft(:,find(abs(normals(1, :)+1) < 1E-3));
GAMMAD = 2;
gf_mesh_set(m, 'region', GAMMAD, fidleft);

pidright = find((abs(pts(1, :)-1.0) < 1E-7) .* (abs(pts(2, :)) < 0.50000001));
fidright = gf_mesh_get(m, 'faces from pid', pidright);
normals = gf_mesh_get(m, 'normal of faces', fidright);
fidright=fidright(:,find(abs(normals(1, :)-1) < 1E-3));
GAMMAN = 3;
gf_mesh_set(m, 'region', GAMMAN, fidright);


% Definition of the finite element methods
ls=gfLevelSet(m, ls_degree);
mls=gfMeshLevelSet(m);
set(mls, 'add', ls);
mf_ls=gfObject(get(ls, 'mf'));
% mimls=gfMeshIm(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
mimls=gfMeshIm(m, gf_integ('IM_TRIANGLE(4)'));

mf_basic=gfMeshFem(m, N);
% gf_mesh_fem_set(mf_basic,'fem',gf_fem('FEM_QK(2,2)'));
gf_mesh_fem_set(mf_basic,'fem',gf_fem('FEM_PK(2,1)'));

% Definition of the initial level-set
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);
ULS=gf_mesh_fem_get(mf_ls, 'eval', { 'x - 10' });
% ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.8-sin(pi*5*x) .* sin(pi*5*y)' });
%ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.9-sin(pi*4*x) .* cos(pi*4*y)' });

F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '-4*(abs(y) < 0.0125)'});

while(1) 

  gf_workspace('push');

  set(ls, 'values', ULS);
  disp('Adapting the mesh');
  set(mls, 'adapt');
  disp('Mesh adapted');
  
  mim = gfMeshIm('levelset',mls,'inside', gf_integ('IM_TRIANGLE(6)'));
  set(mim, 'integ', 4);

  M = gf_asm('mass matrix', mim, mf_basic);
  D = abs(full(diag(M)));
  ind = find(D > DX^N/10000000);
  mf = gf_mesh_fem('partial', mf_basic, ind);

  % Problem definition
  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mf);
  gf_model_set(md, 'add initialized data', 'mu', [mu]);
  gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
               'lambda', 'mu');
  gf_model_set(md, 'add initialized data', 'penalty_param', [1E-8]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', 1, GAMMAD);
  gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
  gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);

  subplot(2,1,1);
  [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off', 'disp_options', 'off');
  set(h2{1},'LineWidth',1);
  set(h2{1},'Color','green');
  
  % Solving the direct problem
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');
  nbd = gf_mesh_fem_get(mf_ls, 'nbdof');
  
%   K = gf_asm('linear elasticity', mim, mf, mf_ls, lambda*ones(1, nbd), mu*ones(1, nbd));
%   disp('Elastic energy');
%   disp(U*K*U');

  %   subplot(1,2,1);
%   gf_plot(mf, U);
%  hold on;
%    [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0, 'pcolor', 'off');
%    set(h2{1},'LineWidth',2);
%    set(h2{1},'Color','green');
%    hold off;
%    colorbar;
  
  mf_g=gfMeshFem(m, 1);
  % gf_mesh_fem_set(mf_g,'fem', ...
  %     gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'));
  % gf_mesh_fem_set(mf_g,'fem', ...
  %   gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,1),FEM_PK_DISCONTINUOUS(1,1))'));
  gf_mesh_fem_set(mf_g,'fem', gf_fem('FEM_PK(2,0)'));
  DU = gf_compute(mf, U, 'gradient', mf_g);
  EPSU = DU + permute(DU, [2 1 3]);
  
  % Computation of shape derivative
  if (N == 2)
    GF1 = (DU(1,1,:) + DU(2,2,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
    GF = reshape(GF1, 1, size(GF1, 3)) - threshold;
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % computation of the topological gradient
  if (N == 2)
    GT = -pi*( (lambda+2*mu) / (2*mu*(lambda+mu)) * 4*mu*GF1 + ...
        2*(lambda-mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:)).^2);
    GT = reshape(GT, 1, size(GT, 3));
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % filtering the gradients
  M = gf_asm('mass matrix', mim, mf_g);
  D = abs(full(diag(M)));
  ind = find(D < DX^N/7);
  GF(ind) = GF(ind) * 0;
  ind = find(D < DX^N/5);
  GT(ind) = GT(ind) * 0 - 20;

  
  
  % LS = gf_compute(mf_ls, ULS, 'interpolate on', mf_g);
  % GT = GT .* (1.-sign(LS))/2;
  [val, i] = max(GT);

  if (val >= -3*threshold)
    point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);
    xc = point(1);
    yc = point(2);
    R = rayon_trous;
    ULS = max(ULS, R^2 - ((x - xc).^2 + (y - yc).^2));
  end;    
  


  
  subplot(2,1,1);
  gf_plot(mf_g, GF, 'disp_options', 'off');
  hold on;
  [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off', 'disp_options', 'off');
  set(h2{1},'LineWidth',1);
  set(h2{1},'Color','green');
  % gf_plot(mf, U);
  hold off;
  colorbar;
  subplot(2,1,2);
  gf_plot(mf_g, GT, 'disp_options', 'off', 'disp_options', 'off');
  colorbar;
  
  pause(0.1);
  disp('drawing done');
  
  % Evolution of the level-set. Computation of v
  
  DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
  NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
  GFF = GF ./ NORMDLS;
  
  if (N == 2)
    V = DLS.*[GFF; GFF];
  else
    disp('Should be adapted ...');
    return;
  end;
  
  disp('convect début');
  mf_vcont=gfMeshFem(m, N);
  gf_mesh_fem_set(mf_vcont,'fem', gf_fem('FEM_PK(2,1)'));
  Mcont = gf_asm('mass matrix', mimls, mf_vcont); % Could be computed only once.
  Fdisc = gf_asm('volumic source', mimls, mf_vcont, mf_g, V);
  Vcont = Mcont \ Fdisc;
  gf_compute(mf_ls, ULS, 'convect', mf_vcont, Vcont, 0.02, 8);
  disp('convect fin');
   
  % Re-initialization of the level-set
  dt = 0.1; NT = 15; ddt = dt / NT;
 
  for t = 0:ddt:dt
  
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    Fdisc = gf_asm('volumic source', mimls, mf_vcont, mf_g, DLS);
    DLScont = Mcont \ Fdisc;
    NORMDLS = sqrt(sum(reshape(DLScont, N, size(DLScont, 1)/N).^2, 1)) + 1e-12;
    SULS = sign(ULS) ./ NORMDLS;
        
    if (N == 2)
      W = DLScont.*reshape([SULS; SULS], N*size(SULS, 2), 1);
    else
      disp('Should be adapted ...');
      return;
    end;
   
    gf_compute(mf_ls, ULS, 'convect', mf_vcont, W, ddt, 1);
    ULS = ULS + ddt * sign(ULS);
  end;

  disp 'norm dls apres : '; AA = sqrt(sum(DLS.^2, 1)); AA(1:4)  
  
  gf_workspace('pop');
  
  
end;
  
  
  
  
  
  
