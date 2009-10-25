
gf_workspace('clear all');

% parameters

ls_degree = 1;
lambda = 1;
mu = 1;
% rayon_trous = 0.05;
% threshold = 0.25; % NX = 10
threshold = 0.4;
NX = 20;
DX = 1./NX;
if (0.2 * NX ~= round(0.2 * NX))
  disp('Bad value for NX');
  return;
end;

% Mesh definition
m=gfMesh('cartesian', -.7:DX:.7, -.7:DX:.7);
% m=gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);

N = gf_mesh_get(m, 'dim');


% Find the region corresponding to Omega
pts =  gf_mesh_get(m, 'pts');
pid = find((abs(pts(1, :)) < 0.50000001) .* (abs(pts(2, :)) < 0.50000001));

OMEGA = 1;
cvid = gf_mesh_get(m, 'cvid from pid', pid);
gf_mesh_set(m, 'region', OMEGA, [cvid ; zeros(1, size(cvid, 2))]);

% Find the boundary GammaD and GammaN (to Omega)
pidleft = find((abs(pts(1, :)+0.5) < 1E-7) .* (abs(pts(2, :)) < 0.50000001));
fidleft = gf_mesh_get(m, 'faces from pid', pidleft);
normals = gf_mesh_get(m, 'normal of faces', fidleft);
fidleft=fidleft(:,find(abs(normals(1, :)+1) < 1E-3));
GAMMAD = 2;
gf_mesh_set(m, 'region', GAMMAD, fidleft);

pidright = find((abs(pts(1, :)-0.5) < 1E-7) .* (abs(pts(2, :)) < 0.50000001));
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
mimls=gfMeshIm(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));

mf_basic=gfMeshFem(m, 2);
gf_mesh_fem_set(mf_basic,'fem',gf_fem('FEM_QK(2,2)'), cvid);


% Definition of the initial level-set
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);
% ULS=gf_mesh_fem_get(mf_ls, 'eval', { 'x - 10' });
ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.2-sin(pi*5*x) .* sin(pi*5*y)' });
% ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.2-sin(pi*2.5*x) .* sin(pi*2.5*y)' });

F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '-2*(abs(y) < 0.025)'});

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
  ind = find(D > DX^N/10000);
  mf = gf_mesh_fem('partial', mf_basic, ind);

  % Problem definition
  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mf);
  gf_model_set(md, 'add initialized data', 'mu', [mu]);
  gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
               'lambda', 'mu', OMEGA);
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', 1, GAMMAD);
  gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
  gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);


  % Solving the direct problem
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');

%   subplot(1,2,1);
%   gf_plot(mf, U);
%   hold on;
%   [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0, 'pcolor', 'off');
%   set(h2{1},'LineWidth',2);
%   set(h2{1},'Color','green');
%   hold off;
%   colorbar;
  
  % Computation of shape derivative
  mf_g=gfMeshFem(m, 1);
  % gf_mesh_fem_set(mf_g,'fem', ...
  %     gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'), cvid);
  gf_mesh_fem_set(mf_g,'fem', ...
     gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,1),FEM_PK_DISCONTINUOUS(1,1))'), cvid);
  DU = gf_compute(mf, U, 'gradient', mf_g);
  
  EPSU = DU + permute(DU, [2 1 3]);
  if (N == 2)
    GF = (DU(1,1,:) + DU(2,2,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2)) - threshold;
    GF = reshape(GF, 1, size(GF, 3));
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % filtering the gradient
  M = gf_asm('mass matrix', mim, mf_g);
  D = abs(full(diag(M)));
  ind = find(D < DX^N/13);
  GF(ind) = GF(ind) * 0;

  
  % subplot(1,2,2);
  gf_plot(mf_g, GF);
  hold on;
  [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off');
  set(h2{1},'LineWidth',2);
  set(h2{1},'Color','green');
  hold off;
  colorbar;
  
  pause(1);
  disp('drawing done');
  
  % Evolution of the level-set. Computation of v
  
  mf_v=gfMeshFem(m, 1);
  gf_mesh_fem_set(mf_v,'fem', ...
      gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'));
  DLS = gf_compute(mf_ls, ULS, 'gradient', mf_v);
  NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
  GFF = gf_compute(mf_g, GF, 'interpolate on', mf_v) ./ NORMDLS;
  

  if (N == 2)
    V = DLS.*[GFF; GFF];
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % Normalisation de V (il ne sert a rien d'aller trop vite ou trop
  % lentement
  
  maxV = sqrt(max(sum(V.^2,1)));
  for i = 1:size(V,2)
      normv = sqrt(sum(V(:,i).^2));
      if (normv > maxV / 100)
        V(:,i) = V(:, i) / normv;
      end;
  end;
  
  
  % Evolution of the level-set. (perturbated) Hamilton-Jacobi equation.
  alpha = 0.001; dt = 0.01; NT = 100; ddt = dt / NT;
  M = gf_asm('mass matrix', mimls, mf_ls);
  K = gf_asm('laplacian', mimls, mf_ls, mf_ls, alpha*ones(gf_mesh_fem_get(mf_ls, 'nbdof'),1));
  B = gf_asm('volumic', ...
             'v=data(mdim(#2),#2);M(#1,#1)+=comp(Base(#1).Grad(#1).Base(#2))(:,:,i,j).v(i,j)', ...
             mimls, mf_ls, mf_v, V);
  
  A = M/ddt + B + K;
  for t = 0:ddt:dt
    ULS = (A \ ((M * ULS') / ddt))';
  end;
  
  
  % Renormalisation de la level-set
  % disp 'norm dls avant : ';
  % sqrt(sum(DLS.^2, 1))
  alpha = 0.005; dt = 0.001; NT = 100; ddt = dt / NT;
  K = gf_asm('laplacian', mimls, mf_ls, mf_ls, alpha*ones(gf_mesh_fem_get(mf_ls, 'nbdof'),1));

  for t = 0:ddt:dt
  
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_v);
    NORMDLS = sqrt(sum(DLS.^2, 1)) + 1e-12;
    SIGNLS = sign(gf_compute(mf_ls, ULS, 'interpolate on', mf_v));
    SISN = SIGNLS ./ NORMDLS;
  
    if (N == 2)
      W = DLS.*[SISN; SISN];
    else
      disp('Should be adapted ...');
      return;
    end;

    B = gf_asm('volumic', ...
               'v=data(mdim(#2),#2);M(#1,#1)+=comp(Base(#1).Grad(#1).Base(#2))(:,:,i,j).v(i,j)', ...
               mimls, mf_ls, mf_v, W);
    L = gf_asm('volumic', 'v=data(#2);V(#1)+=comp(Base(#1).Base(#2))(:,i).v(i)', ...
               mimls, mf_ls, mf_v, SIGNLS);

    A = M/dt + B + K;
    ULS = (A \ ((M * ULS') / dt + L))';
  end;
  % disp 'norm dls apres : ';
  % sqrt(sum(DLS.^2, 1))
  
  % gf_workspace('list static objects')
  % pause;
  
  nb = gf_workspace('nb static objects');
  disp('nb static stored object before pop:');
  disp(nb);
  
  gf_workspace('pop');
  
  disp('nb static stored object after pop:');
  disp(nb);
  
end;
  
  
  
  
  
  
