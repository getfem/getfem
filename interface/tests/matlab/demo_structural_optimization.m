% Matlab GetFEM++ interface
%
% Copyright (C) 2009 Alassane SY, Yves Renard.
%
% This file is a part of GetFEM++
%
% GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 2.1 of the License,  or
% (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
%
%  Shape optimization od a structure with both topological and shape gradient
%   (with a fictitious domain approach).
%
%  This program is used to check that matlab-getfem is working. This is
%  also a good example of use of GetFEM++.
%

clear all;
gf_workspace('clear all');

% parameters
ls_degree = 1;  % Degree of the level-set. Should be one for the moment.
k = 1;          % Degree of the finite element method for u
lambda = 1;     % Lame coefficient
mu = 1;         % Lame coefficient
hole_radius = 0.03;
threshold_shape = 1.;
threshold_topo = 0.5;
NX = 40;
DX = 1./NX;
DEBUG = 0;
if (DEBUG)
  NG = 4;
else
  NG = 2;
end;

% Mesh definition
% m=gfMesh('cartesian', -1:DX:1, -.5:DX:.5);
m=gfMesh('regular simplices', -1:(1/NX):1, -.5:(1/NX):.5);
N = gf_mesh_get(m, 'dim');
pts =  gf_mesh_get(m, 'pts');

% Find the boundary GammaD and GammaN
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
mimls=gfMeshIm(m, gf_integ('IM_TRIANGLE(4)'));
mf_basic=gfMeshFem(m, N);
gf_mesh_fem_set(mf_basic,'fem',gf_fem(sprintf('FEM_PK(2,%d)', k)));

% Definition of the initial level-set
% ULS=gf_mesh_fem_get(mf_ls, 'eval', { 'x - 2' });
% ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.8-sin(pi*5*x) .* sin(pi*5*y)' });
ULS=gf_mesh_fem_get(mf_ls, 'eval', { '-0.6-sin(pi*4*x) .* cos(pi*4*y)' });

% Level-set nodes
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);

% Force on the right part (Neumann condition)
F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '-4*(abs(y) < 0.0125)'});



while(1) % Optimization loop

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
  gf_model_set(md, 'add initialized data', 'penalty_param', [1E-7]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', 1, GAMMAD);
  gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
  gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);


  
  % Solving the direct problem
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');
  nbd = gf_mesh_fem_get(mf_ls, 'nbdof');
  
  % Computation of indicators (computation of K could be avoided)
  K = gf_asm('linear elasticity', mim, mf, mf_ls, lambda*ones(1, nbd), ...
             mu*ones(1, nbd));
  disp(sprintf('Elastic energy: %g', U*K*U'));
  S = gf_asm('volumic','V()+=comp()',mim);
  disp(sprintf('Remaining surface of matieral: %g', S));

%   subplot(1,2,1);
%   gf_plot(mf, U);
%  hold on;
%    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0, 'pcolor', 'off');
%    set(h2{1},'LineWidth',2);
%    set(h2{1},'Color','green');
%    hold off;
%    colorbar;
  
  mf_g=gfMeshFem(m, 1);
  gf_mesh_fem_set(mf_g,'fem', gf_fem(sprintf('FEM_PK_DISCONTINUOUS(2,%d)', k-1)));
  DU = gf_compute(mf, U, 'gradient', mf_g);
  EPSU = DU + permute(DU, [2 1 3]);
  
  % Computation of the shape derivative
  if (N == 2)
    GF1 = (DU(1,1,:) + DU(2,2,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
    GF = reshape(GF1, 1, size(GF1, 3)) - threshold_shape;
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % computation of the topological gradient
  if (N == 2)
    GT = -pi*( (lambda+2*mu) / (2*mu*(lambda+mu)) * 4*mu*GF1 + ...
        2*(lambda-mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:)).^2);
    GT = reshape(GT, 1, size(GT, 3)) + threshold_topo;
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % filtering the gradients
  M = gf_asm('mass matrix', mim, mf_g);
  D = abs(full(diag(M)));
  maxD = max(D);
  ind = find(D < maxD/3);
  GF(ind) = GF(ind) * 0;
  ind = find(D < maxD/1.3);
  GT(ind) = GT(ind) * 0 - 20;

  % Drawing the gradients
  subplot(NG,1,1);
  gf_plot(mf_g, GF, 'disp_options', 'off', 'refine', 1);
  title('Shape derivative');
  hold on;
  [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                  'off', 'disp_options', 'off', 'refine', 3);
  set(h2{1},'LineWidth',1);
  set(h2{1},'Color','green');
  hold off;
  colorbar;
  subplot(NG,1,2);
  % gf_plot(mf_g, GT, 'disp_options', 'off', 'disp_options', 'off', 'refine', 8);
  gf_plot(mf_ls, ULS, 'disp_options', 'off', 'refine', 1);
  title('Level set function');
  colorbar;  
  pause(0.1);

  [val, i] = max(GT);
  disp(sprintf('Max value of the topological gradient: %g', val));

  % Making a new hole (topological optimization)
  if (val > 0)
    point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);
    xc = point(1);
    yc = point(2);
    disp(sprintf('Making a new hole of coordinates (%g, %g)', xc, yc));
    R = hole_radius;
    ULS = max(ULS, (R^2 - (x - xc).^2 - (y - yc).^2)/(2*R));
  end;    
  
  % Evolution of the level-set thank to shape derivative. Computation of v
  DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
  NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
  GFF = GF ./ NORMDLS;
  
  if (N == 2)
    V = DLS.*[GFF; GFF];
  else
    disp('Should be adapted ...');
    return;
  end;
  
  % mf_vcont=gfMeshFem(m, N);
  % gf_mesh_fem_set(mf_vcont,'fem', gf_fem('FEM_PK(2,1)'));
  % Mcont = gf_asm('mass matrix', mimls, mf_vcont); % Could be computed only once.
  % Fdisc = gf_asm('volumic source', mimls, mf_vcont, mf_g, V);
  Mcont = gf_asm('mass matrix', mimls, mf_basic); % Could be computed only once.
  Fdisc = gf_asm('volumic source', mimls, mf_basic, mf_g, V);
  Vcont = Mcont \ Fdisc;

  % Level set convection
%  gf_compute(mf_ls, ULS, 'convect', mf_basic, Vcont, 0.0025, 5);
%  gf_compute(mf_ls, ULS, 'convect', mf_basic, Vcont, 0.0025, 5);
%  gf_compute(mf_ls, ULS, 'convect', mf_basic, Vcont, 0.0025, 5);
%  gf_compute(mf_ls, ULS, 'convect', mf_basic, Vcont, 0.0025, 5);
  gf_compute(mf_ls, ULS, 'convect', mf_basic, Vcont, 0.01, 20);

  if (DEBUG)
    disp('Drawing the level set function after convection');
    subplot(NG,1,3);
    gf_plot(mf_ls, ULS, 'disp_options', 'off');
    colorbar;
    hold on;
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                  'off', 'disp_options', 'off', 'refine', 3);
    set(h2{1},'LineWidth',1);
    set(h2{1},'Color','green');
    hold off;
    pause(0.01);
  end;
   
  % Re-initialization of the level-set
  dt = 0.05; NT = 10; ddt = dt / NT;
 
  for t = 0:ddt:dt
  
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    Fdisc = gf_asm('volumic source', mimls, mf_basic, mf_g, DLS);
    DLScont = Mcont \ Fdisc;
    NORMDLS = sqrt(sum(reshape(DLScont, N, size(DLScont, 1)/N).^2, 1)) + 1e-12;
    SULS = sign(ULS) ./ NORMDLS;
        
    if (N == 2)
      W = DLScont.*reshape([SULS; SULS], N*size(SULS, 2), 1);
    else
      disp('Should be adapted ...');
      return;
    end;
   
    gf_compute(mf_ls, ULS, 'convect', mf_basic, W, ddt, 1);
    ULS = ULS + ddt * sign(ULS);
  end;

  AA = sqrt(sum(DLS.^2, 1));
  disp(sprintf('Norm dls after: %g %g %g %g', AA(1), AA(2), AA(3), AA(4)));
  disp(sprintf('Norm dls after: max = %g, min = %g', max(AA), min(AA)));

  if (DEBUG)
    disp('Drawing the level set function after re-initialization');
    subplot(NG,1,4);
    gf_plot(mf_ls, ULS, 'disp_options', 'off');
    colorbar;
    hold on;
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                  'off', 'disp_options', 'off', 'refine', 3);
    set(h2{1},'LineWidth',1);
    set(h2{1},'Color','green');
    hold off;
  end;

  gf_workspace('pop');
  
  
end;
  
  
  
  
  
  
