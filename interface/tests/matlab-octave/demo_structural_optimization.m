% Copyright (C) 2009-2020 Alassane SY, Yves Renard.
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
%  Shape optimization of a structure with a coupling between topological and
%  shape gradient (with a fictitious domain approach).
%
%  This program is used to check that matlab-getfem is working. This is
%  also a good example of use of GetFEM.
%





clear all;
gf_workspace('clear all');

% parameters

TEST_CASE = 1 % 0 : 2D, initial holes, shape gradient only
              % 1 : 2D, no initial hole, coupling with topological gradient
              % 2 : 3D, initial holes, shape gradient only
              % 3 : 3D, no initial hole, coupling with topological gradient

switch TEST_CASE
  case 0,
    N = 2;
    initial_holes = 1;
  case 1,
    N = 2;
    initial_holes = 0;
  case 2,
    N = 3;
    initial_holes = 1;
  case 3,
    N = 3;
    initial_holes = 0;
end

k = 1;                 % Degree of the finite element method for u
lambda = 1;            % Lame coefficient
mu = 1;                % Lame coefficient

if (N == 2)
  NY = 80              % Number of elements in y direction
  level_set_rate = 0.4 / NY
  reinitialisation_time = 0.005;
  threshold_shape = 0.90;
  if (TEST_CASE == 1)
    threshold_topo  = 1.3;
  else
    threshold_topo  = 0;
  end
  penalty_param = 1E-6;
  nbiter = 400;
  NBDRAW = 20;            % Draw solution each NBDRAW iterations
else
  NY = 30
  level_set_rate = 0.025 / NY;
  reinitialisation_time = 0.0035;
  threshold_shape = 28;
  if (TEST_CASE == 3)
    threshold_topo = 30;
  else
    threshold_topo = 0;
  end
  penalty_param = 1E-6;
  nbiter = 600;
  NBDRAW = 5;            % Draw solution each NBDRAW iterations
end;
hole_radius = max(0.03, 2./NY);  % Hole radius for topological optimization
ls_degree = 1;         % Degree of the level-set. Should be one for the moment.
DEBUG = 0;
if (DEBUG)
  NG = 3;
else
  NG = 2;
end;

% Mesh definition
% m=gfMesh('cartesian', -1:(1/NY):1, -.5:(1/NY):.5);
if (N == 2)
  m=gfMesh('regular simplices', -1:(1/NY):1, -.5:(1/NY):.5);
else
  m=gfMesh('regular simplices', -1:(1/NY):1, -.5:(1/NY):.5, -.5:(1/NY):.5);
end;
pts =  gf_mesh_get(m, 'pts');

% Find the boundary GammaD and GammaN
pidleft = find((abs(pts(1, :)+1.0) < 1E-7));
fidleft = gf_mesh_get(m, 'faces from pid', pidleft);
normals = gf_mesh_get(m, 'normal of faces', fidleft);
fidleft=fidleft(:,find(abs(normals(1, :)+1) < 1E-3));
GAMMAD = 2;
gf_mesh_set(m, 'region', GAMMAD, fidleft);

pidright = find((abs(pts(1, :)-1.0) < 1E-7));
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
if (N == 2)
  mimls=gfMeshIm(m, gf_integ('IM_TRIANGLE(4)'));
else
  mimls=gfMeshIm(m, gf_integ('IM_TETRAHEDRON(5)'));   
end;
mf_basic=gfMeshFem(m, N);
gf_mesh_fem_set(mf_basic,'fem',gf_fem(sprintf('FEM_PK(%d,%d)', N, k)));
mf_g=gfMeshFem(m, 1);
gf_mesh_fem_set(mf_g,'fem', gf_fem(sprintf('FEM_PK_DISCONTINUOUS(%d,%d)', N, k-1)));
mf_cont=gfMeshFem(m, N);
gf_mesh_fem_set(mf_cont,'fem', gf_fem(sprintf('FEM_PK(%d,%d)', N, ls_degree)));
disp(sprintf('There is %d elasticity dofs', gf_mesh_fem_get(mf_basic, 'nbdof')));

disp('Computation of mass matrices and preconditioners');
ti = cputime;
Mcont = gf_asm('mass matrix', mimls, mf_cont);
RMcont = cholinc(Mcont, '0');
Mcontls = gf_asm('mass matrix', mimls, mf_ls);
RMcontls = cholinc(Mcontls, '0');
disp(sprintf('Computation done in %g seconds', cputime-ti));

% Definition of the initial level-set
if (initial_holes)
  if (N == 2)
    ULS=gf_mesh_fem_get(mf_ls, 'eval', ...
        { '(-0.6-sin(pi*4*x).*cos(pi*4*y))/(4*pi)' });
  else
    ULS=gf_mesh_fem_get(mf_ls, 'eval', ...
        { '(-0.6-sin(pi*4*x).*cos(pi*4*y).*cos(pi*4*z))/(4*pi)' });
  end
else
  ULS=gf_mesh_fem_get(mf_ls, 'eval', { 'x - 2' });
end;

% Level-set nodes
P=get(mf_ls, 'basic dof nodes');

% Force on the right part (Neumann condition)
if (N == 2)
  F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '-(abs(y) < 0.05)'});
else
  F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '0', '-20*(abs(y) < 0.05).*(abs(z) < 0.05)'});
end;

% Model definition
set(ls, 'values', ULS);
disp('Adapting the mesh');
set(mls, 'adapt');
if (N == 2)
  mim = gfMeshIm('levelset',mls,'inside', gf_integ('IM_TRIANGLE(6)'));
else
  mim = gfMeshIm('levelset',mls,'inside', gf_integ('IM_TETRAHEDRON(6)'));
end;
disp('Mesh adapted');
set(mim, 'integ', 4);
disp('Integration methods adapted');
mf = gf_mesh_fem('partial', mf_basic, 1:gf_mesh_fem_get(mf_basic, 'nbdof'));

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add initialized data', 'mu', [mu]);
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
             'lambda', 'mu');
% gf_model_set(md, 'add Dirichlet condition with multipliers', ...
%             mim, 'u', 1, GAMMAD);
gf_model_set(md,'add Dirichlet condition with penalization', mim, ...
      'u', 1E5, GAMMAD);
gf_model_set(md, 'add initialized data', 'penalty_param', [penalty_param]);
gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);



% Optimization loop
for niter = 1:nbiter
  ti = cputime;
  gf_workspace('push');

  if (niter > 1)
    set(ls, 'values', ULS);
    disp('Adapting the mesh');
    set(mls, 'adapt');
    disp('Mesh adapted');
    set(mim, 'adapt');
    disp('Integration methods adapted');
  end;

  M = gf_asm('mass matrix', mim, mf_basic);
  D = abs(full(diag(M)));
  ind = find(D > (1/NY)^N/10000000);
  gf_mesh_fem_set(mf, 'set partial', ind); 
  % mf = gf_mesh_fem('partial', mf_basic, ind);

  % Solving the direct problem
  disp('solving the direct problem');
  gf_model_get(md, 'solve', 'max_res', 1E-7); % ,'noisy');
  U = gf_model_get(md, 'variable', 'u');
  nbd = gf_mesh_fem_get(mf_ls, 'nbdof');
  
  % Computation of indicators (computation of K could be avoided)
  K = gf_asm('linear elasticity', mim, mf, mf_ls, lambda*ones(1, nbd), ...
             mu*ones(1, nbd));
  disp(sprintf('Elastic energy at iteration %d: %g', niter, U*K*U'));
  S = gf_asm('volumic','V()+=comp()',mim);
  if (N == 2)
    disp(sprintf('Remaining surface of material: %g', S));
  else
    disp(sprintf('Remaining volume of material: %g', S));
  end

  DU = gf_compute(mf, U, 'gradient', mf_g);
  EPSU = DU + permute(DU, [2 1 3]);
  
  % Computation of the shape derivative
  if (N == 2)
    GF1 = (DU(1,1,:) + DU(2,2,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
  else
    GF1 = (DU(1,1,:) + DU(2,2,:) + DU(3,3,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
  end;
  GF = reshape(GF1, 1, size(GF1, 3)) - threshold_shape;
  
  % computation of the topological gradient
  if (N == 2)
     GT = -pi*( (lambda+2*mu) / (2*mu*(lambda+mu)) * (4*mu*GF1 + ...
          2*(lambda-mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:)).^2));
  else
     GT = -pi*( (lambda+2*mu) / (mu*(9*lambda+14*mu)) * (20*mu*GF1 + ...
        2*(3*lambda-2*mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:) + DU(3,3,:)).^2));
  end;
  GT = reshape(GT, 1, size(GT, 3)) + threshold_topo;
  
  % filtering the gradients
  M = gf_asm('mass matrix', mim, mf_g);
  D = abs(full(diag(M)));
  maxD = max(D);
  ind = find(D < maxD/40);
  % Extension of the gradient into the hole. Too rough ?
  GF(ind) = GF(ind) * 0;
  % Threshold on the gradient
  GF = min(GF, 2*threshold_shape);
  ind = find(D < maxD/1.2);
  GT(ind) = GT(ind) * 0 - 20;

  % Drawing the gradients
  if (mod(niter, NBDRAW) == 0 || niter == 1)
    if (N == 2)
      subplot(NG,1,1);
      gf_plot(mf_g, GF, 'disp_options', 'off', 'refine', 1);
      title('Shape gradient');
      hold on;
      [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                      'off', 'disp_options', 'off', 'refine', 3);
      set(h2{1},'LineWidth',1);
      set(h2{1},'Color','black');
      hold off;
      colorbar;
      if (DEBUG == 0)
        subplot(NG,1,2);
        % gf_plot(mf_g, GT, 'disp_options', 'off', 'disp_options');
        % title('Topological gradient');
        gf_plot(mf_ls, ULS, 'disp_options', 'off', 'refine', 3);
        hold on;
        [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                    'off', 'disp_options', 'off', 'refine', 3);
        set(h2{1},'LineWidth',1);
        set(h2{1},'Color','black');
        hold off;
        title('Level set function');
        colorbar; 
      end; 
      pause(0.3);
    else
      sl=gf_slice({'boundary', {'isovalues', -1, mf_ls, ULS, 0.0}}, m, 5);
      % sl=gf_slice({'isovalues', 0, mf_ls, ULS, 0.0}, m, 5);
      Usl=gf_compute(mf_g, GF,'interpolate on',sl);
      % P=gf_slice_get(sl,'pts'); P=P([1 3 2],:); gf_slice_set(sl,'pts',P);
      gf_plot_slice(sl,'data',Usl,'mesh','on','mesh_slice_edges_color', ...
                    [.7 .7 .7],'mesh_edges_color',[.5 .5 1]);
      colorbar;
      title('Shape gradient on the remaining volume');
      pause(0.3);
    end
    print(gcf,'-dpng','-r450', sprintf('demo_structural_optimi%dD_%d_%d_%d', N, niter, initial_holes, NY));
  end

  [val, i] = max(GT);
  disp(sprintf('Max value of the topological gradient: %g', val));

  % Making a new hole (topological optimization)
  if (val > 0)
    point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);
    if (N == 2) 
      disp(sprintf('Making a new hole whose center is (%g, %g)', point(1), point(2)));
      ULS = max(ULS, (hole_radius^2 - (P(1,:) - point(1)).^2 - ...
                (P(2,:) - point(2)).^2)/(2*hole_radius));
    else
      disp(sprintf('Making a new hole whose center is (%g, %g, %g)', point(1), point(2), point(3)));
      ULS = max(ULS, (hole_radius^2 - (P(1,:) - point(1)).^2 - ...
                (P(2,:) - point(2)).^2 - (P(3,:) - point(3)).^2)/(2*hole_radius));
    end;
  end;    
  

  % Evolution of the level-set thanks to shape derivative. Simple version.
  
  Fdisc = gf_asm('volumic source', mimls, mf_ls, mf_g, GF);
  % Vcont = Mcontls \ Fdisc;
  Vcont = cgs(Mcontls, Fdisc, 1e-8, 100000, RMcontls);
  ULS = ULS - Vcont' * level_set_rate;


  % Evolution of the level-set thank to shape derivative.
  % Hamilton-Jacobi equation. Less stable.
  %  dt = 0.006; NT = 10; ddt = dt / NT;
  %  for t = 0:ddt:dt
  %    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
  %    NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
  %    GFF = GF ./ NORMDLS;
  % 
  %    if (N == 2)
  %      V = DLS.*[GFF; GFF];
  %    else
  %      V = DLS.*[GFF; GFF; GFF];
  %    end;
  %
  %    Fdisc = gf_asm('volumic source', mimls, mf_cont, mf_g, V);
  %    % Vcont = Mcont \ Fdisc;
  %    Vcont = cgs(Mcont, Fdisc, 1e-8, 100000, RMcont);
  %
  %    gf_compute(mf_ls, ULS, 'convect', mf_cont,Vcont,ddt,2, 'extrapolation');
  %  end;

  if (DEBUG && mod(niter, NBDRAW) == 0)
    subplot(NG,1,2);
    gf_plot(mf_ls, ULS, 'disp_options', 'off', 'refine', 3);
    colorbar;
    hold on;
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                  'off', 'disp_options', 'off', 'refine', 3);
    set(h2{1},'LineWidth',1);
    set(h2{1},'Color','black');
    hold off;
    disp('Level set function after convection drawn');
    pause(0.3);
  end;
   
  % Re-initialization of the level-set
  dt = reinitialisation_time; NT = 10; ddt = dt / NT;
  ULS0 = ULS;
  for t = ddt:ddt:dt
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    Fdisc = gf_asm('volumic source', mimls, mf_cont, mf_g, DLS);
    % DLScont = Mcont \ Fdisc;
    DLScont = cgs(Mcont, Fdisc, 1e-8, 100000, RMcont);
    NORMDLS = sqrt(sum(reshape(DLScont, N, size(DLScont, 1)/N).^2,1))+1e-12;
    SULS = sign(ULS0) ./ NORMDLS;
      
    if (N == 2)
      W = DLScont.*reshape([SULS; SULS], N*size(SULS, 2), 1);
    else
      W = DLScont.*reshape([SULS; SULS; SULS], N*size(SULS, 2), 1);
    end;
  
    gf_compute(mf_ls, ULS, 'convect', mf_cont, W, ddt, 1, 'unchanged');

    ULS = ULS + ddt * sign(ULS0);
  end;

  if (DEBUG && mod(niter, 3) == 0)
    AA = sqrt(sum(DLS.^2, 1));
    disp(sprintf('Norm dls after: %g %g %g %g', AA(1), AA(2), AA(3), AA(4)));
    disp(sprintf('Norm dls after: max = %g, min = %g', max(AA), min(AA)));

    subplot(NG,1,3);
    gf_plot(mf_ls, ULS, 'disp_options', 'off');
    colorbar;
    hold on;
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', ...
                  'off', 'disp_options', 'off');
    set(h2{1},'LineWidth',1);
    set(h2{1},'Color','black');
    hold off;
    disp('Level set function after re-initialization drawn');
    pause(0.3);
  end;

  gf_workspace('pop');
  disp(sprintf('this iteration took %g minutes', (cputime-ti)/60));
end;
  
  
  
  
  
  
