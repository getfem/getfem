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
%  Shape optimization of a structure with a coupling between topological and
%  shape gradient (with a fictitious domain approach).
%
%  This program is used to check that matlab-getfem is working. This is
%  also a good example of use of GetFEM++.
%

clear all;
gf_workspace('clear all');

% parameters
initial_holes = 1;     % Pre-existing holes or not.
NY = 40;               % Number of elements in y direction

k = 1;                 % Degree of the finite element method for u
N = 2                 % Dimension of the mesh (2 or 3).
lambda = 1;            % Lame coefficient
mu = 1;                % Lame coefficient
hole_radius = max(0.03, 2./NY);  % Hole radius for topological optimization
if (N == 2)
  CF = k*NY/40.;       % Correction factor. Usefull ?
else
  CF = k*NY/8.;        % Correction factor. Usefull ?
end
threshold_shape = CF * 0.9;
threshold_topo  = CF * 0.2;
NBDRAW = 5;            % Draw solution each NBDRAW iterations
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
  F = gf_mesh_fem_get(mf_basic, 'eval', {'0', '0', '-6*(abs(y) < 0.05).*(abs(z) < 0.05)'});
end;


for niter = 1:100000 % Optimization loop

  gf_workspace('push');

  set(ls, 'values', ULS);
  disp('Adapting the mesh');
  set(mls, 'adapt');
  disp('Mesh adapted');
  
  if (N == 2)
    mim = gfMeshIm('levelset',mls,'inside', gf_integ('IM_TRIANGLE(6)'));
  else
    mim = gfMeshIm('levelset',mls,'inside', gf_integ('IM_TETRAHEDRON(6)'));
  end;
  set(mim, 'integ', 4);

  M = gf_asm('mass matrix', mim, mf_basic);
  D = abs(full(diag(M)));
  ind = find(D > (1/NY)^N/10000000);
  mf = gf_mesh_fem('partial', mf_basic, ind);

  % Problem definition
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
  gf_model_set(md, 'add initialized data', 'penalty_param', [1E-7]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
  gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
  gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);
  
  % Solving the direct problem
  disp('solving the direct problem');
  gf_model_get(md, 'solve', 'max_res', 1E-7,'noisy');
  U = gf_model_get(md, 'variable', 'u');
  nbd = gf_mesh_fem_get(mf_ls, 'nbdof');
  
  % Computation of indicators (computation of K could be avoided)
  K = gf_asm('linear elasticity', mim, mf, mf_ls, lambda*ones(1, nbd), ...
             mu*ones(1, nbd));
  disp(sprintf('Elastic energy at iteration %d: %g', niter, U*K*U'));
  S = gf_asm('volumic','V()+=comp()',mim);
  if (N == 2)
    disp(sprintf('Remaining surface of matieral: %g', S));
  else
    disp(sprintf('Remaining volume of matieral: %g', S));
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
  GF = min(GF, 1.1*threshold_shape);
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
  

  % Evolution of the level-set thank to shape derivative. Simple version.
  Mcontls = gf_asm('mass matrix', mimls, mf_ls); % Could be computed only once.
                                                 % and factorized once !
  Fdisc = gf_asm('volumic source', mimls, mf_ls, mf_g, GF);
  Vcont = Mcontls \ Fdisc;
  ULS = ULS - Vcont' * 0.005;



  % Evolution of the level-set thank to shape derivative.
  % Hamilton-Jacobi equation. Less stable.
  Mcont = gf_asm('mass matrix', mimls, mf_cont); % Could be computed only once.
                                                 % and factorized once !

  if (0)
  dt = 0.006; NT = 10; ddt = dt / NT;
  for t = 0:ddt:dt
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
    GFF = GF ./ NORMDLS;
  
    if (N == 2)
      V = DLS.*[GFF; GFF];
    else
      V = DLS.*[GFF; GFF; GFF];
    end;
  
    Fdisc = gf_asm('volumic source', mimls, mf_cont, mf_g, V);
    Vcont = Mcont \ Fdisc;

    gf_compute(mf_ls, ULS, 'convect', mf_cont, Vcont, ddt, 2, 'extrapolation');
  end;
  end;

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
  dt = 0.01; NT = 4; ddt = dt / NT;
  ULS0 = ULS;
  for t = 0:ddt:dt
    DLS = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    Fdisc = gf_asm('volumic source', mimls, mf_cont, mf_g, DLS);
    DLScont = Mcont \ Fdisc;
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
end;
  
  
  
  
  
  
