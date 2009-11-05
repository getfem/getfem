// Matlab GetFEM++ interface
//
// Copyright (C) 2009 Alassane SY, Yves Renard.
//
// This file is a part of GetFEM++
//
// GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//  Shape optimization of a structure with a coupling between topological and
//  shape gradient (with a fictitious domain approach).
//
//  This program is used to check that matlab-getfem is working. This is
//  also a good example of use of GetFEM++.
//

lines(0);
stacksize('max');
gf_workspace('clear all');

// parameters
ls_degree       = 1;    // Degree of the level-set. Should be one for the moment.
k               = 1;    // Degree of the finite element method for u
lambda          = 1;    // Lame coefficient
mu              = 1;    // Lame coefficient
hole_radius     = 0.03; // 0.03 // Hole radius for topological optimization
initial_holes   = 0;    // Pre-existing holes or not.
threshold_shape = 0.6;
threshold_topo  = 0.6;
NX    = 40; // Number of elements in x direction
NY    = 40; // Number of elements in y direction
NZ    = 40; // Number of elements in z direction
N     = 2;  // Dimension of the mesh (2 or 3).
DEBUG = 0;

if (DEBUG) then
  NG = 4;
else
  NG = 2;
end

// Mesh definition
// m=gf_mesh('cartesian', -1:(1/NX):1, -.5:(1/NY):.5);
if (N == 2) then
  m = gf_mesh('regular simplices', -1:(1/NX):1, -.5:(1/NY):.5);
else
  m = gf_mesh('regular simplices', -1:(1/NX):1, -.5:(1/NY):.5, -.5:(1/NZ):.5);
end
pts =  gf_mesh_get(m, 'pts');

// Find the boundary GammaD and GammaN
pidleft = find((abs(pts(1, :)+1.0) < 1E-7));
fidleft = gf_mesh_get(m, 'faces from pid', pidleft);
normals = gf_mesh_get(m, 'normal of faces', fidleft);
fidleft = fidleft(:,find(abs(normals(1, :)+1) < 1E-3));
GAMMAD  = 2;
gf_mesh_set(m, 'region', GAMMAD, fidleft);

pidright = find((abs(pts(1, :)-1.0) < 1E-7));
fidright = gf_mesh_get(m, 'faces from pid', pidright);
normals  = gf_mesh_get(m, 'normal of faces', fidright);
fidright = fidright(:,find(abs(normals(1, :)-1) < 1E-3));
GAMMAN   = 3;
gf_mesh_set(m, 'region', GAMMAN, fidright);


// Definition of the finite element methods
_ls = gf_levelset(m, ls_degree);
mls = gf_mesh_levelset(m);
gf_mesh_levelset_set(mls, 'add', _ls);
mf_ls = gf_levelset_get(_ls, 'mf');
if (N == 2) then
  mimls = gf_mesh_im(m, gf_integ('IM_TRIANGLE(4)'));
else
  mimls = gf_mesh_im(m, gf_integ('IM_TETRAHEDRON(5)'));   
end
mf_basic = gf_mesh_fem(m, N);
gf_mesh_fem_set(mf_basic,'fem',gf_fem(sprintf('FEM_PK(%d,%d)', N, k)));
mf_g = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mf_g,'fem', gf_fem(sprintf('FEM_PK_DISCONTINUOUS(%d,%d)', N, k-1)));
mf_cont = gf_mesh_fem(m, N);
gf_mesh_fem_set(mf_cont,'fem', gf_fem(sprintf('FEM_PK(%d,%d)', N, ls_degree)));
 
// Definition of the initial level-set
if (initial_holes) then
  if (N == 2) then
    ULS = gf_mesh_fem_get_eval(mf_ls, list('-0.6-sin(%pi*4*x).*cos(%pi*4*y)'));
  else
    ULS = gf_mesh_fem_get_eval(mf_ls, list('-0.6-sin(%pi*4*x).*cos(%pi*4*y).*cos(%pi*4*z)'));
  end
else
  ULS = gf_mesh_fem_get_eval(mf_ls, list('x - 2'));
end

// Level-set nodes
P = gf_mesh_fem_get(mf_ls, 'basic dof nodes');

// Force on the right part (Neumann condition)
if (N == 2) then
  //F = gf_mesh_fem_get_eval(mf_basic, list('0', '-4*(abs(y) < 0.0125)')); // YC: Fill with 0 then with -4 where abs(y) < 0.0125
  F = gf_mesh_fem_get_eval(mf_basic, list('-4*(abs(y) < 0.0125)'));
else
  F = gf_mesh_fem_get_eval(mf_basic, list('0', '0', '-4*(abs(y) < 0.0125).*(abs(z) < 0.0125)'));
end

h = scf();
h.color_map = jetcolormap(255);
title('structural optimization');

while(1) // Optimization loop
  gf_workspace('push');

  gf_levelset_set(_ls, 'values', ULS);
  disp('Adapting the mesh');
  gf_mesh_levelset_set(mls, 'adapt');
  disp('Mesh adapted');
  
  if (N == 2) then
    mim = gf_mesh_im('levelset',mls,'inside', gf_integ('IM_TRIANGLE(6)'));
  else
    mim = gf_mesh_im('levelset',mls,'inside', gf_integ('IM_TETRAHEDRON(6)'));
  end
  gf_mesh_im_set(mim, 'integ', 4);

  M   = gf_asm('mass matrix', mim, mf_basic);
  D   = abs(full(diag(M)));
  if (N == 2) then
    ind = find(D > (1/(NX*NY))/10000000);
  else
    ind = find(D > (1/(NX*NY*NZ))/10000000);
  end
  mf  = gf_mesh_fem('partial', mf_basic, ind);

  // Problem definition
  md = gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mf);
  gf_model_set(md, 'add initialized data', 'mu', [mu]);
  gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
  gf_model_set(md, 'add initialized data', 'penalty_param', [1E-7]);          
  gf_model_set(md, 'add mass brick', mim, 'u', 'penalty_param');
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', 1, GAMMAD);
  gf_model_set(md, 'add initialized fem data', 'Force', mf_basic, F);
  gf_model_set(md, 'add source term brick', mim, 'u', 'Force', GAMMAN);
  
  // Solving the direct problem
  gf_model_get(md, 'solve');
  U   = gf_model_get(md, 'variable', 'u');
  nbd = gf_mesh_fem_get(mf_ls, 'nbdof');
  
  // Computation of indicators (computation of K could be avoided)
  K = gf_asm('linear elasticity', mim, mf, mf_ls, lambda*ones(1, nbd), mu*ones(1, nbd));
  disp(sprintf('Elastic energy: %g', U*K*U'));
  S = gf_asm('volumic','V()+=comp()',mim);
  disp(sprintf('Remaining surface of material: %g', S));

//   subplot(1,2,1);
//   gf_plot(mf, U);
//  hold on;
//    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0, 'pcolor', 'off');
//    set(h2(1),'LineWidth',2);
//    set(h2(1),'Color','green');
//    hold off;
//    colorbar;
  DU   = gf_compute(mf, U, 'gradient', mf_g);

  EPSU = DU + permute(DU, [2 1 3]); // YC: the permute function is really slow

  // Computation of the shape derivative
  if (N == 2) then
    GF1 = (DU(1,1,:) + DU(2,2,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
  else
    GF1 = (DU(1,1,:) + DU(2,2,:) + DU(3,3,:)).^2*lambda + 2*mu*(sum(sum(EPSU.^2, 1), 2));
  end
  GF = matrix(GF1, 1, size(GF1, 3)) - threshold_shape;
  
  // computation of the topological gradient
  if (N == 2) then
     GT = -%pi*( (lambda+2*mu) / (2*mu*(lambda+mu)) * (4*mu*GF1 + ...
          2*(lambda-mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:)).^2));
  else
     GT = -%pi*( (lambda+2*mu) / (mu*(9*lambda+14*mu)) * (20*mu*GF1 + ...
        2*(3*lambda-2*mu)*(lambda+mu)*(DU(1,1,:) + DU(2,2,:) + DU(3,3,:)).^2));
  end
  GT = matrix(GT, 1, size(GT, 3)) + threshold_topo;
  
  // filtering the gradients
  M    = gf_asm('mass matrix', mim, mf_g);
  D    = abs(full(diag(M)));
  maxD = max(D);
  ind  = find(D < maxD/4);
  // Extension of the gradient into the hole. Too rough ?
  GF(ind) = GF(ind) * 0 - threshold_shape/8;
  // Threshold on the gradient
  GF  = min(GF, 4*threshold_shape);
  ind = find(D < maxD/1.3);
  GT(ind) = GT(ind) * 0 - 20;

  // Drawing the gradients
  if (N == 2) then
    clf();
    subplot(NG,1,1);
    drawlater;
    gf_plot(mf_g, GF, 'disp_options', 'off', 'refine', 1);
    title('Shape gradient');
    //[h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', 'off', 'disp_options', 'off', 'refine', 3);
    //set(h2(1),'LineWidth',1);
    //set(h2(1),'Color','green');
    colorbar(min(GF),max(GF));
    subplot(NG,1,2);
    // gf_plot(mf_g, GT, 'disp_options', 'off', 'disp_options', 'off', 'refine', 8);
    gf_plot(mf_ls, ULS, 'disp_options', 'off', 'refine', 1);
    title('Level set function');
    colorbar(min(ULS),max(ULS));
    drawnow
    sleep(100);
  else
    // To be adapted for 3D ...
    clf();
    drawlater;
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', 'off', 'disp_options', 'off', 'refine', 3);
    //set(h2(1),'LineWidth',1);
    //set(h2(1),'Color','green');
    drawnow;
    sleep(100);
  end

  [val, i] = max(GT);
  disp(sprintf('Max value of the topological gradient: %g', val));

  // Making a new hole (topological optimization)
  if (val > 0) then
    point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);
    if (N == 2) then
      disp(sprintf('Making a new hole whose center is (%g, %g)', point(1), point(2)));
      ULS = max(ULS, (hole_radius^2 - (P(1,:) - point(1)).^2 - (P(2,:) - point(2)).^2)/(2*hole_radius));
    else
      disp(sprintf('Making a new hole whose center is (%g, %g, %g)', point(1), point(2), point(3)));
      ULS = max(ULS, (hole_radius^2 - (P(1,:) - point(1)).^2 - (P(2,:) - point(2)).^2 - (P(3,:) - point(3)).^2)/(2*hole_radius));
    end
  end;   
  
  // Evolution of the level-set thank to shape derivative. Computation of v
  DLS     = gf_compute(mf_ls, ULS, 'gradient', mf_g);
  NORMDLS = sqrt(sum(DLS.^2, 1)) + 0.000001;
  GFF     = GF ./ NORMDLS;
  
  if (N == 2) then
    V = DLS.*[GFF; GFF];
  else
    V = DLS.*[GFF; GFF; GFF];
  end
  
  Mcont = gf_asm('mass matrix', mimls, mf_cont); // Could be computed only once.
  Fdisc = gf_asm('volumic source', mimls, mf_cont, mf_g, V);
  Vcont = Mcont \ Fdisc;

  // Level set convection
  gf_compute(mf_ls, ULS, 'convect', mf_cont, Vcont, 0.01, 20,'extrapolation');

  if (DEBUG)
    disp('Drawing the level set function after convection');
    drawlater
    subplot(NG,1,3);
    gf_plot(mf_ls, ULS, 'disp_options', 'off');
    colorbar(min(ULS),max(ULS));
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', 'off', 'disp_options', 'off', 'refine', 3);
    //set(h2(1),'LineWidth',1);
    //set(h2(1),'Color','green');
    sleep(100);
  end
   
  // Re-initialization of the level-set
  dt = 0.05; NT = 20; ddt = dt / NT;
 
  for t = 0:ddt:dt
    DLS     = gf_compute(mf_ls, ULS, 'gradient', mf_g);
    Fdisc   = gf_asm('volumic source', mimls, mf_cont, mf_g, DLS);
    DLScont = Mcont \ Fdisc;
    NORMDLS = sqrt(sum(matrix(DLScont, N, size(DLScont, 1)/N).^2, 1)) + 1e-12;
    SULS    = sign(ULS) ./ NORMDLS;
        
    if (N == 2) then
      W = DLScont.*matrix([SULS; SULS], N*size(SULS, 2), 1);
    else
      W = DLScont.*matrix([SULS; SULS; SULS], N*size(SULS, 2), 1);
    end;
   
    gf_compute(mf_ls, ULS, 'convect', mf_cont, W, ddt, 1);
    ULS = ULS + ddt * sign(ULS);
  end

  if (DEBUG) then
    AA = sqrt(sum(DLS.^2, 1));
    disp(sprintf('Norm dls after: %g %g %g %g', AA(1), AA(2), AA(3), AA(4)));
    disp(sprintf('Norm dls after: max = %g, min = %g', max(AA), min(AA)));

    disp('Drawing the level set function after re-initialization');
    subplot(NG,1,4);
    drawlater;
    gf_plot(mf_ls, ULS, 'disp_options', 'off');
    colorbar(min(ULS),max(ULS));
    [h1,h2]=gf_plot(mf_ls, ULS, 'contour', 0,'pcolor', 'off', 'disp_options', 'off', 'refine', 3);
    //set(h2(1),'LineWidth',1);
    //set(h2(1),'Color','green');
    sleep(100);
  end

  gf_workspace('pop');
end
