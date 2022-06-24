// Scilab GetFEM interface
// Copyright (C) 2013-2020 Tomas Ligursky, Yves Renard.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
// Continuation and bifurcation of contact between a rectangular block and a
// rigid foundation.
//
// This program is used to check that scilab-getfem is working. This is also
// a good example of use of GetFEM.
//

gf_workspace('clear all');
lines(0);
stacksize('max');
if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linux, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level', 1);
gf_util('warning level', 0);

datapath = get_absolute_file_path('demo_continuation_block.sce') + 'data/';

// parameters
plot_mesh = %f;
P = [0; 0];                    // coordinates of the node whose dofs have to
                                // be visualised
lawname = 'Ciarlet Geymonat';
params  = [4000; 120; 30];      // in N/mm^2
friction_coeff = 1;           // coefficient of friction
r = 10;                         // augmentation parameter

// continuation data
// surface tractions on both vertical sides of the block (in N/mm^2)
p_1_left_init = -2; p_1_right_init = -2;
p_2_left_init = -2.4; p_2_right_init = 2.4;
p_1_left_final = 2; p_1_right_final = 2; 
p_2_left_final = 2.4; p_2_right_final = -2.4;

// If the file name sing_point_char is non-empty, the continuation will be
// started from the singular point and the tangent with the index ind_branch
// saved there. Otherwise, the continuation will be initialised according to
// direction, gm0, and eventually X0_char.
sing_point_char = '';
//sing_point_char = 'continuation_step_105_sing.sod';
ind_branch = 4;
direction = 1;
X0_char = '';
gm0 = 0;
nbstep = 500;

niter = 200;   // maximum number of iterations for the initial solver
h_init = 5e-4;
h_max = 5e-1;
h_min = 5e-7;
h_dec = 0.35;
maxit = 5;
thrit = 4;
maxres = 5e-10;
maxdiff = 5e-10;
mincos = 1 - 1e-5;
maxres_solve = 1e-10;
ndir = 20;
nspan = 15;
noisy = 'noisy';

// build a mesh (size in mm)
m = gf_mesh('cartesian Q1', [0:4:40], [0:4:80]);

// selection of the Dirichlet, Neumann and contact boundaries for the block
// Dirichlet on the top, Neumann on the vertical sides, contact at the bottom
GAMMAD = 1; GAMMAC = 2;  GAMMAN = 3;
border  = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);

dirichlet_boundary = border(:, find(normals(2, :) > 1 - 1e-7));
gf_mesh_set(m, 'region', GAMMAD, dirichlet_boundary);
contact_boundary = border(:, find(normals(2, :) < -1 + 1e-7));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);
neumann_boundary = border(:, find(abs(normals(1, :)) > 1 - 1e-7));
gf_mesh_set(m, 'region', GAMMAN, neumann_boundary);

// plot the mesh
if (plot_mesh) then
  scf(0); drawlater; clf();
  //gf_plot_mesh(m, 'regions', [GAMMAC], 'refine', 1);
  gf_plot_slice(gf_slice(list('none'),m,1),'mesh', 'on',...
                'mesh_slice_edges_color', [0 0 0]);
  a = gca(); a.isoview="on"; a.box="hidden_axes";
  a.tight_limits = "on"; a.data_bounds = [0, 0; 40, 80];
  xtitle('Mesh','','');
  drawnow;
end

// finite-element methods
mfu = gf_mesh_fem(m, 2);
gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_QK(2,1)'));
mfvm = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm,'fem',gf_fem('FEM_QK_DISCONTINUOUS(2,1)'));

// integration method
mim = gf_mesh_im(m, 4);

// elasticity model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'params', params);
gf_model_set(md, 'add nonlinear elasticity brick', mim, 'u', lawname,...
             'params');

// zero Dirichlet condition
Ddata = zeros(1, 2);
gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, GAMMAD, 'Ddata');

// parametrised Neumann condition
pdata_init = gf_mesh_fem_get_eval(mfu,...
     list(['(p_1_right_init-p_1_left_init)*x/40 + p_1_left_init',...
           '(p_2_right_init-p_2_left_init)*x/40 + p_2_left_init']));
pdata_final = gf_mesh_fem_get_eval(mfu,...
     list(['(p_1_right_final-p_1_left_final)*x/40 + p_1_left_final',...
           '(p_2_right_final-p_2_left_final)*x/40 + p_2_left_final']));
gf_model_set(md, 'add data', 'gamma', 1);
gf_model_set(md, 'add initialized fem data', 'pdata_init', mfu, pdata_init);
gf_model_set(md, 'add initialized fem data', 'pdata_final', mfu, pdata_final);
gf_model_set(md, 'add fem data', 'pdata_current', mfu, 1);
gf_model_set(md, 'add source term brick', mim, 'u','pdata_current', GAMMAN);

// contact conditions
nbdofu = gf_mesh_fem_get(mfu, 'nbdof');
cdof = gf_mesh_fem_get(mfu, 'dof on region', GAMMAC);
nbc  = size(cdof, 2) / 2;
contact_dof   = cdof(2:2:2*nbc);
contact_nodes = gf_mesh_fem_get(mfu, 'basic dof nodes', contact_dof);

BN   = spzeros(nbc, nbdofu);
ngap = zeros(nbc, 1);
for i = 1:nbc
  BN(i, contact_dof(i)) = -1.0;
  ngap(i) = contact_nodes(2, i);
end

BT = spzeros(nbc, nbdofu);
for i = 1:nbc
  BT(i, contact_dof(i) - 1) = 1.0;
end

gf_model_set(md, 'add variable', 'lambda_n', nbc);
gf_model_set(md, 'add initialized data', 'r', [r]);
gf_model_set(md, 'add variable', 'lambda_t', nbc);
gf_model_set(md, 'add initialized data', 'friction_coeff', [friction_coeff]);
gf_model_set(md, 'add initialized data', 'ngap', ngap);
gf_model_set(md, 'add initialized data', 'alpha', ones(nbc, 1));
gf_model_set(md, 'add basic contact brick', 'u', 'lambda_n', ...
             'lambda_t', 'r', BN, BT, 'friction_coeff', 'ngap', 'alpha', 1);
    
// initialise the continuation
scfac = 1 / (gf_model_get(md, 'nbdof'));
S = gf_cont_struct(md, 'gamma', 'pdata_init', 'pdata_final', ...
  'pdata_current', scfac, 'h_init', h_init, 'h_max', h_max, ...
  'h_min', h_min, 'h_dec', h_dec, 'max_iter', maxit, 'thr_iter', thrit, ...
  'max_res', maxres, 'max_diff', maxdiff, 'min_cos', mincos, ...
  'max_res_solve', maxres_solve, 'singularities', 2, 'non-smooth',...
  'nb_dir', ndir, 'nb_span', nspan, noisy);


if (~isempty(sing_point_char)) then
  load(datapath + sing_point_char);
  X = X_sing; gm = gm_sing;
  T_X = T_X_sing(:, ind_branch);
  T_gm = T_gm_sing(ind_branch);
  h = gf_cont_struct_get(S, 'init step size');
  gf_model_set(md, 'to variables', X);
else
  if (~isempty(X0_char)) then
    load(datapath + X0_char, 'X');
    gf_model_set(md, 'to variables', X);
  end
  gm = gm0;
  pdata_current = (1 - gm) * pdata_init + gm * pdata_final;
  gf_model_set(md, 'variable', 'gamma', [gm]);
  gf_model_set(md, 'variable', 'pdata_current', pdata_current);

  if (~isempty(noisy)) then
    printf('starting computing an initial point\n')
  end
  [iter, converged] = gf_model_get(md, 'solve', 'max_res', maxres, ...
                                   'max_iter', niter, noisy);
  if (converged ~= 1) then
    printf('No initial point found!');
    return
  end

  X = gf_model_get(md, 'from variables');
  [T_X, T_gm, h] = gf_cont_struct_get(S, 'init Moore-Penrose continuation',...
                                      X, gm, direction);
end

U = gf_model_get(md, 'variable', 'u');
lambda_n = gf_model_get(md, 'variable', 'lambda_n');
lambda_t = gf_model_get(md, 'variable', 'lambda_t');
VM = gf_model_get(md, 'compute Von Mises or Tresca', 'u', lawname, ...
                  'params', mfvm);

sl = gf_slice('points', m, P);
U_P = gf_compute(mfu, U, 'interpolate on', sl);
u_nP_hist = zeros(1, nbstep + 1); u_nP_hist(1) = -U_P(2);
u_nP_min = -U_P(2); u_nP_max = -U_P(2);
u_tP_hist = zeros(1, nbstep + 1); u_tP_hist(1) = U_P(1);
u_tP_min = U_P(1); u_tP_max = U_P(1);
indP = find(contact_nodes(1,:) == P(1));
lambda_nP = lambda_n(indP); lambda_tP = lambda_t(indP);
lambda_nP_hist = zeros(1, nbstep + 1); lambda_nP_hist(1) = lambda_nP;
lambda_nP_min = lambda_nP; lambda_nP_max = lambda_nP;
lambda_tP_hist = zeros(1, nbstep + 1); lambda_tP_hist(1) = lambda_tP;
lambda_tP_min = lambda_tP; lambda_tP_max = lambda_tP;
gm_hist = zeros(1, nbstep + 1); gm_hist(1) = gm;
gm_min = gm; gm_max = gm;
//tau = gf_cont_struct_get(S, 'bifurcation test function');
sing_out = [];

fig = scf(1); drawlater; clf();
fig.color_map = jetcolormap(255);
colorbar(min(VM),max(VM));
gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
        'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 1);
xlabel('$x_{1}$'); ylabel('$x_{2}$');
title('Initial deformed configuration');
a = get("current_axes"); a.isoview = "on";
a.tight_limits = "on"; a.data_bounds = [-5, 0; 45, 80];
drawnow;

scf(3); drawlater; clf();
subplot(2,2,1);
if u_nP_max-u_nP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, u_nP_min-1; gm_max+1e-4, u_nP_max+1];
end
plot(gm_hist(1), u_nP_hist(1), 'k.');
xtitle('', 'gamma', 'u_n(P)');
  
subplot(2,2,2);
if u_tP_max-u_tP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, u_tP_min-1; gm_max+1e-4, u_tP_max+1];
end
plot(gm_hist(1), u_tP_hist(1), 'k.');
xtitle('', 'gamma', 'u_t(P)');
  
subplot(2,2,3);
if lambda_nP_max- lambda_nP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, lambda_nP_min-1; ...
                   gm_max+1e-4, lambda_nP_max+1];
end
plot(gm_hist(1), lambda_nP_hist(1), 'k.');
xtitle('', 'gamma', 'lambda_n(P)');
  
subplot(2,2,4);
if lambda_tP_max-lambda_tP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, lambda_tP_min-1; ...
                   gm_max+1e-4, lambda_tP_max+1];
end
plot(gm_hist(1), lambda_tP_hist(1), 'k.');
xtitle('', 'gamma', 'lambda_t(P)');
drawnow;
  
//scf(4); drawlater; clf();
//plot(0, tau, 'k.');
//xtitle('', 'iteration', 'bifurcation test function');
//drawnow;

// continue from the initial point
for step = 1:nbstep
  printf('\nbeginning of step %d\n', step);
  [X, gm, T_X, T_gm, h, h0, sing_label] = gf_cont_struct_get(S, ...
    'Moore-Penrose continuation', X, gm, T_X, T_gm, h);
  if (h == 0) then
    printf('Continuation has failed.\n')
    break
  elseif (~isempty(sing_label)) then
    if (sing_label == 'limit point') then
      s = 'Step ' + sci2exp(step) + ': ' + sing_label;
    elseif (sing_label == 'non-smooth bifurcation point') then
      [X_sing, gm_sing, T_X_sing, T_gm_sing]...
        = gf_cont_struct_get(S, 'sing_data');
      save(datapath + 'continuation_step_' + sci2exp(step) + '_sing.sod',...
           'X_sing', 'gm_sing', 'T_X_sing', 'T_gm_sing');
      s = 'Step ' + sci2exp(step) + ': ' + sing_label + ', '...
            + sci2exp(size(T_X_sing, 2)) + ' branch(es) located';
    end
    sing_out = [sing_out; s];
  end
  
  gf_model_set(md, 'to variables', X);
  U = gf_model_get(md, 'variable', 'u');
  U_P = gf_compute(mfu, U, 'interpolate on', sl);
  lambda_n = gf_model_get(md, 'variable', 'lambda_n');
  lambda_t = gf_model_get(md, 'variable', 'lambda_t');
  VM = gf_model_get(md, 'compute Von Mises or Tresca', 'u', lawname, ...
                    'params', mfvm);
  
  u_nP_hist(step+1) = -U_P(2);
  u_nP_min = min(u_nP_min, -U_P(2)); u_nP_max = max(u_nP_max, -U_P(2));
  u_tP_hist(step+1) = U_P(1);
  u_tP_min = min(u_tP_min, U_P(1)); u_tP_max = max(u_tP_max, U_P(1));
  lambda_nP = lambda_n(indP); lambda_nP_hist(step+1) = lambda_nP;
  lambda_nP_min = min(lambda_nP_min, lambda_nP);
  lambda_nP_max = max(lambda_nP_max, lambda_nP);
  lambda_tP = lambda_t(indP); lambda_tP_hist(step+1) = lambda_tP;
  lambda_tP_min = min(lambda_tP_min, lambda_tP);
  lambda_tP_max = max(lambda_tP_max, lambda_tP);
  gm_hist(step+1) = gm;
  gm_min = min(gm_min, gm); gm_max = max(gm_max, gm);
  [tau, alpha_hist, tau_hist]...
    = gf_cont_struct_get(S, 'bifurcation test function');
  
  fig = scf(2); drawlater; clf();
  fig.color_map = jetcolormap(255);
  colorbar(min(VM),max(VM));
  gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
          'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 1);
  xlabel('$x_{1}$'); ylabel('$x_{2}$');
  title('Current deformed configuration');
  a = get("current_axes"); a.isoview = "on";
  a.tight_limits = "on"; a.data_bounds = [-5, 0; 45, 80];
  drawnow;

  scf(3); drawlater; clf();
  subplot(2,2,1);
  if u_nP_max-u_nP_min < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [gm_min-1e-4, u_nP_min-1; gm_max+1e-4, u_nP_max+1];
  end
  plot(gm_hist(step:step+1), u_nP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), u_nP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), u_nP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'u_n(P)');
  
  subplot(2,2,2);
  if u_tP_max-u_tP_min < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [gm_min-1e-4, u_tP_min-1; gm_max+1e-4, u_tP_max+1];
  end
  plot(gm_hist(step:step+1), u_tP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), u_tP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), u_tP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'u_t(P)');
  
  subplot(2,2,3);
  if lambda_nP_max-lambda_nP_min < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [gm_min-1e-4, lambda_nP_min-1; ...
                     gm_max+1e-4, lambda_nP_max+1];
  end
  plot(gm_hist(step:step+1), lambda_nP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), lambda_nP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), lambda_nP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'lambda_n(P)');
  
  subplot(2,2,4);
  if lambda_tP_max-lambda_tP_min < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [gm_min-1e-4, lambda_tP_min-1; ...
                     gm_max+1e-4, lambda_tP_max+1];
  end
  plot(gm_hist(step:step+1), lambda_tP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), lambda_tP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), lambda_tP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'lambda_t(P)');
  drawnow;
  
//  scf(4); drawlater;
//  l = length(tau_hist);
//  if l > 1 then
//    plot((step - 1) + alpha_hist(1:l), tau_hist(1:l), 'r.-')
//    plot(step-1, tau_hist(1), 'k.');
//  end
//  plot(step, tau, 'k.');
//  drawnow;
  
  printf('end of step n° %d / %d\n', step, nbstep)
  if (abs(gm-0.5) >= 0.5) then
    printf('|gamma - 0.5| >= 0.5, stop\n')
    break
  end
end

nsing = size(sing_out, 1);
if (nsing > 0) then
  printf('\n------------------------------\n')
  printf('   Detected singular points\n')
  printf('------------------------------\n')
  for i = 1:nsing
    printf(sing_out(i) + '\n')
  end
end
