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
// Continuation of contact between a vee-shaped body and a rigid foundation.
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
gf_util('warning level', 3);

datapath = get_absolute_file_path('demo_continuation_vee.sce') + 'data/';

// parameters
plot_mesh = %f;
P = [0; 0];                 // (approximate) coordinates of the node whose
                            // dofs have to be visualised
clambda = 100; cmu = 82;    // Lame coefficients (in GPa)
friction_coeff = 1.8;       // coefficient of friction
r = 10;                     // augmentation parameter

// continuation data
// volume forces (in GN/m^3)
horizontal_force_init = -5; vertical_force_init = 1; 
horizontal_force_final = -5; vertical_force_final = -8;

direction = 1;    // direction of the initial tangent wrt the parameter
X0_char = '';     // initial approximation of the state variable
gm0 = 0;
nbstep = 4000;

niter = 200;   // maximum number of iterations for the initial solver
h_init = 5e-4;
h_max = 1e-1;
h_min = 5e-7;
maxit = 5;
thrit = 4;
maxres = 5e-12;
maxdiff = 5e-12;
mincos = 1 - 1e-5;
maxres_solve = 1e-12;
noisy = 'noisy';

// import a mesh (size in m) and refine it eventually
m = gf_mesh('load', datapath + 'vee_h_0.03.mesh');
for i = 1:0
  gf_mesh_set(m, 'refine');
end

// selection of the Dirichlet and contact boundaries
GAMMAD = 1;  GAMMAC = 2;
border  = gf_mesh_get(m,'outer faces');
normals = gf_mesh_get(m, 'normal of faces', border);

dirichlet_boundary = border(:, find((normals(2, :) > 1 - 1e-7)...
  | ((normals(1, :) - 1 / sqrt(2))^2 + (normals(2, :) - 1 / sqrt(2))^2 ...
    < 1e-7)));
gf_mesh_set(m, 'region', GAMMAD, dirichlet_boundary);
contact_boundary = border(:, find(normals(2, :) < -1 + 1e-7));
gf_mesh_set(m, 'region', GAMMAC, contact_boundary);

// plot the mesh
if (plot_mesh) then
  scf(0); drawlater; clf();
  //gf_plot_mesh(m, 'regions', [GAMMAC], 'refine', 1);
  gf_plot_slice(gf_slice(list('none'),m,1),'mesh', 'on',...
                'mesh_slice_edges_color', [0 0 0]);
  a = gca(); a.isoview="on"; a.box="hidden_axes";
  a.tight_limits = "on"; a.data_bounds = [0, 0; 0.15 + 0.95 / sqrt(2), 0.8];
  xtitle('Mesh','','');
  drawnow;
end

// finite-element methods
mfu = gf_mesh_fem(m, 2);
gf_mesh_fem_set(mfu, 'fem' ,gf_fem('FEM_PK(2,1)'));
mfvm = gf_mesh_fem(m, 1);
gf_mesh_fem_set(mfvm,'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,0)'));

// integration method
mim = gf_mesh_im(m, 4);

// elasticity model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'clambda', [clambda]);
gf_model_set(md, 'add initialized data', 'cmu', [cmu]);
gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', ...
                 'clambda', 'cmu');

// Dirichlet condition
Ddata = zeros(1, 2);
gf_model_set(md, 'add initialized data', 'Ddata', Ddata);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u',...
             mfu, GAMMAD, 'Ddata');

// parametrised volume force
force_init = [horizontal_force_init, vertical_force_init];
force_final = [horizontal_force_final, vertical_force_final];
gf_model_set(md, 'add data', 'gamma', 1);
gf_model_set(md, 'add initialized data', 'force_init', force_init);
gf_model_set(md, 'add initialized data', 'force_final', force_final);
gf_model_set(md, 'add data', 'force_current', 2);
gf_model_set(md, 'add source term brick', mim, 'u', 'force_current');

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
S = gf_cont_struct(md, 'gamma', 'force_init', 'force_final',...
  'force_current', scfac, 'h_init', h_init, 'h_max', h_max, ...
  'h_min', h_min, 'max_iter', maxit, 'thr_iter', thrit, 'max_res', maxres,...
  'max_diff', maxdiff, 'min_cos', mincos, 'max_res_solve', maxres_solve, ...
  'singularities', 1, 'non-smooth', noisy);

if (~isempty(X0_char)) then
  load(datapath + X0_char, 'X');
  gf_model_set(md, 'to variables', X);
end
gm = gm0;
force_current = (1 - gm) * force_init + gm * force_final;
gf_model_set(md, 'variable', 'gamma', [gm]);
gf_model_set(md, 'variable', 'force_current', force_current);

if (~isempty(noisy)) then
  printf('starting computing an initial point\n')
end
[iter, converged] = gf_model_get(md, 'solve', 'max_res', maxres,...
                                 'max_iter', niter, noisy);
if (converged ~= 1) then
  printf('No initial point found!');
  return
end

X = gf_model_get(md, 'from variables');
[T_X, T_gm, h] = gf_cont_struct_get(S, 'init Moore-Penrose continuation',...
                                    X, gm, direction);

U = gf_model_get(md, 'variable', 'u');
lambda_n = gf_model_get(md, 'variable', 'lambda_n');
lambda_t = gf_model_get(md, 'variable', 'lambda_t');
VM = gf_model_get(md, 'compute_isotropic_linearized_Von_Mises_or_Tresca',...
                  'u', 'clambda', 'cmu', mfvm);

pid = gf_mesh_get(m, 'pid from coords', P, 0.002);
P = gf_mesh_get(m, 'pts', pid);
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
sing_out = [];

fig = scf(1); drawlater; clf();
fig.color_map = jetcolormap(255);
colorbar(min(VM),max(VM));
gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
    'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 1);
xlabel('$x_{1}$'); ylabel('$x_{2}$');
title('Initial deformed configuration');
a = get("current_axes"); a.tight_limits = "on"; a.isoview = "on";
a.data_bounds = [-0.1, 0; 0.15 + 0.95 / sqrt(2), 0.8];
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
if lambda_nP_max-lambda_nP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, lambda_nP_min-1;...
                   gm_max+1e-4, lambda_nP_max+1];
end
plot(gm_hist(1), lambda_nP_hist(1), 'k.');
xtitle('', 'gamma', 'lambda_n(P)');

subplot(2,2,4);
if lambda_tP_max-lambda_tP_min < 1e-10 then
  a = gca(); a.tight_limits = "on";
  a.data_bounds = [gm_min-1e-4, lambda_tP_min-1;...
                   gm_max+1e-4, lambda_tP_max+1];
end
plot(gm_hist(1), lambda_tP_hist(1), 'k.');
xtitle('', 'gamma', 'lambda_t(P)');
drawnow;

// continue from the initial point
for step = 1:nbstep
  printf('\nbeginning of step %d\n', step);
  [X, gm, T_X, T_gm, h, h0, sing_label] = gf_cont_struct_get(S, ...
    'Moore-Penrose continuation', X, gm, T_X, T_gm, h);
  if (h == 0) then
    printf('Continuation has failed.\n')
    break
  elseif (sing_label == 'limit point') then
    s = 'Step ' + sci2exp(step) + ': ' + sing_label;
    sing_out = [sing_out; s];
  end
  
  gf_model_set(md, 'to variables', X);
  U = gf_model_get(md, 'variable', 'u');
  U_P = gf_compute(mfu, U, 'interpolate on', sl);
  lambda_n = gf_model_get(md, 'variable', 'lambda_n');
  lambda_t = gf_model_get(md, 'variable', 'lambda_t');
  VM = gf_model_get(md, ...
    'compute_isotropic_linearized_Von_Mises_or_Tresca', 'u', ...
    'clambda', 'cmu', mfvm);
  
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
  
  fig = scf(2); drawlater; clf();
  fig.color_map = jetcolormap(255);
  colorbar(min(VM),max(VM));
  gf_plot(mfvm, VM, 'deformed_mesh', 'on', 'deformation', U, ...
    'deformation_mf', mfu, 'deformation_scale', 1, 'refine', 1);
  xlabel('$x_{1}$'); ylabel('$x_{2}$');
  title('Current deformed configuration');
  a = get("current_axes"); a.tight_limits = "on"; a.isoview = "on";
  a.data_bounds = [-0.1, 0; 0.15 + 0.95 / sqrt(2), 0.8];
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
    a.data_bounds = [gm_min-1e-4, lambda_nP_min-1;...
                     gm_max+1e-4, lambda_nP_max+1];
  end
  plot(gm_hist(step:step+1), lambda_nP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), lambda_nP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), lambda_nP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'lambda_n(P)');
  
  subplot(2,2,4);
  if lambda_tP_max-lambda_tP_min < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [gm_min-1e-4, lambda_tP_min-1;...
                     gm_max+1e-4, lambda_tP_max+1];
  end
  plot(gm_hist(step:step+1), lambda_tP_hist(step:step+1), 'k-');
  plot(gm_hist(1:step), lambda_tP_hist(1:step), 'ko-');
  plot(gm_hist(step+1), lambda_tP_hist(step+1), 'k.');
  xtitle('', 'gamma', 'lambda_t(P)');
  drawnow;
  
  printf('end of step n° %d / %d\n', step, nbstep)
  if (abs(gm) > 1) then
    printf('|gamma| > 1, stop\n')
    break
  end
  //sleep(100);
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

//X = gf_model_get(md, 'from variables');
//save(datapath + "solution.sod", "X");
