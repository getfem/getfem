// Scilab GetFEM interface
// Copyright (C) 2011-2020 Tomas Ligursky, Yves Renard.
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
// Simple example of a bifurcation problem:
// -Delta(u) + u = lambda * exp(u).
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

// continuation data
datapath = get_absolute_file_path('demo_continuation.sce') + 'data/';
// If the file name bp_char is non-empty, the continuation will be started
// from the bifurcation point and the tangent with the index ind_branch
// saved there. Direction of the initial tangent will be determined by
// direction. Otherwise, the continuation will be initialised according to
// direction and lambda0.
bp_char = '';
//bp_char = 'continuation_step_62_bp.mat';
ind_branch = 2;
direction = 1;
lambda0 = 0;
nbstep = 80;

h_init = 2e-2;
h_max = 2e-1;
h_min = 2e-5;
mincos = 0.997;
noisy = 'noisy';

// create a simple cartesian mesh
m = gf_mesh('cartesian', [0:.1:1]);

// create a mesh_fem for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m, 1);
// assign the Q1 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'classical fem', 1);

// integration which will be used
mim = gf_mesh_im(m, 4);

// define the model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add data', 'lambda', 1);
gf_model_set(md, 'add nonlinear generic assembly brick', mim, ...
             '(u-lambda*exp(u))*Test_u');


// initialise the continuation
scfac = 1 / gf_mesh_fem_get(mf, 'nbdof');
S = gf_cont_struct(md, 'lambda', scfac, 'h_init', h_init, 'h_max', h_max, ...
                   'h_min', h_min, 'min_cos', mincos, noisy, ...
                   'singularities', 2);

if (~isempty(bp_char)) then
  load(datapath + bp_char);
  U = U_bp; lambda = lambda_bp;
  T_U = direction * T_U_bp(:, ind_branch);
  T_lambda = direction * T_lambda_bp(ind_branch);
  h = gf_cont_struct_get(S, 'init step size');
else
  lambda = lambda0;
  gf_model_set(md, 'variable', 'lambda', [lambda]);
  
  if (~isempty(noisy)) then
    printf('starting computing an initial point\n');
  end
  gf_model_get(md, 'solve', noisy, 'max_iter', 100);
  U = gf_model_get(md, 'variable', 'u');
  [T_U, T_lambda, h] = ...
    gf_cont_struct_get(S, 'init Moore-Penrose continuation', ...
                       U, lambda, direction);
end

U_hist = zeros(1, nbstep + 1); lambda_hist = zeros(1, nbstep + 1);
U_hist(1) = U(1); lambda_hist(1) = lambda;
//tau = gf_cont_struct_get(S, 'bifurcation test function');
sing_out = [];

scf(0); drawlater; clf();
subplot(2,1,1);
plot(lambda_hist(1), U_hist(1), 'k.');
xtitle('', 'lambda', 'u(0)');
subplot(2,1,2)
  if max(U)-min(U) < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [0, min(U)-1; 1, max(U)+1];
  end
gf_plot_1D(mf, U, 'style', 'k.-');
xtitle('', 'x', 'u');
drawnow;

//scf(1); drawlater; clf();
//plot(0, tau, 'k.');
//xtitle('', 'iteration', 'bifurcation test function');
//drawnow;

// continue from the initial point
for step = 1:nbstep
  //sleep(1000);
  printf('\nbeginning of step %d\n', step);
  [U, lambda, T_U, T_lambda, h, h0, sing_label] = ...
    gf_cont_struct_get(S, 'Moore-Penrose continuation',...
                       U, lambda, T_U, T_lambda, h);
  if (h == 0) then
    return
  elseif (~isempty(sing_label)) then
    if (sing_label == 'limit point') then
      s = 'Step ' + sci2exp(step) + ': ' + sing_label;
    elseif (sing_label == 'smooth bifurcation point') then
      [U_bp, lambda_bp, T_U_bp, T_lambda_bp]...
        = gf_cont_struct_get(S, 'sing_data');
      s = 'Step ' + sci2exp(step) + ': ' + sing_label + ', '...
         + sci2exp(size(T_U_bp, 2)) + ' branch(es) located';
      save(datapath + 'continuation_step_' + sci2exp(step) + '_bp.mat',...
       'U_bp', 'lambda_bp', 'T_U_bp', 'T_lambda_bp');
    end
    sing_out = [sing_out; s];
  end
  
  U_hist(step+1) = U(1); lambda_hist(step+1) = lambda;
//  tau = gf_cont_struct_get(S, 'bifurcation test function');

  scf(0); drawlater; clf();
  subplot(2,1,1);
  plot(lambda_hist(1:step+1), U_hist(1:step+1), 'k-');
  plot(lambda_hist(1:step), U_hist(1:step), 'ko');
  plot(lambda_hist(step+1), U_hist(step+1), 'k.');
  xtitle('', 'lambda', 'u(0)');
  subplot(2,1,2)
  if max(U)-min(U) < 1e-10 then
    a = gca(); a.tight_limits = "on";
    a.data_bounds = [0, min(U)-1; 1, max(U)+1];
  end
  gf_plot_1D(mf, U, 'style', 'k.-');
  xtitle('', 'x', 'u');
  drawnow;

//  scf(1); drawlater;
//  plot(step, tau, 'k.');
//  drawnow;
  
  printf('end of step n° %d / %d\n', step, nbstep)
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
