// Scilab GetFEM++ interface
// Copyright (C) 2011-2011 Tomas Ligursky, Yves Renard.
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
// Simple example of the bifurcation problem: -Delta(u) + u = lambda exp(u)
//
// This program is used to check that scilab-getfem is working. This is also
// a good example of use of GetFEM++.
//

lines(0);
stacksize('max');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');
lambda = 0;
direction = 1;
nbstep = 80;

maxit = 5;
thrit = 4;
minang = 0.993;
maxres_solve = 1.e-7;
noisy = 'very_noisy';

h_init = 1e-3;
h_max = 2e-1;
h_min = 1e-5;

// create a simple cartesian mesh
m = gf_mesh('cartesian', [0:.1:1]);

// create a mesh_fem for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m, 1);
// assign the Q1 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'classical fem', 1);

// Integration which will be used
mim = gf_mesh_im(m, 4);


md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add basic nonlinear brick', mim, 'u', 'u-lambda*exp(u)', '1-lambda*exp(u)', 'lambda');

scfac = 1 / gf_mesh_fem_get(mf, 'nbdof');


if (~isempty(noisy)) then
    printf('computing initial point\n');
end
gf_model_get(md, 'solve', noisy, 'max_iter', 100, 'max_res', maxres_solve);
[T_U, T_lambda, h, tau1, tau2, tau3] = gf_model_get(md, 'init Moore-Penrose continuation', 'lambda', scfac, direction, noisy, 'h_init', h_init);
U = gf_model_get(md, 'variable', 'u');
//printf('U = '); disp(U); printf('lambda = %e\n', lambda);
//printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

U_hist = zeros(1, nbstep + 1); lambda_hist = zeros(1, nbstep + 1);
U_hist(1) = U(1); lambda_hist(1) = lambda;
//tau_hist = zeros(1, nbstep + 1); tau_hist(1) = tau3;

scf(0); drawlater; clf();
subplot(2,1,1);
plot(lambda_hist(1), U_hist(1), 'k.');
xtitle('', 'lambda', 'u(0)');
subplot(2,1,2)
gf_plot_1D(mf, U, 'style', 'k.-');
xtitle('', 'x', 'u');
drawnow;

//scf(1); drawlater; clf();
//plot(0, tau_hist(1), 'k.');
//xtitle('', 'iteration', 'tau');
//drawnow;


for step = 1:nbstep
  sleep(1000);
  printf('\nbeginning of step %d\n', step);
  [T_U, T_lambda, h, tau1, tau2, tau3] = gf_model_get(md, 'Moore-Penrose continuation', 'lambda', scfac, T_U, T_lambda, h, noisy, 'max_iter', maxit, 'thr_iter', thrit, 'min_ang', minang, 'h_init', h_init, 'h_max', h_max, 'h_min', h_min, 'tau1', tau1, 'tau2', tau2, 'tau3', tau3);
  
  if (h == 0) then
    printf('Continuation has failed');
    break;
  end
  
  U = gf_model_get(md, 'variable', 'u');
  lambda = gf_model_get(md, 'variable', 'lambda');
  
//  printf('U = '); disp(U); printf('lambda = %e\n', lambda);
//  printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

  U_hist(step+1) = U(1); lambda_hist(step+1) = lambda;
//  tau_hist(step + 1) = tau3;

  scf(0); drawlater; clf();
  subplot(2,1,1);
  plot(lambda_hist(1:step+1), U_hist(1:step+1), 'k-');
  plot(lambda_hist(1:step), U_hist(1:step), 'ko');
  plot(lambda_hist(step+1), U_hist(step+1), 'k.');
  xtitle('', 'lambda', 'u(0)');
  subplot(2,1,2)
  gf_plot_1D(mf, U, 'style', 'k.-');
  xtitle('', 'x', 'u');
  drawnow;

//  scf(1); drawlater; clf();
//  plot(0:step, tau_hist(1:step + 1), 'k.-');
//  xtitle('', 'iteration', 'tau');
//  drawnow;
  
  // calculate the determinant of the augmented Jacobian directly
//  lambda = lambda + 1e-8; gf_model_set(md, 'variable', 'lambda', [lambda]);
//  gf_model_get(md, 'assembly', 'build_rhs');
//  F1 = gf_model_get(md, 'rhs');
//  lambda = lambda - 1e-8; gf_model_set(md, 'variable', 'lambda', [lambda]);
//  gf_model_get(md, 'assembly', 'build_all');
//  F0 = gf_model_get(md, 'rhs');
//  J(1:11,1:11) = gf_model_get(md, 'tangent_matrix');
//  J(1:11,12) = ((1 / 1e-8) * (F0 - F1))';
//  J(12,1:11) = T_U; J(12,12) = T_lambda; detJ = det(J);
//  printf('J = '); disp(J); printf('det(J) = %e\n', detJ);
//  scf(2); drawlater;
//  plot(step, detJ, 'k.');
//  xtitle('', 'iteration', 'tau');
//  drawnow;
  
  printf('end of step n° %d', step); printf(' / %d\n', nbstep);
end