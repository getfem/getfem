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
// Simple example of a bifurcation problem: Delta(u) + u = lambda exp(u)
//
// This program is used to check that scilab-getfem is working. This is also
// a good example of use of GetFEM++.
//

lines(0);
stacksize('max');

gf_workspace('clear all');
lambda = 0;
direction = 1;
nbstep = 70;

maxit = 5;
thrit = 4;
maxres_solve = 1.e-7;
noisy = 'very_noisy';

h_init = 1e-3;
h_max = 1e-2;
h_min = 1e-6;

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


nb_dof = gf_mesh_fem_get(mf, 'nbdof') + 1;
h_init = h_init * nb_dof;
h_max = h_max * nb_dof;
h_min = h_min * nb_dof;

if (~isempty(noisy)) then
    printf('computing initial point\n');
end
gf_model_get(md, 'solve', noisy, 'max_iter', 100, 'max_res', maxres_solve);
[T_U, T_lambda, h] = gf_model_get(md, 'init Moore-Penrose continuation', 'lambda', direction, noisy, 'h_init', h_init);
U = gf_model_get(md, 'variable', 'u');
printf('U = '); disp(U); printf('lambda = %e\n', lambda);
printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

U_hist = zeros(1, nbstep); lambda_hist = zeros(1, nbstep);
U_hist(1) = U(1); lambda_hist(1) = lambda;

clf();
drawlater;
subplot(2,1,1);
plot(lambda_hist(1), U_hist(1), 'k.');
xtitle('', 'lambda', 'u(0)');
subplot(2,1,2)
gf_plot_1D(mf, U, 'style', 'k.-');
xtitle('', 'x', 'u');
drawnow;


for step = 1:nbstep
  sleep(1000);
  printf('\nbeginning of step %d\n', step);
  [T_U, T_lambda, h] = gf_model_get(md, 'Moore-Penrose continuation', 'lambda', T_U, T_lambda, h, noisy, 'max_iter', maxit, 'thr_iter', thrit, 'h_init', h_init, 'h_max', h_max, 'h_min', h_min);
  U = gf_model_get(md, 'variable', 'u');
  lambda = gf_model_get(md, 'variable', 'lambda');
  printf('U = '); disp(U); printf('lambda = %e\n', lambda);
  printf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1)));

  U_hist(step+1) = U(1); lambda_hist(step+1) = lambda;

  clf();
  drawlater;
  subplot(2,1,1);
  plot(lambda_hist(1:step+1), U_hist(1:step+1), 'k-');
  plot(lambda_hist(1:step), U_hist(1:step), 'ko');
  plot(lambda_hist(step+1), U_hist(step+1), 'k.');
  xtitle('', 'lambda', 'u(0)');
  subplot(2,1,2)
  gf_plot_1D(mf, U, 'style', 'k.-');
  xtitle('', 'x', 'u');
  drawnow;
  printf('end of step n° %d', step); printf(' / %d\n', nbstep);
end