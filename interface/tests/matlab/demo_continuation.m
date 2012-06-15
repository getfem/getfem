% Copyright (C) 2011-2012 Tomas Ligursky, Yves Renard.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
% Simple exemple othe f bifurcation problem: -Delta(u) + u = lambda exp(u)
%
% This program is used to check that matlab-getfem is working. This is also
% a good example of use of GetFEM++.
%

gf_workspace('clear all');

lambda0 = 0;
direction = 1;
nbstep = 70;

maxit = 5;
thrit = 4;
maxres_solve = 1.e-7;
noisy = 'very_noisy'

h_init = 1e-3;
h_max = 2e-1;
h_min = 1e-6;

with_dirichlet = true;

% create a simple cartesian mesh
m = gf_mesh('cartesian', [0:.1:1]);

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% assign the P1 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'classical fem', 1);

% Integration which will be used
mim = gf_mesh_im(m, 4);

% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add initialized data', 'lambda', [lambda0]);
gf_model_set(md, 'add basic nonlinear brick', mim, 'u', 'u-lambda*exp(u)', '1-lambda*exp(u)', 'lambda');
if (with_dirichlet)
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 1);
end;

scfac = 1 / gf_mesh_fem_get(mf, 'nbdof');

if (noisy) disp('computing initial point\n'); end
gf_model_get(md, 'solve', noisy, 'max iter', 100, 'max_res', maxres_solve);
[T_U, T_lambda, h] = gf_model_get(md, 'init Moore-Penrose continuation', 'lambda', scfac, direction, noisy);
U = gf_model_get(md, 'variable', 'u');
lambda = gf_model_get(md, 'variable', 'lambda');
disp('U = '); disp(U); disp(sprintf('lambda = %e\n', lambda));
disp(sprintf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1))));

U_hist = zeros(1, nbstep + 1); lambda_hist = zeros(1, nbstep + 1);
U_hist(1) = max(U); lambda_hist(1) = lambda;

figure(1);
subplot(2,1,1);
plot(lambda_hist(1), U_hist(1), 'k.');
xlabel('lambda'); ylabel('max(u)');
if (with_dirichlet) axis([0 4 0 10]); else axis([0 0.4 0 11]); end
subplot(2,1,2)
gf_plot_1D(mf, U, 'style', 'k.-');
if (with_dirichlet) axis([0 1 0 10]); else axis([0 1 0 11]); end  
xlabel('x'); ylabel('u');
pause(1);


for step = 1:nbstep
  disp(sprintf('\nbeginning of step %d\n', step));
  [T_U, T_lambda, h] = gf_model_get(md, 'Moore-Penrose continuation', 'lambda', scfac, T_U, T_lambda, h, ...
      noisy, 'max_iter', maxit, 'thr_iter', thrit, 'h_init', h_init, 'h_max', h_max, 'h_min', h_min);
  U = gf_model_get(md, 'variable', 'u');
  lambda = gf_model_get(md, 'variable', 'lambda');
  % disp('U = '); disp(U);
  disp(sprintf('lambda = %e\n', lambda));
  % disp(sprintf('lambda - U(1) * exp(-U(1)) = %e\n', lambda - U(1) * exp(-U(1))));
   
  U_hist(step+1) = max(U); lambda_hist(step+1) = lambda;
    
  subplot(2,1,1);
  plot(lambda_hist(1:step+1), U_hist(1:step+1), 'k-');
  hold on;
  plot(lambda_hist(1:step), U_hist(1:step), 'ko');
  plot(lambda_hist(step+1), U_hist(step+1), 'k.');
  hold off;
  if (with_dirichlet) axis([0 4 0 10]); else axis([0 0.4 0 11]); end
  xlabel('lambda'); ylabel('max(u)');
  subplot(2,1,2)
  gf_plot_1D(mf, U, 'style', 'k.-');
  if (with_dirichlet) axis([0 1 0 10]); else axis([0 1 0 11]); end
  xlabel('x'); ylabel('u');
  pause(0.25);
  disp(sprintf('end of step n° %d', step)); disp(sprintf(' / %d\n', nbstep));
end



% gf_plot(mf,U,'mesh','on','contour',.01:.01:.1); 
% colorbar; title('computed solution');
