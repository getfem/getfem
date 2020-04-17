% Copyright (C) 2011-2020 Tomas Ligursky, Yves Renard.
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
% Simple example of the bifurcation problem: -Delta(u) + u = lambda * exp(u)
%
% This program is used to check that matlab-getfem is working. This is also
% a good example of use of GetFEM.
%

gf_workspace('clear all');
gf_util('trace level', 1);
gf_util('warning level', 3);

% continuation data
datapath = 'data/';
% If the file name bp_char is non-empty, the continuation will be started
% from the bifurcation point and the tangent with the index ind_tangent
% saved there, direction of that tangent will be determined by direction.
% Otherwise, the continuation will be initialised according to direction and
% lambda0.
bp_char = '';
%bp_char = 'continuation_step_62_bp.mat';
ind_tangent = 2;
direction = 1;
lambda0 = 0;
nbstep = 80;

h_init = 2e-2;
h_max = 2e-1;
h_min = 2e-5;
mincos = 0.997;
noisy = 'noisy';

with_dirichlet = false;

% create a simple cartesian mesh
m = gf_mesh('cartesian', [0:.1:1]);

% create a mesh_fem for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% assign the P1 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf, 'classical fem', 1);

% integration which will be used
mim = gf_mesh_im(m, 4);

% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);

% define the model
md = gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add data', 'lambda', 1);
gf_model_set(md, 'add nonlinear term', mim, '(u-lambda*exp(u))*Test_u');
        
if (with_dirichlet)
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', mf, 1);
end;

% initialise the continuation
scfac = 1 / gf_mesh_fem_get(mf, 'nbdof');
S = gf_cont_struct(md, 'lambda', scfac, 'h_init', h_init, 'h_max', h_max,...
                   'h_min', h_min, 'min_cos', mincos, noisy,...
                   'singularities', 2);

if (bp_char)
  load([datapath bp_char]);
  U = U_bp; lambda = lambda_bp;
  T_U = direction * T_U_bp(:, ind_tangent);
  T_lambda = direction * T_lambda_bp(ind_tangent);
  h = gf_cont_struct_get(S, 'init step size');
else
  lambda = lambda0;
  gf_model_set(md, 'variable', 'lambda', [lambda]);
  
  if (noisy) disp('starting computing an initial point'); end
  gf_model_get(md, 'solve', noisy, 'max iter', 100);
  U = gf_model_get(md, 'variable', 'u');
  [T_U, T_lambda, h] = ...
    gf_cont_struct_get(S, 'init Moore-Penrose continuation', ...
                       U, lambda, direction);
end

U_hist = zeros(1, nbstep + 1); lambda_hist = zeros(1, nbstep + 1);
U_hist(1) = U(1); lambda_hist(1) = lambda;

figure(1);
subplot(2,1,1);
plot(lambda_hist(1), U_hist(1), 'k.');
xlabel('lambda'); ylabel('U(1)');
if (with_dirichlet) axis([0 4 0 15]); else axis([0 0.4 0 15]); end
subplot(2,1,2)
gf_plot_1D(mf, U, 'style', 'k.-');
if (with_dirichlet) axis([0 1 0 15]); else axis([0 1 0 15]); end  
xlabel('x'); ylabel('u');
pause(1);

sing_out = {};
% continue from the initial point
for step = 1:nbstep
  disp(sprintf('\nbeginning of step %d', step));
  [U, lambda, T_U, T_lambda, h, h0, sing_label] = ...
    gf_cont_struct_get(S, 'Moore-Penrose continuation', ...
                       U, lambda, T_U, T_lambda, h);
                   
  % gf_model_get(md, 'test tangent matrix', 1E-8, 20, 0.0001);
                       
  if (h ==0) return
  elseif (sing_label)
    if (strcmp(sing_label, 'limit point'))
      s = ['step ' sprintf('%d', step) ': limit point' ];
    elseif (strcmp(sing_label, 'smooth bifurcation point'))
     [U_bp, lambda_bp, T_U_bp, T_lambda_bp]...
       = gf_cont_struct_get(S, 'sing_data');
     % save([datapath 'continuation_step_' sprintf('%d', step) '_bp.mat'], ...
     %     'U_bp', 'lambda_bp', 'T_U_bp', 'T_lambda_bp');
     s = ['step ' sprintf('%d', step) ': smooth bifurcation point, '...
          sprintf('%d', size(T_U_bp, 2)) ' branch(es) located'];
    end
    sing_out(size(sing_out, 1)+1,1) = {s};
  end
  
  U_hist(step+1) = U(1); lambda_hist(step+1) = lambda;
    
  subplot(2,1,1);
  plot(lambda_hist(1:step+1), U_hist(1:step+1), 'k-');
  hold on;
  plot(lambda_hist(1:step), U_hist(1:step), 'ko');
  plot(lambda_hist(step+1), U_hist(step+1), 'k.');
  hold off;
  if (with_dirichlet) axis([0 4 0 15]); else axis([0 0.4 0 15]); end
  xlabel('lambda'); ylabel('U(1)');
  subplot(2,1,2)
  gf_plot_1D(mf, U, 'style', 'k.-');
  if (with_dirichlet) axis([0 1 0 15]); else axis([0 1 0 15]); end
  xlabel('x'); ylabel('u');
  pause(0.25);
  disp(sprintf('end of step %d / %d', step, nbstep));
end

nsing = size(sing_out, 1);
if (nsing)
  disp('')
  disp('------------------------------')
  disp('   Detected singular points   ')
  disp('------------------------------')
  for i = 1:nsing
    disp(sing_out(i,:))
  end
end


% gf_plot(mf,U,'mesh','on','contour',.01:.01:.1); 
% colorbar; title('computed solution');
