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
%  Shape optimization with topological gradient
%   (with a fictitious domain approach).
%
%  This program is used to check that matlab-getfem is working. This is
%  also a good example of use of GetFEM.
%

gf_workspace('clear all');

% parameters
NX=30;
ls_degree = 1;
alpha = 1;
beta = 1;
rayon_trous = 0.2;

% Finite element and integration methods definition
m=gfMesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
% m=gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
% mim=gfMeshIm(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));

mf_basic=gfMeshFem(m, 1);
gf_mesh_fem_set(mf_basic,'fem',gf_fem('FEM_QK(2,2)'));
mls=gfMeshLevelSet(m);
ls=gfLevelSet(m, ls_degree);
set(mls, 'add', ls);
mf_ls=gfObject(get(ls, 'mf'));
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:);
ULS=1000*ones(1,numel(x));


% Loop on the topological optimization
while(1) 

  gf_workspace('push');

  set(ls, 'values', ULS);

  set(mls, 'adapt');

  mim_bound = gfMeshIm('levelset',mls,'boundary', gf_integ('IM_TRIANGLE(6)'));
  mim = gfMeshIm('levelset',mls,'outside', gf_integ('IM_TRIANGLE(6)'));
  set(mim, 'integ', 4);
  mf_mult=gfMeshFem(m); set(mf_mult, 'fem', gf_fem('FEM_QK(2,1)'));

  M = gf_asm('mass matrix', mim, mf_basic);
  D = abs(full(diag(M)));
  ind = find(D > 1E-8);

  mf = gf_mesh_fem('partial', mf_basic, ind);


  S = gf_asm('volumic','V()+=comp()',mim);
  disp('remaining surface :'); disp(S);


  % Problem definition (Laplace(u) + u = f)
  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mf);
  gf_model_set(md, 'add Laplacian brick', mim, 'u');
  gf_model_set(md, 'add fem data', 'VolumicData', mf_basic);
  gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
  gf_model_set(md, 'add initialized data', 'rho', [1.]);
  gf_model_set(md, 'add mass brick', mim, 'u', 'rho');
  gf_model_set(md, 'add multiplier', 'mult_dir', mf_mult, 'u');
  % To be completely robust, a stabilization should be used on the Dirichlet
  % boundary to ensure the inf-sup condition (Nitsche or Barbosa-Hughes)
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
	       mim_bound, 'u', 'mult_dir', -1);

  % Solving the direct problem.
  U0 = gf_mesh_fem_get(mf_basic, 'eval', ...
                       { '0.4*(3.*sin(pi*(x+y)) + ((x-0.5).^10 + (y-0.5).^10 + (x+0.5).^10 + (y+0.5).^10))' });
  gf_model_set(md, 'variable', 'VolumicData', U0);
  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');

  subplot(2,1,1);
  gf_plot(mf, U);
  hold on;
  [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off');
  set(h2{1},'LineWidth',2);
  set(h2{1},'Color','green');
  colorbar;
  title('u');
  hold off;


  % Solving the adjoint problem.
  UBASIC = gf_compute(mf, U, 'interpolate on', mf_basic);
  F = 2*(UBASIC-U0);
  gf_model_set(md, 'variable', 'VolumicData', F);
  gf_model_get(md, 'solve');
  W = gf_model_get(md, 'variable', 'u');


  % Computation of the topological gradient
  mf_g=gfMeshFem(m, 1);
  gf_mesh_fem_set(mf_g,'fem', ...
   gf_fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'));
  DU = gf_compute(mf, U, 'gradient', mf_g);
  DW = gf_compute(mf, W, 'gradient', mf_g);
  nbdof = gf_mesh_fem_get(mf_g, 'nbdof');
  DU = reshape(DU, 2, nbdof);
  DW = reshape(DW, 2, nbdof);

  UU = gf_compute(mf, U, 'interpolate on', mf_g);
  UU0 = gf_compute(mf_basic, U0, 'interpolate on', mf_g);
  LS = gf_compute(mf_ls, ULS, 'interpolate on', mf_g);

  G = (-4*pi*( alpha*(DU(1,:).^2 + DU(2,:).^2 + DU(1,:).*DW(1,:) + ...
       DU(2,:).*DW(2,:)) + beta*(UU-UU0).^2)) .* (sign(LS)+1.)/2;
  
  subplot(2,1,2);
  gf_plot(mf_g, G);
  title('Topological gradient');
  colorbar;
  pause(0.01);

  % Find the point where the topological gradient is minimum
  [val, i] = min(G);
  if (val >= -12)
    disp('Topological optimization finished.');
    return;
  end;

  point = gf_mesh_fem_get(mf_g, 'basic dof nodes', [i]);


  gf_workspace('pop');

  % Updating the level set to add the hole
  R = -(val+7) / 200;
  xc = point(1);
  yc = point(2);
  ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);

end;

