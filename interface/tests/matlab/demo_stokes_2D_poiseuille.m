% Copyright (C) 2005-2012 Julien Pommier.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

% this example uses the "old" gf_solve instead of the bricks
% framework..

gf_workspace('clear all');
disp('validation for 2D stokes with rectangular elements : Poiseuille flow with cartesian mesh');
clear pde; clear bc;
pde.type = 'stokes';
pde.viscos=1.0;
pde.bound{1}.type = 'Dirichlet';
pde.bound{1}.R  = {'y.*(y-1)',0}; pde.bound{1}.H={1 0; 0 1};

m=gf_mesh('cartesian',[0:.3:5],[0:.2:1]);
pde.mf_u=gf_mesh_fem(m,2); % U mesh_fem (vector field -> qdim=2)
pde.mf_p=gf_mesh_fem(m,1); % Pression mesh_fem
pde.mf_d=gf_mesh_fem(m,1); % Data mesh_fem (boundary conditions, source terms etc)
pde.mim=gf_mesh_im(m,  gf_integ('IM_EXACT_PARALLELEPIPED(2)'));
gf_mesh_fem_set(pde.mf_u,'fem',gf_fem('FEM_QK(2,2)'));
gf_mesh_fem_set(pde.mf_d,'fem',gf_fem('FEM_QK(2,1)'));
gf_mesh_fem_set(pde.mf_p,'fem',gf_fem('FEM_QK(2,1)'));
s=gf_mesh_fem_get(pde.mf_u,'char');pde.mf_u=gf_mesh_fem('from string',s,m);
all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
gf_mesh_set(m, 'boundary', 1, all_faces);
[U,P,pde]=gf_solve(pde);
pde
subplot(3,1,1); gf_plot(pde.mf_u,U(:)','dir','x','deformation',U,'deformation_scale','10%','deformed_mesh','on'); colorbar;
%subplot(3,1,2); gf_pdeplot(pde.mf_u,U(:,2)); colorbar; 
subplot(3,1,2); gf_plot(pde.mf_p,P(:)','deformation',U,'deformation_mf',pde.mf_u); colorbar;
subplot(3,1,3); gf_plot(pde.mf_u,U(:)','mesh','on'); hold on; gf_plot(pde.mf_p,P(:)','refine',1); hold off; colorbar; 
disp('Note that the dirichlet condition was described on a P1 fem');
disp('(visible on the deformed mesh: on boundaries, the deformation');
disp('is linear, hence there is a small error on the computed solution.');
