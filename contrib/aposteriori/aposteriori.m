% Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
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

% addpath ~/source++/getfem++/contrib/aposteriori/

gf_workspace('clear all');
mesh = gf_mesh('load', 'aposteriori.meshfem2');
mf = gf_mesh_fem('load', 'aposteriori.meshfem2', mesh);
mf_vm = gf_mesh_fem('load', 'aposteriori.meshfem_vm2', mesh);
U = load('aposteriori.U2')';
VM = load('aposteriori.VM2')';

% gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, ...
%     'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', '5%');



figure(2);
gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 2, 'deformation', U, ...
	'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', 1.0);
colorbar;
pause(0.001);

meshh = gf_mesh('load', 'aposteriori.meshh');
figure(1); gf_plot_mesh(meshh);

%figure(1); gf_plot_mesh(mesh);
% hold on;  gf_plot_mesh(meshh, 'convexes', 'on'); hold off

% caxis([1E3 2e8]);
% a = 1e-4; axis([-a a -a a]);

