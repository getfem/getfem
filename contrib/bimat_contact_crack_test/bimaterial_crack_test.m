% Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
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


% addpath ~/source++/getfem_toolbox
% addpath ~/source++/getfem++/contrib/bimat_contact_crack_test/

gf_workspace('clear all');
mesh = gf_mesh('load', 'bimaterial_crack.meshfem');
mf = gf_mesh_fem('load', 'bimaterial_crack.meshfem', mesh);
mf_vm = gf_mesh_fem('load', 'bimaterial_crack.meshfem_vm', mesh);
U = load('bimaterial_crack.U')';
VM = load('bimaterial_crack.VM')';
% VM = max(1, VM);
% VM = log(VM);
% VM = max(1E6, VM);
% VM = min(1E8, VM);
% for i = 1:size(VM, 2),
% if (VM(i) > 1E5)
%   VM(i) = 1E5;
% end;
% end;
% clear VM2; VM2(1:2:2*size(VM, 2)) = VM; VM2(2*size(VM, 2)) = 0;
gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', '5%');
caxis([1E3 2e8]);
a = 1e-4; axis([-a a -a a]);
colorbar;
