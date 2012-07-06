% Copyright (C) 2012-2012 Yves Renard, Julien Pommier.
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
mesh=gf_mesh('load','xfem_stab_unilat_contact_friction.meshfem');    
mf=gf_mesh_fem('load','xfem_stab_unilat_contact_friction.meshfem', mesh);
U=load('xfem_stab_unilat_contact_friction.U');   
gf_plot(mf,U','mesh','off','norm','on','deformed_mesh','on','deformation_scale',1, 'deformation_mf', mf, 'deformation', U')
%caxis([0 0.009])