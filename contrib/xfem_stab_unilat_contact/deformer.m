% Copyright (C) 2010-2012 Saber Amdouni.
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

mesh=gf_mesh('load','xfem_stab_unilat_contact.meshfem');    
mf=gf_mesh_fem('load','xfem_stab_unilat_contact.meshfem', mesh);
U=load('xfem_stab_unilat_contact.U');   
gf_plot(mf,U','mesh','off','norm','on','deformed_mesh','on','deformation_scale',1, 'deformation_mf', mf, 'deformation', U')
