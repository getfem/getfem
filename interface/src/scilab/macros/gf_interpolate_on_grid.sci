// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function [G,varargout]=gf_interpolate_on_grid(mf,U,varargin)
//  function G=gf_interpolate_on_grid(mf,U,X,Y,...)
//  interpolates a field defined on mesh_fem 'mf' on
//  a cartesian grid [X(1),X(2),...] x [Y(1),Y(2),...] x ...

dim=gf_get_mesh_dim(mf);

if (length(varargin) ~= dim) then
  error('wrong number of arguments');
end

if (gf_nb_dof(mf) ~= length(U(:,1))) then 
  error(sprintf('wrong dimensions for U, should be %d instead of %d',gf_nb_dof(mf),size(U,1)));
end;

// creates the cartesian mesh
mc = new_mesh;
gf_cartesian_mesh(mc, varargin(:));

// use basic Q1 interpolation on this mesh
fem_c=QK_fem(dim,1);lst=new_intset; 

// count the total number of elements
nb_elt=1;
npts = [];
for i=1:dim
  npts(i)=length(varargin(i));
  nb_elt = nb_elt*(npts(i)-1);
end;

// builds the integration method on a paralellepipedic cell
// YC: function ID too long
pfi=gf_intmethod_approx_simplex(1,3);
for i=1:dim, 
  pfi=gf_intmethod_approx_product(pfi, pfi);
end

add_to_intset(lst,1,nb_elt);
mf_c = new_mesh_fem(mc);
set_finite_element(mf_c, lst,fem_c, pfi);

Uc = gf_interpolate_on_other_mesh(mf, mf_c, U');
Uc=Uc';

xy = gf_get_interpolation_pts(mf_c); xy=xy';
[XY,I]=sortrows(xy);

Uc=Uc(I,:);
G=matrix(Uc,[npts size(Uc,2)]);
if (length(varargout)==1) then
  varargout(1)=I;
end
endfunction
