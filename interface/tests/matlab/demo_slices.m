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

% not working, not part of the getfem-interface distrib

[mf]=gfMeshFem('load','signorini_cou.mesh_fem'); m=mf.linked_mesh;
load signorini_cou.data; U=signorini_cou';

mfdu=gf_mesh_fem(m,1);
% the P2 fem is not derivable across elements, hence we use a discontinuous
% fem for the derivative of U.
gf_mesh_fem_set(mfdu,'fem',gf_fem('FEM_PRODUCT(FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,1),FEM_PK_DISCONTINUOUS(1,1)),FEM_PK_DISCONTINUOUS(1,1))'));

% on output size(DU)=[3,3,nbdof(mfdu)]
DU=gf_compute(mf,U,'gradient',mfdu);

% from the derivative, we compute the von mises stress
VM=zeros(1,gf_mesh_fem_get(mfdu,'nbdof'));
N=gf_mesh_get(m,'dim');
for i=1:size(DU,3),
  t=DU(:,:,i);
  E=(t+t')/2;
  VM(i) = sum(E(:).^2) - (1./N)*sum(diag(E))^2;
end;
lambda=1;
VM = 4*lambda^2*VM;



nrefine=6;
sl1=gf_slice({'boundary',{'none'}},m,nrefine);
c=[0.1;0.1;20.1];x=[1;0;0];y=[0;1;0];z=[0;0;1];
sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y},{'planar',+1,c,z}}},m,nrefine);
%sl2=gf_slice({'boundary',{'union',{'planar',+1,c,x},{'planar',+1,c,y}}},m,nrefine);



P=gf_slice_get(sl2,'pts'); dP=gf_compute(mf,U,'interpolate on',sl2); gf_slice_set(sl2, 'pts', P+dP);

VMsl=gf_compute(mfdu,VM,'interpolate on',sl2);
figure(1); h=gf_plot_slice(sl2,'mesh','on','data',VMsl); view(-80,-15); axis off; camlight;
figure(2); h=gf_plot_slice(sl1,'mesh_faces','on','mesh','on'); view(-85,-15); axis off; camlight; set(h,'facecolor',[.8 0 0]);

