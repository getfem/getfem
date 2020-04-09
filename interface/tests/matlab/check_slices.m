% Copyright (C) 2005-2020 Julien Pommier.
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


function check_slices(iverbose,idebug)
  global gverbose;
  global gdebug;  
  if (nargin >= 1),
    gverbose = iverbose;
    if (nargin == 2),
      gdebug = idebug;
    else 
      gdebug = 0; end;
  else 
    gverbose = 0;
  end;
  gf_workspace('clear all');
  m=gf_mesh('triangles grid',[-5:1:5],[-4:.8:4]);
%  m=gf_mesh('triangles grid',[-1 1],[-1 1]);
  gf_mesh_get(m,'cvid');
  gf_mesh_set(m,'del convex',[1]);
  mf=gf_mesh_fem(m,1);
  gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,2)'))
  U=gf_mesh_fem_get(mf,'eval', {'x.*x + y.*y'});
  sl=gf_slice({'planar',0,[.5;0],[1;0]},m,3);
  pp=gf_slice_get(sl,'pts');
  gfassert('abs(pp(1,:)-.5)<1e-15');
  sl2=gf_slice('points',m,pp(:,1:3));
  pp2=gf_slice_get(sl2,'pts');
  gfassert('abs(pp2(1,:)-.5)<1e-15');
  
%  n=8;sl=gf_slice(m,{'isovalues',-1,mf,U,0.25},n);
  sl=gf_slice({'isovalues',-1,mf,U,16.0},m,4);
%  gf_plot_slice(sl,'mesh','on','data',gf_compute(mf,U,'interpolate on',sl)); colorbar;
  pp=gf_slice_get(sl,'pts');
  gfassert('max(sqrt(sum(pp.^2,1)))<4.0000001');
  
  sl=gf_slice({'isovalues',0,mf,U,9.0},m,7);
  pp=gf_slice_get(sl,'pts');
  gfassert('max(abs(3-sqrt(sum(pp.^2,1))))<0.0015');

  

  N=1;m=gf_mesh('triangles grid',[-N:(2*N/3):N],[-N:(N/5):N]);
  m2=gf_mesh('cartesian',[-N:(N/5):N]+.1,[-N:(N/7):N]+.1);
  sl=gf_slice({'mesh',m2},m,3);%gf_plot_slice(sl,'mesh_faces','on');
  a=gf_slice_get(sl,'area') - 1.9*1.9;
  gfassert('a < 1e-10');
