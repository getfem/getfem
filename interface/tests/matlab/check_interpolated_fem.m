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


gf_workspace('clear all');
clf;
m1=gf_mesh('regular_simplices', 0:.5:2, 0:.4:1, 'degree', 2, 'noised');
%m1=gf_mesh('regular_simplices', 0:1:2, 0:.5:1, 'degree', 2, 'noised');
gf_plot_mesh(m1, 'refine' ,5, 'curved','on'); hold on;
mf1 = gfMeshFem(m1); 
mim1 = gfMeshIm(m1, gfInteg('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),4)'));
set(mf1, 'fem', gfFem('FEM_PK(2, 1)'));


m2=gfMesh('regular_simplices', 0:.3:3, -.2:.4:1.2, 'degree', 1,'noised');
%m2=gf_mesh('regular_simplices', [0 3], [0 1], 'degree', 1, 'noised');
%gf_plot_mesh(m2, 'refine' ,5, 'curved','on'); hold on;
mf2 = gfMeshFem(m2); 
mim2 = gfMeshIm(m2,gfInteg('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),4)'));
%mim2 = gfMeshIm(m2, gfInteg('IM_TRIANGLE(6)'));
set(mf2, 'fem', gfFem('FEM_PK(2, 1)'));

f = gfFem('interpolated fem', mf1, mim2)

set(mf2, 'fem', f);
gf_workspace('stats');



mf3=gfMeshFem(m2);
set(mf3, 'fem', gfFem('FEM_PK(2,1)'));
set(mf3, 'fem', gfFem('FEM_PK(2, 0)'), [1 2 3 5]);
mf4=gfMeshFem('sum', mf2, mf3);

set(m2, 'del convex', 4);


mf = mf4; nbd = get(mf, 'nbdof');
gf_plot(mf, rand(1, nbd), 'refine', 16);
%for i=1:nbd, 
%  U=zeros(1,nbd); U(i)=1;
%  disp(sprintf('dof %d/%d', i, nbd));
%  gf_plot(mf,U,'refine',16, 'mesh','on');
%  pause
%end;

gf_workspace('stats');
gf_delete(f);

gf_fem_get(f, 'char')
