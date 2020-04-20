% Copyright (C) 2006-2020 Julien Pommier.
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
m=gf_mesh('regular_simplices', -1:.2:1, -1:.2:1, 'degree', 2, 'noised');
%m=gf_mesh('cartesian', -1:.33:1, -1:.33:1);
ls=gfLevelSet(m, 2, 'x*x + y*y - 0.7*0.7', 'x-.4')
%ls=gfLevelSet(m, 2, 'x + y - 0.2'); %, 'x-5')
%ls=gfLevelSet(m, 2, 'x + y - 0.2', 'x-5')
ls2=gf_levelset(m, 2, '0.6*x*x + (y-0.1)*(y-0.1) - 0.6*0.6');
ls3=gf_levelset(m, 4, 'x*x + (y+.08)*(y+.08) - 0.05*0.05');

mls=gfMeshLevelSet(m)
set(mls, 'add', ls);
if 1,
  set(mls, 'sup', ls);
  set(mls, 'add', ls);
  set(mls, 'add', ls2);
  set(mls, 'add', ls2);
  set(mls, 'add', ls2);
  set(mls, 'add', ls3);
end;
set(mls, 'adapt');

gfObject(get(mls, 'linked_mesh'))

lls = gf_mesh_levelset_get(mls, 'levelsets')

cm = gfObject(get(mls, 'cut_mesh'))

ctip = get(mls, 'crack_tip_convexes')


mf=gfMeshFem(m); set(mf, 'classical_fem', 1);
mfls=gfMeshFem('levelset',mls,mf);

%gf_workspace('stats');

nbd = get(mfls,'nbdof')
if 1,
  sl=gfSlice({'none'}, mls, 2);
  %for i=1:nbd,
%  U=zeros(1,nbd);U(i)=1;
    U=rand(1,nbd);
    gf_plot(mfls,U,'refine',4,'zplot','on');
    %pause;
%end;
  hold on;
  %gf_plot_mesh(cm, 'curved', 'on','refine',8,'edges_width',2);
  gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
  hold off;
  %caxis([0 2]); colorbar;
else
  for i=1:nbd,
    U=zeros(1,nbd); U(i)=1;
    gf_plot(mfls,U,'refine',16);
    hold on;
    gf_plot_mesh(cm, 'curved', 'on','refine',8);
    hold on;
    gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
    hold off;
    pause
  end;
end;

