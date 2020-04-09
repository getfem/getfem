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

disp('this file should be launched after demo_tripod.m as it assumes the tripod mesh and solutions are in memory')

cylinder = true; % plane or cylnder slice

gfObject(m)
%m=gf_mesh('from string',sm);
%mfu=gf_mesh_fem('from string',smfu,m);
%mfdu=gf_mesh_fem('from string',smfdu,m);

disp('plotting ... can also take some minutes!');
%r=[0.7 .7 .7]; l = r(end,:); s=63; s1=20; s2=25; s3=48;s4=55; for i=1:s, c1 = max(min((i-s1)/(s2-s1),1),0);c2 = max(min((i-s3)/(s4-s3),1),0); r(end+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; end; colormap(r); colorbar;
c=[0 0 .5; 0 .2 1; 0 .5 .8; 0 .9 0; .4 1 0; .9 .7 0; .8 0 0; 1 0 0];  colormap(c);
cnt=1;

if (cylinder)
    I = 46.1:-.1:4;
else
    I = -10.3:+.1:12;
end

for r=I
  clf
  if (cylinder)
    sl=gf_slice({'boundary',{'cylinder',-1,[0;0;0],[0;1;0],r}},mfu,U*10, 5);
  else
    sl=gf_slice({'boundary',{'planar',-1,[0;r;0],[0;1;0]}},mfu,U*10,5);
  end
  Usl=gf_compute(mfdu,VM,'interpolate on',sl);
  P=gf_slice_get(sl,'pts'); P=P([1 3 2],:); gf_slice_set(sl,'pts',P);
  gf_plot_slice(sl,'data',Usl,'mesh','on','mesh_slice_edges_color',[.7 .7 .7],'mesh_edges_color',[.5 .5 1]);
  view(0,30);
  caxis([0 7]);
  axis([-48 48 -5 15 -48 48]); axis off; camzoom(1.4*1.3); campan(+.8,-.2);
  camlight;
  print('-dpng','-r120',sprintf('tripod_slice_p%03d',cnt));
  cnt=cnt+1;
  gfObject(sl)
  pause(1)
  gf_delete(sl);
end

