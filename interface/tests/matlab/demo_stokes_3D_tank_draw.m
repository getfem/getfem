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

if (exist('U')~=1 | exist('P') ~= 1),
  error('run demo_stokes_3D_tank2 first');
end;
clf;

% slice the mesh with two half spaces
sl=gf_slice({'boundary',{'intersection',{'planar',+1,[0;0;0],[0;1;0]},{'planar',+1,[0;0;0],[1;0;0]}}},m,6);
Usl=gf_compute(mfu,U,'interpolate on', sl);
Psl=gf_compute(mfp,P,'interpolate on', sl);
gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');

sl=gf_slice({'boundary',{'intersection',{'planar',+1,[0;0;6],[0;0;-1]},{'planar',+1,[0;0;0],[0;1;0]}}},m,6);
Usl=gf_compute(mfu,U,'interpolate on', sl);
Psl=gf_compute(mfp,P,'interpolate on', sl);
hold on;gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off'); hold off;
  
  
sl2=gf_slice({'boundary',{'planar',+1,[0;0;0],[0;1;0]}},m,6,setdiff(all_faces',TOPfaces','rows')');
hold on; gf_plot_slice(sl2,'mesh_faces','off','mesh','on','pcolor','off'); hold off;

% streamline "starting" points
hh=[1 5 9 12.5 16 19.5];
H=[zeros(2,numel(hh));hh];

% compute the streamlines
tsl=gf_slice('streamlines',mfu,U,H);
Utsl=gf_compute(mfu,U,'interpolate on', tsl);
% render them with "tube plot"
hold on; [a,h]=gf_plot_slice(tsl,'mesh','off','tube_radius',.2,'tube_color','white'); hold off;

% a nice colormap
caxis([0 .7]);
c=[0 0 1; 0 .5 1; 0 1 .5; 0 1 0; .5 1 0; 1 .5 0; 1 .4 0; 1 0 0; 1 .2 0; 1 .4 0; 1 .6 0; 1 .8 0];
colormap(c); camlight;
