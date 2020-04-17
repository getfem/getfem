// Copyright (C) 2010-2020 Yann COLLETTE.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 3 of the License,  or
// (at your option) any later version along with the GCC Runtime Library
// Exception either version 3.1 or (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License and GCC Runtime Library Exception for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_stokes_3D_tank_draw.sce');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

if (exists('U')~=1 | exists('P') ~= 1) then
  error('run demo_stokes_3D_tank2 first');
end

// a nice colormap
c = [0.0 0.0 1.0; 
     0.0 0.5 1.0; 
     0.0 1.0 0.5; 
     0.0 1.0 0.0; 
     0.5 1.0 0.0; 
     1.0 0.5 0.0; 
     1.0 0.4 0.0; 
     1.0 0.0 0.0; 
     1.0 0.2 0.0; 
     1.0 0.4 0.0; 
     1.0 0.6 0.0; 
     1.0 0.8 0.0];

h = scf();
h.color_map = c;

// slice the mesh with two half spaces
sl  = gf_slice(list('boundary',list('intersection',list('planar',+1,[0;0;0],[0;1;0]),list('planar',+1,[0;0;0],[1;0;0]))),m,6);
Usl = gf_compute(mfu,U,'interpolate on', sl);
Psl = gf_compute(mfp,P,'interpolate on', sl);

drawlater;
gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');
drawnow;

sl  = gf_slice(list('boundary',list('intersection',list('planar',+1,[0;0;6],[0;0;-1]),list('planar',+1,[0;0;0],[0;1;0]))),m,6);
Usl = gf_compute(mfu,U,'interpolate on', sl);
Psl = gf_compute(mfp,P,'interpolate on', sl);

drawlater;
gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');
drawnow;

sl2 = gf_slice(list('boundary',list('planar',+1,[0;0;0],[0;1;0])),m,6,_setdiff(all_faces',TOPfaces','rows')');
drawlater;
gf_plot_slice(sl2,'mesh_faces','off','mesh','on','pcolor','off');
drawnow;

// streamline "starting" points
hh = [1 5 9 12.5 16 19.5];
H  = [zeros(2,length(hh));hh];

// compute the streamlines
tsl  = gf_slice('streamlines',mfu,U,H);
Utsl = gf_compute(mfu,U,'interpolate on', tsl);

// render them with "tube plot"
drawlater;
[a,h] = gf_plot_slice(tsl,'mesh','off','tube_radius',.2,'tube_color','red');
title('Demo Stokes Tank 3D');
drawnow;

h.color_map = c;

printf('demo stokes_3D_tank_draw terminated\n');
