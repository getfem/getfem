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

path = get_absolute_file_path('demo_tripod_anim.sce');

printf('demo tripod_anim started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

// You should run demo_tripod first ...
//m    = gf_mesh('import','gid', path + '/data/tripod.GiD.msh');
//mfu  = gf_mesh_fem('from string', smfu, m);
//mfdu = gf_mesh_fem('from string', smfdu, m);

//drawlater;
//gf_plot_mesh(m,'cvlst',gf_mesh_get(m,'outer faces'),'curved','on','edges_color',[1 0 0]);
//drawnow;
pr   = 1;
haut = 0;

h = scf();
h.color_map = jetcolormap(255);

Index = 0;

//for r=[6 14 26]
for r=6:4:60
  printf('slicing...\n'); tic;
  //sl = gf_slice(list('cylinder', 0, [0;0;0], [0;1;0], r), m, 4);
  //sl = gf_slice(m,list('boundary',list('cylinder', [0;0;0], [0;1;0], 15)),4);

  //sl = gf_slice(m,list('union',list('cylinderb', [0;0;0], [0;1;0], 15),list('cylinderb', [0;0;0], [0;1;0], 7)),4);
  sl = gf_slice(list('boundary',list('diff',list('cylinder', -1, [0;0;0], [0;1;0], r),list('cylinder', -1, [0;0;0], [0;1;0], r-4))),m, 6);
  //sl = gf_slice(m,list('none'),4, gf_mesh_get(m,'outer faces'));
  //sl = gf_slice(m,list('boundary', list('none')),4);
  //sl = gf_slice(m,list('ballb', [0;0;0], 10),2);
  //sl = gf_slice(m,list('planarb',[0;0;0],[0;0;1]),1);
  printf('..........done in %3.2f sec\n',toc());

  P      = gf_slice_get(sl,'pts'); 
  P(2,:) = P(2,:) - haut;
  sl
  D = gf_compute(mfdu,VM,'interpolate on',sl);
  
  drawlater;
  gf_plot_slice(sl, 'mesh','on','data',D,'pcolor','on','mesh_edges_color',[1 1 .7]);
  h.color_map = jetcolormap(255);

  a = gca();
  a.view = '3d';
  a.data_bounds = [-30 -15 -50;
                    50  15  50];

  xlabel('');
  ylabel('');
  zlabel('');

  a.axes_visible = ['off','off','off'];
  a.box = 'off';
  drawnow;

  // use:
  // convert -delay 50 -loop 0 wave*.png animatewave.gif
  // To produce the animated gif image.
  // Convert is an ImageMagick tool.
  xs2png(h.figure_id, path + sprintf('/tripod%02d.png',Index));
  
  Index = Index + 1;
  pr   = r;
  haut = haut + 24;
end

printf('demo tripod_anim terminated\n');
