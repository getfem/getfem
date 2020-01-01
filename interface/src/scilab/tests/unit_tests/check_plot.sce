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

gf_workspace('clear all');
lines(0);

m2   = gf_mesh('triangles grid',[0:.1:1],[0:.1:1]);
m22  = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);
m3   = gf_mesh('cartesian',[0:.1:1],[0:.15:1],[0:.2:1]);
mf2  = gf_mesh_fem(m2,1);
mf22 = gf_mesh_fem(m22,1);
mf3  = gf_mesh_fem(m3,1);
mf2v = gf_mesh_fem(m2,2);
mf3v = gf_mesh_fem(m3,3);
gf_mesh_fem_set(mf2 ,'fem',gf_fem('FEM_PK(2,2)'));//,gf_integ('IM_TRIANGLE(5)'));
gf_mesh_fem_set(mf22,'fem',gf_fem('FEM_QK(2,2)'));//,gf_integ('IM_EXACT_PARALLELEPIPED(2)'));
gf_mesh_fem_set(mf2v,'fem',gf_fem('FEM_PK(2,1)'));//,gf_integ('IM_TRIANGLE(5)'));
gf_mesh_fem_set(mf3 ,'fem',gf_fem('FEM_QK(3,1)'));//,gf_integ('IM_EXACT_PARALLELEPIPED(3)'));
gf_mesh_fem_set(mf3v,'fem',gf_fem('FEM_QK(3,1)'));//,gf_integ('IM_NC_PARALLELEPIPED(3,2)'));
U2  = gf_mesh_fem_get_eval(mf2,list(list('x.*y')));
U2v = gf_mesh_fem_get_eval(mf2v,list(list('x.*y','1-x+y.*y')));
U3  = gf_mesh_fem_get_eval(mf3,list(list('z')));
U3v = gf_mesh_fem_get_eval(mf3v,list(list('z',0,'x+y')));
gf_workspace('push');
sl2 = gf_slice(list('none'),m2,2);
sl3 = gf_slice(list('none'),m3,2,gf_mesh_get(m3,'outer faces'));

scf();
drawlater;
title('plot 1');
subplot(2,1,1);
gf_plot_slice(sl2);
subplot(2,1,2);
gf_plot_slice(sl3); // Plot nothing
drawnow;

gf_workspace('pop');
gf_workspace('push');
sl2 = gf_slice(list('boundary'),m2,3);
sl3 = gf_slice(list('boundary',list('none')),m3,3);

scf();
title('plot 2');
drawlater;
subplot(3,1,1);
//gf_plot_slice(sl2);
subplot(3,1,2);
//gf_plot_slice(sl2, 'mesh','on','mesh_edges_color', [0 0 1], 'mesh_edges_width', 2, 'mesh_faces','on', 'mesh_faces_color', [1 0 0]);
subplot(3,1,3);
gf_plot_slice(sl3, 'tube','on','mesh_edges_color', [0 0 1], 'mesh_edges_width', 2, 'mesh_faces','on', 'mesh_faces_color', [1 0 0]); // Error
drawnow;

sl4 = gf_slice(list('planar',0,[.5;.5;.5],[0;0;1]),sl3);
P   = gf_slice_get(sl4,'pts');

scf()
title('plot 3');
drawlater;
//gf_plot_slice(sl4,'tube','on','tube_radius',0.05*abs(sin(P(2,:)*10))+0.01);
drawnow;

P = gf_slice_get(sl3,'pts');

scf();
title('plot 4');
drawlater;
gf_plot_slice(sl3, 'data',0.05*abs(sin(P(2,:)*10))+0.01);
drawnow;

gf_workspace('pop');
