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
m = gf_mesh('triangles grid',[-5:1:5],[-4:.8:4]);
//  m = gf_mesh('triangles grid',[-1 1],[-1 1]);
gf_mesh_get(m,'cvid');
gf_mesh_set(m,'del convex',[1]);
mf = gf_mesh_fem(m,1);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,2)'))
//U  = gf_mesh_fem_get(mf,'eval', list('x.*x + y.*y'));
U  = gf_mesh_fem_get_eval(mf, list(list('x.*x + y.*y')));
sl = gf_slice(list('planar',0,[.5;0],[1;0]),m,3);
pp = gf_slice_get(sl,'pts');
gfassert('abs(pp(1,:)-.5)<1e-15');
sl2 = gf_slice('points',m,pp(:,1:3));
pp2 = gf_slice_get(sl2,'pts');
gfassert('abs(pp2(1,:)-.5)<1e-15');
//  n=8;sl = gf_slice(m,list('isovalues',-1,mf,U,0.25),n);
sl = gf_slice(list('isovalues',-1,mf,U,16.0),m,4);
h = scf();
h.color_map = jetcolormap(255);
title('plot 1');
drawlater;
gf_plot_slice(sl,'mesh','on','data',gf_compute(mf,U,'interpolate on',sl)); 
colorbar(min(U),max(U));
drawnow;
pp = gf_slice_get(sl,'pts');
gfassert('max(sqrt(sum(pp.^2,1)))<4.0000001');
sl = gf_slice(list('isovalues',0,mf,U,9.0),m,7);
pp = gf_slice_get(sl,'pts');
gfassert('max(abs(3-sqrt(sum(pp.^2,1))))<0.0015');
N=1;
m  = gf_mesh('triangles grid',[-N:(2*N/3):N],[-N:(N/5):N]);
m2 = gf_mesh('cartesian',[-N:(N/5):N]+.1,[-N:(N/7):N]+.1);
sl = gf_slice(list('mesh',m2),m,3); 
h = scf();
h.color_map = jetcolormap(255);
title('plot 1');
drawlater;
gf_plot_slice(sl,'mesh_faces','on');
drawnow;
a  = gf_slice_get(sl,'area') - 1.9*1.9;
gfassert('a < 1e-10');
