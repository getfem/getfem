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

// Need to compile getfem with qhull first.
lines(0);
gf_workspace('clear all');

h = scf();
h.color_map = jetcolormap(255);

m = gf_mesh('regular_simplices', -1:.2:1, -1:.2:1, 'degree', 2, 'noised');
ls1 = gf_levelset(m, 2, 'sqr(x) + sqr(y) - sqr(0.7)', 'x-.4');
ls2 = gf_levelset(m, 2, '0.6*sqr(x) + sqr(y-0.1) - sqr(0.6)');
ls3 = gf_levelset(m, 4, 'sqr(x) + sqr(y+.08) - sqr(0.05)');
mls = gf_mesh_levelset(m);
gf_mesh_levelset_set(mls, 'add', ls1);
if 1 then
  gf_mesh_levelset_set(mls, 'sup', ls1);
  gf_mesh_levelset_set(mls, 'add', ls1);
  gf_mesh_levelset_set(mls, 'add', ls2);
  gf_mesh_levelset_set(mls, 'add', ls2);
  gf_mesh_levelset_set(mls, 'add', ls2);
  gf_mesh_levelset_set(mls, 'add', ls3);
end
gf_mesh_levelset_set(mls, 'adapt');
gf_mesh_levelset_get(mls, 'linked_mesh');
lls = gf_mesh_levelset_get(mls, 'levelsets');
cm = gf_mesh_levelset_get(mls, 'cut_mesh');
ctip = gf_mesh_levelset_get(mls, 'crack_tip_convexes');
mf   = gf_mesh_fem(m); gf_mesh_fem_set(mf, 'classical_fem', 1);
mfls = gf_mesh_fem('levelset',mls,mf);
//gf_workspace('stats');
nbd = gf_mesh_fem_get(mfls,'nbdof');
if 1 then
  sl = gf_slice(list('none'), mls, 2);
  U  = rand(1,nbd);
  drawlater;
  gf_plot(mfls,U,'refine',4,'zplot','on');
  gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
  drawnow;
  colorbar(min(U),max(U));
else
  for i=1:nbd
    U = zeros(1,nbd); U(i)=1;
    drawlater;
    gf_plot(mfls,U,'refine',16);
    gf_plot_mesh(cm, 'curved', 'on','refine',8);
    gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
    drawnow;
    pause
  end
end
