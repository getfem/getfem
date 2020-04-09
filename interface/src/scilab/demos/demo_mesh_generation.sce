// Scilab GetFEM interface
//
// Copyright (C) 2011-2020 Yves Renard, Yann Collette.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// Mesh generation with the experimental meshing procedure of Getfem which
// uses simple primitives to describe the mesh geometry. 
//
// This program is used to check that matlab-getfem is working. This is also
// a good example of use of GetFEM.
//

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all'); clear all;

N = 3;   // dimension of the mesh
K = 2;   // degree of the mesh (for curved boundaries)
if (N == 1) then
  mo = gf_mesher_object('ball', [0], 2);
  fixed_vertices = [0];
  h = 0.5;
elseif (N == 2) then
  mo = gf_mesher_object('ball', [0 4], 2);
  fixed_vertices = [0; 4];
  h = 0.5;
elseif (N == 3) then
  if (0) then
    mo1 = gf_mesher_object('ball', [0 0 1], 2);
    mo2 = gf_mesher_object('ball', [0 0 -1], 2);
    mo3 = gf_mesher_object('intersect', mo1, mo2);
    mo4 = gf_mesher_object('ball', [0 0 0], 1.3);
    mo5 = gf_mesher_object('union', mo4, mo3);
    mo6 = gf_mesher_object('ball', [-1 0 0], 1);
    mo  = gf_mesher_object('set minus', mo5, mo6);
    fixed_vertices = []; h = 0.3; 
  else
    alpha = %pi/5;
    L   = 20;
    R   = L * tan(alpha) * 0.7;
    mo1 = gf_mesher_object('cone', [0 0 0], [0 0 1], L, alpha);
    mo2 = gf_mesher_object('cylinder', [0 0 L], [0, 0, 1], L, R);
    mo  = gf_mesher_object('union', mo1, mo2);
    fixed_vertices = []; h = 2;
  end
elseif (N == 4) then
  mo = gf_mesher_object('ball', [0 0 0 4], 2);
  fixed_vertices = [0; 0; 0; 4];
  h = 1;
else
  error('It is not very reasonable to build a mesh in dimension greater than 4 !');
end
 
m = gf_mesh('generate', mo, h, K, fixed_vertices);

hh = scf();
hh.color_map = gf_colormap('chouette');
if (N <= 2) then
  gf_plot_mesh(m);
elseif (N == 3) then
  mf = gf_mesh_fem(m, 1);
  gf_mesh_fem_set(mf, 'classical fem', K);
  VAL = gf_mesh_fem_get_eval(mf, list('x+y+z'));
  gf_plot(mf, VAL, 'mesh', 'on', 'cvlst', gf_mesh_get(mf,'outer faces'), 'refine', 4);
  // axis on; camlight;
end
