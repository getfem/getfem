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

path = get_absolute_file_path('demo_stokes_3D_tank.sce');

printf('demo stokes_3D_tank started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

disp('3D stokes demonstration on a quadratic mesh');

compute = input('  1:compute the solution\n  0:load a previously computed solution\n ? ');

viscosity = 10;

R1  = list(list('9-(y.^2+(z-6.0).^2)'),list(0),list(0));
R2  = list(list('9-(y.^2+(z-6.0).^2)'),list(0),list(0));
R4  = list(list(0),list(0),list(0));

m = gf_mesh('import','GiD',path + 'data/tank_quadratic_2500.GiD.msh');
mfu = gf_mesh_fem(m,3);
mfp = gf_mesh_fem(m,1);
mfd = gf_mesh_fem(m,1);
mim = gf_mesh_im(m, gf_integ('IM_TETRAHEDRON(5)'));

gf_mesh_fem_set(mfu,'fem',gf_fem('FEM_PK(3,2)'));
gf_mesh_fem_set(mfd,'fem',gf_fem('FEM_PK(3,2)'));
gf_mesh_fem_set(mfp,'fem',gf_fem('FEM_PK(3,1)'));

all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));

P = gf_mesh_get(m,'pts');
INpid    = find(abs(P(1,:)+25) < 1e-4);
OUTpid   = find(abs(P(1,:)-25) < 1e-4);
TOPpid   = find(abs(P(3,:)-20) < 1e-4);
INfaces  = gf_mesh_get(m, 'faces from pid', INpid);
OUTfaces = gf_mesh_get(m, 'faces from pid', OUTpid);
TOPfaces = gf_mesh_get(m, 'faces from pid', TOPpid);

gf_mesh_set(m, 'region', 1, INfaces);
gf_mesh_set(m, 'region', 2, OUTfaces);
gf_mesh_set(m, 'region', 3, TOPfaces);
gf_mesh_set(m, 'region', 4, _setdiff(all_faces',union(union(INfaces',OUTfaces','r'),TOPfaces','r'),'rows')');

disp(sprintf('nbdof: mfu=%d, mfp=%d',gf_mesh_fem_get(mfu,'nbdof'),gf_mesh_fem_get(mfp,'nbdof')));

if (compute) then
  md = gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mfu);
  gf_model_set(md, 'add initialized data', 'lambda', [0]);
  gf_model_set(md, 'add initialized data', 'mu', [viscosity]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
  gf_model_set(md, 'add fem variable', 'p', mfp);
  gf_model_set(md, 'add linear incompressibility brick', mim, 'u', 'p');
  gf_model_set(md, 'add variable', 'mult_spec', 1);
  
  gf_model_set(md, 'add constraint with multipliers', 'p', 'mult_spec', sparse(ones(1, gf_mesh_fem_get(mfp, 'nbdof'))), [0]);
  gf_model_set(md, 'add initialized data', 'NeumannData', [0 -10 0]);
  gf_model_set(md, 'add source term brick', mim, 'u', 'NeumannData', 1);
  gf_model_set(md, 'add initialized fem data', 'Dir1data', mfd, gf_mesh_fem_get_eval(mfd, R1));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 1, 'Dir1data');
  gf_model_set(md, 'add initialized fem data', 'Dir2data',  mfd, gf_mesh_fem_get_eval(mfd, R2));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 2, 'Dir2data');
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 3);
  gf_model_set(md, 'add initialized fem data', 'Dir3data', mfd, gf_mesh_fem_get_eval(mfd, R4));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 4, 'Dir3data');

  disp('running solve... can take some minutes and needs ~600MB of memory');
  t0 = timer(); 

  gf_model_get(md, 'solve', 'lsolver', 'superlu', 'noisy');
  disp(sprintf('solve done in %.2f sec', timer()-t0));

  U = gf_model_get(md, 'variable', 'u');
  P = gf_model_get(md, 'variable', 'p');
  
  save(path + '/demo_stokes_3D_tank_UP.mat',U,P);
  disp('[the solution has been saved in ''demo_stokes_3D_tank_UP.mat'']');
else
  load(path + '/demo_stokes_3D_tank_UP.mat');
end

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

printf('demo stokes_3D_tank terminated\n');
