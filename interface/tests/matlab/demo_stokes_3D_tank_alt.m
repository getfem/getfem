% Copyright (C) 2005-2012 Julien Pommier.
%
% This file is a part of GETFEM++
%
% Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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


% old example, which uses the deprecated gf_solve function
% see demo_stokes_3D_tank.m for a "modern" version

disp('3D stokes demonstration on a quadratic mesh -- 512MB of memory needed for the solve!!');
compute=input('  1:compute the solution\n  0:load a previously computed solution\n ? ');
gf_workspace('clear all');
global verbosity; verbosity=1;
clear pde; 
pde.type = 'stokes';
pde.viscos=1.0;
pde.F = {0,0,0};
pde.bound{1}.R  = {'9-(y.^2+(z-6.0).^2)',0,0};
pde.bound{2}.R  = {'9-(y.^2+(z-6.0).^2)',0,0};
pde.bound{3}.R  = {0,0,0};pde.bound{3}.H = {0,0,0;0,0,0;0,0,1};pde.bound{3}.G  = {0,0,0};
pde.bound{4}.R  = {0,0,0};
pde.bound{1}.type = 'Dirichlet';
pde.bound{2}.type = 'Dirichlet';
%pde.bound{3}.type = 'Mixed'; 
pde.bound{3}.type = 'Dirichlet';
pde.bound{4}.type = 'Dirichlet';

m=gf_mesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh');
pde.mf_u=gf_mesh_fem(m,3);
mfulag=gf_mesh_fem(m,3);
pde.mf_p=gf_mesh_fem(m,1);
pde.mf_d=gf_mesh_fem(m,1);
pde.mim=gf_mesh_im(m, gf_integ('IM_TETRAHEDRON(5)'));
% this is a good example of the usefullness of the cubic bubble
% -> if not used, the pression has strange values
%gf_mesh_fem_set(pde.mf_u,'fem',gf_fem('FEM_PK_WITH_CUBIC_BUBBLE(3,2)')
gf_mesh_fem_set(pde.mf_u,'fem',gf_fem('FEM_PK(3,2)'));
gf_mesh_fem_set(pde.mf_d,'fem',gf_fem('FEM_PK(3,2)'));
%gf_mesh_fem_set(pde.mf_p,'fem',gf_fem('FEM_PK_DISCONTINUOUS(3,1)'))
gf_mesh_fem_set(pde.mf_p,'fem',gf_fem('FEM_PK(3,1)'));

% we use a P3 mesh fem for interpolation of the U field, since
% because of its cubic bubble function, the pde.mf_u is not lagrangian 
gf_mesh_fem_set(mfulag,'fem',gf_fem('FEM_PK(3,1)'));
all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));

P=gf_mesh_get(m,'pts');
INpid=find(abs(P(1,:)+25) < 1e-4);
OUTpid=find(abs(P(1,:)-25) < 1e-4);
TOPpid=find(abs(P(3,:)-20) < 1e-4);
INfaces=gf_mesh_get(m, 'faces from pid', INpid);
OUTfaces=gf_mesh_get(m, 'faces from pid', OUTpid);
TOPfaces=gf_mesh_get(m, 'faces from pid', TOPpid);
gf_mesh_set(m, 'boundary', 1, INfaces);
gf_mesh_set(m, 'boundary', 2, OUTfaces);
gf_mesh_set(m, 'boundary', 3, TOPfaces);
gf_mesh_set(m, 'boundary', 4, setdiff(all_faces',union(union(INfaces',OUTfaces','rows'),TOPfaces','rows'),'rows')');

disp(sprintf('nbdof: mf_u=%d, mf_p=%d',gf_mesh_fem_get(pde.mf_u,'nbdof'),gf_mesh_fem_get(pde.mf_p,'nbdof')));
if (compute),
  % unfortunately, the basic stokes solver is very slow...
  % on this problem, the fastest way is to reduce to a (full) linear system on the pression...
  % drawback: matlab will be killed if you don't have 512MB of memory
  pde.solver = 'brute_stokes'; 
  tic; [U,P]=gf_solve(pde); disp(sprintf('solve done in %.2f sec', toc));
  save demo_stokes_3D_tank_UP U P;
  disp('[the solution has been saved in "demo_stokes_3D_tank_UP.mat"]');
else
  load demo_stokes_3D_tank_UP;
end;

mfu=pde.mf_u;
mfp=pde.mf_p;
mfd=pde.mf_d;

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

