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


disp('3D stokes demonstration on a quadratic mesh');
compute=input('  1:compute the solution\n  0:load a previously computed solution\n ? ');
gf_workspace('clear all');
global verbosity; verbosity=1;


viscosity=10;

R1  = {'9-(y.^2+(z-6.0).^2)';0;0};
R2  = {'9-(y.^2+(z-6.0).^2)';0;0};
R4  = {0;0;0};


m=gfMesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh');
mfu=gfMeshFem(m,3);
mfp=gfMeshFem(m,1);
mfd=gfMeshFem(m,1);
mim=gfMeshIm(m, gfInteg('IM_TETRAHEDRON(5)'));

set(mfu,'fem',gfFem('FEM_PK(3,2)'));
set(mfd,'fem',gfFem('FEM_PK(3,2)'));
set(mfp,'fem',gfFem('FEM_PK(3,1)'));

all_faces = get(m, 'outer faces', get(m, 'cvid'));

P=get(m,'pts');
INpid=find(abs(P(1,:)+25) < 1e-4);
OUTpid=find(abs(P(1,:)-25) < 1e-4);
TOPpid=find(abs(P(3,:)-20) < 1e-4);
INfaces=get(m, 'faces from pid', INpid);
OUTfaces=get(m, 'faces from pid', OUTpid);
TOPfaces=get(m, 'faces from pid', TOPpid);
set(m, 'region', 1, INfaces);
set(m, 'region', 2, OUTfaces);
set(m, 'region', 3, TOPfaces);
set(m, 'region', 4, setdiff(all_faces',union(union(INfaces',OUTfaces','rows'),TOPfaces','rows'),'rows')');

disp(sprintf('nbdof: mfu=%d, mfp=%d',get(mfu,'nbdof'),get(mfp,'nbdof')));

if (compute),

  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mfu);
  gf_model_set(md, 'add initialized data', 'lambda', [0]);
  gf_model_set(md, 'add initialized data', 'mu', [viscosity]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', ...
               mim, 'u', 'lambda', 'mu');
  gf_model_set(md, 'add fem variable', 'p', mfp);
  gf_model_set(md, 'add linear incompressibility brick', mim, 'u', 'p');
  gf_model_set(md, 'add variable', 'mult_spec', 1);
  
  gf_model_set(md, 'add constraint with multipliers', 'p', 'mult_spec', ...
               sparse([ones(1, get(mfp, 'nbdof'))]), [0]);
  gf_model_set(md, 'add initialized data', 'NeumannData', [0 -10 0]);
  gf_model_set(md, 'add source term brick', mim, 'u', 'NeumannData', 1);
  gf_model_set(md, 'add initialized fem data', 'Dir1data', mfd, ...
               get(mfd, 'eval', R1));
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', mfu, 1, 'Dir1data');
  gf_model_set(md, 'add initialized fem data', 'Dir2data',  mfd, ...
               get(mfd, 'eval', R2));
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', mfu, 2, 'Dir2data');
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', mfu, 3);
  gf_model_set(md, 'add initialized fem data', 'Dir3data',  mfd, ...
               get(mfd, 'eval', R4));
  gf_model_set(md, 'add Dirichlet condition with multipliers', ...
               mim, 'u', mfu, 4, 'Dir3data');


  disp('running solve... can take some minutes and needs ~600MB of memory');
  
  t0=cputime; 

  gf_model_get(md, 'solve', 'lsolver', 'superlu', 'noisy');
  disp(sprintf('solve done in %.2f sec', cputime-t0));

  U = gf_model_get(md, 'variable', 'u');
  P = gf_model_get(md, 'variable', 'p');
  
  save demo_stokes_3D_tank_UP U P;
  disp('[the solution has been saved in "demo_stokes_3D_tank_UP.mat"]');
else
  load demo_stokes_3D_tank_UP;
end;

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

