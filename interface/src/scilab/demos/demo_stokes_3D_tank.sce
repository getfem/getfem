// YC: Object oriented example

disp('3D stokes demonstration on a quadratic mesh');

compute=input('  1:compute the solution\n  0:load a previously computed solution\n ? ');

gf_workspace('clear all');

global verbosity; verbosity=1;

viscosity = 10;

R1  = list('9-(y.^2+(z-6.0).^2)';0;0);
R2  = list('9-(y.^2+(z-6.0).^2)';0;0);
R4  = list(0;0;0);

m=gfMesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh');
mfu = gfMeshFem(m,3);
mfp = gfMeshFem(m,1);
mfd = gfMeshFem(m,1);
mim = gfMeshIm(m, gfInteg('IM_TETRAHEDRON(5)'));

set(mfu,'fem',gfFem('FEM_PK(3,2)'));
set(mfd,'fem',gfFem('FEM_PK(3,2)'));
set(mfp,'fem',gfFem('FEM_PK(3,1)'));

all_faces = get(m, 'outer faces', get(m, 'cvid'));

P=get(m,'pts');
INpid    = find(abs(P(1,:)+25) < 1e-4);
OUTpid   = find(abs(P(1,:)-25) < 1e-4);
TOPpid   = find(abs(P(3,:)-20) < 1e-4);
INfaces  = get(m, 'faces from pid', INpid);
OUTfaces = get(m, 'faces from pid', OUTpid);
TOPfaces = get(m, 'faces from pid', TOPpid);

set(m, 'region', 1, INfaces);
set(m, 'region', 2, OUTfaces);
set(m, 'region', 3, TOPfaces);
//set(m, 'region', 4, setdiff(all_faces',union(union(INfaces',OUTfaces','r'),TOPfaces','r'),'rows')'); // YC:
set(m, 'region', 4, setdiff(all_faces',union(union(INfaces',OUTfaces','r'),TOPfaces','r'))');

disp(sprintf('nbdof: mfu=%d, mfp=%d',get(mfu,'nbdof'),get(mfp,'nbdof')));

if (compute) then
  md=gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mfu);
  gf_model_set(md, 'add initialized data', 'lambda', [0]);
  gf_model_set(md, 'add initialized data', 'mu', [viscosity]);
  gf_model_set(md, 'add isotropic linearized elasticity brick', mim, 'u', 'lambda', 'mu');
  gf_model_set(md, 'add fem variable', 'p', mfp);
  gf_model_set(md, 'add linear incompressibility brick', mim, 'u', 'p');
  gf_model_set(md, 'add variable', 'mult_spec', 1);
  
  gf_model_set(md, 'add constraint with multipliers', 'p', 'mult_spec', sparse([ones(1, get(mfp, 'nbdof'))]), [0]);
  gf_model_set(md, 'add initialized data', 'NeumannData', [0 -10 0]);
  gf_model_set(md, 'add source term brick', mim, 'u', 'NeumannData', 1);
  gf_model_set(md, 'add initialized fem data', 'Dir1data', mfd, get(mfd, 'eval', R1));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 1, 'Dir1data');
  gf_model_set(md, 'add initialized fem data', 'Dir2data',  mfd, get(mfd, 'eval', R2));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 2, 'Dir2data');
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 3);
  gf_model_set(md, 'add initialized fem data', 'Dir3data', mfd, get(mfd, 'eval', R4));
  gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mfu, 4, 'Dir3data');

  disp('running solve... can take some minutes and needs ~600MB of memory');
  
  t0=cputime; 

  gf_model_get(md, 'solve', 'lsolver', 'superlu', 'noisy');
  disp(sprintf('solve done in %.2f sec', cputime-t0));

  U = gf_model_get(md, 'variable', 'u');
  P = gf_model_get(md, 'variable', 'p');
  
  save('demo_stokes_3D_tank_UP.mat',U,P);
  disp('[the solution has been saved in ''demo_stokes_3D_tank_UP.mat'']');
else
  load('demo_stokes_3D_tank_UP.mat');
end;

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

