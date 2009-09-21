disp('3D stokes demonstration on a quadratic mesh');

compute = input('  1:compute the solution\n  0:load a previously computed solution\n ? ');

gf_workspace('clear all');

global verbosity; verbosity = 1;

viscosity = 10;

R1  = list('9-(y.^2+(z-6.0).^2)';0;0);
R2  = list('9-(y.^2+(z-6.0).^2)';0;0);
R4  = list(0;0;0);

m = gf_mesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh');
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
//gf_mesh_set(m, 'region', 4, _setdiff(all_faces',union(union(INfaces',OUTfaces','r'),TOPfaces','r'),'rows')'); // YC:
gf_mesh_set(m, 'region', 4, _setdiff(all_faces',union(union(INfaces',OUTfaces','r'),TOPfaces','r'))');

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
  
  t0=timer(); 

  gf_model_get(md, 'solve', 'lsolver', 'superlu', 'noisy');
  disp(sprintf('solve done in %.2f sec', timer()-t0));

  U = gf_model_get(md, 'variable', 'u');
  P = gf_model_get(md, 'variable', 'p');
  
  save('demo_stokes_3D_tank_UP.mat',U,P);
  disp('[the solution has been saved in ''demo_stokes_3D_tank_UP.mat'']');
else
  load('demo_stokes_3D_tank_UP.mat');
end;

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

