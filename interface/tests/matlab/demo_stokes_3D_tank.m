disp('3D stokes demonstration on a quadratic mesh -- 512MB of memory needed for the solve!!');
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
  % unfortunately, the basic stokes solver is very slow...
  % on this problem, the fastest way is to reduce to a (full) linear system on the pression...
  % drawback: matlab will be killed if you don't have 512MB of
  % memory

  if (1),
    b0 = gfMdBrick('isotropic_linearized_elasticity',mim,mfu)
    set(b0, 'param','lambda', 0);
    set(b0, 'param','mu', viscosity);
    b1 = gfMdBrick('linear incompressibility term', b0, mfp);
  else
    % does not work..
    b1 = gfMdBrick('navier stokes', mim, mfu, mfp, viscosity);
  end;
  b2 = gfMdBrick('source term', b1, 1);
  set(b2, 'param', 'source_term', mfd, get(mfd, 'eval', {0;-10;0}));
  DIRICHLET_TYPE = 'penalized';
  bconst = gfMdBrick('constraint', b2, 'augmented', 1);
  set(bconst, 'constraints', [ones(1, get(mfp, 'nbdof'))], 0);
  bdir1 = gfMdBrick('dirichlet', bconst   , 1, mfu, DIRICHLET_TYPE);
  bdir2 = gfMdBrick('dirichlet', bdir1, 2, mfu, DIRICHLET_TYPE);
  bdir3 = gfMdBrick('dirichlet on normal component', bdir2, 3, mfd, DIRICHLET_TYPE);
  bdir4 = gfMdBrick('dirichlet', bdir3, 4, mfu, DIRICHLET_TYPE);
  set(bdir1, 'param', 'R', mfd, get(mfd, 'eval', R1));
  set(bdir2, 'param', 'R', mfd, get(mfd, 'eval', R2));
  set(bdir3, 'param', 'R', 0);
  set(bdir4, 'param', 'R', mfd, get(mfd, 'eval', R4));

  mds=gfMdState(bdir4)

  disp('running solve... can take some minutes and needs ~500MB of memory');
  
  t0=cputime; 
  get(bdir4, 'solve', mds, 'noisy');
  disp(sprintf('solve done in %.2f sec', cputime-t0));

  S=get(mds, 'state');
  U=S(1:get(mfu,'nbdof'));
  P=S(numel(U) + (1:get(mfp,'nbdof')));
  
  save demo_stokes_3D_tank_UP U P;
  disp('[the solution has been saved in "demo_stokes_3D_tank_UP.mat"]');
else
  load demo_stokes_3D_tank_UP;
end;

disp('Got a solution, now you can call demo_stokes_3D_tank_draw to generate graphics');

