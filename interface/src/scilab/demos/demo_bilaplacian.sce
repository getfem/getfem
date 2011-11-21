lines(0);
stacksize('max');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

printf("demo bilaplacian started\n");

N  = 2;
NX = 10;
NY = 14;

m=gf_mesh('regular simplices',0:0.4/NX:0.4, 0:1.2/NY:1.2);
// m=gf_mesh('cartesian',0:1/NX:1, 0:1/NY:1);
// m=gf_mesh('cartesian',0:0.4/NX:0.4, 0:1.2/NY:1.2);


useKL = 1; // use the Kirchhoff-Love plate model, or just a pure
           // bilaplacian problem

D = 1.0;   // Flexion modulus
newbricks = 1; // Use new bricks or old ones.

if useKL then NU=0.3; end; // poisson ratio (0 <= NU <= 1)

mim = gf_mesh_im(m); 
mfu = gf_mesh_fem(m); 
mfd = gf_mesh_fem(m);

gf_mesh_im_set(mim, 'integ',gf_integ('IM_TRIANGLE(13)'));
gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_ARGYRIS'));
gf_mesh_fem_set(mfd, 'fem',gf_fem('FEM_PK(2,5)'));

//gf_mesh_im_set(mim, 'integ',gf_integ('IM_GAUSS_PARALLELEPIPED(2,10)'));
//gf_mesh_fem_set(mfu, 'fem',gf_fem('FEM_REDUCED_QUADC1_COMPOSITE'));
//gf_mesh_fem_set(mfd, 'fem',gf_fem('FEM_QK(2,3)'));

flst = gf_mesh_get(m, 'outer_faces');
n    = gf_mesh_get(m, 'normal of faces', flst);

ftop     = flst(:,find(abs(n(1,:)-1) < 1e-5));
fbottom  = flst(:,find(abs(n(1,:)+1) < 1e-5));
fleft    = flst(:,find(abs(n(2,:)+1) < 1e-5));
fright   = flst(:,find(abs(n(2,:)-1) < 1e-5));

FORCE_BOUNDARY          = 1;
MOMENTUM_BOUNDARY       = 2;
SIMPLE_SUPPORT_BOUNDARY = 3;
CLAMPED_BOUNDARY        = 4;

gf_mesh_set(m, 'region', FORCE_BOUNDARY, fright);
gf_mesh_set(m, 'region', SIMPLE_SUPPORT_BOUNDARY, [fleft ftop fbottom]);
gf_mesh_set(m, 'region', CLAMPED_BOUNDARY, [fleft ftop fbottom]);
gf_mesh_set(m, 'region', MOMENTUM_BOUNDARY, [ftop fbottom]);

FT = 2.;
sol_u = gf_mesh_fem_get_eval(mfd,list(list(sprintf('sin(%g*(x+y))',FT))));
sol_f = sol_u*FT*FT*FT*FT*N*N*D;
sol_lapl_u = -FT*FT*sol_u*N;

if (newbricks) then // uses new bricks
  md = gf_model('real');
  gf_model_set(md, 'add fem variable', 'u', mfu);

  if useKL
    gf_model_set(md, 'add initialized data', 'D', [D]);
    gf_model_set(md, 'add initialized data', 'nu', [NU]);
    gf_model_set(md, 'add Kirchhoff-Love plate brick', mim, 'u', 'D', 'nu');
    M = zeros(N,N, gf_mesh_fem_get(mfd,'nbdof'));
  else
    gf_model_set(md, 'add initialized data', 'D', [D]);
    gf_model_set(md, 'add bilaplacian brick', mim, 'u', 'D');
    M = zeros(1, gf_mesh_fem_get(mfd, 'nbdof'));
  end;

  gf_model_set(md, 'add initialized fem data', 'VolumicData', mfd, ...
    	       gf_mesh_fem_get_eval(mfd, list(list('1-(x-y).^2'))));
  gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
  
  gf_model_set(md, 'add initialized fem data', 'M', mfd, M);
  gf_model_set(md, 'add normal derivative source term brick', mim, 'u', ...
    	       'M', MOMENTUM_BOUNDARY);

  if (useKL) then
    H = zeros(N, N, gf_mesh_fem_get(mfd, 'nbdof'));
    F = zeros(N, gf_mesh_fem_get(mfd, 'nbdof'));
    gf_model_set(md, 'add initialized fem data', 'H', mfd, H);
    gf_model_set(md, 'add initialized fem data', 'F', mfd, F);
    gf_model_set(md, 'add Kirchhoff-Love Neumann term brick', mim, 'u', ...
		 'H', 'F', FORCE_BOUNDARY);
  else
    F = zeros(1, N, gf_mesh_fem_get(mfd, 'nbdof'));
    gf_model_set(md, 'add initialized fem data', 'F', mfd, F);
    gf_model_set(md, 'add normal source term brick', mim, 'u', 'F', ...
                 FORCE_BOUNDARY);
  end;
 
  gf_model_set(md, ...
 	     'add normal derivative Dirichlet condition with penalization', ...
 	       mim, 'u', 1e10, CLAMPED_BOUNDARY);
 
  gf_model_set(md, 'add Dirichlet condition with penalization', ...
    	       mim, 'u', 1e10, SIMPLE_SUPPORT_BOUNDARY);

  t0 = timer(); 
  gf_model_get(md, 'solve', 'noisy');
  U = gf_model_get(md, 'variable', 'u');
  disp(sprintf('solve done in %.2f sec', timer()-t0));
else // uses old bricks
  if useKL then
    b0 = gf_mdbrick('bilaplacian', mim, mfu, 'Kirchhoff-Love')
    gf_mdbrick_set(b0, 'param','D', D);
    gf_mdbrick_set(b0, 'param','nu', NU);
    M = zeros(1, gf_mesh_fem_get(mfd,'nbdof'));
  else
    b0 = gf_mdbrick('bilaplacian', mim, mfu);
    gf_mdbrick_set(b0, 'param','D', D);
    M = zeros(1, gf_mesh_fem_get(mfd, 'nbdof'));
  end

  b1 = gf_mdbrick('source term', b0);
  gf_mdbrick_set(b1, 'param', 'source_term', mfd, gf_mesh_fem_get_eval(mfd, list(list('1-(x-y).^2'))));

  b2 = gf_mdbrick('normal derivative source term',b1,MOMENTUM_BOUNDARY);
  gf_mdbrick_set(b2, 'param', 'source_term', mfd,M);

  if (useKL) then
    H = zeros(N*N, gf_mesh_fem_get(mfd, 'nbdof'));
    F = zeros(N, gf_mesh_fem_get(mfd, 'nbdof'));
    b3 = gf_mdbrick('neumann Kirchhoff-Love source term',b2,FORCE_BOUNDARY);
    gf_mdbrick_set(b3, 'param', 'M', mfd, H);
    gf_mdbrick_set(b3, 'param', 'divM', mfd, F);
  else
    F = zeros(1, N, gf_mesh_fem_get(mfd, 'nbdof'));
    b3 = gf_mdbrick('normal source term', b2, FORCE_BOUNDARY);
    gf_mdbrick_set(b3, 'param', 'normal_source_term', mfd, F);
  end

  b4 = gf_mdbrick('dirichlet on normal derivative', b3, mfd, CLAMPED_BOUNDARY, 'penalized');

  b5 = gf_mdbrick('dirichlet', b4, SIMPLE_SUPPORT_BOUNDARY, mfd, 'penalized');

  mds = gf_mdstate(b5);
  disp('running solve... ');
  t0  = timer(); 

  gf_mdbrick_get(b5, 'solve', mds, 'noisy');
  disp(sprintf('solve done in %.2f sec', timer()-t0));

  U = gf_mdstate_get(mds, 'state');
end

h = scf();
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mfu,U,'mesh','on');
colorbar(min(U),max(U));
h.color_map = jetcolormap(255);
drawnow;

disp(sprintf('H2 norm of the solution: %g', gf_compute(mfu,U,'H2 norm', mim)));

//err=gf_compute(mfu,U,'interpolate on',mfd) - sol_u;

//disp(sprintf('H1 norm of the error: %g', gf_compute(mfd,err,'H1 norm', mim)));
//disp(sprintf('H2 norm of the error: %g', gf_compute(mfd,err,'H2 norm', mim)));

printf("demo bilaplacian terminated\n");
