disp('This demo use levelset to impose (weakly) a Dirichlet condition on a part of an ');
disp('implicit boundary defined by the zero of the levelset and a Neumann condition on ');
disp('the remaining part of that boundary. A Poisson problem');

clear all;
gf_workspace('clear all');
NX=5;
ls_degree = 1;
R = 0.4;


m=gf_mesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5, -.5:(1/NX):.5);
%m=gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
ls=gf_levelset(m, ls_degree);
ls2=gf_LevelSet(m, ls_degree, 'with_secondary');

mf_ls=gfObject(gf_levelset_get(ls, 'mf'));
P=get(mf_ls, 'basic dof nodes');
x = P(1,:); y = P(2,:); z = P(3,:);

ULS = ((x.^2 + y.^2 + z.^2).^1.5 - R^3);
ULS2 = x;
gf_levelset_set(ls, 'values', ULS);
gf_levelset_set(ls2, 'values', ULS, ULS2);


mlsint=gfMeshLevelSet(m); set(mlsint, 'add', ls); set(mlsint, 'adapt');
mlsbound=gfMeshLevelSet(m); set(mlsbound, 'add', ls2); set(mlsbound, 'adapt');


mim_bound = gfMeshIm('levelset', mlsbound, 'boundary', gf_integ('IM_TETRAHEDRON(6)'));
set(mim_bound, 'integ', 4);

mim_bound2 = gfMeshIm('levelset', mlsint, 'boundary', gf_integ('IM_TETRAHEDRON(6)'));
set(mim_bound2, 'integ', 4);

mim_int = gfMeshIm('levelset', mlsint, 'inside', gf_integ('IM_TETRAHEDRON(6)'));
set(mim_int, 'integ', 4);

mfu0=gfMeshFem(m,1); set(mfu0, 'fem', gf_fem('FEM_QK(3,2)'));
mf_mult=gfMeshFem(m,1); set(mf_mult, 'fem', gf_fem('FEM_QK(3,1)'));


% Some verifications
A=gf_asm('volumic','V()+=comp()',mim_bound);
disp(sprintf('area : %g should be %g', A, 4*pi*R^2/2));
A=gf_asm('volumic','V()+=comp()',mim_bound2);
disp(sprintf('area : %g should be %g', A, 4*pi*R^2));
V=gf_asm('volumic','V()+=comp()',mim_int);
disp(sprintf('volume : %g should be %g', V, 4*pi*R^3/3.));


% partial mesh fem
dof_out = get(mfu0, 'dof from im', mim_int);
cv_out = get(mim_int, 'convex_index');
cv_in = setdiff(gf_mesh_get(m, 'cvid'), cv_out);
mfu = gfMeshFem('partial', mfu0, dof_out, cv_in);

md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add Laplacian brick', mim_int, 'u');

Volumic_data = gf_mesh_fem_get(mfu0, 'eval', { '60*sqrt(x.^2+y.^2+z.^2)' });
gf_model_set(md, 'add initialized fem data', 'VolumicData', mfu0, Volumic_data);
gf_model_set(md, 'add source term brick', mim_int, 'u', 'VolumicData');

surface_data = gf_mesh_fem_get(mfu0, 'eval', { '-15*(x.^2+y.^2+z.^2)' });
gf_model_set(md, 'add initialized fem data', 'SurfaceData', mfu0, surface_data);
gf_model_set(md, 'add source term brick', mim_bound2, 'u', 'SurfaceData');
gf_model_set(md, 'add initialized fem data', 'SurfaceData2', mfu0, -surface_data);
gf_model_set(md, 'add source term brick', mim_bound, 'u', 'SurfaceData2');

gf_model_set(md, 'add multiplier', 'mult_dir', mf_mult, 'u');
gf_model_set(md, 'add Dirichlet condition with multipliers', ...
	     mim_bound, 'u', 'mult_dir', -1);

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');





Sol_U = gf_mesh_fem_get(mfu0, 'eval', { sprintf('5*((%g)^3-(x.^2+y.^2+z.^2).^1.5)', R) });
ERRL2 = gf_compute(mfu, U, 'L2 dist', mim_int, mfu0, Sol_U);
ERRH1 = gf_compute(mfu, U, 'H1 semi dist', mim_int, mfu0, Sol_U);
disp(sprintf(' erreur L2 = %g\n erreur semi H1 = %g', ERRL2, ERRH1));





gf_plot(mfu0, U, 'mesh','on', 'cvlst', get(m, 'outer faces'), 'refine', 2);
gf_colormap('chouette');
