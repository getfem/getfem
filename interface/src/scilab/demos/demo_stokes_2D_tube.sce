// this example uses the "old" gf_solve instead of the bricks
// framework..

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

path = get_absolute_file_path('demo_stokes_2D_tube.sce');

printf('demo stokes_2D_tube started\n');

disp('2D stokes demonstration on a quadratic mesh');

pde = init_pde();

pde('type')   = 'stokes';
pde('viscos') = 1.0;
pde('bound')  = list();

pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list('-y.*(y-5)',0));
pde('bound')($)('H')    = [];

pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list(0,'+(x-20).*(x-25)'));
pde('bound')($)('H')    = [];

pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list(0,0));
pde('bound')($)('H')    = [];

m = gf_mesh('import','GiD', path + '/data/tube_2D_spline.GiD.msh');

mfulag = gf_mesh_fem(m,2);

pde('mf_u') = [];
pde('mf_p') = [];
pde('mf_d') = [];
pde('mim')  = [];

pde('mf_u') = gf_mesh_fem(m,2);
pde('mf_p') = gf_mesh_fem(m,1);
pde('mf_d') = gf_mesh_fem(m,1);
pde('mim')  = gf_mesh_im(m,gf_integ('IM_TRIANGLE(5)'));

// this is a good example of the usefullness of the cubic bubble
// -> if not used, the pressure has strange values
gf_mesh_fem_set(pde('mf_u'),'fem',gf_fem('FEM_PK_WITH_CUBIC_BUBBLE(2,2)'));
gf_mesh_fem_set(pde('mf_d'),'fem',gf_fem('FEM_PK(2,2)'));
gf_mesh_fem_set(pde('mf_p'),'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

// we use a P3 mesh fem for interpolation of the U field, since
// because of its cubic bubble function, the pde.mf_u is not lagrangian 
gf_mesh_fem_set(mfulag,'fem',gf_fem('FEM_PK(2,3)'));

all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
P         = gf_mesh_get(m,'pts');
INpid     = find(abs(P(1,:)) < 1e-4);
OUTpid    = find(abs(P(2,:)+20) < 1e-4);
INfaces   = gf_mesh_get(m, 'faces from pid', INpid);
OUTfaces  = gf_mesh_get(m, 'faces from pid', OUTpid);

gf_mesh_set(m, 'boundary', 1, INfaces);
gf_mesh_set(m, 'boundary', 2, OUTfaces);
gf_mesh_set(m, 'boundary', 3, _setdiff(all_faces',union(INfaces',OUTfaces','r'),'rows')');

tic;
[U,P] = gf_solve(pde);
disp(sprintf('solve done in %.2f sec', toc));

Ul = gf_compute(pde('mf_u'),U,'interpolate on',mfulag);

h = scf();
h.color_map = jetcolormap(255);

drawlater;
subplot(2,2,1); 
gf_plot(mfulag,Ul,'norm','on','deformation',Ul,'deformation_scale',0.1,	'deformed_mesh','on');
colorbar(min(Ul),max(Ul));
title('|U| plotted on the deformed mesh');

subplot(2,2,2); 
gf_plot(pde.mf_p,P(:)','deformation',U,'deformation_mf',pde('mf_u')); 
colorbar(min(P),max(P)); //YC: P or U ?
title('Pression on the deformed mesh');

subplot(2,2,3); 
gf_plot(mfulag,Ul(:)','mesh','on','meshopts',list());
gf_plot(pde('mf_p'),P(:)','refine',1);
colorbar(min(Ul),max(Ul)); 
title('Quiver plot of U, with color plot of the pression');

subplot(2,2,4); 
gf_plot(mfulag,Ul(:)','mesh','on','meshopts',list(), 'quiver_density',100,'quiver_scale',0.4); 
gf_plot(pde('mf_p'),P(:)'); 
colorbar(min(Ul),max(Ul));
title('Quiver plot zoomed');
drawnow;

printf('demo stokes_2D_tube terminated\n');
