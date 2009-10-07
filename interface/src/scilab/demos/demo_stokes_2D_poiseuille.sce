// this example uses the "old" gf_solve instead of the bricks framework..

lines(0);
stacksize('max');

gf_workspace('clear all');

disp('validation for 2D stokes with rectangular elements : Poiseuille flow with cartesian mesh');

clear pde; clear bc;

pde = init_pde();

pde('type') = 'stokes'; // YC:
pde('asm')('viscos') = 1.0;
pde('viscos') = 1.0;
pde('bound') = list();
pde('bound')(1) = mlist(['bound','type','R','H']);
pde('bound')(1)('type') = 'Dirichlet';
pde('bound')(1)('R')    = list('y.*(y-1)',0);
pde('bound')(1)('H') = list([1 0; 0 1]);

m = gf_mesh('cartesian',[0:.3:5],[0:.2:1]);

pde('mf_u') = [];
pde('mf_p') = [];
pde('mf_d') = [];
pde('mim')  = [];

pde('mf_u') = gf_mesh_fem(m,2); // U mesh_fem (vector field -> qdim=2)
pde('mf_p') = gf_mesh_fem(m,1); // Pression mesh_fem
pde('mf_d') = gf_mesh_fem(m,1); // Data mesh_fem (boundary conditions, source terms etc)
pde('mim')  = gf_mesh_im(m,  gf_integ('IM_EXACT_PARALLELEPIPED(2)'));

gf_mesh_fem_set(pde('mf_u'),'fem',gf_fem('FEM_QK(2,2)'));
gf_mesh_fem_set(pde('mf_d'),'fem',gf_fem('FEM_QK(2,1)'));
gf_mesh_fem_set(pde('mf_p'),'fem',gf_fem('FEM_QK(2,1)'));

s = gf_mesh_fem_get(pde('mf_u'),'char');
pde('mf_u')  = gf_mesh_fem('from string',s,m);
all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
gf_mesh_set(m, 'boundary', 1, all_faces);

[U,P,pde] = gf_solve(pde);

pde

h = scf();
h.color_map = jetcolormap(255);
drawlater;
subplot(3,1,1); 
gf_plot(pde('mf_u'),U(:)','dir','x','deformation',U,'deformation_scale',0.1,'deformed_mesh','on'); 
colorbar(min(U),max(U));

subplot(3,1,2); 
gf_plot(pde('mf_p'),P(:)','deformation',U,'deformation_mf',pde('mf_u')); 
colorbar(min(P),max(P));

subplot(3,1,3); 
gf_plot(pde('mf_u'),U(:)','mesh','on'); 
gf_plot(pde('mf_p'),P(:)','refine',1); 
colorbar(min(P),max(P)); //YC: U or P ?
h.color_map = jetcolormap(255);
drawnow;

disp('Note that the dirichlet condition was described on a P1 fem');
disp('(visible on the deformed mesh: on boundaries, the deformation');
disp('is linear, hence there is a small error on the computed solution.');
