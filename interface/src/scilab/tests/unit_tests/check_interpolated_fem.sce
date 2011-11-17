gf_workspace('clear all');
lines(0);

h = scf();
h.color_map = jetcolormap(255);

m1 = gf_mesh('regular_simplices', 0:.5:2, 0:.4:1, 'degree', 2, 'noised');
//m1 = gf_mesh('regular_simplices', 0:1:2, 0:.5:1, 'degree', 2, 'noised');
drawlater;
gf_plot_mesh(m1, 'refine' ,5, 'curved','on');
drawnow;
mf1  = gf_mesh_fem(m1); 
mim1 = gf_mesh_im(m1, gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),4)'));
gf_mesh_fem_set(mf1, 'fem', gf_fem('FEM_PK(2, 1)'));
m2 = gf_mesh('regular_simplices', 0:.3:3, -.2:.4:1.2, 'degree', 1,'noised');
//m2 = gf_mesh('regular_simplices', [0 3], [0 1], 'degree', 1, 'noised');
//drawlater;
//gf_plot_mesh(m2, 'refine' ,5, 'curved','on');
//drawnow;
mf2  = gf_mesh_fem(m2); 
mim2 = gf_mesh_im(m2,gf_integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),4)'));
//mim2 = gf_meshI_im(m2, gf_integ('IM_TRIANGLE(6)'));
gf_mesh_fem_set(mf2, 'fem', gf_fem('FEM_PK(2, 1)'));
f = gf_fem('interpolated fem', mf1, mim2)
gf_mesh_fem_set(mf2, 'fem', f);
gf_workspace('stats');
mf3 = gf_mesh_fem(m2);
gf_mesh_fem_set(mf3, 'fem', gf_fem('FEM_PK(2,1)'));
gf_mesh_fem_set(mf3, 'fem', gf_fem('FEM_PK(2, 0)'), [1 2 3 5]);
mf4 = gf_mesh_fem('sum', mf2, mf3);
gf_mesh_set(m2, 'del convex', 4);
mf  = mf4; 
nbd = gf_mesh_fem_get(mf, 'nbdof');
drawlater;
gf_plot(mf, rand(1, nbd), 'refine', 16); // YC: There is a little plot bug here ...
drawnow;
//for i=1:nbd, 
//  U=zeros(1,nbd); U(i)=1;
//  disp(sprintf('dof %d/%d', i, nbd));
//  drawlater;
//  gf_plot(mf,U,'refine',16, 'mesh','on');
//  drawnow;
//  pause
//end;
gf_workspace('stats');
gf_delete(f);
gf_fem_get(f, 'char'); // YC: logic error here: f not found anymore
