m  = gf_mesh('cartesian',[0:.1:1],[0:.1:1]); 
printf('SCIGETFEM: m'); disp(m);
mf = gf_mesh_fem(m,1); // create a meshfem of for a field of dimension 1 (i.e. a scalar field)
printf('SCIGETFEM: mf'); disp(mf);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));
sl = gf_slice(list('planar',0,[.5;0],[1;0]),m,3);
[P,T1,T2] = gf_slice_get(sl,'edges');
printf('SCIGETFEM: sl - P'); disp(P)
printf('SCIGETFEM: sl - T1'); disp(T1)
printf('SCIGETFEM: sl - T2'); disp(T2)

pstr = gf_fem_get(gf_fem('FEM_QK(2,2)'), 'poly_str');
printf('SCIGETFEM: pstr'); disp(pstr);

mim = gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));
printf('SCIGETFEM: mim'); disp(mim);
border = gf_mesh_get(m,'outer faces');
printf('SCIGETFEM: border'); disp(border);
gf_mesh_set(m, 'region', 42, border); // create the region

// the boundary edges appears in red
h = scf();
drawlater;
gf_plot_mesh(m, 'regions', [42], 'vertices','on','convexes','on'); 
drawnow;

md = gf_model('real');
printf('SCIGETFEM: md'); disp(md);
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');

//R = gf_mesh_fem_get(mf, 'eval', list('(x-.5).^2 + (y-.5).^2 + x/5 - y/3'));
R = gf_mesh_fem_get_eval(mf, list('(x-.5).^2 + (y-.5).^2 + x/5 - y/3'));
printf('SCIGETFEM: R'); disp(R);
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, R);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');
gf_model_get(md, 'listvar');
gf_model_get(md, 'solve');

U = gf_model_get(md, 'variable', 'u');
printf('SCIGETFEM: U'); disp(U);

h = scf();
h.color_map = jetcolormap(255);
drawlater;
gf_plot(mf, U, 'mesh','on');
h.color_map = jetcolormap(255);
drawnow;
