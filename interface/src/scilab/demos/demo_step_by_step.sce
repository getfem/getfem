lines(0);
stacksize('max');

path = get_absolute_file_path('demo_step_by_step.sce');

printf('demo step_by_step started\n');

gf_workspace('clear all');

// creation of a simple cartesian mesh
m = gf_mesh('cartesian', 0:0.1:1.1, 0:0.1:1.1);

// create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m, 1);

// assign the Q2 fem to all convexes of the MeshFem
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

// view the expression of its basis functions on the reference convex
printf('The expression of its basis functions on the reference convex\n');
disp(gf_fem_get(gf_fem('FEM_QK(2,2)'),'poly_str')); 

// an exact integration will be used
mim = gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));

// detect the border of the mesh
border = gf_mesh_get(m,'outer_faces');

// mark it as boundary #42
gf_mesh_set(m,'region',42, border);

// generic elliptic brick
b0 = gf_mdbrick('generic elliptic', mim, mf);

// list of parameters of the brick
printf('List of parameters of the brick\n');
disp(gf_mdbrick_get(b0,'param_list'));

// add a Dirichlet condition on the domain boundary

b1 = gf_mdbrick('dirichlet', b0, 42, mf, 'penalized');

// change Dirichlet condition

R = gf_mesh_fem_get_eval(mf,list(list('x(1)*(x(1)-1) - x(2)*(x(2)-1)')));
gf_mdbrick_set(b1,'param','R', mf, R);

// created model state
mds = gf_mdstate('real');
gf_mdbrick_get(b1,'solve',mds);

// extracted solution
sol = gf_mdstate_get(mds,'state');

// export computed solution
gf_mesh_fem_get(mf,'export_to_pos','sol.pos',sol,'Computed solution');

printf('demo step_by_step terminated\n');
