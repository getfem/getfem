lines(0);
stacksize('max');

path = get_absolute_file_path('demo_finallaplacian.sce');

printf('demo finallaplacian started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

// boundary names
up    = 101; // m.region(101)
down  = 102; // m.region(102)
left  = 103; // m.region(103)
right = 104; // m.region(104)

// importing the mesh and boundary
m = gf_mesh('import','gmsh',path + 'data/quad.msh');

// create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);

// assign the P1 fem to all convexes of the MeshFem
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,1)'));

// an exact integration will be used
mim = gf_mesh_im(m,gf_integ('IM_EXACT_SIMPLEX(2)'));

// interpolate the exact solution
U = gf_mesh_fem_get_eval(mf, list(list('x(1)*(x(1)-1)*x(2)*(x(2)-1)')));

// its Normal derivative on the domain boundary right (left is -DU)
DU = gf_mesh_fem_get_eval(mf, list(list('(2*x(1)-1)*x(2)*(x(2)-1)')));

// its laplacian
F = gf_mesh_fem_get_eval(mf, list(list('-(2*(x(1)**2+x(2)**2)-2*(x(1)+x(2)))')));

// generic elliptic brick
b0 = gf_mdbrick('generic_elliptic',mim,mf);

// add a Dirichlet condition on the domain boundary up
b1 = gf_mdbrick('dirichlet',b0,up,mf,'penalized');
gf_mdbrick_set(b1,'param','R',mf,U);

// add a Dirichlet condition on the domain boundary down
b2 = gf_mdbrick('dirichlet',b1,down,mf,'penalized');
gf_mdbrick_set(b1,'param','R',mf,U);

// add a source term
b3 = gf_mdbrick('source_term',b2);
gf_mdbrick_set(b3,'param','source_term',mf,F);

// add a Neumann condition on the domain boundary left
b4 = gf_mdbrick('source_term',b3,left);
gf_mdbrick_set(b4,'param','source_term',mf,-DU);

// add a Neumann condition on the domain boundary right
b5 = gf_mdbrick('source_term',b4,right);
gf_mdbrick_set(b5,'param','source_term',mf,DU);

// model state
mds = gf_mdstate('real');
gf_mdbrick_get(b5,'solve',mds);

// extracted solution
csol = gf_mdstate_get(mds,'state');
vsol = [0*csol,0*csol,csol];
vsol = matrix(vsol,3,length(vsol)/3);

// export
gf_mesh_fem_get(mf,'export_to_pos', path + '/sol.pos',U,'Exact solution', ...
                csol,'Computed solution', ...
                vsol,'Displacement', ...
                abs(csol-U),'abs difference');

printf('you can visualize the result using gmsh %s/sol.pos\n', path);

printf('demo finallaplacian terminated\n');
