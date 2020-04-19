% creation of a simple cartesian mesh
m = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);

% we enable vertices and convexes labels
gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

% assign the same integration method on all convexes
mim=gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));

% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #42
gf_mesh_set(m, 'region', 42, border);
gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red

% empty real model
md = gf_model('real');

% declare that "u" is an unknown of the system
% on the finite element method `mf`
gf_model_set(md, 'add fem variable', 'u', mf);

% add generic elliptic brick on "u"
gf_model_set(md, 'add Laplacian brick', mim, 'u');

% add Dirichlet condition
Uexact = gf_mesh_fem_get(mf, 'eval', {'(x-.5).^2 + (y-.5).^2 + x/5 - y/3'});
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');

% add source term
f = gf_mesh_fem_get(mf, 'eval', { '2(x^2+y^2)-2(x+y)+20x^3' });
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, f);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

% solve the linear system
gf_model_get(md, 'solve');

% extracted solution
u = gf_model_get(md, 'variable', 'u');
% display
gf_plot(mf, u, 'mesh','on');
