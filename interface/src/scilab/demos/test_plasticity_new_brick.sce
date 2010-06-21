gf_workspace('clear all');

//

// We try to compute a plasticity problem with a Von Mises crierion
// For convenience we consider an homogenous Dirichlet condition on the left
// of the domain and an easy computed Neumann Condition on the right

//

// Initialize used data
L = 100;
H = 20;
lambda = 121150;
mu = 80769;
von_mises_threshold = 8000;
//sigma_n = 1;
//u_n = 0;
f = [0 -200; 0 -300; 0 0];

// Create the mesh
m = gf_mesh('triangles grid', [0:4:L], [0:2:H]);

// Plotting
h_graph = scf();
h_graph.color_map = jetcolormap(256);
drawlater;
gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');
drawnow;

// Define used MeshIm
mim = gf_mesh_im(m);  gf_mesh_im_set(mim, 'integ', gf_integ('IM_TRIANGLE(6)')); // Gauss methods on triangles

// Define used MeshFem
mf_u = gf_mesh_fem(m,2); gf_mesh_fem_set(mf_u, 'fem', gf_fem('FEM_PK(2,2)'));
mf_data = gf_mesh_fem(m); gf_mesh_fem_set(mf_data, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,2)'));
mf_sigma = gf_mesh_fem(m); gf_mesh_fem_set(mf_sigma, 'fem', gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

// Find the border of the domain
P = gf_mesh_get(m, 'pts');
pidleft  = find(abs(P(1,:))<1e-6);
pidright = find(abs(P(1,:) - L)<1e-6);
fleft    = gf_mesh_get(m, 'faces from pid', pidleft); 
fright   = gf_mesh_get(m, 'faces from pid', pidright);

// Assign boundary numbers
gf_mesh_set(m,'boundary',1,fleft);  // for Dirichlet condition
gf_mesh_set(m,'boundary',2,fright); // for Neumann condition

// Plotting
//gf_plot_mesh(m, 'region', [1]);
//gf_plot_mesh(m, 'region', [2]);

// Create the model
md = gf_model('real');
//nbd_sigma = gf_model_get(mf_sigma, 'nbdof');
//nbd_data = gf_model_get(mf_data, 'nbdof');

// Declare that u is the unknown of the system on mf_u
// 2 is the number of version of the data stored, for the time integration scheme 
gf_model_set(md, 'add fem variable', 'u', mf_u, 2); 
//gf_model_set(md, 'variable', 'u', u_n, 1);

// Declare that lambda is a data of the system on mf_data
gf_model_set(md, 'add initialized data', 'lambda', lambda);

// Declare that mu is a data of the system on mf_data
gf_model_set(md, 'add initialized data', 'mu', mu);

// Declare that von_mises_threshold is a data of the system on mf_data
gf_model_set(md, 'add initialized data', 'von_mises_threshold', von_mises_threshold);

// Declare that sigma is a data of the system on mf_sigma
// 2 is the number of version of the data stored, for the time integration scheme
gf_model_set(md, 'add fem data', 'sigma', mf_sigma, 4, 2);
// Try to set an initial value to sigma
//gf_model_set(md, 'variable', 'sigma', sigma_n*ones(1, nbd_sigma), 1);

// Add plasticity brick on u
gf_model_set(md, 'add plasticity brick', mim, 'u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');

// Add homogeneous Dirichlet condition to u on the the left hand side of the domain
//u_Dir = gf_mesh_fem_get(mf_u, 'eval', 0);
//u_Dir = gf_mesh_fem_get_eval(mf_u, 0);
//gf_model_set(md, 'add initialized fem data', 'DirichletData', mf_u, u_Dir);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);

// Add a source term to the system
//gf_model_set(md,'add initialized fem data', 'VolumicData', mf_data, gf_mesh_fem_get(mf_data, 'eval',list(list(f(1,1)),list(f(2,1)))));
gf_model_set(md,'add initialized fem data', 'VolumicData', mf_data, gf_mesh_fem_get_eval(mf_data, list(list(f(1,1)),list(f(2,1)))));
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData', 2);

// Solve the system
gf_model_get(md, 'solve', 'very noisy', 'max_iter', 1000, 'max_res', 1e-6);

// Extract the solution
u = gf_model_get(md, 'variable', 'u', 1);

// Extract final sigma_np1
//sigma = gf_model_get(md, 'variable', 'sigma', 1);

// Display the solution u
drawlater;
gf_plot(mf_u, u, 'mesh', 'on');
drawnow;
