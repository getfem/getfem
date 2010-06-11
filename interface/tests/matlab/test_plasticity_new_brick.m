clear all;
gf_workspace('clear all');

%%


% We try to compute a plasticity problem with a Von Mises crierion
% For convenience we consider an homogenous Dirichlet condition on the left
% of the domain and an easy computed Neumann Condition on the right


%%



% Initialize used data
L = 100;
H = 20;
lambda = 121150;
mu = 80769;
von_mises_threshold = 8000;
%sigma_n = 1;
%u_n = 0;
f = [0 -200; 0 -300; 0 0]


% Create the mesh
m = gfMesh('triangles grid', [0:4:L], [0:2:H]);

% Plotting
gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

% Define used MeshIm
mim=gfMeshIm(m);  set(mim, 'integ',gfInteg('IM_TRIANGLE(6)')); % Gauss methods on triangles

% Define used MeshFem
mf_u=gfMeshFem(m,2); set(mf_u, 'fem',gfFem('FEM_PK(2,1)'));
mf_data=gfMeshFem(m); set(mf_data, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,0)'));
mf_sigma=gfMeshFem(m); set(mf_sigma, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,0)'));


% Find the border of the domain
P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6);
pidright=find(abs(P(1,:) - L)<1e-6);
fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);

% Assign boundary numbers
set(m,'boundary',1,fleft); % for Dirichlet condition
set(m,'boundary',2,fright); % for Neumann condition

% Plotting
%gf_plot_mesh(m, 'region', [1]);
%gf_plot_mesh(m, 'region', [2]);

% Create the model
md = gfModel('real');
nbd_sigma=get(mf_sigma, 'nbdof');
nbd_data=get(mf_data, 'nbdof');

% Declare that u is the unknown of the system on mf_u
% 2 is the number of version of the data stored, for the time integration scheme 
set(md, 'add fem variable', 'u', mf_u, 2); 
%set(md, 'variable', 'u', u_n, 1);

% Declare that lambda is a data of the system on mf_data
set(md, 'add initialized data', 'lambda', lambda*ones(1,nbd_data));

% Declare that mu is a data of the system on mf_data
set(md, 'add initialized data', 'mu', mu*ones(1,nbd_data));

% Declare that von_mises_threshold is a data of the system on mf_data
set(md, 'add initialized data', 'von_mises_threshold', von_mises_threshold*ones(1,nbd_data));

% Declare that sigma is a data of the system on mf_sigma
% 2 is the number of version of the data stored, for the time integration scheme
set(md, 'add fem data', 'sigma', mf_sigma, 1, 2);
% Try to set an initial value to sigma
%set(md, 'variable', 'sigma', sigma_n*ones(1, nbd_sigma), 1);

% Add plasticity brick on u
set(md, 'add plasticity brick', mim, 'u', 'lambda', 'mu', 'von_mises_threshold', 'sigma');

% Add homogeneous Dirichlet condition to u on the the left hand side of the domain
%u_Dir = get(mf_u, 'eval', 0);
%set(md, 'add initialized fem data', 'DirichletData', mf_u, u_Dir);
set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf_u, 1);

% Add a source term to the system
set(md,'add initialized fem data', 'VolumicData', mf_data, get(mf_data, 'eval',{f(1,1);f(2,1)}));
set(md, 'add source term brick', mim, 'u', 'VolumicData',2);


% Solve the system
%get(md,'solve');% 'solve', 'very noisy', 'max_iter', 1000, 'max_res', 1e-6);

% Extract the solution
%u = get(md, 'variable', 'u', 1);

% Extract final sigma_np1
%sigma = get(md, 'variable', 'sigma', 1);

% Display the solution u
%gf_plot(mf_u, u, 'mesh', 'on');











