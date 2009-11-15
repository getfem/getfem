% trace on;
gf_workspace('clear all');
m = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);
%m=gf_mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')

% create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf_mesh_fem(m,1);
% assign the Q2 fem to all convexes of the mesh_fem,
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

% Integration which will be used
mim = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,4)'));
%mim = gf_mesh_im(m, gf_integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,5),4)'));
% detect the border of the mesh
border = gf_mesh_get(m,'outer faces');
% mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);
gf_plot_mesh(m, 'regions', [1]); % the boundary edges appears in red
pause(1);

% interpolate the exact solution
Uexact = gf_mesh_fem_get(mf, 'eval', { 'y.*(y-1).*x.*(x-1)+x.^5' });
% its second derivative
F      = gf_mesh_fem_get(mf, 'eval', { '-(2*(x.^2+y.^2)-2*x-2*y+20*x.^3)' });


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mf);
gf_model_set(md, 'add Laplacian brick', mim, 'u');
gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 1, 'DirichletData');

gf_model_get(md, 'solve');
U = gf_model_get(md, 'variable', 'u');

% Version with old bricks
% b0=gf_mdbrick('generic elliptic',mim,mf);
% b1=gf_mdbrick('dirichlet', b0, 1, mf, 'penalized');
% gf_mdbrick_set(b1, 'param', 'R', mf, Uexact); 
% b2=gf_mdbrick('source term',b1);
% gf_mdbrick_set(b2, 'param', 'source_term', mf, F);
% mds=gf_mdstate(b1);
% gf_mdbrick_get(b2, 'solve', mds)
% U=gf_mdstate_get(mds, 'state');

disp(sprintf('H1 norm of error: %g', gf_compute(mf,U-Uexact,'H1 norm',mim)));

subplot(2,1,1); gf_plot(mf,U,'mesh','on','contour',.01:.01:.1); 
colorbar; title('computed solution');

subplot(2,1,2); gf_plot(mf,U-Uexact,'mesh','on'); 
colorbar;title('difference with exact solution');
