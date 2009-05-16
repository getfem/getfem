% Example of automatic refinement of the mesh
% In this example, the refinement will focus on the
% transition between the Dirichlet and the Neumann boundary.

gf_workspace('clear all'); 
%clear all; clf;
L=100; H=22;
N=2;
if (N == 2), % 2D beam
  m=gfMesh('regular simplices',0:10:L, 0:11:H);
  mim=gfMeshIm(m);    set(mim, 'integ',gfInteg('IM_TRIANGLE(6)'));
  mfu=gfMeshFem(m,N); set(mfu, 'fem',gfFem('FEM_PK(2,2)'));
  mfd=gfMeshFem(m);   set(mfd, 'fem',gfFem('FEM_PK(2,1)'));
  mf0=gfMeshFem(m);   set(mf0, 'fem',gfFem('FEM_PK(2,0)'));
  mfdu=gfMeshFem(m);  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(2,2)'));
else         % 3D beam
  m=gfMesh('regular simplices',0:10:L, 0:11:H, 0:11:H);
  mim=gfMeshIm(m);    set(mim, 'integ',gfInteg('IM_TETRAHEDRON(5)'));
  mfu=gfMeshFem(m,N); set(mfu, 'fem',gfFem('FEM_PK(3,2)'));
  mfd=gfMeshFem(m);   set(mfd, 'fem',gfFem('FEM_PK(3,1)'));
  mf0=gfMeshFem(m);   set(mf0, 'fem',gfFem('FEM_PK(3,0)'));
  mfdu=gfMeshFem(m);  set(mfdu,'fem',gfFem('FEM_PK_DISCONTINUOUS(3,1)'));
end;

lambda=121150; mu=80769;

P=get(m,'pts');
fleft =gf_mesh_get(m,'faces from pid',find(abs(P(1,:))<1e-6));
fright=gf_mesh_get(m,'faces from pid',find(abs(P(1,:) - L)<1e-6));
% assign boundary numbers
gf_mesh_set(m,'boundary',1,fleft);
gf_mesh_set(m,'boundary',2,fright);

F=zeros(N,1); F(2) = -20; % the external force


md=gf_model('real');
gf_model_set(md, 'add fem variable', 'u', mfu);
gf_model_set(md, 'add initialized data', 'lambda', [lambda]);
gf_model_set(md, 'add initialized data', 'mu', [mu]);
gf_model_set(md, 'add isotropic linearized elasticity brick', ...
	     mim, 'u', 'lambda', 'mu');
gf_model_set(md, 'add initialized data', 'VolumicData', F);
gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');
gf_model_set(md, 'add Dirichlet condition with multipliers', ...
	     mim, 'u', mfu, 1);














% b0=gfMdBrick('isotropic_linearized_elasticity',mim, mfu);
% b1=gfMdBrick('dirichlet',b0,1,mfu,'penalized');
% b2=gfMdBrick('source term',b1,2);

% set(b0, 'param', 'lambda', lambda);
% set(b0, 'param', 'mu', mu);

% mds=gfMdState(b2)

for step=1:8,
  dd=get(mf0, 'basic dof from cvid');
  
  % set(b2, 'param','source_term', mfd, F);

  gf_model_get(md, 'solve');
  U = gf_model_get(md, 'variable', 'u');

  % get(b2, 'solve', mds, 'very noisy'); %, 'lsolver', 'superlu');
  
  % U=get(mds, 'state'); U=U(1:get(mfu, 'nbdof'));
  
  VM = gf_model_get(md, 'compute isotropic linearized Von Mises or Tresca', 'u', 'lambda', 'mu', mfdu);

  % VM = get(b0, 'von mises', mds, mfdu);

  subplot(2,1,1);
  if (N==3) opt = {'cvlst', get(m,'outer_faces')}; 
  else opt = {}; end;
  gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,...
	  'deformation_mf',mfu,'refine', 4, 'deformation_scale',1, opt{:}); 
  gf_colormap('chouette');
  caxis([0 1e7]); colorbar; 
  title('Von Mises stress');
  
  ERR=gf_compute(mfu,U,'error estimate', mim);
  E=ERR; E(dd)=ERR;
  subplot(2,1,2);
  gf_plot(mf0, E, 'mesh','on', 'refine', 1, opt{:}); colorbar;
  title('Error estimate')
  disp('press a key..'); pause;
  set(m, 'refine', find(ERR > 1e-3));
  set(m, 'optimize structure');
end;
