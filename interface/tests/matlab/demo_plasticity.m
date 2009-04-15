L=100; H=20;
m=gfMesh('triangles grid',0:4:L, 0:2:H);


mim=gfMeshIm(m);  set(mim, 'integ',gfInteg('IM_TRIANGLE(6)'));
mfu=gfMeshFem(m,2); set(mfu, 'fem',gfFem('FEM_PK(2,1)'));
mfd=gfMeshFem(m); set(mfd, 'fem',gfFem('FEM_PK(2,1)'));
mf0=gfMeshFem(m); set(mf0, 'fem',gfFem('FEM_PK(2,0)'));
mfdu=gfMeshFem(m); set(mfdu, 'fem',gfFem('FEM_PK_DISCONTINUOUS(2,1)'));
lambda=121150;
mu=80769;
von_mises_threshold=8000;

P=get(m, 'pts');
pidleft=find(abs(P(1,:))<1e-6);
pidright=find(abs(P(1,:) - L)<1e-6);

fleft =get(m,'faces from pid',pidleft); 
fright=get(m,'faces from pid',pidright);
% assign boundary numbers
set(m,'boundary',1,fleft);
set(m,'boundary',2,fright);

b0=gfMdBrick('small deformations plasticity',mim,mfu, von_mises_threshold);
set(b0, 'param','lambda',lambda);
set(b0, 'param','mu',mu);
b1=gfMdBrick('generalized dirichlet',b0,1);
b2=gfMdBrick('source term',b1,2);

mds=gfMdState(b2)


VM=zeros(1,get(mfdu, 'nbdof'));

F=[0 -200; 0 -300; 0 0]';
nbstep = size(F,2);

dd=get(mf0, 'basic dof from cvid');

for step=1:nbstep,

  
  %b2.set('param','source_term', mfd, get(mfd, 'eval',{0;-400*sin(step*pi/2)}));
  set(b2, 'param','source_term', mfd, get(mfd, 'eval',{F(1,step);F(2,step)}));
  get(b2, 'solve', mds, 'very noisy', 'max_iter', 1000, 'max_res', 1e-6);
  
  
  U=get(mds, 'state'); U=U(1:get(mfu, 'nbdof'));
  VM = get(b0, 'von mises', mds, mfdu);
  max(abs(VM))
  
  subplot(2,1,1);
  gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,'deformation_mf',mfu,'refine', 4, 'deformation_scale',1); 
  colorbar;
  caxis([0 10000]);
  
  ERR=gf_compute(mfu,U,'error estimate', mim);
  E=ERR; E(dd)=ERR;
  subplot(2,1,2);
  gf_plot(mf0, E, 'mesh','on', 'refine', 1); colorbar;
  pause;
end;
