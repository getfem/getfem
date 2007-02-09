function res=check_solve(iverbose,idebug)
% check: gf_mesh('cartesian_mesh'), gf_mesh('triangles grid'), 
% gf_compute('interpolate on Q1 grid')
% gf_compute('gradient'), gf_mesh_fem_set('boundary'), gf_solve(stokes)
% gf_compute('interpolate on'), discontinuous fems, L2&H1 norms
  
  global gverbose;
  global gdebug;  
  if (nargin >= 1),
    gverbose = iverbose;
    if (nargin == 2),
      gdebug = idebug;
    else 
      gdebug = 0; end;
  else 
    gverbose = 0;
  end;
  tic;
  gf_workspace('clear all');
  m1=gf_mesh('empty',1);
  pde.type = 'stokes';
  pde.viscos=1.0;
  pde.bound{1}.type = 'Dirichlet';
  pde.bound{1}.R  = {'-y.*(y-1)',0};
  m=gf_mesh('cartesian',[0:.3:4.8],[0:.2:1]);
  pde.mf_u=gf_mesh_fem(m,2);
  pde.mf_p=gf_mesh_fem(m,1);
  pde.mf_d=gf_mesh_fem(m,1);
  pde.mim=gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));
  mf_comp =gf_mesh_fem(m,2);
  gf_mesh_fem_set(pde.mf_u,'fem',gf_fem('FEM_QK(2,2)'));
  gf_mesh_fem_set(mf_comp ,'fem',gf_fem('FEM_QK(2,2)'));
  % the piecewise linear mf_d will induce a small error on the solution
  % since the dirichlet condition is parabolic
  gf_mesh_fem_set(pde.mf_d,'fem',gf_fem('FEM_QK(2,1)'));
  gf_mesh_fem_set(pde.mf_p,'fem',gf_fem('FEM_QK(2,1)'));
  all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
  for mf=[pde.mf_u pde.mf_p pde.mf_d],
    gf_mesh_set(m, 'boundary', 1, all_faces);
  end;
  [U,P]=gf_solve(pde);  
  Uco=gf_compute(pde.mf_u, U, 'interpolate on', mf_comp); % tests interpolation on same mesh
  dof=gf_mesh_fem_get(mf_comp, 'dof nodes'); Xc=dof(1,1:2:end); Yc=dof(2,1:2:end);
  Uex=[-Yc.*(Yc-1); zeros(1,numel(Xc))]; Uex=Uex(:)';

  % norm tests
  l2m = gf_compute(mf_comp,Uex,'L2 norm',pde.mim);
  assert('abs(l2m-.4)<1e-10');
  disp(sprintf('L2 norm %f',l2m));
  l2err = gf_compute(mf_comp,Uex-Uco,'L2 norm',pde.mim);
  disp(sprintf('L2 err %f', l2err));
  h1err = gf_compute(mf_comp,Uex-Uco,'H1 norm',pde.mim);
  disp(sprintf('H1 err %f', h1err));
  disp(sprintf('done in %.2f sec.',toc));
  %gf_plot(mf_comp,Uex-Uco,'norm'); colorbar;
  assert('abs(l2err)<0.016'); % 0.015926
  assert('abs(h1err)<0.0665'); % 0.066215
  
  [Uq,Iq,mfq]=gf_compute(pde.mf_u, U, 'interpolate on Q1 grid', 'regular h', [.05, .05]);
  [XX,YY]=meshgrid(0:.05:4.8,0.:0.05:1);  XX=XX'; YY=YY';
  UU=-YY.*(YY-1);
  assert('max(max(abs(UU-squeeze(Uq(1,:,:))))) < 0.01001');
  assert('max(max(abs(Uq(2,:,:)))) < 0.0030');
  
  %gradient check
  mtri=gf_mesh('triangles grid',[0 .2 .4 .8 1:.3:4.8],[0:.2:.6 .9 1]);
  mf_DU =gf_mesh_fem(mtri,2);
  gf_mesh_fem_set(mf_DU,'fem',gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));
  mim2=gf_mesh_im(mtri,2);
  asserterr('gf_compute(pde.mf_u, U, ''gradient'', mf_DU)'); % cause different meshes;
  UU=gf_compute(pde.mf_u, U, 'interpolate on', mf_DU); % tests interpolation on other mesh
  DU = gf_compute(mf_DU, UU, 'gradient', mf_DU);
  dof=gf_mesh_fem_get(mf_DU, 'dof nodes'); Xc=dof(1,1:2:end); Yc=dof(2,1:2:end);
  DUex=[1-2*Yc; zeros(1,numel(Xc))]; DUex=DUex(:)';
  diff=norm(DUex-DU(2,:));
  assert('diff>4.62 & diff<4.64');
  
  % yes the error on the derivative is quite big. This is because we interpolated
  % U on mf_DU which is piecewise linear
  % the 3 plots below illustrate this
  %subplot(3,1,1); gf_plot(mf_DU, DU(2,:),'mesh','x'); colorbar; %dUx/dy
  %subplot(3,1,2); gf_plot(mf_DU, DUex,'mesh','x'); colorbar;
  %subplot(3,1,3); gf_plot(mf_DU, DUex-DU(2,:),'mesh','x'); colorbar;    
  d2=gf_compute(mf_DU,DUex-DU(2,:),'L2 norm',mim2);
  assert('d2>0.28 & d2 < 0.29'); % 0.2866
  
  %gradient of vector fields
  DU2 = gf_compute(mf_DU, [UU;2*UU;3*UU], 'gradient', mf_DU);
  d=permute(cat(3,DU,2*DU,3*DU),[1,3,2]);
  d(1:10)
  DU2(1:10)
  
  assert('max(abs(DU2(:)-d(:))) < 2e-15');
  