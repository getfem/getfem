disp('This demo use levelset to impose (weakly) a Dirichlet condition on an');
disp('implicit boundary defined by the zero of the levelset');

%clear all;
gf_workspace('clear all');
NX=18;
ls_degree = 2;



m=gfMesh('cartesian', -.5:(1/NX):.5, -.5:(1/NX):.5);
%m=gfMesh('triangles grid', -.5:(1/NX):.5, -.5:(1/NX):.5);
ls=gfLevelSet(m, ls_degree);
ls2=gfLevelSet(m, ls_degree, 'with_secondary');

mf_ls=gfObject(get(ls, 'mf'))
mf_ls2=gfObject(get(ls2, 'mf'))

gf_workspace('stats')

gf_workspace('push');
gf_workspace('pop');

P=get(mf_ls, 'dof nodes');
x = P(1,:); y = P(2,:);



%ULS = ((x + 0.25).^2 + (y - 0.4).^2) - 0.05^2;
%ULS = min(ULS, ((x - 0.25).^2 + (y - 0.4).^2) - 0.05^2);

ULS=1000*ones(1,numel(x));
rand('state',1);

if 0,
  for ix=1:5,
    for iy=1:5,
      xc = ((ix-1)/4) * 0.8 - 0.4;
      yc = ((iy-1)/4) * 0.8 - 0.4;
      if (mod(iy,2)==0),
	xc=xc + 0.05;
      else
	xc=xc - 0.05;
      end;
      R = 0.03 + 0.005*(iy-1);
      ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);
    end;
  end;
else
  for i=1:8,
    xc = rand() - 0.5;
    yc = rand() - 0.5;
    R = rand() * 0.09 + 0.02;
    ULS = min(ULS, ((x - xc).^2 + (y - yc).^2) - R^2);
  end;
end;

set(ls, 'values', ULS);
%set(ls, 'values', 'x-.04');
%x = P(1,:); y = P(2,:);
%P2=get(mf_ls2, 'dof nodes', 
ULS2=1000*ones(1,numel(x));
ULS2s=1000*ones(1,numel(x));
for i=1:1,
  xc = 0; %rand() - 0.5;
  yc = 0.0; %rand() - 0.5;
  theta = pi/3; %pi*rand();
  n = [-sin(theta) cos(theta)];
  
  R = 0.19; %rand() * 0.09 + 0.02;
  ULS2 = min(ULS2, ((x-xc)*n(1) + (y-yc)*n(2)));
  %ULS2s = min(ULS2s, ((x - xc).^2 + (y - yc).^2) - R^2);
  ULS2s = min(ULS2s, (abs(y - yc)+abs(x-xc) - R));
end;

set(ls2, 'values', ULS2, ULS2s); %'-y-x+.2'); %, '(y-.2)^2 - 0.04');

%set(ls, 'simplify');

mls=gfMeshLevelSet(m);
set(mls, 'add', ls);
set(mls, 'add', ls2);
set(mls, 'adapt');

mim_bound = gfMeshIm('levelset',mls,'boundary(a+b)', gf_integ('IM_TRIANGLE(6)')); %, gf_integ('IM_QUAD(5)'));
mim = gfMeshIm('levelset',mls,'all(a+b)', gf_integ('IM_TRIANGLE(6)'));
set(mim, 'integ', 4);

mfu0=gfMeshFem(m,1); set(mfu0, 'fem', gf_fem('FEM_QK(2,3)'));
mfdu=gfMeshFem(m,1); set(mfdu, 'fem', gf_fem('FEM_QK_DISCONTINUOUS(2,2)'));
mf_mult0=gfMeshFem(m); set(mf_mult0, 'fem', gf_fem('FEM_QK(2,0)'));

%A=gf_asm('volumic','U=data(#1);V()+=comp(Base(#1))(i).U(i)',mim_bound,mfu,U)
A=gf_asm('volumic','V()+=comp()',mim_bound)

%clf; gf_plot_mesh(get(mls,'cut mesh'));
%return

%gf_plot_mesh(get(mls, 'cut_mesh'), 'curved', 'on');

%hold on; gf_plot(mf_ls, ULS);

dof_mult = get(mf_mult0, 'dof from im', mim_bound);

dof_out = get(mfu0, 'dof from im', mim);
cv_out = get(mim, 'convex_index');
cv_in = setdiff(get(m, 'cvid'), cv_out);
%[imlst,cv2im]=get(mim, 'integ');
%im_none=[]
%for i=1:numel(imlst),
%  if strcmp(gf_integ_get(imlst(i),'char'), 'IM_NONE'),
%    im_none=i;
%  end;
%end;
%if numel(im_none),
%  cv_in = find(cv2im==im_none);
%  cv_out = find(cv2im ~= im_none & cv2im > 0);  
%end;

mf_mult = gfMeshFem('partial', mf_mult0, dof_mult);
mfu = gfMeshFem('partial', mfu0, dof_out, cv_in);
set(mf_mult, 'qdim', 2);
set(mfu,'qdim',2);

B = gf_asm('mass matrix', mim_bound, mf_mult, mfu);


nbd=get(mfu, 'nbdof');

b0 = gfMdBrick('isotropic_linearized_elasticity',mim, mfu);
b1 = gfMdBrick('source term', b0);
set(b1, 'param', 'source term', [0; -9.81]);
b2 = gfMdBrick('constraint', b1, 'augmented');
set(b2, 'constraints', B, zeros(1,size(B,1)));

mds=gfMdState(b2)
get(b2, 'solve', mds, 'very noisy'); %, 'lsolver', 'superlu');
U=get(mds, 'state'); U=U(1:nbd);

VM = get(b0, 'von mises', mds, mfdu);
gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,...
	'deformation_mf',mfu,'refine', 8, 'cvlst', cv_out); 
%gf_plot(mfu,U,'norm','on','deformed_mesh','on', 'deformation',U,...
%	'deformation_mf',mfu,'refine', 8, 'cvlst', cv_out); 
hold on;
set(mfu,'qdim',1); Unorm=sqrt(U(1:2:end).^2 + U(2:2:end).^2);
[h1,h2]=gf_plot(mfu, Unorm,'contour',0.00001,'pcolor','off');
set(h2{1},'LineWidth',2);
set(h2{1},'Color','white');

[h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off');
set(h2{1},'LineWidth',1);
set(h2{1},'Color','black');
%[h1,h2]=gf_plot(mf_ls2, get(ls2,'values'), 'contour',
%0,'pcolor','off');

h2=line([xc + R*n(2); xc - R*n(2)],[yc - R*n(1), yc + R*n(1)]);
set(h2,'LineWidth',1);
set(h2,'Color','blue');



hold off; caxis([0 10]);
gf_colormap('chouette');