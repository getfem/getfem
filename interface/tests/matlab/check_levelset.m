gf_workspace('clear all');
clf;
m=gf_mesh('regular_simplices', -1:.2:1, -1:.2:1, 'degree', 2, 'noised');
%m=gf_mesh('cartesian', -1:.33:1, -1:.33:1);
ls=gfLevelSet(m, 2, 'x^2 + y^2 - 0.7^2', 'x-.4')
%ls=gfLevelSet(m, 2, 'x + y - 0.2'); %, 'x-5')
%ls=gfLevelSet(m, 2, 'x + y - 0.2', 'x-5')
ls2=gf_levelset(m, 2, '0.6*x^2 + (y-0.1)^2 - 0.6^2');
ls3=gf_levelset(m, 4, 'x^2 + (y+.08)^2 - 0.05^2');

mls=gfMeshLevelset(m)
set(mls, 'add', ls);
if 1,
  set(mls, 'del', ls);
  set(mls, 'add', ls);
  set(mls, 'add', ls2);
  set(mls, 'add', ls2);
  set(mls, 'add', ls2);
  set(mls, 'add', ls3);
end;
set(mls, 'adapt');

gfObject(get(mls, 'linked_mesh'))

lls = gf_mesh_levelset_get(mls, 'levelsets')

cm = gfObject(get(mls, 'cut_mesh'))

ctip = get(mls, 'crack_tip_convexes')


mf=gfMeshFem(m); set(mf, 'classical_fem', 1);
mfls=gfMeshFem('levelset',mf,mls);

%gf_workspace('stats');

nbd = get(mfls,'nbdof')
if 1,
  sl=gfSlice({'none'}, mls, 2);
  %for i=1:nbd,
%  U=zeros(1,nbd);U(i)=1;
    U=rand(1,nbd);
    gf_plot(mfls,U,'refine',4,'zplot','on');
    %pause;
%end;
  hold on;
  %gf_plot_mesh(cm, 'curved', 'on','refine',8,'edges_width',2);
  gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
  hold off;
  %caxis([0 2]); colorbar;
else
  for i=1:nbd,
    U=zeros(1,nbd); U(i)=1;
    gf_plot(mfls,U,'refine',16);
    hold on;
    gf_plot_mesh(cm, 'curved', 'on','refine',8);
    hold on;
    gf_plot_mesh(m, 'curved', 'on','refine',8, 'edges_color', [0 0 0]);
    hold off;
    pause
  end;
end;

