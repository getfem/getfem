// You should run demo_tripod first ...
//m    = gf_mesh('import','gid','../meshes/tripod.GiD.msh');
//mfu  = gf_mesh_fem('from string', smfu, m);
//mfdu = gf_mesh_fem('from string', smfdu, m);

//gf_plot_mesh(m,'cvlst',gf_mesh_get(m,'outer faces'),'curved','on','edges_color',[1 0 0]);
pr   = 1;
haut = 0;
for r=[6 14 26]
  disp('slicing...'); tic;
  //sl = gf_slice(list('cylinder', 0, [0;0;0], [0;1;0], r), m, 4);
  //sl = gf_slice(m,list('boundary',list('cylinder', [0;0;0], [0;1;0], 15)),4);

  //sl = gf_slice(m,list('union',list('cylinderb', [0;0;0], [0;1;0], 15),list('cylinderb', [0;0;0], [0;1;0], 7)),4);
  sl = gf_slice(list('boundary',list('diff',list('cylinder', -1, [0;0;0], [0;1;0], r),list('cylinder', -1, [0;0;0], [0;1;0], r-4))),m, 6);
  //sl = gf_slice(m,list('none'),4, gf_mesh_get(m,'outer faces'));
  //sl = gf_slice(m,list('boundary', list('none')),4);
  //sl = gf_slice(m,list('ballb', [0;0;0], 10),2);
  //sl = gf_slice(m,list('planarb',[0;0;0],[0;0;1]),1);
  disp(sprintf('..........done in %3.2f sec',toc));

  P      = gf_slice_get(sl,'pts'); 
  P(2,:) = P(2,:)-haut;
  sl
  D    = gf_compute(mfdu,VM,'interpolate on',sl);
  gf_plot_slice(sl, 'mesh','on','data',D,'pcolor','on','mesh_edges_color',[1 1 .7]);
  pr   = r;
  haut = haut+24;
end

gf_workspace('stats');
gf_colormap('earth');

