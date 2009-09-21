disp('this file should be launched after demo_tripod.m as it assumes the tripod mesh and solutions are in memory')

//m    = gf_mesh('from string',sm);
//mfu  = gf_mesh_fem('from string',smfu,m);
//mfdu = gf_mesh_fem('from string',smfdu,m);

disp('plotting ... can also take some minutes!');

f = scf();
c = [0 0 .5; 0 .2 1; 0 .5 .8; 0 .9 0; .4 1 0; .9 .7 0; .8 0 0; 1 0 0];  
f.color_map = c;
cnt = 1;
for r=-10.3:+.1:12 //46.1:-.1:4,
  clf;
  //sl=gf_slice(mfu,U*10,{'boundary',{'cylinder',-1,[0;0;0],[0;1;0],r}},5);
  sl  = gf_slice(list('boundary',list('planar',-1,[0;r;0],[0;1;0])),mfu,U*10,5);
  Usl = gf_compute(mfdu,VM,'interpolate on',sl);
  P   = gf_slice_get(sl,'pts'); P=P([1 3 2],:); gf_slice_set(sl,'pts',P);
  drawlater;
  gf_plot_slice(sl,'data',Usl,'mesh','on','mesh_slice_edges_color',[.7 .7 .7],'mesh_edges_color',[.5 .5 1]);
  drawnow;
  //view(0,30);
  //caxis([0 7]);
  //axis([-48 48 -5 15 -48 48]); axis off; camzoom(1.4*1.3); campan(+.8,-.2);
  //camlight;
  xs2png(sprintf('tripod_slice_p%03d',cnt));
  cnt = cnt+1;
  sleep(1000)
  gf_delete(sl);
end

