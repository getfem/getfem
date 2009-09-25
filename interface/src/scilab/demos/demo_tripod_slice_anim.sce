disp('this file should be launched after demo_tripod.sce as it assumes the tripod mesh and solutions are in memory')

//m    = gf_mesh('from string',sm);
//mfu  = gf_mesh_fem('from string',smfu,m);
//mfdu = gf_mesh_fem('from string',smfdu,m);

disp('plotting ... can also take some minutes!');

h = scf();
c = [0.0 0.0 0.5;
     0.0 0.2 1.0;
     0.0 0.5 0.8;
     0.0 0.9 0.0;
     0.4 1.0 0.0;
     0.9 0.7 0.0;
     0.8 0.0 0.0;
     1.0 0.0 0.0];  
h.color_map = c;

cnt = 1;
for r=-10.3:+.1:12 //46.1:-.1:4,
  clf;
  //sl = gf_slice(mfu,U*10,list('boundary',list('cylinder',-1,[0;0;0],[0;1;0],r)),5);
  sl  = gf_slice(list('boundary',list('planar',-1,[0;r;0],[0;1;0])),mfu,U*10,5);
  Usl = gf_compute(mfdu,VM,'interpolate on',sl);
  P   = gf_slice_get(sl,'pts'); 
  P   = P([1 3 2],:); 
  gf_slice_set(sl,'pts',P);
  
  drawlater;
  gf_plot_slice(sl,'data',Usl,'mesh','on','mesh_slice_edges_color',[.7 .7 .7],'mesh_edges_color',[.5 .5 1]);
  h.color_map = c;
  drawnow;

  xs2png(sprintf('tripod_slice_p%03d',cnt));
  
  cnt = cnt+1;
  sleep(1000)
  gf_delete(sl);
end

