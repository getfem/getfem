lines(0);
stacksize('max');

path = get_absolute_file_path('demo_tripod_anim.sce');

printf('demo tripod_anim started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

// You should run demo_tripod first ...
//m    = gf_mesh('import','gid', path + '/data/tripod.GiD.msh');
//mfu  = gf_mesh_fem('from string', smfu, m);
//mfdu = gf_mesh_fem('from string', smfdu, m);

//drawlater;
//gf_plot_mesh(m,'cvlst',gf_mesh_get(m,'outer faces'),'curved','on','edges_color',[1 0 0]);
//drawnow;
pr   = 1;
haut = 0;

h = scf();
h.color_map = jetcolormap(255);

Index = 0;

//for r=[6 14 26]
for r=6:4:60
  disp('slicing...'); tic;
  //sl = gf_slice(list('cylinder', 0, [0;0;0], [0;1;0], r), m, 4);
  //sl = gf_slice(m,list('boundary',list('cylinder', [0;0;0], [0;1;0], 15)),4);

  //sl = gf_slice(m,list('union',list('cylinderb', [0;0;0], [0;1;0], 15),list('cylinderb', [0;0;0], [0;1;0], 7)),4);
  sl = gf_slice(list('boundary',list('diff',list('cylinder', -1, [0;0;0], [0;1;0], r),list('cylinder', -1, [0;0;0], [0;1;0], r-4))),m, 6);
  //sl = gf_slice(m,list('none'),4, gf_mesh_get(m,'outer faces'));
  //sl = gf_slice(m,list('boundary', list('none')),4);
  //sl = gf_slice(m,list('ballb', [0;0;0], 10),2);
  //sl = gf_slice(m,list('planarb',[0;0;0],[0;0;1]),1);
  printf('..........done in %3.2f sec',toc());

  P      = gf_slice_get(sl,'pts'); 
  P(2,:) = P(2,:) - haut;
  sl
  D = gf_compute(mfdu,VM,'interpolate on',sl);
  
  drawlater;
  gf_plot_slice(sl, 'mesh','on','data',D,'pcolor','on','mesh_edges_color',[1 1 .7]);
  h.color_map = jetcolormap(255);

  a = gca();
  a.view = '3d';
  a.data_bounds = [-30 -15 -50;
                    50  15  50];

  xlabel('');
  ylabel('');
  zlabel('');

  a.axes_visible = ['off','off','off'];
  a.box = 'off';
  drawnow;

  // use:
  // convert -delay 50 -loop 0 wave*.png animatewave.gif
  // To produce the animated gif image.
  // Convert is an ImageMagick tool.
  xs2png(h.figure_id, path + sprintf('/tripod%02d.png',Index));
  
  Index = Index + 1;
  pr   = r;
  haut = haut + 24;
end

printf('demo tripod_anim terminated\n');
