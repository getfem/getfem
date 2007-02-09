if (exist('U')~=1 | exist('P') ~= 1),
  error('run demo_stokes_3D_tank2 first');
end;
clf;

% slice the mesh with two half spaces
sl=gf_slice({'boundary',{'intersection',{'planar',+1,[0;0;0],[0;1;0]},{'planar',+1,[0;0;0],[1;0;0]}}},m,6);
Usl=gf_compute(mfu,U,'interpolate on', sl);
Psl=gf_compute(mfp,P,'interpolate on', sl);
gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');

sl=gf_slice({'boundary',{'intersection',{'planar',+1,[0;0;6],[0;0;-1]},{'planar',+1,[0;0;0],[0;1;0]}}},m,6);
Usl=gf_compute(mfu,U,'interpolate on', sl);
Psl=gf_compute(mfp,P,'interpolate on', sl);
hold on;gf_plot_slice(sl,'mesh_faces','on','mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off'); hold off;
  
  
sl2=gf_slice({'boundary',{'planar',+1,[0;0;0],[0;1;0]}},m,6,setdiff(all_faces',TOPfaces','rows')');
hold on; gf_plot_slice(sl2,'mesh_faces','off','mesh','on','pcolor','off'); hold off;

% streamline "starting" points
hh=[1 5 9 12.5 16 19.5];
H=[zeros(2,numel(hh));hh];

% compute the streamlines
tsl=gf_slice('streamlines',mfu,U,H);
Utsl=gf_compute(mfu,U,'interpolate on', tsl);
% render them with "tube plot"
hold on; [a,h]=gf_plot_slice(tsl,'mesh','off','tube_radius',.2,'tube_color','white'); hold off;

% a nice colormap
caxis([0 .7]);
c=[0 0 1; 0 .5 1; 0 1 .5; 0 1 0; .5 1 0; 1 .5 0; 1 .4 0; 1 0 0; 1 .2 0; 1 .4 0; 1 .6 0; 1 .8 0];
colormap(c); camlight;