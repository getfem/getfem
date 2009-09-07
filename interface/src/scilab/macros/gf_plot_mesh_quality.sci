function gf_plot_mesh_quality(m)
try 
  gf_workspace('push');
  gf_plot_mesh_quality_aux(m);
catch
  disp('error in gf_plot_mesh_quality : ' + lasterror());
end;
gf_workspace('pop');
endfunction

function gf_plot_mesh_quality_aux(m)
clf
q     = gf_mesh_get(m,'quality');
cvlst = find(q < .1);
sl0   = gf_slice(list('boundary'),m,2); // gf_plot_slice(sl0,'mesh_edges_color',[.5 1 .5]);
//hold on;
sl = gf_slice(list('none'), m, 10, gf_mesh_get(m,'outer faces',cvlst));
mf = gf_mesh_fem(m); gf_mesh_fem_set(mf,'classical fem',0,1);
U  = zeros(1,gf_mesh_fem_get(mf,'nbdof'));
for cv=cvlst
  U(gf_mesh_fem_get(mf, 'dof from cv', cv)) = q(cv);
end
Q = gf_compute(mf,U,'interpolate on',sl);
gf_plot_slice(sl,'data',Q,'mesh_edges_color',[0 0 0],'mesh_edges','on');
gf_plot_slice(sl,'mesh_edges_color',[0 0 0],'mesh_edges','on','mesh_edges_width',2);
endfunction

