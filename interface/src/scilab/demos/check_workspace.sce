gf_workspace('clear all');
gf_workspace('stats');
gf_workspace('push');
 
m  = gf_mesh('empty',1);
mf = gf_mesh_fem(m);
gf_workspace('stats');
gf_workspace('pop');
gf_workspace('push','foo');

m  = gf_mesh('empty',2);
mf = gf_mesh_fem(m);
gf_workspace('keep',mf);
gf_workspace('pop');
gf_workspace('stats');
gf_delete(mf);
asserterr('gf_delete(mf)');

gf_workspace('stats');
gf_workspace('clear all');

