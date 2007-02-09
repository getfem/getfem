function check_workspace(iverbose,idebug)
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
  gf_workspace('clear all');
  gf_workspace('stats');
  gf_workspace('push');
  m=gf_mesh('empty',1);
  mf=gf_mesh_fem(m);
  gf_workspace('stats');
  gf_workspace('pop');
  gf_workspace('push','foo');
  m=gf_mesh('empty',2);
  mf=gf_mesh_fem(m);
  gf_workspace('keep',mf);
  gf_workspace('pop');
  gf_workspace('stats');
  gf_delete(mf);
  asserterr('gf_delete(mf)');
  gf_workspace('stats');
  gf_workspace('clear all');
