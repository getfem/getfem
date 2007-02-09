function check_oo(iverbose,idebug)
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
  m1=gfMesh('empty',1);
  assert('m1.nbpts==0');
  assert('m1.dim==1');
  p=[0 1 0 1.1; 0 0 1 1]; t = [1 2 3 0; 2 3 4 0]'; 
  m2=gf_mesh('pt2D',p,t); 
  m2=gfMesh(m2);
  assert('isa(m2,''gfMesh'')');
  m3=gfMesh('empty',3);
  set(m3,'add convex',gfGeoTrans('GT_QK(3,1)'),...
	 [0 1 0 1 0 1 0 1;...
	  0 0 1 1 0 0 1 1;...
	  0 0 0 0 1 1 1 1]);
  assert('m3.nbpts==8');
  assert('m3.pts(8)==[1;1;1]');
  assert('m3.pts([3 5])==[0 0; 1 0; 0 1]');
  asserterr('m3.pts(9)');
  asserterr('m3.pts(-1)');
  asserterr('m3.pts(0)');
  asserterr('m3.pts(''kjk'')');
  assert('numel(m3.pid_from_cvid(1))==8');
  assert('m2.nbcvs==2');
  gf_delete(m1);
  m1=gfMesh('cartesian',1:.1:5); 
  mf1=gfMeshFem(m1,2);
  mim=gfMeshIm(m1);
  assert('class(mf1)==''gfMeshFem''');
  assert('mf1.qdim==2');
  assert('mf1.mesh.dim==1');
  assert('mf1.mesh.pts(2)==1.1');
  asserterr('set(m1,''fem'',gfFem(''FEM_PK(1,2)''),gfInteg(''IM_EXACT_SIMPLEX(1)''))');
  set(mf1,'fem',gfFem('FEM_PK(1,2)')); 
  set(mim,'integ',gfInteg('IM_EXACT_SIMPLEX(1)'));
  assert('mf1.nbdof==162');
  e=get(mf1.mesh,'outer faces'); assert('e(1,:)==[1 40]');
  e=get(mf1.mesh,'outer faces', 1:mf1.mesh.nbcvs-2); assert('e(1,:)==[1 38]');
  p10=m1.pts(mf1.mesh.pid_from_cvid(10));
  set(mf1.mesh,'del convex',10);
  assert('m1.nbcvs==39');
  assert('mf1.nbdof==160'); % check the mesh_fem was correctly updated
  assert('isempty(m1.pid_from_cvid(10))');
  n=get(m1,'normal of face',9,1); 
  assert('abs(n-1)<1e-15');
  n=get(m1,'normal of face',9,2); assert('abs(n+1)<1e-15');
  asserterr('get(m1,''normal of face'', 10,1)');  
  set(mf1.mesh,'add convex',gfGeoTrans('GT_QK(1,1)'), p10);
  assert('mf1.mesh.nbcvs==40');
  assert('mf1.nbdof==160');
  s=char(mf1);
  assert('numel(s)>2000');
  gt=mf1.mesh.geotrans(3);
  cvs=mf1.mesh.cvstruct(1);
  fem=mf1.fem(2:3);
  integ=mim.integ(1);
