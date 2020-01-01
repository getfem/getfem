// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function check_oo(iverbose,idebug)
  [nargout,nargin] = argn();

  global gverbose;
  global gdebug;  
  if (nargin >= 1) then
    gverbose = iverbose;
    if (nargin == 2) then
      gdebug = idebug;
    else 
      gdebug = 0; end;
  else 
    gverbose = 0;
  end

  gf_workspace('clear all');

  m1 = gfMesh('empty',1);
  gfassert('m1.nbpts==0');
  gfassert('m1.dim==1');
  p  = [0 1 0 1.1; 0 0 1 1]; t = [1 2 3 0; 2 3 4 0]'; 
  m2 = gf_mesh('pt2D',p,t);
  m2 = gfMesh(m2);
  gfassert('gf_typeof(m2)==''gfMesh''');
  m3 = gfMesh('empty',3);
  set(m3,'add convex',gfGeoTrans('GT_QK(3,1)'),...
	 [0 1 0 1 0 1 0 1;...
	  0 0 1 1 0 0 1 1;...
	  0 0 0 0 1 1 1 1]);
  gfassert('m3.nbpts==8');
  gfassert('m3.pts(8)==[1;1;1]');
  gfassert('m3.pts([3 5])==[0 0; 1 0; 0 1]');
  asserterr('m3.pts(9)');
  asserterr('m3.pts(-1)');
  asserterr('m3.pts(0)');
  asserterr('m3.pts(''kjk'')');
  gfassert('length(m3.pid_from_cvid(1))==8');
  gfassert('m2.nbcvs==2');
  gf_delete(m1);
  m1  = gfMesh('cartesian',1:.1:5); 
  mf1 = gfMeshFem(m1,2);
  mim = gfMeshIm(m1);
  gfassert('class(mf1)==''gfMeshFem''');
  gfassert('mf1.qdim==2');
  gfassert('mf1.mesh.dim==1');
  gfassert('mf1.mesh.pts(2)==1.1');
  asserterr('set(m1,''fem'',gfFem(''FEM_PK(1,2)''),gfInteg(''IM_EXACT_SIMPLEX(1)''))');
  set(mf1,'fem',gfFem('FEM_PK(1,2)')); 
  set(mim,'integ',gfInteg('IM_EXACT_SIMPLEX(1)'));
  gfassert('mf1.nbdof==162');
  e   = get(mf1.mesh,'outer faces'); gfassert('e(1,:)==[1 40]');
  e   = get(mf1.mesh,'outer faces', 1:mf1.mesh.nbcvs-2); gfassert('e(1,:)==[1 38]');
  p10 = m1.pts(mf1.mesh.pid_from_cvid(10));
  set(mf1.mesh,'del convex',10);
  gfassert('m1.nbcvs==39');
  gfassert('mf1.nbdof==160'); // check the mesh_fem was correctly updated
  gfassert('isempty(m1.pid_from_cvid(10))');
  n = get(m1,'normal of face',9,1); 
  gfassert('abs(n-1)<1e-15');
  n = get(m1,'normal of face',9,2); gfassert('abs(n+1)<1e-15');
  asserterr('get(m1,''normal of face'', 10,1)');  
  set(mf1.mesh,'add convex',gfGeoTrans('GT_QK(1,1)'), p10);
  gfassert('mf1.mesh.nbcvs==40');
  gfassert('mf1.nbdof==162');
  s = char(mf1);
  gfassert('numel(s)>1700');
  gt    = mf1.mesh.geotrans(3);
  cvs   = mf1.mesh.cvstruct(1);
  fem   = mf1.fem(2:3);
  integ = mim.integ(1);
