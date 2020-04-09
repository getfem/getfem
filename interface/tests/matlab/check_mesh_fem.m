% Copyright (C) 2005-2020 Julien Pommier.
%
% This file is a part of GetFEM
%
% GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
% under  the  terms  of the  GNU  Lesser General Public License as published
% by  the  Free Software Foundation;  either version 3 of the License,  or
% (at your option) any later version along with the GCC Runtime Library
% Exception either version 3.1 or (at your option) any later version.
% This program  is  distributed  in  the  hope  that it will be useful,  but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% License and GCC Runtime Library Exception for more details.
% You  should  have received a copy of the GNU Lesser General Public License
% along  with  this program;  if not, write to the Free Software Foundation,
% Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.


function check_mesh_fem(iverbose,idebug)
  global gverbose;
  global gdebug;  
  if (nargin >= 1),
    gverbose = iverbose;
    if (nargin == 2),
      gdebug = idebug;
    else 
      gdebug = 0;
    end;
  else 
    gverbose = 0;
  end;
  s=['BEGIN POINTS LIST' 10,...
     'POINT  0  0 0 0' 10,...
     'POINT  1  -4 6 2' 10,...
     'POINT  2  0 6 0' 10,...
     'POINT  3  0 2 0' 10,...
     'POINT  4  -2 6 2' 10,...
     'POINT  5  0 4 0' 10,...
     'POINT  6  -1.5 4.5 .5' 10,...
     'POINT  7  1 2 0' 10,...
     'POINT  8  1.5 1.5 0' 10,...
     'POINT  9   5 5 0' 10,...
     'POINT  10  2 1 0' 10,...
     'POINT  11  6 3 0' 10,...
     'POINT  12  2 0 0' 10,...
     'POINT  13  6 0 0' 10,...
     'POINT  14  2 4 0' 10,...
     'POINT  15  4 2 0' 10,...
     'POINT  46  4 4 0' 10,...
     'POINT  17  3 6 0' 10,...
     'POINT  18  2 -2 2' 10,...
     'POINT  19  2 -2 -2' 10,...
     'POINT  20  6 -2 2' 10,...
     'POINT  21  6 -2 -2' 10,...
     'POINT  22  2 -1 1' 10,...
     'POINT  23  2 -2.5 0' 10,...
     'POINT  24  2 -1 -1' 10,...
     'POINT  25  6 -1 1' 10,...
     'POINT  26  6 -2.5 0' 10,...
     'POINT  27  6 -1 -1    ' 10,...
     'POINT  28  -1 6 -1 ' 10,...
     'POINT  29  -1 2 -1' 10,...
     'POINT  30  +1 6 -2' 10,...
     'POINT  31  +1 2 -2' 10,...
     'POINT  32  0 6 -3' 10,...
     'POINT  33  0 2 -3' 10,...
     'POINT  34  2 -5 -2' 10,...
     'POINT  35  2 -4 0' 10,...
     'POINT  36  4 -5 2' 10,...
     'POINT  37  6 -5 -2' 10,...
     'POINT  38  6 -5 0' 10,...
     'POINT  49  6 -5 2' 10,...
     'END POINTS LIST' 10,...
     'BEGIN MESH STRUCTURE DESCRIPTION' 10,...
     'CONVEX 0    GT_PK(2,2)      1 4 2 6 5 3' 10,...
     'CONVEX 1    GT_QK(2,1)      2 17 3 7' 10,...
     'CONVEX 2    GT_QK(2,2)      7 8 10 14 46 15 17 9 11' 10,...
     'CONVEX 3    GT_QK(2,1)      10 12 11 13' 10,...
     'CONVEX 4    GT_PRODUCT(GT_PK(2,2),GT_PK(1,1)) 12 22 18 24 23 19 13 25 20 27 26 21' 10,...
     'CONVEX 5    GT_PRODUCT(GT_PK(1,1),GT_PK(1,3)) 2 3 28 29 30 31 32 33' 10,...
     'CONVEX 8    GT_PRODUCT(GT_QK(2,1),GT_PK(1,2)) 19 21 34 37 23 26 35 38 18 20 36 49' 10,...
     'END MESH STRUCTURE DESCRIPTION' 10,...
     'BEGIN MESH_FEM' 10,...
     ' CONVEX 0 FEM_PK(2,2)' 10,...
     ' CONVEX 1 FEM_QK(2,2)' 10,...
     ' CONVEX 2 FEM_QK(2,3)' 10,...
     ' CONVEX 3 FEM_QK(2,2)' 10,...
     ' CONVEX 4 FEM_PRODUCT(FEM_PK(2,2),FEM_PK(1,2))' 10,...
     ' CONVEX 5 FEM_PRODUCT(FEM_PK(1,2),FEM_PK(1,3))' 10,...
     ' CONVEX 8 FEM_PRODUCT(FEM_QK(2,3),FEM_PK(1,2))' 10,...
     'END MESH_FEM\n'];
  m=gf_mesh('from string',s);
  mf=gf_mesh_fem('from string',s,m);
  s2=gf_mesh_fem_get(mf,'char');
  mf2=gf_mesh_fem('from string',s2,m);
  gf_mesh_fem_get(mf,'nbdof');
  gf_mesh_fem_get(mf2,'nbdof');
  
  N=gf_mesh_get(m,'dim');
  npt=gf_mesh_get(m,'nbpts');
  gfassert('N==3 & npt==40');
  ncv=gf_mesh_get(m,'nbcvs');
  gfassert('ncv==7');
  lastcv=gf_mesh_get(m, 'max cvid');
  gfassert('lastcv==9');
  lastpid=gf_mesh_get(m, 'max pid');
  gfassert('lastpid==50');
  [d,c]=gf_mesh_get(mf, 'pid from cvid',[2 6]);
  gfassert('c==[1 5 13]');
  gfassert('d==[3 18 4 8 3 4 29 30 31 32 33 34]');
  [d,c]=gf_mesh_get(mf, 'pid from cvid',1:gf_mesh_get(m,'max cvid'));
  [d,c]=gf_mesh_get(mf, 'pid from cvid');
  for i=[-1 0 -10],
    asserterr('gf_mesh_get(m, ''pid from cvid'',i)');
  end;
  P=gf_mesh_get(m,'pts');
  V=gf_mesh_get(m, 'pid from coords', P);
  pid=gf_mesh_get(m,'pid');
  find(V~=-1)
  pid
  P
  gf_mesh_get(m, 'char')
  gfassert('find(V~=-1)==pid');
  a=gf_mesh_get(m, 'faces from pid', pid);
  b=[1 1 1 2 1 3 2 1 2 2 2 3 2 4 6 1 6 2 6 3 6 4 3 1 3 2 3 3 3 4 4 ...
     1 4 2 4 3 4 4 5 1 5 2 5 3 5 4 5 5 9 1 9 2 9 3 9 4 9 5 9 6];
  gfassert('a(:)==b(:)');  
  
  for i=[-1 0 48 49]
    asserterr('gf_mesh_get(m, ''faces from pid'', i)');
  end;
  a=gf_mesh_get(m, 'outer faces');
  b = [5 2 5 3 5 4 5 5 9 1 9 2 9 3 9 5 9 6];
  gfassert('a(:)==b(:)');
  a=gf_mesh_get(m, 'outer faces',3, [4 5]);
  gfassert('a(:)==[5 1 5 2 5 3 5 4 5 5]''');
  % gf_mesh_get(m, 'outer faces',[4 6 7 8])
  asserterr('gf_mesh_get(m, ''outer faces'',[4 6 7 8])');
  asserterr('gf_mesh_get(m, ''outer faces'',3, [0])');
  E=gf_mesh_get(m, 'edges');
  asserterr('gf_mesh_get(m, ''edges'', [0])');
  E=gf_mesh_get(m, 'curved edges',10);
  E=gf_mesh_get(m, 'curved edges',8);
  asserterr('gf_mesh_get(m, ''curved edges'',-1)');
  gfassert('abs(sum(sum(sum(E)))-1.872e3) < 2');
  asserterr('gf_mesh_get(m, ''triangulated surface'', 3)');
  Z=gf_mesh_get(m, 'triangulated surface', 4,gf_mesh_get(m, 'outer faces',[4 5]));
  % gfassert('size(Z)==[9 160]');
  gfassert('size(Z)==[9 128]');
  Z=gf_mesh_get(m, 'curved edges', 4, gf_mesh_get(m, 'outer faces',[4 5]));
  ZZ=gf_mesh_get(m, 'curved edges', 4, [4 5]);
  for i=0:7
    if (i > 0 & i < 7),
      n=gf_mesh_get(m, 'normal of face', 5, 3, i);
      gfassert('norm(n-[0    0.7071    0.7071]) < 1e-3');
      nn(i,:)=gf_mesh_get(m, 'normal of face', 5, 1, i);
    else
      asserterr('gf_mesh_get(m, ''normal of face'', 5, 3, i)');
    end;
  end;
  zz=[0 0 0 0 0 0 -0.894427 -1 -0.894427 -0.894427 -1 -0.894427 0.447214 0 -0.447214 0.447214 0 -0.447214];  
  gfassert('norm(nn(:)''-zz)<1e-5'); %8.9465e-07
  asserterr('gf_mesh_get(m, ''normal of faces'', [1 -1])');
  N=gf_mesh_get(m, 'normal of faces', gf_mesh_get(m, 'outer faces',[5 9]));
  s2=gf_mesh_get(m,'char');
  gfassert('length(s2)>500');
  m2=gf_mesh('from string',s);

  gf_mesh_fem_get(mf,'nbdof') % should be 99 or 100 (element 0 and 1 are not really neigbhor but can be viewed as such)
  d=gf_mesh_fem_get(mf,'basic dof from cv',[1 5])
  gfassert(['d==[1 2 3 4 5 6 29 32 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49]']);
  d=gf_mesh_fem_get(mf,'basic dof from cv',[1 5;1 2])
  gfassert('d==[3 5 6 29 32 34 37 39 42 44 47 49]');
  d=gf_mesh_fem_get(mf,'basic dof from cvid',5)
  gfassert('d==[29 35 36 37 38 39 32 40 41 42 43 44 34 45 46 47 48 49]');
  
  s2=gf_mesh_get(mf,'char');
  gfassert('length(s2)>500');
  m2=gf_mesh('from string',s);
  mf2=gf_mesh_fem('from string',s);
  mf3=gf_mesh_fem('from string',s,m2);
  gf_mesh_fem_set(mf2,'qdim',2);
  s2=gf_mesh_fem_get(mf2,'char')
  s3=gf_mesh_get(mf2,'char')
  % ~bug here: doesn't work if s2 and s3 are reversed
  mf2=gf_mesh_fem('from string',[s3 s2]); 
  
  
  d=gf_mesh_fem_get(mf2,'basic dof from cv',[1 5]);
%  dd=[1 2 3 4 5 6 7 8 9 10 11 12 73 74 79 80 83 84 85 86 87 88 89 90 91 92 ...
%      93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 ...
%      112 113 114];
%  gfassert('d==dd');
  
  d=gf_mesh_fem_get(mf2,'basic dof from cv',[1 5;1 2]);
%  dd=[5 6 9 10 11 12 73 74 79 80 83 84 89 90 93 94 99 100 103 104 109 110 113 ...
%      114];  
%  gfassert('d==dd');
  d=gf_mesh_fem_get(mf2,'basic dof from cvid',5);
%  dd=[73 74 85 86 87 88 89 90 91 92 93 94 79 80 95 96 97 98 99 100 101 102 ...
%      103 104 83 84 105 106 107 108 109 110 111 112 113 114];  
%  gfassert('d==dd');

  [f,c]=gf_mesh_get(mf2, 'geotrans');
  gfassert('c(2)==c(4)');

  fs1=gf_geotrans_get(f(c(6)),'char');
  gfassert('fs1==''GT_PRODUCT(GT_PK(1,1),GT_PK(1,3))''');
  [f,c]=gf_mesh_get(mf2, 'cvstruct');
  gfassert('c(2)==c(4)');
  [f,c]=gf_mesh_fem_get(mf2, 'fem');
  gfassert('c(2)==c(4)');
  fs1=gf_fem_get(f(c(5)),'char');
  gfassert('fs1==''FEM_PRODUCT(FEM_PK(2,2),FEM_PK(1,2))''');
%  [f,c]=gf_mesh_fem_get(mf2, 'integ');
%  gfassert('c(2)==c(3)');
%  fs1=gf_integ_get(f(c(3)),'char');
%  gfassert('fs1==''IM_QUAD(5)''');
  
  %test for non conformal dof
  m = gf_mesh('triangles grid',[0:.5:1], [0:.5:1]);
  mf_u = gf_mesh_fem(m,2);
  mf_d = gf_mesh_fem(m,1);
  gf_mesh_fem_set(mf_u,'fem',gf_fem('FEM_PK(2,1)'),[1:5 7]);
  gf_mesh_fem_set(mf_u,'fem',gf_fem('FEM_PK(2,3)'),8);
  gf_mesh_fem_set(mf_d,'fem',gf_fem('FEM_PK(2,1)'));
  gf_mesh_set(m, 'boundary', 3, [1 1 1; 1 2 3]);
  gf_mesh_set(m, 'boundary', 7, [3 4; 3 2]);
  cl=[1:5 7 8];
  asserterr('gf_mesh_fem_get(mf_u, ''non conformal basic dof'')');
  d=gf_mesh_fem_get(mf_u, 'non conformal basic dof',cl)
  %gf_plot_mesh(mf_u, 'dof', 'on');
  gfassert('d==[19 20 21 22 23 24 29 30]');

  f=gf_mesh_fem_get(mf2, 'fem');
  f5=gf_mesh_fem_get(mf2, 'fem',5);
  asserterr('gf_mesh_fem_get(mf_u, ''fem of cvs'')');
  asserterr('gf_mesh_fem_get(mf_u, ''fem of cvs'',7)');
  asserterr('gf_mesh_fem_get(mf_u, ''fem of cvs'',6)');
  gf_mesh_fem_get(mf_u, 'is_lagrangian',cl);
  gf_mesh_fem_get(mf_u, 'is equivalent',cl);
  gf_mesh_fem_get(mf_u, 'is_polynomial',cl);
  
  me=gf_eltm('base', f5);
%  ME=gf_mesh_fem_get(mf2,'eltm',me,5);
%  MME=[-0.0444444 1.15556 0.0222222 1.15556 1.24444 0.0222222 -0.177778 4.62222,...
%       0.0888889 4.62222 4.97778 0.0888889 -0.0444444 1.15556 0.0222222 1.15556,...
%       1.24444 0.0222222]';
%  gfassert('norm(ME-MME)<1e-4');
  m=gf_mesh_fem_get(mf2,'linked_mesh');

  oo=gf_mesh_get(mf2,'outer faces');
  oo=oo(:,find(oo(2,:)~=0));
  gf_mesh_set(m,'boundary',51,oo);

  o=gf_mesh_get(m,'boundary',51);
  gfassert('size(o)==size(oo) & sum(sum(o))==sum(sum(oo))');

  o=gf_mesh_get(mf2,'boundary',1);
  gfassert('isempty(o)');
  gf_mesh_set(gf_mesh_fem_get(mf2,'linked mesh'),'boundary',1,oo(:,1));
  o=gf_mesh_get(mf2,'boundary',1);
  gfassert('o==oo(:,1)');

  o=gf_mesh_get(mf2,'boundaries');
  gfassert('o==[1 51]');
  
  gf_mesh_set(gf_mesh_fem_get(mf2,'linked mesh'),'delete boundary',1);

  o=gf_mesh_get(mf2,'boundary',1);
  gfassert('isempty(o)');
  o=gf_mesh_get(mf2,'boundaries');
  gfassert('o==51');
  
  % test region intersect/merge/setdiff
  R1 = [1 2 5 6 6 6; 1 2 1 2 3 2];
  gf_mesh_set(m,'region', 10, R1);
  r1 = gf_mesh_get(m, 'region', 10);
  R2 = [5 6 3 4 6; 1 3 1 2 1];
  gf_mesh_set(m,'region', 11, R2);
  r2 = gf_mesh_get(m, 'region', 11);
  gf_mesh_set(m,'region merge', 11, 10);
  rr=gf_mesh_get(m,'region',11);
  RR=union(R1',R2','rows')';
  gfassert('rr==RR');
  
  gf_mesh_set(m,'region', 11, R2);
  gf_mesh_set(m,'region subtract', 11, 10);
  rr=gf_mesh_get(m,'region',11);
  RR=setdiff(R2',R1','rows')';
  gfassert('rr==RR');

  gf_mesh_set(m,'region', 11, R2);
  gf_mesh_set(m,'region intersect', 11, 10);
  rr=gf_mesh_get(m,'region',11);
  RR=intersect(R2',R1','rows')';
  gfassert('rr==RR');

  
  
  asserterr('gf_mesh_set(m, ''del point'', [3])');
  o=gf_mesh_get(m,'pid from cvid', 3);
  gfassert('o==[8 9 11 15 47 16 18 10 12]');

  gf_mesh_set(m,'del convex',3);
  gf_mesh_set(m,'del convex',2);
  c=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 25,...
     26 27 28 29 30 31 32 33 34 35 36 37 38 39 47 50];
  d=gf_mesh_get(m,'pid');
  gfassert('d(:)==c(:)');
  
  d=gf_mesh_fem_get(mf2,'basic dof on region', 0:100);
  
  % test for optimize_structure
  maxpid=gf_mesh_get(m,'max pid');
  maxcvid=gf_mesh_get(m,'max cvid');
  np=gf_mesh_get(m,'nbpts');
  ncv=gf_mesh_get(m,'nbcvs');
  gfassert('np < maxpid');
  gfassert('ncv < maxcvid');

  
  gf_mesh_set(m,'optimize structure', false);

  maxpid=gf_mesh_get(m,'max pid');
  maxcvid=gf_mesh_get(m,'max cvid');
  gfassert('np == maxpid');
  gfassert('ncv == maxcvid');

  gf_mesh_set(m,'del convex',2);
  disp('-----------------------------PLOP---------------------------------');
  
  disp('-----------------------------PLOP---------------------------------');
  
  % test gradient/hessian
  m=gfMesh('empty',2); 
  %gf_mesh_set(m,'add convex',gfGeoTrans('GT_PK(2,2)'),...
  %[0 0; .6 0; 1.2 0; 0 .4; .6 .4; 0 0.8]');
  gf_mesh_set(m,'add convex',gfGeoTrans('GT_PK(2,1)'),[1 1; 1.1 0;0.9 1.3]');
  mf=gfMeshFem(m);
  gf_mesh_fem_set(mf, 'classical fem', 4); %'gfFem('FEM_PK(2,3)'));
  U=rand(3, mf.nbdof());
  %U(1) = 1;
  
  DU=gf_compute(mf, U, 'gradient', mf);
  D2U=gf_compute(mf, DU, 'gradient', mf);
  
  D2U2=gf_compute(mf, U, 'hessian', mf);
  gfassert('max(max(abs(D2U(:)-D2U2(:)))) < 1e-9');
