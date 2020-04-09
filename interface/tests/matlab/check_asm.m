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

function check_asm(iverbose,idebug)
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
    gverbose = 0; gdebug = 0;
  end;

  gf_workspace('clear all');
  p=[0 1 0 1.5; 0 0 1 1];
  t=[1 2 3 0; 2 3 4 0]';
  m=gf_mesh('pt2D',p,t);
  z=rand(200,22);
  zz=rand(200,22);
  z=z+zz;
  clear z;
  zzz=rand(200,22);
  zz=zz+zz;
  
  
  mf=gf_mesh_fem(m,1);
  mim=gf_mesh_im(m);
  mf3=gf_mesh_fem(m,3);
  gf_mesh_fem_set(mf,'fem', gf_fem('FEM_PK(2,1)'));
  gf_mesh_fem_set(mf3,'fem', gf_fem('FEM_PK(2,2)'));
  gf_mesh_im_set(mim,'integ',gf_integ('IM_TRIANGLE(3)'));
  v=gf_asm('volumic','V(#1)+=comp(Base(#1).Base(#1)(i))',mim,mf);
  asserterr('gf_asm(''volumic'',''V(#1)+=comp(Base(#2))'',mf)');
  a = gf_compute(mf,v','l2 norm',mim);
  b = gf_compute(mf,1i*v','l2 norm',mim);
  gfassert('a==b');
  a = gf_compute(mf,v','h1 norm',mim);
  b = gf_compute(mf,1i*v','h1 norm',mim);
  gfassert('a==b');
  
  X=gf_asm('volumic','V(#1,#2)+=comp(Base(#1).Base(#1))',mim,mf,mf);
  gfassert('max(abs((X-X'')))<1e-15');
  X=gf_asm('volumic','V(#1,#1,#1,#1)+=comp(Base(#1).Base(#1).Base(#1).Base(#1))',mim,mf);
  gfassert('size(X)==[4 4 4 4]');
  X=gf_asm('volumic','M(#1,#2)+=comp(Grad(#1).vBase(#2))(:,z,:,i)',mim,mf,mf3);
  gfassert('size(X)==[4 27]');
  gfassert('abs(sum(sum(abs(X)))-10.5) < 8e-15');
  asserterr('gf_asm(''volumic'',''V(#1)+=comp(Base(#1))'',mim,mf3)');
  X=gf_asm('volumic','V(qdim(#1),#1)+=comp(vBase(#1)){2,1}',mim,mf3);
  for i=1:size(X,1)
    for  j=1:size(X,2)
       if (abs(X(i,j)) < 1E-10)
           X(i,j) = 0;
       end
    end
  end
  gfassert('nnz(X)==15');
  xnnz=find(X);
  zz=[10 14 18 28 32 36 37 41 45 55 59 63 64 68 72];
  gfassert('xnnz(:)==zz(:)');
  X2=gf_asm('volumic','V(3,#1)+=comp(vBase(#1)){2,1}',mim,mf3);
  for i=1:size(X2,1)
    for  j=1:size(X2,2)
       if (abs(X2(i,j)) < 1E-10)
           X2(i,j) = 0;
       end
    end
  end
  gfassert('X2==X');
  X=gf_asm('volumic','V(#1,mdim(#1),mdim(#1))+=comp(Hess(#1))',mim,mf);
  gfassert('X==0');
  X=gf_asm('volumic','V(#1,qdim(#1),mdim(#1),mdim(#1))+=comp(vHess(#1))',mim,mf3);
  gfassert('abs(sum(sum(sum(sum(X))))) < 1e-14');
  asserterr('gf_asm(''volumic'',''V(#1)+=1'')');
  
  H=[.1 .1 0 0; 0 0 0 0; 0 0 0 1]; R=[4 0 1];
  [HH,RR]=gf_spmat_get(sparse(H),'dirichlet nullspace',R);
  gfassert('max(max(abs(full(HH)-[0 -sqrt(2)/2; 0 sqrt(2)/2; 1 0; 0 0]))) < 1e-15');

  gfassert('max(abs(RR-[20 20 0 1]))<1e-14');
  
  % Test on high level generic assembly
  V = ones(1, gf_mesh_fem_get(mf, 'nbdof'));
  K = gf_asm('generic', mim, 2, 'Grad_u.Grad_u/2', -1, 'u', 1, mf, V);
  K2 = gf_asm('laplacian', mim, mf, mf, V);
  gfassert('norm(K-K2, inf) < 1E-12');
  
  gf_asm('undefine function', 'myf');
  gf_asm('define function', 'myf', 1, '6*sin(t)+2*t*t');
  a = gf_asm('generic', mim, 0, 'myf(X(1))', -1);
  gfassert('abs(a-5.48) < 2E-4');
  
  
  m_1d= gf_mesh('cartesian', 1:1:10);
  mf_1d=gf_mesh_fem(m_1d,3);
  gf_mesh_fem_set(mf_1d,'fem', gf_fem('FEM_PK(1,1)'));
  mim=gf_mesh_im(m_1d, 2);
  U = ones(1, gf_mesh_fem_get(mf_1d, 'nbdof'));
  P = ones(1, gf_mesh_fem_get(mf_1d, 'nbdof'));
  K=gf_asm('generic', mim, 2, 'Grad_u.Test_p + p.Grad_Test_u - p.Test_u', -1, 'u', 1, mf_1d, U, 'p', 1, mf_1d, P);
  K1=gf_asm('generic', mim, 2, 'Grad_u.Test_p', -1, 'u', 1, mf_1d, U, 'p', 1, mf_1d, P);
  K2=gf_asm('generic', mim, 2, 'p.Grad_Test_u', -1, 'u', 1, mf_1d, U, 'p', 1, mf_1d, P);
  K3=gf_asm('generic', mim, 2, 'p.Test_u', -1, 'u', 1, mf_1d, U, 'p', 1, mf_1d, P);
  err = norm(full(K-K1-K2+K3))
  gfassert('err < 1E-10');
  
  
  
  
  
  
