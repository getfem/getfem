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

gf_workspace('clear all');
if (0) then
  asserterr('gf_spmat(''empty'',-5)');
  asserterr('gf_spmat(''empty'',2:3)');
  asserterr('gf_spmat(''empty'',0)');
end;
// TEST EMPTY COPY FULL
A  = gf_spmat('empty', 5,6);
B  = gf_spmat('empty', 11111);
C  = gf_spmat('copy', A);
C  = sprand(50,50,.1); C(2,2)=1+2*%i; I = 1:40; J = [6 7 8 3 10];
D  = gf_spmat('copy', C, I, J);
DD = gf_spmat_get(D,'full');
// gfassert('and(DD==C(I,J))'); // Bug scilab 4783
gfassert('and(full(DD)==full(C(I,J)))');
asserterr('gf_spmat(D,''full'',100)');
asserterr('gf_spmat(D,''full'',10,-1)');
// TEST MULT
A  = gf_spmat('identity', 11111);
C  = gf_spmat('mult',A,B);
n  = gf_spmat_get(C,'nnz'); gfassert('n==0');
C  = gf_spmat('mult',A,A);
n  = gf_spmat_get(C,'nnz'); gfassert('n==11111');
M1 = sprand(50,43,.1);
M2 = sprand(43,14,.3);
C  = gf_spmat('mult',M1,M2);
C  = gf_spmat_get(C, 'full');
P  = full(M1*M2);
gfassert('max(max(abs(C-P)))<1e-13');
asserterr('gf_spmat(''mult'',M2,M1);'); 
//TEST ADD
d  = rand(1,size(P,1));
D  = gf_spmat('diag', d');
M1 = sprand(50,50,.1); C(2,2)=1+2*%i;
M2 = sprand(50,50,.1); C(2,2)=1+2*%i;  
C  = gf_spmat('add',M1, M2);
C  = gf_spmat_get(C, 'full');
gfassert('max(max(abs(C-full(M1+M2))))<1e-13');
C = gf_spmat('add',M1, real(M2));
C = gf_spmat_get(C, 'full');
gfassert('max(max(abs(C-full(M1+real(M2)))))<1e-13');
// TEST DIAG
K = gf_spmat('diag', [1 1; 2 3; 4 5; 6 7],[0 -2],6,9);
// NNZ
gf_spmat_get(K,'full');
//gfassert('gf_spmat_get(K,''nnz'')==8'); // YC: 6 - a voir
//-->gf_spmat_get(K,'full')
// ans  =
// 
//    1.    0.    0.    0.    0.    0.    0.    0.    0.  
//    0.    2.    0.    0.    0.    0.    0.    0.    0.  
//    5.    0.    4.    0.    0.    0.    0.    0.    0.  
//    0.    7.    0.    6.    0.    0.    0.    0.    0.  
//    0.    0.    0.    0.    0.    0.    0.    0.    0.  
//    0.    0.    0.    0.    0.    0.    0.    0.    0.  
cK=gf_spmat('diag', [1 1*%i; 2 3*%i; 4 5; 6*%i 7; 5 5; 6 -2],[0 -1],6,9);
gfassert('gf_spmat_get(cK,''nnz'')==11');  
C = gf_spmat('add',K,cK);
gfassert('gf_spmat_get(C,''is_complex'')');
// MULT VECTOR
fK  = gf_spmat_get(K,'full');
fcK = gf_spmat_get(cK,'full');
V6  = rand(6,1); 
V9  = rand(9,1);
W6  = gf_spmat_get(K,'mult',V9);
gf_spmat_get(K, 'full');
W9  = gf_spmat_get(K,'tmult',V6);
gfassert('max(abs(W6(:)-fK*V9))<1e-13');
gfassert('max(abs(W9(:)-fK''*V6))<1e-13');
W6  = gf_spmat_get(cK,'mult',V9);
W9  = gf_spmat_get(cK,'tmult',V6);
gfassert('max(abs(W6(:)-fcK*V9))<1e-13');
gfassert('max(abs(W9(:)-fcK''*V6))<1e-13');
V6  = rand(6,1) + 1*%i*rand(6,1); V9=rand(9,1) + 1*%i*rand(9,1);
asserterr('gf_spmat_get(K,''mult'',V9)');
W6  = gf_spmat_get(cK,'mult',V9);
W9  = gf_spmat_get(cK,'tmult',V6);
gfassert('max(abs(W6(:)-fcK*V9))<1e-13');
gfassert('max(abs(W9(:)-fcK''*V6))<1e-13');
// STORAGE, SIZE, IS_COMPLEX, CSC_IND CSC_VAL
gf_spmat_get(cK, 'storage');
gfassert('gf_spmat_get(cK, ''size'')==[6 9]');
gfassert('gf_spmat_get(sparse(fcK), ''size'')==[6 9]');
[jc,ir]=gf_spmat_get(cK, 'csc_ind');
v=gf_spmat_get(cK, 'csc_val');
// CLEAR
gf_spmat_set(K, 'to_wsc'); gf_spmat_set(cK, 'to_wsc');  
KK = gf_spmat('copy',K); gf_spmat_set(KK,'clear');
gfassert('gf_spmat_get(KK,''nnz'')==0');
KK = gf_spmat('copy',cK); gf_spmat_set(KK,'clear');
gfassert('gf_spmat_get(KK,''nnz'')==0');
for i=1:20, 
  if (modulo(i,2)==0) then
    gf_spmat_set(K, 'to_wsc'); 
  else 
    gf_spmat_set(K, 'to_csc'); 
  end
  KK = gf_spmat('copy',cK); gf_spmat_set(KK,'scale',int32(-1));
  C  = gf_spmat('add',cK,KK);
  gfassert('gf_spmat_get(C,''nnz'')==0');
  C = gf_spmat('copy',cK); gf_spmat_set(C,'transpose');
  gfassert('and(gf_spmat_get(C,''full'')==fcK.'')');
  C = gf_spmat('copy',cK); gf_spmat_set(C,'transconj');
  gfassert('and(gf_spmat_get(C,''full'')==fcK'')');
  C = gf_spmat('copy',cK); gf_spmat_set(C,'conjugate');
  gfassert('and(gf_spmat_get(C,''full'')==conj(fcK))');
end
gf_spmat_set(cK,'to_complex');
C = gf_spmat('copy',K); gf_spmat_set(C,'to_complex');
gfassert('gf_spmat_get(C,''is_complex'')');
gf_spmat_set(C,'clear');
B = [1 1 1 1 1 2; ...
6 5 4 3 2 1; ...
7 8 5 3 2 1]';
gf_spmat_set(C,'diag', B(:,1));
gf_spmat_set(C,'diag', B(:,2:3), [-2 +2]);
CC = full(spdiags(B, [0 -2 2], 6, 9));
P  = gf_spmat_get(C,'full');
gfassert('and(CC==P)');
//L1 = gf_spmat_get(C,'diag', [0 -2 2]);
//L2 = spdiags(sparse(CC),[0 -2 2]);
//gfassert('L1==L2');
K  = sprand(50,50,.1) + 20*speye(50,50); K(2,3)=.4;
gK = gf_spmat('copy',K);
gf_spmat_set(gK, 'to_csc');
asserterr('gf_spmat_set(gK, ''assign'', 1, 1, 1)');
gf_spmat_set(gK, 'to_wsc');
gf_spmat_set(gK, 'assign', 2, 2, 1+2*%i);
gf_spmat_set(gK, 'add', 2, 2:4, 2*%i*ones(1,3));
A = gf_spmat_get(gK, 'full', 2, 2:4);
B = full(K(2,2:4)); B(1)=1+2*%i; B=B+2*%i;
gfassert('max(abs(A-B))<1e-13');
gf_workspace('clear all')
// luinc not yet defined under Scilab
if 0 then
  m   = gf_mesh('cartesian',[1:10],[1:10]);
  mf  = gf_mesh_fem(m,1);
  gf_mesh_fem_set(mf,'classical fem', 1);
  mim = gf_mesh_im(m, 0); // integration of degree 0
  A   = gf_asm('laplacian',mim,mf,mf,ones(1,gf_mesh_fem_get(mf,'nbdof')));
  A   = A +.1*speye(size(A,1),size(A,1));
  B   = rand(gf_mesh_fem_get(mf,'nbdof'),1);
  //[L,U] = luinc(A,'0');
  [L,U] = lu(full(A));
  X1 = gf_linsolve('cg',A,B);
  mm = gf_spmat('copy',inv(L));
  p  = gf_precond('spmat',mm);
  gf_workspace('stats')
  X2 = gf_linsolve('cg',A,B,gf_precond('spmat',speye(size(A))));
  gfassert('norm(X1-X2)<1e-13');
end
