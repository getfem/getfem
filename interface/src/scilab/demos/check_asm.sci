function check_asm(iverbose,idebug)
global gverbose;
global gdebug;  

[nargout,nargin] = argn();

if (nargin >= 1) then
  gverbose = iverbose;
  if (nargin == 2) then
    gdebug = idebug;
  else 
    gdebug = 0; 
  end
else 
  gverbose = 0;
end

gf_workspace('clear all');

p = [0 1 0 1.5;
     0 0 1 1];
t = [1 2 3 0; 
     2 3 4 0]';
m = gf_mesh('pt2D',p,t);

mf  = gf_mesh_fem(m,1);
mim = gf_mesh_im(m,gf_integ('IM_EXACT_SIMPLEX(2)'));
asserterr('gf_asm(''volumic'',''V(#1)+=comp(Base(#1))'',mim,mf)');

mf3 = gf_mesh_fem(m,3);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,1)'));
gf_mesh_fem_set(mf3,'fem',gf_fem('FEM_PK(2,2)'));
gf_mesh_im_set(mim,'integ',gf_integ('IM_TRIANGLE(3)'));
v = gf_asm('volumic','V(#1)+=comp(Base(#1).Base(#1)(i))',mim,mf)
asserterr('gf_asm(''volumic'',''V(#1)+=comp(Base(#2))'',mf)');

a = gf_compute(mf,v','l2 norm',mim);
b = gf_compute(mf,1*%i*v','l2 norm',mim);
assert('a==b');

a = gf_compute(mf,v','h1 norm',mim);
b = gf_compute(mf,1*%i*v','h1 norm',mim);
assert('a==b');

X=gf_asm('volumic','V(#1,#2)+=comp(Base(#1).Base(#1))',mim,mf,mf);
assert('max(abs((X-X'')))<1e-15');

X=gf_asm('volumic','V(#1,#1,#1,#1)+=comp(Base(#1).Base(#1).Base(#1).Base(#1))',mim,mf);
assert('size(X)==[4 4 4 4]');

// YC: bug here
//X=gf_asm('volumic','M(#1,#2)+=comp(Grad(#1).vBase(#2))(:,z,:,i)',mim,mf,mf3);
//
//assert('size(X)==[4 27]');
//assert('abs(sum(sum(abs(X)))-10.5) < 8e-15');
//asserterr('gf_asm(''volumic'',''V(#1)+=comp(Base(#1))'',mim,mf3)');

X=gf_asm('volumic','V(qdim(#1),#1)+=comp(vBase(#1)){2,1}',mim,mf3);
assert('nnz(X)==27');

xnnz=find(X);
zz=[1 5 9 10 14 18 19 23 27 28 32 36 37 41 45 46 50 54 55 59 63 64 68 72 73 77 81];
assert('xnnz(:)==zz(:)');

X2=gf_asm('volumic','V(3,#1)+=comp(vBase(#1)){2,1}',mim,mf3);
assert('X2==X');

X=gf_asm('volumic','V(#1,mdim(#1),mdim(#1))+=comp(Hess(#1))',mim,mf);
assert('X==0');

X=gf_asm('volumic','V(#1,qdim(#1),mdim(#1),mdim(#1))+=comp(vHess(#1))',mim,mf3);
assert('abs(sum(sum(sum(sum(X))))) < 1e-14');
asserterr('gf_asm(''volumic'',''V(#1)+=1'')');

H = [0.1 0.1 0 0; 
     0   0   0 0; 
     0   0   0 1]; 
R = [4 0 1];

[HH,RR]=gf_spmat_get(sparse(H),'dirichlet nullspace',R);
disp(full(HH))
disp(full(RR))
assert('max(max(abs(full(HH)-[0 -sqrt(2)/2; 0 sqrt(2)/2; 1 0; 0 0]))) < 1e-15');
assert('max(abs(RR-[20 20 0 1]))<1e-14');
endfunction

