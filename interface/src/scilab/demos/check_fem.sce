gf_workspace('clear all');
lines(0);

f   = gf_fem('FEM_PK(3,4)');
dim = gf_fem_get(f,'dim');
assert('dim==3');
tdim = gf_fem_get(f,'target_dim');  
assert('tdim==1');
nbd = gf_fem_get(f,'nbdof');
assert('nbd==35');
is_pol = gf_fem_get(f,'is_polynomial');
assert('is_pol');
is_lag = gf_fem_get(f,'is_lagrange');
assert('is_lag');
is_equ = gf_fem_get(f,'is_equivalent');
assert('is_equ');
p = gf_fem_get(f,'pts');
assert('size(p)==[3 35]');
ed = gf_fem_get(f,'estimated_degree');
assert('ed==4');
Z = [5 -8 3 0 0 -8 12 -4 0 3 -4 1 0 0 0 -8 12 -4 0 12 -16 4 -4 4 0 3 -4 1 -4 ...
4 1 0 0 0 0]';
z = gf_fem_get(f,'base_value',[.5;.5;.5]);
assert('norm(Z-z) < 1e-13');
asserterr('gf_fem_get(f,''base_value'',[.5;.5])');
DZ = [77 -152 84 -8 -1 -104 192 -96 8 30 -48 18 0 0 0 -104 192 -96 8 120 -192 ...
72 -24 24 0 30 -48 18 -24 24 0 0 0 0 0 77 -104 30 0 0 -152 192 -48 0 ...
84 -96 18 -8 8 -1 -104 120 -24 0 192 -192 24 -96 72 8 30 -24 0 -48 24 ...
18 0 0 0 0 77 -104 30 0 0 -104 120 -24 0 30 -24 0 0 0 0 -152 192 -48 0 ...
192 -192 24 -48 24 0 84 -96 18 -96 72 18 -8 8 8 -1];  
dz = gf_fem_get(f,'grad_base_value',[.5;.5;.5]);
assert('norm(DZ(:)-dz(:)*3) < 1e-12'); // 2.8432e-13 on sgi O2K / CC debug mode
assert('size(dz)==[35 1 3]');
HZ = [284 -704 552 -128 -4 -288 672 -480 96 48 -96 48 0 0 0 -288 672 -480 96 ...
192 -384 192 0 0 0 48 -96 48 0 0 0 0 0 0 0 284 -496 228 -16 0 -496 816 ...
-336 16 228 -336 108 -16 16 0 -288 432 -144 0 432 -576 144 -144 144 0 ...
48 -48 0 -48 48 0 0 0 0 0 284 -496 228 -16 0 -288 432 -144 0 48 -48 0 ...
0 0 0 -496 816 -336 16 432 -576 144 -48 48 0 228 -336 108 -144 144 0 ...
-16 16 0 0 284 -496 228 -16 0 -496 816 -336 16 228 -336 108 -16 16 0 ...
-288 432 -144 0 432 -576 144 -144 144 0 48 -48 0 -48 48 0 0 0 0 0 284 ...
-288 48 0 0 -704 672 -96 0 552 -480 48 -128 96 -4 -288 192 0 0 672 -384 ...
0 -480 192 96 48 0 0 -96 0 48 0 0 0 0 284 -288 48 0 0 -496 432 -48 0 ...
228 -144 0 -16 0 0 -496 432 -48 0 816 -576 48 -336 144 16 228 -144 0 ...
-336 144 108 -16 0 16 0 284 -496 228 -16 0 -288 432 -144 0 48 -48 0 0 ...
0 0 -496 816 -336 16 432 -576 144 -48 48 0 228 -336 108 -144 144 0 -16 ...
16 0 0 284 -288 48 0 0 -496 432 -48 0 228 -144 0 -16 0 0 -496 432 -48 ...
0 816 -576 48 -336 144 16 228 -144 0 -336 144 108 -16 0 16 0 284 -288 ...
48 0 0 -288 192 0 0 48 0 0 0 0 0 -704 672 -96 0 672 -384 0 -96 0 0 552 ...
-480 48 -480 192 48 -128 96 96 -4];
hz = gf_fem_get(f,'hess_base_value',[.5;.5;.5]);
assert('norm(HZ(:)-hz(:)*3) < 1e-12'); // 7.9986e-13 on sgi O2K / CC debug mode
assert('size(hz)==[35 1 3 3]');
f = gf_fem('FEM_HERMITE(1)');
f = gf_fem('FEM_HERMITE(3)');
f = gf_fem('FEM_PK_DISCONTINUOUS(2,1)');
f = gf_fem('FEM_P1_NONCONFORMING');
f = gf_fem('FEM_PK_WITH_CUBIC_BUBBLE(2,1)');
ed = gf_fem_get(f,'estimated_degree');
assert('ed==3');
asserterr('gf_fem(''FEM_PK_WITH_CUBIC_BUBBLE(2,4)'')'); // YC: logic error here
f = gf_fem('FEM_PK_PRISM_HIERARCHICAL(3,3)');
nbd = gf_fem_get(f,'nbdof');
assert('nbd==40');
is_pol = gf_fem_get(f,'is_polynomial');
assert('is_pol');
is_lag = gf_fem_get(f,'is_lagrange');
assert('~is_lag');
is_equ = gf_fem_get(f,'is_equivalent');
assert('is_equ');
P = [0 0 0 3 0 0 0 3 0 1 0 0 2 0 0 0 1 0 1 1 0 2 1 0 0 2 0 1 2 0 0 0 3 3 0 ...
3 0 3 3 1 0 3 2 0 3 0 1 3 1 1 3 2 1 3 0 2 3 1 2 3 0 0 1 3 0 1 0 3 1 1 ...
0 1 2 0 1 0 1 1 1 1 1 2 1 1 0 2 1 1 2 1 0 0 2 3 0 2 0 3 2 1 0 2 2 0 2 ...
0 1 2 1 1 2 2 1 2 0 2 2 1 2 2];
p = gf_fem_get(f,'pts');
assert('norm(P(:)-p(:)*3)<1e-15'); // exactly 0 on sgi O2K
ed = gf_fem_get(f,'estimated_degree');
assert('ed==6');
