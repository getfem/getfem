function check_slices(iverbose,idebug)
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

m = gf_mesh('triangles grid',[-5:1:5],[-4:.8:4]);
//  m = gf_mesh('triangles grid',[-1 1],[-1 1]);
gf_mesh_get(m,'cvid');
gf_mesh_set(m,'del convex',[1]);
mf = gf_mesh_fem(m,1);
gf_mesh_fem_set(mf,'fem',gf_fem('FEM_PK(2,2)'))
//U  = gf_mesh_fem_get(mf,'eval', list('x.*x + y.*y'));
U  = gf_mesh_fem_get_eval(mf, list('x.*x + y.*y'));
sl = gf_slice(list('planar',0,[.5;0],[1;0]),m,3);
pp = gf_slice_get(sl,'pts');
assert('abs(pp(1,:)-.5)<1e-15');

sl2 = gf_slice('points',m,pp(:,1:3));
pp2 = gf_slice_get(sl2,'pts');
assert('abs(pp2(1,:)-.5)<1e-15');

//  n=8;sl = gf_slice(m,list('isovalues',-1,mf,U,0.25),n);
sl = gf_slice(list('isovalues',-1,mf,U,16.0),m,4);
//  gf_plot_slice(sl,'mesh','on','data',gf_compute(mf,U,'interpolate on',sl)); colorbar; // YC
pp = gf_slice_get(sl,'pts');
assert('max(sqrt(sum(pp.^2,1)))<4.0000001');

sl = gf_slice(list('isovalues',0,mf,U,9.0),m,7);
pp = gf_slice_get(sl,'pts');
assert('max(abs(3-sqrt(sum(pp.^2,1))))<0.0015');

N=1;
m  = gf_mesh('triangles grid',[-N:(2*N/3):N],[-N:(N/5):N]);
m2 = gf_mesh('cartesian',[-N:(N/5):N]+.1,[-N:(N/7):N]+.1);
sl = gf_slice(list('mesh',m2),m,3); //gf_plot_slice(sl,'mesh_faces','on');
a  = gf_slice_get(sl,'area') - 1.9*1.9;
assert('a < 1e-10');
endfunction
