lines(0);
stacksize('max');

path = get_absolute_file_path('demo_plasticity.sce');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

L = 100;
H = 20;

m = gf_mesh('triangles grid',0:4:L, 0:2:H);

mim  = gf_mesh_im(m);    gf_mesh_im_set(mim,   'integ', gf_integ('IM_TRIANGLE(6)'));
mfu  = gf_mesh_fem(m,2); gf_mesh_fem_set(mfu,  'fem',   gf_fem('FEM_PK(2,1)'));
mfd  = gf_mesh_fem(m);   gf_mesh_fem_set(mfd,  'fem',   gf_fem('FEM_PK(2,1)'));
mf0  = gf_mesh_fem(m);   gf_mesh_fem_set(mf0,  'fem',   gf_fem('FEM_PK(2,0)'));
mfdu = gf_mesh_fem(m);   gf_mesh_fem_set(mfdu, 'fem',   gf_fem('FEM_PK_DISCONTINUOUS(2,1)'));

lambda = 121150;
mu     = 80769;
von_mises_threshold = 8000;

P = gf_mesh_get(m, 'pts');
pidleft  = find(abs(P(1,:))<1e-6);
pidright = find(abs(P(1,:) - L)<1e-6);

fleft  = gf_mesh_get(m,'faces from pid',pidleft); 
fright = gf_mesh_get(m,'faces from pid',pidright);

// assign boundary numbers
gf_mesh_set(m,'boundary',1,fleft);
gf_mesh_set(m,'boundary',2,fright);

b0 = gf_mdbrick('small deformations plasticity',mim,mfu, von_mises_threshold);
gf_mdbrick_set(b0, 'param','lambda',lambda);
gf_mdbrick_set(b0, 'param','mu',mu);
b1 = gf_mdbrick('generalized dirichlet',b0,1);
b2 = gf_mdbrick('source term',b1,2);

mds = gf_mdstate(b2);

VM = zeros(1,gf_mesh_fem_get(mfdu, 'nbdof'));

//F=[0 -200; 0 -300; 0 0]';
F = [0 -330]';
t = [0 0.9032 1 1.1 1.3 1.5 1.7 1.74 1.7 1.5 1.3 1.1 1 0.9032 0.7 0.5 0.3 0.1 0];

//nbstep = size(F,2);
nbstep = size(t,2);

dd = gf_mesh_fem_get(mf0, 'basic dof from cvid');

h = scf();
h.color_map = jetcolormap(255);

for step=1:nbstep
  gf_mdbrick_set(b2, 'param','source_term', mfd, gf_mesh_fem_get_eval(mfd,list(list(F(1,1)),list(F(2,1)*t(step)))));
  gf_mdbrick_get(b2, 'solve', mds, 'very noisy', 'max_iter', 1000, 'max_res', 1e-6);
    
  T  = gf_mdstate_get(mds, 'tangent_matrix');
  U  = gf_mdstate_get(mds, 'state'); 
  U  = U(1:gf_mesh_fem_get(mfu, 'nbdof'));
  VM = gf_mdbrick_get(b0, 'von mises', mds, mfdu);
  printf('max(abs(VM)) = %f\n', max(abs(VM)));

  drawlater;
  clf();
  subplot(2,1,1);
  gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,'deformation_mf',mfu,'refine', 4, 'deformation_scale',1); 
  colorbar(min(U),max(U));
  
  n = t(step);
  title(['Von Mises criterion for t = ', string(n)]);

  ERR   = gf_compute(mfu,U,'error estimate', mim);
  E     = ERR;
  E(dd) = ERR;
  subplot(2,1,2);
  gf_plot(mf0, E, 'mesh','on', 'refine', 1); 
  colorbar(min(E),max(E));
  title('Error estimate');
  h.color_map = jetcolormap(255);
  drawnow;
  sleep(1000);
end

printf('demo plasticity terminated\n');
