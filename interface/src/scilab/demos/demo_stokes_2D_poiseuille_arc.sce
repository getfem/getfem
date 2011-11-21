// this example uses the "old" gf_solve instead of the bricks framework..

lines(0);
stacksize('max');

path = get_absolute_file_path('demo_stokes_2D_poiseuille_arc.sce');

printf('demo stokes_2D_poiseuille_arc started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

disp('2D stokes demonstration on a quadratic mesh');

pde = init_pde();

pde('type')   = 'stokes';
pde('asm')('viscos') = 1.0;
pde('F') = list(list(0,0));

R0 = 9; R1 = 10;

pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list(sprintf('(y-%d).*(y-%d)',-R0,-R1),0));
pde('bound')($)('H')    = [];
pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list(0,sprintf('(x-%d).*(x-%d)',R0,R1)));
pde('bound')($)('H')    = [];
pde = add_empty_bound(pde);
pde('bound')($)('type') = 'Dirichlet';
pde('bound')($)('R')    = list(list(0,0));
pde('bound')($)('H')    = [];

//m = gf_mesh('import','GiD', path + '/data/tube_2D_spline.GiD.msh');

Ku=3; Kp=1;
Nt=8; Nr=4;

m      = gf_mesh('empty',2);
dtheta = 3*%pi/2*1/Nt; 
R      = R0+(R1-R0)*(0:Nr-1)/(Nr-1);
gt_order = min(Ku,3);
gt       = gf_geotrans(sprintf('GT_PRODUCT(GT_PK(1,%d),GT_PK(1,1))',gt_order));
ddtheta  = dtheta/gt_order;
for i=1:Nt
  for j=1:Nr-1
    ti=(i-1)*dtheta:ddtheta:i*dtheta;
    X = [R(j)*cos(ti) R(j+1)*cos(ti)];
    Y = [R(j)*sin(ti) R(j+1)*sin(ti)];
    gf_mesh_set(m,'add convex',gt,[X;Y]);
  end
end

fem_u = gf_fem(sprintf('FEM_QK(2,%d)',Ku));
fem_p = gf_fem(sprintf('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,%d),FEM_PK_DISCONTINUOUS(1,%d))',Kp,Kp));
fem_d = gf_fem(sprintf('FEM_QK(2,2)'));

mfulag   = gf_mesh_fem(m,2);

pde('mf_u') = [];
pde('mf_p') = [];
pde('mf_d') = [];
pde('mim')  = [];

pde('mf_u') = gf_mesh_fem(m,2);
pde('mf_p') = gf_mesh_fem(m,1);
pde('mf_d') = gf_mesh_fem(m,1);
pde('mim')  = gf_mesh_im(m, gf_integ('IM_GAUSS_PARALLELEPIPED(2,10)'));

// this is a good example of the usefullness of the cubic bubble
// -> if not used, the pression has strange values
gf_mesh_fem_set(pde('mf_u'),'fem',fem_u);
gf_mesh_fem_set(pde('mf_p'),'fem',fem_u);
gf_mesh_fem_set(pde('mf_d'),'fem',fem_u);

mfulag = pde('mf_u');

// we use a P3 mesh fem for interpolation of the U field, since
// because of its cubic bubble function, the pde.mf_u is not lagrangian 
//gf_mesh_fem_set(mfulag,'fem',gf_fem('FEM_PK(2,3)'), gf_integ('IM_TRIANGLE(5)'));
//   ____
//  /    \
//  |   OUT
//  \__IN

all_faces = gf_mesh_get(m, 'outer faces', gf_mesh_get(m, 'cvid'));
P         = gf_mesh_get(m,'pts');
INpid     = find(abs(P(1,:)) < 1e-4 & P(2,:)<0);
OUTpid    = find(abs(P(2,:)) < 1e-4 & P(1,:)>0);
INfaces   = gf_mesh_get(m, 'faces from pid', INpid);
OUTfaces  = gf_mesh_get(m, 'faces from pid', OUTpid);

gf_mesh_set(m, 'boundary', 1, INfaces);
gf_mesh_set(m, 'boundary', 2, OUTfaces);
gf_mesh_set(m, 'boundary', 3, _setdiff(all_faces',union(INfaces',OUTfaces','r'),'rows')');

tic; 
[U,P,pde2] = gf_solve(pde); 
disp(sprintf('solve done in %.2f sec', toc()));

Ul = gf_compute(pde('mf_u'),U,'interpolate on',mfulag);

h = scf();
h.color_map =  jetcolormap(255);
drawlater;
subplot(2,2,1); 
gf_plot(mfulag,Ul,'norm','on','deformation',Ul,'deformation_scale', 0.05, 'deformed_mesh','on');
colorbar(min(Ul),max(Ul));
title('|U| plotted on the deformed mesh');

subplot(2,2,2); 
gf_plot(pde('mf_p'),P(:)','deformation',U,'deformation_mf',pde('mf_u'),'deformation_scale',0.05); 
colorbar(min(U),max(U));
title('Pression on the deformed mesh');

subplot(2,2,3); 
gf_plot(mfulag,Ul(:)','mesh','on','quiver_density',300,'quiver_scale',0.4); //,'meshopts',ilist('regions',1,'curved','on'));
gf_plot(pde('mf_p'),P(:)');
colorbar(min(Ul),max(Ul));
title('Quiver plot of U, with color plot of the pression');

subplot(2,2,4); 
gf_plot(mfulag,Ul(:)','mesh','on','quiver_density',300,'quiver_scale',0.4);//...'meshopts',list('regions',1,'curved','on'));
gf_plot(pde('mf_p'),P(:)');
title('Quiver plot zoomed');
h.color_map = jetcolormap(255);
drawnow;

printf('demo stokes_2D_poiseuille_arc terminated\n');
