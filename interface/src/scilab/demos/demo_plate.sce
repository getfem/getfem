lines(0);
stacksize('max');

path = get_absolute_file_path('demo_plate.sce');

printf('demo plate started\n');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

// Plate problem test.

NX        = 10.0;
mixed     = %t;
thickness = 0.01;

m    = gf_mesh('regular simplices', 0:1/NX:1.01, 0:1/NX:1.01);
mfut = gf_mesh_fem(m,2);
mfu3 = gf_mesh_fem(m,1);
mfth = gf_mesh_fem(m,2);
mfd  = gf_mesh_fem(m,1);

gf_mesh_fem_set(mfut,'fem',gf_fem('FEM_PK(2,2)'));
gf_mesh_fem_set(mfu3,'fem',gf_fem('FEM_PK(2,1)'));
gf_mesh_fem_set(mfth,'fem',gf_fem('FEM_PK(2,2)'));
gf_mesh_fem_set(mfd, 'fem',gf_fem('FEM_PK(2,2)'));
mim = gf_mesh_im(m, gf_integ('IM_TRIANGLE(5)'));

// get the list of faces whose normal is [-1,0]
flst = gf_mesh_get(m,'outer_faces');
fnor = gf_mesh_get(m,'normal_of_faces',flst);

fleft = flst(:,find(abs(fnor(1,:)+1) < 1e-14));
fright= flst(:,find(abs(fnor(1,:)-1) < 1e-14));
//fleft = compress(abs(fnor[1,:]+1) < 1e-14, flst, axis=1);
//fright= compress(abs(fnor[1,:]-1) < 1e-14, flst, axis=1);

CLAMPED_BOUNDARY = 1;
gf_mesh_set(m,'region',CLAMPED_BOUNDARY, fleft);
SIMPLE_SUPPORT_BOUNDARY = 2
gf_mesh_set(m,'region',SIMPLE_SUPPORT_BOUNDARY, fright);
E  = 1e3;
Nu = 0.3;
Lambda = E*Nu/((1+Nu)*(1-2*Nu));
Mu     = E/(2*(1+Nu));

if ~mixed then
  b0 = gf_mdbrick('isotropic_linearized_plate',mim,mim,mfut,mfu3,mfth,thickness);
else
  b0 = gf_mdbrick('mixed_isotropic_linearized_plate',mim,mfut,mfu3,mfth,thickness);
end

b1 = gf_mdbrick('plate_source_term', b0);
gf_mdbrick_set(b1,'param', 'M', mfd, [ones(1,441); gf_mesh_fem_get_eval(mfd,list(list('x(2)*x(2)/1000')))]); // YC: pb here
b2 = gf_mdbrick('plate clamped support', b1, CLAMPED_BOUNDARY, 'augmented');
b3 = gf_mdbrick('plate simple support', b2, SIMPLE_SUPPORT_BOUNDARY, 'augmented');
b4 = b3;

if mixed then
  b4 = gf_mdbrick('plate closing', b3);
end

mds = gf_mdstate(b4);

printf('running solve...\n');
gf_mdbrick_get(b4,'solve',mds, 'noisy', 'lsolver','superlu');
printf('solve done!\n');

U   = gf_mdstate_get(mds,'state');
nut = gf_mesh_fem_get(mfut,'nbdof');
nu3 = gf_mesh_fem_get(mfu3,'nbdof');
nth = gf_mesh_fem_get(mfth,'nbdof');

ut = U(1:nut); // YC: nut+1 ?
u3 = U(nut+1:(nut+nu3));
th = U((nut+nu3)+1:(nut+nu3+nth));
sl = gf_slice(list('none'), mfu3, 4);
gf_slice_get(sl,'export_to_vtk', path + '/plate.vtk', mfu3, u3, 'Displacement');
gf_slice_get(sl,'export_to_pos', path + '/plate.pos', mfu3, u3,'Displacement');

printf('You can view the solution with (for example):\n');
printf('mayavi -d %s/plate.vtk -f WarpScalar -m BandedSurfaceMap\n', path);
printf('or\n');
printf('mayavi2 -d %s/plate.vtk -f WarpScalar -m Surface\n', path);
printf('or\n');
printf('gmsh %s/plate.pos\n', path);

printf('demo plate terminated\n');
