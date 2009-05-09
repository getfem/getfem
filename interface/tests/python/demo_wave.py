## 2D scalar wave equation (helmholtz) demonstration. 
from numpy import *
from getfem import *
import os

make_check=os.environ.has_key('srcdir');

filename='../meshes/holed_disc_with_quadratic_2D_triangles.msh';
if (make_check):
    filename=os.environ['srcdir']+'/'+filename

## Parameters
PK=3; k = 1.;

## Mesh and MeshFems
m=Mesh('import','gid',filename);
mfu=MeshFem(m,1);
mfu.set_fem(Fem('FEM_PK(2,%d)' % (PK,)));
mfd=MeshFem(m,1);
mfd.set_fem(Fem('FEM_PK(2,%d)' % (PK,)));
mim=MeshIm(m, Integ('IM_TRIANGLE(13)'));

## Boundary selection
P=m.pts(); # get list of mesh points coordinates
Psqr=sum(P*P, 0);
cobj=(Psqr < 1*1+1e-6);
cout=(Psqr > 10*10-1e-2);
pidobj=compress(cobj, range(0, m.nbpts()))
pidout=compress(cout, range(0, m.nbpts()))
fobj=m.faces_from_pid(pidobj)
fout=m.faces_from_pid(pidout)
ROBIN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
m.set_region(DIRICHLET_BOUNDARY,fobj)
m.set_region(ROBIN_BOUNDARY,fout)

## Interpolate the exact solution on mfd (assuming it is a Lagrange fem)
wave_expr = ('cos(%f*x[1]+.2)+complex(0.,1.)*sin(%f*x[1]+.2)' % (k,k));
Uinc=mfd.eval(wave_expr);

## Model Bricks
md = Model('complex')
md.add_fem_variable('u', mfu)
md.add_initialized_data('k', [k]);
md.add_Helmholtz_brick(mim, 'u', 'k');
md.add_initialized_data('Q', [complex(0.,1.)*k]);
md.add_Fourier_Robin_brick(mim, 'u', 'Q', ROBIN_BOUNDARY);
md.add_initialized_fem_data('DirichletData', mfd, Uinc);
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, DIRICHLET_BOUNDARY,
                                            'DirichletData');
#md.add_Dirichlet_condition_with_penalization(mim, 'u', 1e12,
#                                             DIRICHLET_BOUNDARY,
#                                             'DirichletData');

## Solving the problem
md.solve();
U = md.variable('u');

if (not(make_check)):

    sl=Slice(('none',), mfu, 8)
    sl.export_to_vtk('wave.vtk', mfu, real(U)*1, 'rWave',
                     mfu, imag(U)*1, 'iWave')

    print 'You can view the solution with (for example):'
    print 'mayavi2 -d wave.vtk -f WarpScalar -m Surface'
