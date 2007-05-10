from getfem import *
from numarray import *

print '3D stokes demonstration on a quadratic mesh -- 512MB of memory needed for the solve!!'

viscosity = 10


m=Mesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh')
print 'mesh loaded!'
mfu=MeshFem(m,3) # velocity
mfulag=MeshFem(m,3)
mfp=MeshFem(m,1) # pressure
mfd=MeshFem(m,1) # data
mfe=MeshFem(m,1) 
mim=MeshIm(m, Integ('IM_TETRAHEDRON(5)'))

mfu.set_fem(Fem('FEM_PK(3,2)'))
mfd.set_fem(Fem('FEM_PK(3,2)'))
mfp.set_fem(Fem('FEM_PK(3,1)'))
mfe.set_fem(Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))

print 'nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof())


P=m.pts()
r = range(0, m.nbpts());
INpid=compress(abs(P[0,:]+25) < 1e-4, r)
OUTpid=compress(abs(P[0,:]-25) < 1e-4, r)
TOPpid=compress(abs(P[2,:]-20) < 1e-4, r)
INfaces=m.faces_from_pid(INpid)
OUTfaces=m.faces_from_pid(OUTpid)
TOPfaces=m.faces_from_pid(TOPpid)

m.set_region(1, INfaces)
m.set_region(2, OUTfaces)
m.set_region(3, TOPfaces)
m.set_region(4, m.outer_faces())
m.region_substract(4, 1)
m.region_substract(4, 2)
m.region_substract(4, 3)

b0 = MdBrick('isotropic_linearized_elasticity',mim,mfu)
b0.set_param('lambda', 0);
b0.set_param('mu', viscosity);
b1 = MdBrick('linear incompressibility term', b0, mfp)
b2 = MdBrick('source term', b1, 1)
b2.set_param('source_term', [0,-10,0])

DIRICHLET_TYPE = 'penalized'
bconst = MdBrick('constraint', b2, 'augmented', 1)
bconst.set_constraints(ones((1,mfp.get('nbdof'))), [0]);
bdir1 = MdBrick('dirichlet', bconst, 1, mfu, DIRICHLET_TYPE);
bdir2 = MdBrick('dirichlet', bdir1, 2, mfu, DIRICHLET_TYPE);
bdir3 = MdBrick('dirichlet on normal component',
                  bdir2, 3, mfd, DIRICHLET_TYPE);
bdir4 = MdBrick('dirichlet', bdir3, 4, mfu, DIRICHLET_TYPE);


D = mfd.dof_nodes();
x = D[0,:]
y = D[1,:]
z = D[2,:]

bdir1.set_param('R', mfd, [9-(y*y+(z-6)*(z-6)),0*x,0*x])
bdir2.set_param('R', mfd, [9-(y*y+(z-6)*(z-6)),0*x,0*x])
bdir3.set_param('R', mfd, [0])
bdir4.set_param('R', mfd, [0,0,0])

mds=MdState(bdir4)
print 'running solve...'
bdir4.solve(mds, 'noisy', 'lsolver','superlu')
print 'solve done!'


VM=b0.von_mises(mds, mfe)
S=mds.state()
U=S[0:mfu.nbdof()] # extract the velocity
P=S[mfu.nbdof():(mfu.nbdof()+mfp.nbdof())] # and the pressure

mfu.save('tank_3D.mfu','with_mesh')
mfp.save('tank_3D.mfp','with_mesh')
U.tofile('tank_3D.U')
P.tofile('tank_3D.P')

mfe.save('tank_3D.mfe')
VM.tofile('tank_3D.VM')
#memstats()
