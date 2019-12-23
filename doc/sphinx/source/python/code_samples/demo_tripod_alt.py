# This is the "old" tripod demo, which uses the low level approach:
# building the linear system by hand, handling Dirichlet, calling the solver etc...

import numpy as np

import getfem as gf

# parameters
file_msh = 'tripod.GiD.msh'
degree = 1
linear = False
incompressible = False # ensure that degree > 1 when incompressible is on..
E = 1e3
Nu = 0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))

# create a Mesh object (importing)
m = gf.Mesh('import','gid',file_msh)
m.set('optimize structure')

# create a MeshFem object
mfu = gf.MeshFem(m,3) # displacement
mfd = gf.MeshFem(m,1) # data
mfe = gf.MeshFem(m,1) # for plot von-mises
# assign the FEM
mfu.set_fem(gf.Fem('FEM_PK(3,%d)' % (degree,)))
mfd.set_fem(gf.Fem('FEM_PK(3,0)'))
mfe.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,%d,0.01)' % (degree,)))

# build a MeshIm object
mim = gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(5)'))

print 'nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof())

# detect some boundary of the mesh
P = m.pts()
ctop = (abs(P[1,:] - 13) < 1e-6)
cbot = (abs(P[1,:] + 10) < 1e-6)
pidtop = np.compress(ctop, range(0, m.nbpts()))
pidbot = np.compress(cbot, range(0, m.nbpts()))
ftop = m.faces_from_pid(pidtop)
fbot = m.faces_from_pid(pidbot)
# create boundary region
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
m.set_region(NEUMANN_BOUNDARY,ftop)
m.set_region(DIRICHLET_BOUNDARY,fbot)

# assembly
nbd = mfd.nbdof()
print "nbd: ",nbd
F = gf.asm_boundary_source(NEUMANN_BOUNDARY, mim, mfu, mfd, np.repeat([[0],[-100],[0]],nbd,1))
print "F.shape: ",F.shape
print "mfu.nbdof(): ",mfu.nbdof()
print "np.repeat([[0],[-100],[0]],nbd,1).shape:",np.repeat([[0],[-100],[0]],nbd,1).shape

K = gf.asm_linear_elasticity(mim, mfu, mfd, np.repeat([Lambda], nbd), np.repeat([Mu], nbd))
print "K.info: ",K.info # Spmat instance
print "np.repeat([Lambda], nbd).shape:",np.repeat([Lambda], nbd).shape
print "np.repeat([Mu], nbd).shape:",np.repeat([Mu], nbd).shape

# handle Dirichlet condition
(H,R) = gf.asm_dirichlet(DIRICHLET_BOUNDARY, mim, mfu, mfd, mfd.eval('numpy.identity(3)'), mfd.eval('[0,0,0]'))
print "H.info: ",H.info # Spmat instance
print "R.shape: ",R.shape
print "mfd.eval('numpy.identity(3)').shape: ",mfd.eval('numpy.identity(3)').shape
print "mfd.eval('[0,0,0]').shape: ",mfd.eval('[0,0,0]').shape

(N,U0) = H.dirichlet_nullspace(R)
print "N.info: ",N.info # Spmat instance
print "U0.shape: ",U0.shape

Nt = gf.Spmat('copy',N)
Nt.transpose()
KK = Nt*K*N
FF = Nt*F # FF = Nt*(F-K*U0)

# solve ...
P = gf.Precond('ildlt',KK)
UU = gf.linsolve_cg(KK,FF,P)
print "UU.shape:",UU.shape
U = N*UU+U0
print "U.shape:",U.shape

# post-processing
sl = gf.Slice(('boundary',), mfu, degree)

# compute the Von Mises Stress
DU = gf.compute_gradient(mfu,U,mfe)
VM = np.zeros((DU.shape[2],),'d')
Sigma = DU

for i in range(DU.shape[2]):
  d = np.array(DU[:,:,i])
  E = (d+d.T)*0.5
  Sigma[:,:,i]=E
  VM[i] = np.sum(E**2) - (1./3.)*np.sum(np.diagonal(E))**2

print 'Von Mises range: ', VM.min(), VM.max()

# export results to VTK (you can use http://mayavi.sourceforge.net/ to view these results )
# i.e. with  "mayavi -d tripod.vtk -m BandedSurfaceMap -f WarpVector"
#sl.export_to_vtk('tripod.vtk', 'ascii', mfe, VM,'Von Mises Stress', mfu, U, 'Displacement')
#sl.export_to_vtk('tripod_edges.vtk','edges')

# export to OpenDX
#sl.export_to_dx('tripod.dx', 'ascii', mfe, VM,'Von Mises Stress')


# export the displacement and the stress tensor field
# can be viewed with mayavi -d ./tripod_ev.vtk -f WarpVector -m TensorGlyphs
SigmaSL = gf.compute_interpolate_on(mfe,Sigma,sl)
#sl.export_to_vtk('tripod_ev.vtk', mfu, U, 'Displacement', SigmaSL, 'stress')

#print 'You can view the tripod with (for example) mayavi:'
#print 'mayavi -d ./tripod.vtk -f WarpVector -m BandedSurfaceMap'

# export to Gmsh
sl.export_to_pos('tripod.pos', mfe, VM,'Von Mises Stress', mfu, U, 'Displacement')
sl.export_to_pos('tripod_ev.pos', mfu, U, 'Displacement', SigmaSL, 'stress')
