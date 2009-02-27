# This is the "old" tripod demo, which uses the low level approach:
# building the linear system by hand, handling Dirichlet, calling the solver etc...

from getfem import *
from numpy import *

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
m = Mesh('import','gid',file_msh)
m.set('optimize structure')

# create a MeshFem object
mfu = MeshFem(m,3) # displacement
mfd = MeshFem(m,1) # data
mfe = MeshFem(m,1) # for plot von-mises
# assign the FEM
mfu.set_fem(Fem('FEM_PK(3,%d)' % (degree,)))
mfd.set_fem(Fem('FEM_PK(3,0)'))
mfe.set_fem(Fem('FEM_PK_DISCONTINUOUS(3,%d,0.01)' % (degree,)))

# build a MeshIm object
mim = MeshIm(m, Integ('IM_TETRAHEDRON(5)'))

print 'nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof())

# detect some boundary of the mesh
P = m.pts()
ctop = (abs(P[1,:] - 13) < 1e-6)
cbot = (abs(P[1,:] + 10) < 1e-6)
pidtop = compress(ctop, range(0, m.nbpts()))
pidbot = compress(cbot, range(0, m.nbpts()))
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
F = asm_boundary_source(NEUMANN_BOUNDARY, mim, mfu, mfd, repeat([[0],[-100],[0]],nbd,1))
print "F.shape: ",F.shape
print "mfu.nbdof(): ",mfu.nbdof()
print "repeat([[0],[-100],[0]],nbd,1).shape:",repeat([[0],[-100],[0]],nbd,1).shape

K = asm_linear_elasticity(mim, mfu, mfd, repeat([Lambda], nbd), repeat([Mu], nbd))
print "K.info: ",K.info # Spmat instance
print "repeat([Lambda], nbd).shape:",repeat([Lambda], nbd).shape
print "repeat([Mu], nbd).shape:",repeat([Mu], nbd).shape

# handle Dirichlet condition
(H,R) = asm_dirichlet(DIRICHLET_BOUNDARY, mim, mfu, mfd, mfd.eval('numpy.identity(3)'), mfd.eval('[0,0,0]'))
print "H.info: ",H.info # Spmat instance
print "R.shape: ",R.shape
print "mfd.eval('numpy.identity(3)').shape: ",mfd.eval('numpy.identity(3)').shape
print "mfd.eval('[0,0,0]').shape: ",mfd.eval('[0,0,0]').shape

(N,U0) = H.dirichlet_nullspace(R)
print "N.info: ",N.info # Spmat instance
print "U0.shape: ",U0.shape

Nt = Spmat('copy',N)
Nt.transpose()
KK = Nt*K*N
FF = Nt*F # FF = Nt*(F-K*U0)

# solve ...
P = Precond('ildlt',KK)
UU = linsolve_cg(KK,FF,P)
print "UU.shape:",UU.shape
U = N*UU+U0
print "U.shape:",U.shape

u = U*1.0
u.shape=(-1,3)
mfu.export_to_pos("ver.pos",u.transpose(),"ver")

# post-processing
sl=Slice(('boundary',), mfu, degree)

# compute the Von Mises Stress
DU = compute_gradient(mfu,U,mfe)
VM = zeros((DU.shape[2],),'d')
Sigma = DU

for i in range(0, DU.shape[2]):
  d = array(DU[:,:,i])
  E = (d+transpose(d))*0.5
  Sigma[:,:,i]=E
  VM[i] = sum(E**2) - (1./3.)*sum(diagonal(E))**2

print 'Von Mises range: ', VM.min(), VM.max()

# export results to VTK (you can use http://mayavi.sourceforge.net/ to view these results )
# i.e. with  "mayavi -d tripod.vtk -m BandedSurfaceMap -f WarpVector"
sl.export_to_vtk('tripod.vtk', 'ascii',mfe,  VM,'Von Mises Stress', mfu, U, 'Displacement')
sl.export_to_vtk('tripod_edges.vtk','edges')

# export to OpenDX
sl.export_to_dx('tripod.dx', 'ascii', mfe, VM,'Von Mises Stress')

# export the displacement and the stress tensor field
# can be viewed with mayavi -d ./tripod_ev.vtk -f WarpVector -m TensorGlyphs
SigmaSL = compute_interpolate_on(mfe,Sigma,sl)
sl.export_to_vtk('tripod_ev.vtk', mfu, U, 'Displacement', SigmaSL, 'stress')

print 'You can view the tripod with (for example) mayavi:'
print 'mayavi -d ./tripod.vtk -f WarpVector -m BandedSurfaceMap'
