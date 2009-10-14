#This is the "modern" tripod demo, which uses the getfem model bricks

import getfem as gf
import numpy as np

# parameters
file_msh = 'tripod.GiD.msh'
degree = 2
linear = False
incompressible = False # ensure that degree > 1 when incompressible is on..
E = 1e3
Nu = 0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))

# create a Mesh object (importing)
m = gf.Mesh('import','gid',file_msh)
m.set('optimize_structure')

# create a MeshFem object
mfu = gf.MeshFem(m,3) # displacement
mfp = gf.MeshFem(m,1) # pressure
mfe = gf.MeshFem(m,3) # for plot displacement
mff = gf.MeshFem(m,1) # for plot von-mises
# assign the FEM
mfu.set_fem(gf.Fem('FEM_PK(3,%d)' % (degree,)))
mfp.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,0)'))
mfe.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))
mff.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))

# build a MeshIm object
mim = gf.MeshIm(m,gf.Integ('IM_TETRAHEDRON(5)'))

print 'nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(),m.nbpts(),mfu.qdim(),mfu.fem()[0].char(),mfu.nbdof())

# detect some boundary of the mesh
P = m.pts()
ctop = (abs(P[1,:] - 13) < 1e-6)
cbot = (abs(P[1,:] + 10) < 1e-6)
pidtop = np.compress(ctop,range(0,m.nbpts()))
pidbot = np.compress(cbot,range(0,m.nbpts()))
ftop = m.faces_from_pid(pidtop)
fbot = m.faces_from_pid(pidbot)
# create boundary region
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
m.set_region(NEUMANN_BOUNDARY,ftop)
m.set_region(DIRICHLET_BOUNDARY,fbot)

# the model bricks
if linear:
  b0 = gf.MdBrick('isotropic_linearized_elasticity',mim,mfu)
  b0.set_param('lambda',Lambda)
  b0.set_param('mu',Mu)
  if (incompressible):
    b1 = gf.MdBrick('linear_incompressibility term',b0,mfp)
  else:
    b1 = b0
else:
  # large deformation with a linearized material law.. not a very good choice!
  if (incompressible):
    b0 = gf.MdBrick('nonlinear_elasticity',mim,mfu,'Mooney_Rivlin')
    b0.set_param('params',[Lambda,Mu])
    b1 = gf.MdBrick('nonlinear_elasticity_incompressibility_term',b0,mfp)
  else:
    b0 = gf.MdBrick('nonlinear_elasticity',mim,mfu,'SaintVenant_Kirchhoff')
    #b0 = gf.MdBrick('nonlinear_elasticity',mim,mfu,'Ciarlet_Geymonat')
    b0.set_param('params',[Lambda,Mu])
    b1 = b0

b2 = gf.MdBrick('source_term',b1,NEUMANN_BOUNDARY)
b2.set_param('source_term',[0,-10,0])
b3 = gf.MdBrick('dirichlet',b2,DIRICHLET_BOUNDARY,mfu,'penalized')

# create model state
mds = gf.MdState(b3)
# running solve...
b3.solve(mds,'noisy','lsolver','superlu')

# extracted solution
U = mds.state()

# post-processing
VM = b0.von_mises(mds,mff)

# export U and VM in a pos file
sl = gf.Slice(('boundary',),mfu,1)
sl.export_to_pos('sol.pos', mfu, U, 'Displacement', mff, VM, 'Von Mises Stress')

# save solution
#mfu.save('tripod.mf','with_mesh')
#U.tofile('tripod.U')
#mff.save('tripod.mff')
#VM.tofile('tripod.VM')
#gf.memstats()
