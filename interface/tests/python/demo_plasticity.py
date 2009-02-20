# demonstration for small deformations plasticty
# with optional graphical vizualisation (requires tvtk)

from getfem import *
from numpy import *

import locale # (bug #13014)
with_graphics=True
try:
    import getfem_tvtk
except:
    print "\n** Could NOT import getfem_tvtk -- graphical output disabled **\n"
    import time
    time.sleep(2)
    with_graphics=False
locale.setlocale(locale.LC_NUMERIC,('en_US','UTF8')) # (bug #13014)

L=100
H=20

m=Mesh('triangles grid', arange(0, L + 0.01, 4), arange(0, H + 0.01, 2))

mim=MeshIm(m, Integ('IM_TRIANGLE(6)'))
mfu=MeshFem(m,2)
mfd=MeshFem(m)
mf0=MeshFem(m)
mfdu=MeshFem(m)

mfu.set_fem(Fem('FEM_PK(2,1)'))
mfd.set_fem(Fem('FEM_PK(2,1)'))
mf0.set_fem(Fem('FEM_PK(2,0)'))
mfdu.set_fem(Fem('FEM_PK_DISCONTINUOUS(2,1)'));

Lambda=121150
Mu=80769
von_mises_threshold=8000

P=m.pts()
pidleft=compress((abs(P[0,:])<1e-6), range(0, m.nbpts()))
pidright=compress((abs(P[0,:] - L)<1e-6), range(0, m.nbpts()))

fleft  = m.faces_from_pid(pidleft)
fright = m.faces_from_pid(pidright)

# assign boundary numbers
m.set_region(1,fleft)
m.set_region(2,fright)

b0=MdBrick('small deformations plasticity',mim,mfu, von_mises_threshold)
b0.set_param('lambda',Lambda)
b0.set_param('mu',Mu)
b1=MdBrick('generalized dirichlet',b0,1)
b2=MdBrick('source term',b1,2)

mds=MdState(b2)

F=array([[0,-200],[0, -300], [0, 200], [0, 0]])
nbstep = F.shape[0]

dd=mf0.dof_from_cvid()

print 'nbstep:', nbstep
for step in range(0, nbstep):
    print 'step %d' % (step,)
    b2.set_param('source_term', mfd, [F[step,0],F[step,1]])
    b2.solve(mds, 'very noisy', 'max_iter', 1000, 'max_res', 1e-6)

    U=mds.state()[0:mfu.nbdof()]
    VM = b0.von_mises(mds, mfdu)
  
    #subplot(2,1,1);
    #gf_plot(mfdu,VM,'deformed_mesh','on', 'deformation',U,'deformation_mf',mfu,'refine', 4, 'deformation_scale',1); 
    #colorbar;
    #caxis([0 10000]);
    
    ERR=compute_error_estimate(mfu,U,mim)
    #E=ERR; E(dd)=ERR;
    #subplot(2,1,2);
    #gf_plot(mf0, E, 'mesh','on', 'refine', 1); colorbar;

    if with_graphics:
        fig = getfem_tvtk.Figure()
        fig.show(mfu, deformation=U, deformation_scale=1, data=(mfdu,VM))
        print "Press Q to continue.."
        fig.loop()
