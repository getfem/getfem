#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2004-2020 Yves Renard, Julien Pommier.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#
############################################################################
"""  This is the "old" tripod demo, which uses the low level approach:
  building the linear system by hand, handling Dirichlet, calling the
  solver etc...

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
from numpy import *

from getfem import *

print('importing the mesh..',)
m=Mesh('import','gid','../meshes/tripod.GiD.msh')
print('done!')
mfu=MeshFem(m,3)
mfd=MeshFem(m,1)
mfe=MeshFem(m,1)
mim=MeshIm(m, Integ('IM_TETRAHEDRON(5)'))
degree = 1	
mfu.set_fem(Fem('FEM_PK(3,%d)' % (degree,)));
mfe.set_fem(Fem('FEM_PK_DISCONTINUOUS(3,%d,0.01)' % (degree,)))
mfd.set_fem(Fem('FEM_PK(3,0)'))

print('nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof()))

P=m.pts()

ctop=(abs(P[1,:] - 13) < 1e-6);
cbot=(abs(P[1,:] + 10) < 1e-6);
pidtop=compress(ctop, list(range(0, m.nbpts())))
pidbot=compress(cbot, list(range(0, m.nbpts())))

ftop=m.faces_from_pid(pidtop)
fbot=m.faces_from_pid(pidbot)
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2

m.set_region(NEUMANN_BOUNDARY,ftop)
m.set_region(DIRICHLET_BOUNDARY,fbot)

E=1e3
Nu=0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu =E/(2*(1+Nu))

F = asm_boundary_source(NEUMANN_BOUNDARY, mim, mfu,mfd,repeat([[0],[-100],[0]], mfd.nbdof(),1))
K = asm_linear_elasticity(mim, mfu, mfd, repeat([Lambda], mfd.nbdof()),
                          repeat([Mu], mfd.nbdof()));

# handle Dirichlet condition
(H,R)=asm_dirichlet(DIRICHLET_BOUNDARY, mim, mfu, mfd,
                    mfd.eval('identity(3)',globals(),locals()),
                    mfd.eval('[0,0,0]'));
(N,U0)=H.dirichlet_nullspace(R)
Nt=Spmat('copy',N); Nt.transpose()
KK=Nt*K*N
FF=Nt*F

# solve ...
print("preconditioner..")
P=Precond('ildlt',KK)
print("solving..."),
UU=linsolve_cg(KK,FF,P)
print("done!")
U=N*UU+U0

# post-processing
sl=Slice(('boundary',), mfu, degree)

# compute the Von Mises Stress
DU=compute_gradient(mfu,U,mfe)
VM=zeros((DU.shape[2],),'d')
Sigma=DU

for i in range(0, DU.shape[2]):
  d = array(DU[:,:,i]);
  E = (d+transpose(d))*0.5
  Sigma[:,:,i]=E
  VM[i] = sum(E.ravel()**2) - (1./3.)*sum(diagonal(E))**2

print('Von Mises range: ', VM.min(), VM.max())

# export results to VTK you can use
# i.e. with  "mayavi2 -d tripod.vtk -f WarpScalar -m Surface"
sl.export_to_vtk('tripod.vtk', 'ascii',mfe,  VM,'Von Mises Stress', mfu, U, 'Displacement')

sl.export_to_vtk('tripod_edges.vtk','edges')

# export to OpenDX
sl.export_to_dx('tripod.dx', 'ascii', mfe, VM, 'Von Mises Stress')
# export the displacement and the stress tensor field
# can be viewed with mayavi -d ./tripod_ev.vtk -f WarpVector -m TensorGlyphs
SigmaSL = compute_interpolate_on(mfe, Sigma, sl);
sl.export_to_vtk('tripod_ev.vtk', mfu, U, 'Displacement', SigmaSL, 'stress')
# export to Gmsh POS
sl.export_to_pos('tripod.pos', mfe, VM, 'Von Mises Stress', mfu, U, 'Displacement')

print('You can view the tripod with (for example) mayavi:')
print('mayavi2 -d tripod.vtk -f WarpScalar -m Surface')
print('or')
print('gmsh tripod.pos')
