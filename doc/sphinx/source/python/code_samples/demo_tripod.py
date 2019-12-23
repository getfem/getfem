#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2004-2017 Yves Renard, Julien Pommier.
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
"""  This is the "modern" tripod demo, which uses the getfem model bricks
  importing the mesh.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

  $Id$
"""

import numpy as np

import getfem as gf

with_graphics=True
try:
    import getfem_tvtk
except:
    print("\n** Could NOT import getfem_tvtk -- graphical output disabled **\n")
    import time
    time.sleep(2)
    with_graphics=False


m=gf.Mesh('import','gid','../meshes/tripod.GiD.msh')
print('done!')
mfu=gf.MeshFem(m,3) # displacement
mfp=gf.MeshFem(m,1) # pressure
mfd=gf.MeshFem(m,1) # data
mim=gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(5)'))
degree = 2
linear = False
incompressible = False # ensure that degree > 1 when incompressible is on..

mfu.set_fem(gf.Fem('FEM_PK(3,%d)' % (degree,)))
mfd.set_fem(gf.Fem('FEM_PK(3,0)'))
mfp.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,0)'))

print('nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof()))

P=m.pts()
print('test', P[1,:])
ctop=(abs(P[1,:] - 13) < 1e-6)
cbot=(abs(P[1,:] + 10) < 1e-6)
pidtop=np.compress(ctop, list(range(0, m.nbpts())))
pidbot=np.compress(cbot, list(range(0, m.nbpts())))

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


md = gf.Model('real')
md.add_fem_variable('u', mfu)
if linear:
    md.add_initialized_data('cmu', Mu)
    md.add_initialized_data('clambda', Lambda)
    md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambda', 'cmu')
    if incompressible:
        md.add_fem_variable('p', mfp)
        md.add_linear_incompressibility_brick(mim, 'u', 'p')
else:
    md.add_initialized_data('params', [Lambda, Mu]);
    if incompressible:
        lawname = 'Incompressible Mooney Rivlin';
        md.add_finite_strain_elasticity_brick(mim, lawname, 'u', 'params')
        md.add_fem_variable('p', mfp);
        md.add_finite_strain_incompressibility_brick(mim, 'u', 'p');
    else:
        lawname = 'SaintVenant Kirchhoff';
        md.add_finite_strain_elasticity_brick(mim, lawname, 'u', 'params');
  

md.add_initialized_data('VolumicData', [0,-1,0]);
md.add_source_term_brick(mim, 'u', 'VolumicData');

# Attach the tripod to the ground
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 2);

print('running solve...')
md.solve('noisy', 'max iter', 1);
U = md.variable('u');
print('solve done!')


mfdu=gf.MeshFem(m,1)
mfdu.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(3,1)'))
if linear:
  VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u','clambda','cmu', mfdu);
else:
  VM = md.compute_finite_strain_elasticity_Von_Mises(lawname, 'u', 'params', mfdu);

# post-processing
sl=gf.Slice(('boundary',), mfu, degree)

print('Von Mises range: ', VM.min(), VM.max())

# export results to VTK
sl.export_to_vtk('tripod.vtk', 'ascii', mfdu,  VM, 'Von Mises Stress', mfu, U, 'Displacement')
sl.export_to_pos('tripod.pos', mfdu, VM, 'Von Mises Stress', mfu, U, 'Displacement')

gf.memstats()

print('You can view the tripod with (for example) mayavi:')
print('mayavi2 -d tripod.vtk -f WarpVector -m Surface')
print('or')
print('gmsh tripod.pos')

# mfu.save('tripod.mf', 'with_mesh')
# U.tofile('tripod.U')
# mfdu.save('tripod.mfe')
# VM.tofile('tripod.VM')

if with_graphics:
  fig = getfem_tvtk.Figure()
  fig.show(mfu, deformation=U, data=(mfdu,VM), deformation_scale='20%')
  print("Press Q to continue..")
  fig.set_colormap('tripod')
  fig.loop()
