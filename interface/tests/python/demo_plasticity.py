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
"""  Demonstration for small deformations plasticty (without hardening),
  with optional graphical vizualisation (requires tvtk).

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""

import numpy as np

import getfem as gf

with_graphics=False
try:
    import getfem_tvtk
except:
    print("\n** Could NOT import getfem_tvtk -- graphical output disabled **\n")
    import time
    time.sleep(2)
    with_graphics=False

L=100
H=20

gf.util('trace_level', 1);

m=gf.Mesh('triangles grid', np.arange(0, L + 0.01, 4), np.arange(0, H + 0.01, 2))

mim=gf.MeshIm(m, gf.Integ('IM_TRIANGLE(6)'))
mim_data = gf.MeshImData(mim, -1, [2, 2])
mfu=gf.MeshFem(m,2)
mfxi=gf.MeshFem(m,1)
mfd=gf.MeshFem(m)
mfdu=gf.MeshFem(m)

mfu.set_fem(gf.Fem('FEM_PK(2,2)'))
mfxi.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,2)'))
mfd.set_fem(gf.Fem('FEM_PK(2,2)'))
mfdu.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,1)'))

Lambda=121150
Mu=80769
sigma_y=4000

P=m.pts()
pidleft=np.compress((abs(P[0,:])<1e-6), list(range(0, m.nbpts())))
pidright=np.compress((abs(P[0,:] - L)<1e-6), list(range(0, m.nbpts())))

fleft  = m.faces_from_pid(pidleft)
fright = m.faces_from_pid(pidright)

# assign boundary numbers
m.set_region(1,fleft)
m.set_region(2,fright)

md = gf.Model('real')
md.add_fem_variable('u', mfu)
md.add_fem_data('Previous_u', mfu)
md.add_im_data('Previous_Ep', mim_data)
md.add_fem_data('xi', mfxi)
md.add_fem_data('Previous_xi', mfxi)
md.add_initialized_data('lambda', Lambda)
md.add_initialized_data('mu', Mu)
md.add_initialized_data('sigma_y', sigma_y)
md.add_small_strain_elastoplasticity_brick(mim, 'plane strain Prandtl Reuss',
                                           'displacement only', 'u', 'xi', 'Previous_Ep', 'lambda', 'mu', 'sigma_y')
md.add_initialized_data('VolumicData', [0,0])
md.add_source_term_brick(mim, 'u', 'VolumicData')
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 1)


F=np.array([[0,-4.],[0, -5.], [0, -4.], [0, 2.], [0, 0]])
nbstep = F.shape[0]

print('nbstep:', nbstep)
for step in range(0, nbstep):
    print('step %d' % (step,))
    md.set_variable('VolumicData', [F[step,0],F[step,1]])
    md.solve('noisy', 'lsearch', 'simplest',  'alpha min', 0.8, 'max_iter', 100, 'max_res', 1e-6)
    U = md.variable('u')
    md.small_strain_elastoplasticity_next_iter\
      (mim, 'plane strain Prandtl Reuss', 'displacement only', 'u', 'xi', 'Previous_Ep', 'lambda', 'mu', 'sigma_y')
    
    VM = md.small_strain_elastoplasticity_Von_Mises\
         (mim, mfdu, 'plane strain Prandtl Reuss', 'displacement only', 'u', 'xi', 'Previous_Ep', 'lambda', 'mu', 'sigma_y')

    ERR = gf.compute_error_estimate(mfu, U, mim)

    if with_graphics:
        fig = getfem_tvtk.Figure()
        fig.show(mfu, deformation=U, deformation_scale=1, data=(mfdu,VM))
        print("Press Q to continue..")
        fig.loop()
