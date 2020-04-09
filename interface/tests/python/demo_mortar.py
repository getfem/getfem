#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2011-2020 Yves Renard.
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

""" Elasticity mortar problem test.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  It does a partition of the mesh into two disjoint areas, and then
  solves the linear elasticity problem with a mortar join on  the
  interface between the two areas
"""

import numpy as np

# import basic modules
import getfem as gf

# Parameters
NX = 9
radius = 0.3
xc = 0.5
yc = 0.5

# creation of a simple cartesian mesh
m = gf.Mesh('cartesian', np.arange(0,1+0.5/NX,1./NX), np.arange(0,1+0.5/NX,1./NX))

(pid,idx) = m.pid_from_cvid()

P = m.pts()

is_in_circle = (P[0,:]-xc)**2+(P[1,:]-yc)**2 <= radius**2

areap = np.zeros(idx.size-1)
for cv in range(idx.size-1):
   if all(is_in_circle[pid[idx[cv]:idx[cv+1]]]):
      areap[cv] = 1

mfu = gf.MeshFem(m,2)
mfd = gf.MeshFem(m,1)
mfm = gf.MeshFem(m,2)
mfdu= gf.MeshFem(m)

mim = gf.MeshIm(m,5)

mfu.set_fem(gf.Fem('FEM_QK(2,2)'))
mfd.set_fem(gf.Fem('FEM_QK(2,1)'))
mfm.set_fem(gf.Fem('FEM_QK(2,2)'))
mfdu.set_fem(gf.Fem('FEM_QK_DISCONTINUOUS(2,2)'))

mfu.set_dof_partition(areap)

b_in  = m.outer_faces(np.nonzero(areap==1))
b_out = m.outer_faces(np.nonzero(areap==0))
b_border = m.outer_faces()
b_out = np.array(tuple(set(tuple(r) for r in b_out.transpose())-
                       set(tuple(r) for r in b_border.transpose()))).transpose()

fleft = m.faces_from_pid(np.nonzero(abs(P[1])<1e-6))
fright = m.faces_from_pid(np.nonzero(abs(P[1]-1.)<1e-6))

# assign boundary numbers
m.set_region(1,fleft)
m.set_region(2,fright)

MORTAR_BOUNDARY_IN = 40
MORTAR_BOUNDARY_OUT = 41
m.set_region(MORTAR_BOUNDARY_IN, b_in)
m.set_region(MORTAR_BOUNDARY_OUT, b_out)

indm = mfm.basic_dof_on_region(MORTAR_BOUNDARY_OUT)
expr = 'M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,i)'
M =   gf.asm_boundary(MORTAR_BOUNDARY_IN, expr, mim, mfm, mfu)
M = M-gf.asm_boundary(MORTAR_BOUNDARY_OUT, expr, mim, mfm, mfu)
M = gf.Spmat('copy', M, indm, list(range(M.size()[1])))

md = gf.Model('real')
md.add_fem_variable('u', mfu);
md.add_initialized_data('lambda', [1])
md.add_initialized_data('mu', [1])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'lambda', 'mu')

md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 1)

F = mfd.eval('[0,y+2]')

md.add_initialized_fem_data('VolumicData', mfd, F)
md.add_source_term_brick(mim, 'u', 'VolumicData')
md.add_variable('mult_spec', indm.size)
md.add_constraint_with_multipliers('u', 'mult_spec', M, np.zeros(indm.size))


md.solve();
U = md.get('variable', 'u')

VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'lambda', 'mu', mfdu)

mfd.export_to_vtk('mortar.vtk', 'ascii', mfdu,  VM, 'Von Mises Stress', mfu, U, 'Displacement')
