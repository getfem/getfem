#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2009-2020 Yves Renard, Luis Saavedra.
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
"""  This demo use levelset to impose (weakly) a Dirichlet condition on an
  implicit boundary defined by the zero of the levelset.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np
from scipy import rand, setdiff1d

import getfem as gf

eps = 1.0/40
ls_degree = 2

m = gf.Mesh('cartesian', np.arange(-.5,.5+eps,eps), np.arange(-.5,.5+eps,eps))
#m = gf.Mesh('triangles grid', np.arange(-.5,.5+eps,eps), np.arange(-.5,.5+eps,eps))
ls = gf.LevelSet(m,ls_degree)
ls2 = gf.LevelSet(m, ls_degree, 'with_secondary')

mf_ls = ls.mf()
mf_ls2 = ls2.mf()

P = mf_ls.basic_dof_nodes()
x = P[0,:]
y = P[1,:]
#ULS = ((x + 0.25)**2 + (y - 0.4)**2) - 0.05**2
#ULS = min(ULS, ((x - 0.25)**2 + (y - 0.4)**2) - 0.05**2)

ULS = 1000*np.ones(x.shape)

if True:
  for ix in range(0,5):
    for iy in range(0,5):
      xc = (ix/4) * 0.8 - 0.4
      yc = (iy/4) * 0.8 - 0.4
      if (iy%2)==1:
        xc = xc + 0.05
      else:
        xc = xc - 0.05
      R = 0.03 + 0.005*iy
      ULS = np.minimum(ULS, ((x - xc)**2 + (y - yc)**2) - R**2);
else:
  for i in range(8):
    xc = rand() - 0.5
    yc = rand() - 0.5
    R = rand() * 0.09 + 0.02
    ULS = np.minimum(ULS, ((x - xc)**2 + (y - yc)**2) - R**2);
ls.set_values(ULS)

ULS2 = 1000*np.ones(x.shape);
ULS2s = 1000*np.ones(x.shape);
for i in range(1):
  xc = 0.0 # rand() - 0.5
  yc = 0.0 # rand() - 0.5
  theta = np.pi/3 # np.pi*rand()
  n = [-np.sin(theta), np.cos(theta)]

  R = 0.19 #rand() * 0.09 + 0.02
  ULS2 = np.minimum(ULS2, ((x-xc)*n[0] + (y-yc)*n[1]))
  #ULS2s = np.minimum(ULS2s, ((x - xc)**2 + (y - yc)**2) - R**2)
  ULS2s = np.minimum(ULS2s, (abs(y - yc)+abs(x-xc) - R))

ls2.set_values(ULS2, ULS2s) # '-y-x+.2') # '(y-.2)^2 - 0.04')

mls = gf.MeshLevelSet(m)
mls.add(ls)
mls.add(ls2)
mls.adapt()
mls.cut_mesh().export_to_pos('ver.pos')

mim_bound = gf.MeshIm('levelset',mls,'boundary(a+b)', gf.Integ('IM_TRIANGLE(6)')) #, gf.Integ('IM_QUAD(5)'))
mim = gf.MeshIm('levelset',mls,'all(a+b)', gf.Integ('IM_TRIANGLE(6)'))
mim.set_integ(4)

mfu0 = gf.MeshFem(m,2)
mfu0.set_fem(gf.Fem('FEM_QK(2,3)'))

mfdu = gf.MeshFem(m,1)
mfdu.set_fem(gf.Fem('FEM_QK_DISCONTINUOUS(2,2)'))

mf_mult = gf.MeshFem(m,2)
mf_mult.set_fem(gf.Fem('FEM_QK(2,1)'))

A = gf.asm('volumic','V()+=comp()',mim_bound)

#mls.cut_mesh().export_to_pos('mls.pos','cut mesh')
#mf_ls.export_to_pos('mf_ls.pos',ULS,'ULS')

dof_out = mfu0.dof_from_im(mim)
cv_out = mim.convex_index()
cv_in = setdiff1d(m.cvid(),cv_out)

# mfu = gf.MeshFem('partial', mfu0, dof_out, cv_in)

md = gf.Model('real')
md.add_fem_variable('u', mfu0)
md.add_initialized_data('lambda', [1])
md.add_initialized_data('mu', [1])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'lambda', 'mu')
md.add_initialized_data('VolumicData', [0, 10])
md.add_source_term_brick( mim, 'u', 'VolumicData')
md.add_multiplier('mult_dir', mf_mult, 'u')
md.add_Dirichlet_condition_with_multipliers(mim_bound, 'u', 'mult_dir', -1)
md.solve()

U = md.variable('u')
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'lambda', 'mu', mfdu)

mfdu.export_to_pos('vm.pos', VM, 'Von Mises', mfu0, U, 'deformation')

mf_ls.export_to_pos('ls.pos',ls.values(0),'ls values 0')

print('You can view the solution with (for instance):')
print('gmsh vm.pos ls.pos')
