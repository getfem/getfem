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
"""  test levelset.

  This program is used to check that python-getfem is working. This is
  also a good example of use of python-getfem..

  $Id$
"""
import numpy as np
from scipy import rand

import getfem as gf

eps = 1.0/10

m = gf.Mesh('regular_simplices', np.arange(-1,1+eps,eps), np.arange(-1,1+eps,eps), 'degree', 2, 'noised')
#m = gf.Mesh('cartesian', np.arange(-1,1+eps,eps), np.arange(-1,1+eps,eps))

ls1 = gf.LevelSet(m, 2, 'y', 'x')
#ls1 = gf.LevelSet(m, 2, 'sqr(x) + sqr(y) - sqr(0.7)', 'x-.4')
#ls1 = gf.LevelSet(m, 2, 'x + y - 0.2') #, 'x-5')
#ls1 = gf.LevelSet(m, 2, 'x + y - 0.2', 'x-5')
#ls2 = gf.LevelSet(m, 2, '0.6*sqr(x) + sqr(y-0.1) - sqr(0.6)');
#ls3 = gf.LevelSet(m, 4, 'sqr(x) + sqr(y+.08) - sqr(0.05)');
ls2 = gf.LevelSet(m, 2, 'y+0.1', 'x')
ls3 = gf.LevelSet(m, 2, 'y-0.1', 'x')

mls = gf.MeshLevelSet(m)

mls.add(ls1)
if True:
  mls.sup(ls1)
  mls.add(ls1)
  mls.add(ls2)
  mls.add(ls2)
  mls.add(ls2)
  mls.add(ls3)
mls.adapt()

#print(mls.linked_mesh())

lls = mls.levelsets()

cm = mls.cut_mesh()

ctip = mls.crack_tip_convexes()

mf = gf.MeshFem(m)
mf.set_classical_fem(1)

mfls = gf.MeshFem('levelset',mls,mf)

gf.memstats()

nbd = mfls.nbdof()

if True:
  sl = gf.Slice(('none',), mls, 2);
  U = rand(1,nbd);
  sl.export_to_pos('slU.pos',mfls,U,'U')
  mfls.export_to_pos('U.pos',U,'U')
  cm.export_to_pos('cm.pos')
  m.export_to_pos('m.pos')
else:
  sl = gf.Slice(('none',), mls, 1);
  for i in range(nbd):
    U = np.zeros(nbd)
    U[i] = 1
    sl.export_to_pos('slU'+str(i)+'.pos',mfls,U,'U'+str(i))
    mfls.export_to_pos('U'+str(i)+'.pos',U,'U'+str(i))
  cm.export_to_pos('cm.pos')
  m.export_to_pos('m.pos')
