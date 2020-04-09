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
from numpy import *

import getfem
import getfem_tvtk

mfu=getfem.MeshFem('load','tank_3D.mfu')
m=mfu.linked_mesh()
mfp=getfem.MeshFem('load','tank_3D.mfp',m)
U = fromfile('tank_3D.U', 'd')
P = fromfile('tank_3D.P', 'd')

sl=getfem.Slice(('boundary',('intersection',('planar',+1,[0,0,0],[0,1,0]),('planar',+1,[0,0,0],[1,0,0]))),m,3);

print("importing tvtk..")
print("import done")

fig = getfem_tvtk.Figure(gui='tvtk')

fig.show(sl, data=(mfp, P), vdata=(mfu,U), edges=False)

fig.show(sl, data=(mfp, P), edges=False)


old=fig.scalar_range()

sl=getfem.Slice(('boundary',('intersection',('planar',+1,[0,0,6],[0,0,-1]),
                             ('planar',+1,[0,0,0],[0,1,0]))),m,3);
fig.show(sl, data=(mfp, P), scalar_bar=True, edges=False)
fig.scalar_range((-40,40));

#print(fig.scalar_range())

m.set_region(42, m.outer_faces());
m.region_subtract(42, 3);

sl2=getfem.Slice(('boundary',('planar',+1,[0,0,0],[0,1,0])),m,3, m.region(42))
fig.show(sl2, faces=False, edges_color=(1,1,1))


hh=array([[1, 5, 9, 12.5, 16, 19.5]]);
H=concatenate((zeros((2,hh.size)), hh));

tsl=getfem.Slice('streamlines', mfu, U, H);

fig.show(tsl, tube_color=(1,1,1));

fig.loop()
