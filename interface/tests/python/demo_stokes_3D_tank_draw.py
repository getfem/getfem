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
"""  Display of the results of demo_stokes_3D_tank.py (to be run first),
  with the visualisation module getfem.plot (requires pyvista).

  Several slices are drawn in the same interactive scene: pressure on two
  cuts of the tank, velocity arrows, edges of the mesh and streamlines.

  $Id$
"""
import numpy as np

import getfem as gf
try:
    import getfem.getfem_plot as gfp
except ImportError:  # in the build tree, getfem is the module getfem.py, not the package
    import getfem_plot as gfp

mfu=gf.MeshFem('load','tank_3D.mfu')
m=mfu.linked_mesh()
mfp=gf.MeshFem('load','tank_3D.mfp',m)
U = np.fromfile('tank_3D.U', 'd')
P = np.fromfile('tank_3D.P', 'd')
prange = (-40, 40)

# pressure and velocity arrows on the skin of the quarter x > 0, y > 0
sl=gf.Slice(('boundary',('intersection',('planar',+1,[0,0,0],[0,1,0]),
                                        ('planar',+1,[0,0,0],[1,0,0]))),m,3)
pl = gfp.plot_slice(sl, mfp, P, name="p", clim=prange, backend="pyvista",
                    show_axes=False, show=False)
gfp.plot_slice(sl, mfu, U, name="u", faces=False, quiver=True, colorbar=False,
               show_axes=False, plotter=pl)

# pressure on the skin of the half y > 0, below z = 6
sl=gf.Slice(('boundary',('intersection',('planar',+1,[0,0,6],[0,0,-1]),
                                        ('planar',+1,[0,0,0],[0,1,0]))),m,3)
gfp.plot_slice(sl, mfp, P, name="p", clim=prange, colorbar=False,
               show_axes=False, plotter=pl)

# white edges of the mesh on the outer faces, except the top ones (region 3)
m.set_region(42, m.outer_faces())
m.region_subtract(42, 3)
sl2=gf.Slice(('boundary',('planar',+1,[0,0,0],[0,1,0])),m,3, m.region(42))
gfp.plot_slice(sl2, mesh_edges=True, faces=False, edge_color="white",
               show_axes=False, plotter=pl)

# streamlines starting on the axis of the tank
hh=np.array([[1, 5, 9, 12.5, 16, 19.5]])
H=np.concatenate((np.zeros((2,hh.size)), hh))
tsl=gf.Slice('streamlines', mfu, U, H)
gfp.plot_slice(tsl, line_color="white", plotter=pl,
               title="Stokes flow in a tank: pressure, velocity, streamlines")

pl.show()
