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
"""  test export.

  This program is used to check that python-getfem is working. This is
  also a good example of use of python-getfem..

  $Id$
"""
import getfem as gf
import numpy as np
import sys
try:
    import pyvista as pv
except:
    print("\n\n** Could not load pyvista. Did you install it ?\n")
    print("   ( https://docs.pyvista.org/getting-started/installation.html ) **\n\n")
    sys.exit()

convex_connectivity = np.array([0, 1, 1, 2])
mesh = gf.Mesh("cartesian", [0.0, 1.0, 2.0])
pts = mesh.pts()[0]


filenames = [
    "check_mesh_ascii.vtk",
    "check_mesh_binary.vtk",
    "check_mesh_ascii.vtu",
    "check_mesh_binary.vtu",
]

mesh.export_to_vtk(filenames[0], "ascii")
mesh.export_to_vtk(filenames[1])
mesh.export_to_vtu(filenames[2], "ascii")
mesh.export_to_vtu(filenames[3])

for filename in filenames:
    print(filename)
    unstructured_grid = pv.read(filename)

    expected = pts
    actual = unstructured_grid.points[:, 0]
    np.testing.assert_equal(expected, actual, "export of mesh pts is not correct.")

    expected = convex_connectivity
    actual = unstructured_grid.cell_connectivity
    np.testing.assert_equal(expected, actual, "export of mesh convex is not correct.")
