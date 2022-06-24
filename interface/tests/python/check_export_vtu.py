#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2020-2020 Tetsuo Koyama.
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


file_name = "check_mesh_ascii.vtk"
mesh.export_to_vtk(file_name, "ascii")
unstructured_grid = pv.read(file_name)
expected = pts
actual = unstructured_grid.points[:, 0]
np.testing.assert_equal(expected, actual, "export of mesh pts is not correct.")
expected = convex_connectivity
actual = unstructured_grid.cell_connectivity
np.testing.assert_equal(expected, actual, "export of mesh convex is not correct.")

file_name = "check_mesh_binary.vtk"
mesh.export_to_vtk(file_name)
unstructured_grid = pv.read(file_name)
expected = pts
actual = unstructured_grid.points[:, 0]
np.testing.assert_equal(expected, actual, "export of mesh pts is not correct.")
expected = convex_connectivity
actual = unstructured_grid.cell_connectivity
np.testing.assert_equal(expected, actual, "export of mesh convex is not correct.")

file_name = "check_mesh_ascii.vtu"
mesh.export_to_vtu(file_name, "ascii")
unstructured_grid = pv.read(file_name)
expected = pts
actual = unstructured_grid.points[:, 0]
np.testing.assert_equal(expected, actual, "export of mesh pts is not correct.")
expected = convex_connectivity
actual = unstructured_grid.cell_connectivity
np.testing.assert_equal(expected, actual, "export of mesh convex is not correct.")

file_name = "check_mesh_binary.vtu"
mesh.export_to_vtu(file_name)
unstructured_grid = pv.read(file_name)
expected = pts
actual = unstructured_grid.points[:, 0]
np.testing.assert_equal(expected, actual, "export of mesh pts is not correct.")
expected = convex_connectivity
actual = unstructured_grid.cell_connectivity
np.testing.assert_equal(expected, actual, "export of mesh convex is not correct.")

mfu = gf.MeshFem(mesh, 1)
mfu.set_classical_fem(1)
U1 = np.array([2.0, 1.0, 0.0])

file_name = "check_meshfem_ascii.vtk"
mfu.export_to_vtk(file_name, "ascii", U1, "U1")
unstructured_grid = pv.read(file_name)
expected = U1
actual = unstructured_grid.point_arrays["U1"]
np.testing.assert_equal(expected, actual, "export of U1 is not correct.")

file_name = "check_meshfem_binary.vtk"
mfu.export_to_vtk(file_name, U1, "U1")
unstructured_grid = pv.read(file_name)
expected = U1
actual = unstructured_grid.point_arrays["U1"]
np.testing.assert_equal(expected, actual, "export of U1 is not correct.")

file_name = "check_meshfem_ascii.vtu"
mfu.export_to_vtu(file_name, "ascii", U1, "U1")
unstructured_grid = pv.read(file_name)
expected = U1
actual = unstructured_grid.point_arrays["U1"]
np.testing.assert_equal(expected, actual, "export of U1 is not correct.")

file_name = "check_meshfem_binary.vtu"
mfu.export_to_vtu(file_name, U1, "U1")
unstructured_grid = pv.read(file_name)
expected = U1
actual = unstructured_grid.point_arrays["U1"]
np.testing.assert_equal(expected, actual, "export of U1 is not correct.")

sl = gf.Slice(("boundary",), mesh, 1)
U2 = np.array([3.0, 2.0, 1.0, 0.0])

file_name = "check_slice_ascii.vtk"
sl.export_to_vtk(file_name, "ascii", U2, "U2")
unstructured_grid = pv.read(file_name)
expected = U2
actual = unstructured_grid.point_arrays["U2"]
np.testing.assert_equal(expected, actual, "export of U2 is not correct.")

file_name = "check_slice_binary.vtk"
sl.export_to_vtk(file_name, U2, "U2")
unstructured_grid = pv.read(file_name)
expected = U2
actual = unstructured_grid.point_arrays["U2"]
np.testing.assert_equal(expected, actual, "export of U2 is not correct.")

file_name = "check_slice_ascii.vtu"
sl.export_to_vtu(file_name, "ascii", U2, "U2")
unstructured_grid = pv.read(file_name)
expected = U2
actual = unstructured_grid.point_arrays["U2"]
np.testing.assert_equal(expected, actual, "export of U2 is not correct.")

file_name = "check_slice_binary.vtu"
sl.export_to_vtu(file_name, U2, "U2")
unstructured_grid = pv.read(file_name)
expected = U2
actual = unstructured_grid.point_arrays["U2"]
np.testing.assert_equal(expected, actual, "export of U2 is not correct.")
