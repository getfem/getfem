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
import vtk
import getfem as gf

m0 = gf.Mesh("cartesian", [0, 1])

filenames = ["check_m0_ascii.vtu", "check_m0_binary.vtu"]

m0.export_to_vtu(filenames[0], "ascii")
filenames.append(filenames[0])

m0.export_to_vtu(filenames[1])
filenames.append(filenames[1])

for filename in filenames:
    print(filename)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()
    cell_data = output.GetCellData()
    nbpts = output.GetNumberOfPoints()
    nbcvs = output.GetNumberOfCells()
    array_name = cell_data.GetArrayName(0)
    assert nbpts == m0.nbpts(), "Number of points is not correct."
    assert nbcvs == m0.nbcvs(), "Number of cells is not correct."
