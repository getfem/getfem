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
import numpy as np

import getfem as gf

m0 = gf.Mesh('cartesian',[0,1,2,3],[0,1,2],[-3,-2])

m0.add_convex(gf.GeoTrans('GT_QK(2,2)'),[[0,0,0,.4,.6,.5,1.2,1,1],
                                         [0,.3,1,0,.5,1,-.1,.5,1.1],
                                         [0,0,0,0,0,0,0,0,0]])

m0.add_convex(gf.GeoTrans('GT_PK(2,2)'),[[2, 2,2,2.6,2.5,3.2],
                                         [0,.3,1,  0, .5,  0],
                                         [0, 0,0,  0,  0,  0]])

m0.add_convex(gf.GeoTrans('GT_PK(1,3)'),[[3.1, 2.8, 3.2, 3.7],
                                         [  0,  .3,  .7, 1.3],
                                         [  0,   0,   0,   0]])

m0.add_convex(gf.GeoTrans('GT_PRISM(3,1)'),[[0, 1, 0, 0, 1, 0],
                                            [0, 0, 1, 0, 0, 1],
                                            [0, 0, 0, 1, 1, 1]]);


m0.add_convex(gf.GeoTrans('GT_PK(3,2)'),[np.array([0, .5, 1, 0, .5, 0, 0, .5, 0, 0])-1.5,
                                         np.array([0, 0, 0, .5, .5, 1, 0, 0, .5, 0])-1,
                                         np.array([0, 0, 0, 0, 0, 0, .5, .5, .5, 1])-1.1])

m0.add_convex(gf.GeoTrans('GT_QK(3,2)'),[np.array([0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1])-1.5,
                                         np.array([0, 0, 0, .5, .5, .5, 1, 1, 1, 0, 0, 0, .5, .5, .5, 1, 1, 1, 0, 0, 0, .5, .5, .5, 1, 1, 1])-1,
                                         np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1, 1, 1, 1, 1, 1, 1, 1, 1])])

m0.add_convex(gf.GeoTrans('GT_PRISM_INCOMPLETE_P2'),[[1,.5,0,.5,0,0,  1,  0,  0,1,.5,0,.5,0,0],
                                                   [0, 1,2, 0,1,0,  0,  2,  0,0, 1,2, 0,1,0],
                                                   [2, 2,2, 2,2,2,3.5,3.5,3.5,4,4, 4,4, 4,4]]);


m1 = gf.Mesh('cartesian',[0,1,2,3],[0,1,2],[-3,-2])

mf0 = gf.MeshFem(m0); mf0.set_classical_fem(1)
mf1 = gf.MeshFem(m1); mf1.set_classical_fem(1)

sl = gf.Slice(('boundary',),m0,6)
U = np.random.standard_normal(mf0.nbdof())

# VTK:
m0.export_to_vtk('check_export0.vtk','quality')
m1.export_to_vtk('check_export1.vtk','quality')
mf0.export_to_vtk('check_export2.vtk','ascii')
mf1.export_to_vtk('check_export3.vtk','ascii')



# DX:
try:
    m0.export_to_dx('check_export0.dx')
except RuntimeError as detail:
    print(detail)

m1.export_to_dx('check_export0.dx','ascii','edges')
m1.export_to_dx('check_export0.dx','ascii','append')
sl.export_to_dx('check_export1.dx','ascii')
sl.export_to_dx('check_export1.dx','append','edges')
sl.export_to_dx('check_export1.dx','append','serie','rndserie',mf0,U)
sl.export_to_dx('check_export1.dx','append','serie','rndserie',mf0,U)

# POS:
m0.export_to_pos('check_export0.pos','m0')
m1.export_to_pos('check_export1.pos','m1')
mf0.export_to_pos('check_export2.pos','mf0')
mf1.export_to_pos('check_export3.pos','mf1')
sl.export_to_pos('check_export4.pos','sl',mf0,U,'U')
