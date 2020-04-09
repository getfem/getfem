#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2004-2020 Julien Pommier.
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

import numpy as np

# import basic modules
import getfem as gf

# creation of a simple cartesian mesh
m = gf.Mesh('cartesian', np.arange(0,1.1,0.1), np.arange(0,1.1,0.1))

# create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = gf.MeshFem(m, 1)
# assign the Q2 fem to all convexes of the MeshFem
mf.set_fem(gf.Fem('FEM_QK(2,2)'))

# view the expression of its basis functions on the reference convex
print(gf.Fem('FEM_QK(2,2)').poly_str())

# an exact integration will be used
mim = gf.MeshIm(m, gf.Integ('IM_GAUSS_PARALLELEPIPED(2,4)'))

# detect the border of the mesh
border = m.outer_faces()
# mark it as boundary #42
m.set_region(42, border)

# empty real model
md = gf.Model('real')

# declare that "u" is an unknown of the system
# on the finite element method `mf`
md.add_fem_variable('u', mf)

# add generic elliptic brick on "u"
md.add_Laplacian_brick(mim, 'u');

# add Dirichlet condition
g = mf.eval('x*(x-1) - y*(y-1)')
md.add_initialized_fem_data('DirichletData', mf, g)
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf, 42, 'DirichletData')

# add source term
#f = mf.eval('0')
#md.add_initialized_fem_data('VolumicData', mf, f)
#md.add_source_term_brick(mim, 'u', 'VolumicData')

# solve the linear system
md.solve()

# extracted solution
u = md.variable('u')

# export computed solution
mf.export_to_pos('u.pos',u,'Computed solution')
