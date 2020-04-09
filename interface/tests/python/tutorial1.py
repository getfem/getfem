#!/usr/bin/env python
# Python GetFEM interface
#
# Copyright (C)  2015-2020 Julien Pommier.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#
############################################################################

import numpy as np

import getfem as gf

NX = 10
m = gf.Mesh('cartesian', np.arange(0,1+1./NX,1./NX), np.arange(0,1+1./NX,1./NX))
mf = gf.MeshFem(m,1) # create a meshfem of for a field of dimension 1
mf.set('fem',gf.Fem('FEM_QK(2,2)'))

print (gf.Fem('FEM_QK(2,2)').poly_str())

# mim=gf.MeshIm(m, gf.Integ('IM_EXACT_PARALLELEPIPED(2)')); // not allowed
mim=gf.MeshIm(m, gf.Integ('IM_GAUSS_PARALLELEPIPED(2, 4)'));


border = m.outer_faces()
m.set_region(42, border)  # create the region B42 (:-


md=gf.Model('real')
md.add_fem_variable('u', mf)
md.add_Laplacian_brick(mim, 'u')
R = mf.eval('(x-.5)*(x-.5) + (y-.5)*(y-.5) + x/5 - y/3')
md.add_initialized_fem_data('DirichletData', mf, R)
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf, 42, 'DirichletData')

md.variable_list()

md.solve()
U = md.variable('u')
