#!/usr/bin/env python
# coding: utf-8
# Python GetFEM interface
#
# Copyright (C) 2019-2020 Tetsuo Koyama.
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
"""  2D Poisson problem test in unit disk \Omega.
  -\Delta u=1 in \Omega, u=0 on \delta\Omega

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
# Import basic modules
import getfem as gf

## Parameters
h = 0.1 # approximate diameter of the elements.

# Create a unit disk mesh
mo = gf.MesherObject('ball', [1.0, 1.0], 1.0)
mesh = gf.Mesh('generate', mo, h, 2)
mesh.translate([-1.0, -1.0])

# Create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu = gf.MeshFem(mesh, 1)
mfrhs = gf.MeshFem(mesh, 1)
# assign the Classical Fem
elements_degree = 2
mfu.set_classical_fem(elements_degree)
mfrhs.set_classical_fem(elements_degree)

#  Integration method used
mim = gf.MeshIm(mesh, pow(elements_degree,2))

# Boundary selection
flst = mesh.outer_faces()

# Mark it as boundary
DIRICHLET_BOUNDARY = 1
mesh.set_region(DIRICHLET_BOUNDARY, flst)

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfu.eval('(1-x*x-y*y)/4')

# Interpolate the source term
F = mfrhs.eval('1')

# Model
md = gf.Model('real')

# Main unknown
md.add_fem_variable('u', mfu)

# Laplacian term on u
md.add_Laplacian_brick(mim, 'u')

# Volumic source term
md.add_initialized_fem_data('VolumicData', mfrhs, F)
md.add_source_term_brick(mim, 'u', 'VolumicData')

# Dirichlet condition on the boundary.
md.add_Dirichlet_condition_with_multipliers(mim, 'u', elements_degree - 1, DIRICHLET_BOUNDARY)

# Assembly of the linear system and solve.
md.solve()

# Main unknown
U = md.variable('u')
L2error = gf.compute(mfu, U-Ue, 'L2 norm', mim)
H1error = gf.compute(mfu, U-Ue, 'H1 norm', mim)
print('Error in L2 norm : ', L2error)
print('Error in H1 norm : ', H1error)

# Export data
mfu.export_to_pos('unit_disk.pos', Ue,'Exact solution',
                                    U,'Computed solution')
print('You can view the solution with (for example):')
print('gmsh unit_disk.pos')


if (H1error > 1e-3):
    print('Error too large !')
    exit(1)
