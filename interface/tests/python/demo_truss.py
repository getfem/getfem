#!/usr/bin/env python
# coding: utf-8
# Python GetFEM++ interface
#
# Copyright (C) 2020-2020 Tetsuo Koyama.
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
""" 1D Truss problem test

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

  |-> x
//|        EA         |-> u0
//|---------+---------+=====> P
//|
  |<---L--->|<---L--->|

  Exact solution of displacement u0 is u0 = (P*x)/(E*A).

  $Id$
"""
# Import basic modules
import getfem as gf
import numpy as np

# Parameters
E = 21E6               # Young Modulus (N/cm^2)
A = 900.0              # Area (cm)
L = 200.0              # Length of element (cm)
P = 2.0                # Force

# Create a 1D mesh
mesh = gf.Mesh("cartesian", np.arange(0.0, 2*L, L))

# Create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu = gf.MeshFem(mesh, 1)
mfrhs = gf.MeshFem(mesh, 1)
# assign the Classical Fem
elements_degree = 1
mfu.set_classical_fem(elements_degree)
mfrhs.set_classical_fem(elements_degree)

#  Integration method used
mim = gf.MeshIm(mesh, pow(elements_degree,2))

# Boundary selection
fleft = mesh.outer_faces_with_direction(v=[-1.0], angle=0.01)
fright = mesh.outer_faces_with_direction(v=[+1.0], angle=0.01)

# Mark it as boundary
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
mesh.set_region(NEUMANN_BOUNDARY, fright)
mesh.set_region(DIRICHLET_BOUNDARY, fleft)

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfu.eval("(" + str(P) + "*x)/(" + str(E) +"*" + str(A) + ")")

# Interpolate the source term
F = mfrhs.eval(str(P))

# Model
md = gf.Model("real")

# Main unknown
md.add_fem_variable("u", mfu)

# Truss term on u
md.add_Truss_brick(mim, "u", "E", "A")

# Force term
md.add_initialized_fem_data("ForceData", mfrhs, F)
md.add_source_term_brick(mim, "u", "ForceData")

# Dirichlet condition on the boundary.
md.add_Dirichlet_condition_with_multipliers(mim, "u", elements_degree - 1, DIRICHLET_BOUNDARY)

# Assembly of the linear system and solve.
md.solve()

# Main unknown
U = md.variable("u")
L2error = gf.compute(mfu, U-Ue, "L2 norm", mim)
H1error = gf.compute(mfu, U-Ue, "H1 norm", mim)
print("Error in L2 norm : ", L2error)
print("Error in H1 norm : ", H1error)

# Export data
mfu.export_to_vtk("truss.vtk", Ue,"Exact solution",
                               U,"Computed solution")
print("You can view the solution with (for example):")
print("You can view solutions with for instance:\nmayavi2 -d truss.vtk -f WarpVector -m Surface")


if (H1error > 1e-3):
    print("Error too large !")
    exit(1)
