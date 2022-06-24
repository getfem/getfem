#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2021-2021 Tetsuo Koyama.
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
""" 1D Truss problem test

  This program is used to check that python-getfem is working. This is
  also a good example of use of python-getfem..

  |-> x
//|        EA         |-> u0
//|---------+---------+=====> P
//|
  |<---L--->|<---L--->|

  Exact solution of displacement u0 is u0 = (P*x)/(E*A).

  $Id$
"""
import getfem as gf
import numpy as np

#
# Physical parameters
#
E = 21E6               # Young Modulus (N/cm^2)
A = 9000.0             # Section Area (cm^2)
L = 2000.0             # Length of element (cm)
P = 200.0              # Force (N)

#
# Numerical parameters
#

elements_degree = 1       # Degree of the finite element methods

#
# Mesh generation.
#
mesh = gf.Mesh("cartesian", np.array([0, L, 2*L]))

#
# Boundary selection
#

fleft = mesh.outer_faces_with_direction(v=[-1.0], angle=0.01)
fright = mesh.outer_faces_with_direction(v=[+1.0], angle=0.01)
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
mesh.set_region(NEUMANN_BOUNDARY, fright)
mesh.set_region(DIRICHLET_BOUNDARY, fleft)
#
# Definition of finite elements methods and integration method
#

mfu = gf.MeshFem(mesh, 1)  # Finite element for the elastic displacement
mfu.set_classical_fem(elements_degree)
mim = gf.MeshIm(mesh, elements_degree*2)   # Integration method

#
# Model definition
#

md = gf.Model("real")
md.add_fem_variable("u", mfu)       # Displacement of the structure

# Truss problem
md.add_initialized_data("D", [E*A])
md.add_generic_elliptic_brick(mim, "u", "D")

F = mfu.eval(str(P))
md.add_initialized_fem_data("ForceData", mfu, F)
md.add_source_term_brick(mim, "u", "ForceData", NEUMANN_BOUNDARY)

md.add_Dirichlet_condition_with_multipliers(mim, "u", elements_degree - 1, DIRICHLET_BOUNDARY)

#
# Model solve
#
md.solve()


#
# Solution export
#
U = md.variable("u")

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfu.eval("(" + str(P) + "*x)/(" + str(E) +"*" + str(A) + ")")
assert np.allclose(U, Ue)

