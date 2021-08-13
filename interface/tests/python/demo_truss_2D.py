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
""" 2D Truss problem test

  This program is used to check that python-getfem is working. This is
  also a good example of use of python-getfem.

  $Id$
"""
import getfem as gf
import numpy as np

#
# Physical parameters
#
E = 21e6  # Young Modulus (N/cm^2)
A = 9000.0  # Section Area (cm^2)
L = 2000.0  # Length of element (cm)
P = 200.0  # Force (N)
Q = 300.0  # Force (N)

#
# Numerical parameters
#
theta = np.pi / 6.0
# theta = 0.0
elements_degree = 1  # Degree of the finite element methods

x1 = 0.0
y1 = L / 2.0
x2 = np.sqrt(3) * L / 2.0
y2 = L / 2.0
x3 = 0.0
y3 = 0.0

#
# Mesh generation.
#
mesh = gf.Mesh("empty", 2)
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x1, x2], [y1, y2]])
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x2, x3], [y2, y3]])

#
# Boundary selection
#
pid1 = mesh.pid_from_coords([x1, y1])
pid2 = mesh.pid_from_coords([x2, y2])
pid3 = mesh.pid_from_coords([x3, y3])
f1 = mesh.faces_from_pid(pid1)
f2 = mesh.faces_from_pid(pid2)
f3 = mesh.faces_from_pid(pid3)
BOUNDARY1 = 1
BOUNDARY2 = 2
BOUNDARY3 = 3
mesh.set_region(BOUNDARY1, f1)
mesh.set_region(BOUNDARY2, f2)
mesh.set_region(BOUNDARY3, f3)

#
# Definition of finite elements methods and integration method
#
mfu = gf.MeshFem(mesh, 2)  # Finite element for the elastic displacement
mfu.set_classical_fem(elements_degree)
mim = gf.MeshIm(mesh, elements_degree * 2)  # Integration method

#
# Model definition
#
md = gf.Model("real")
md.add_fem_variable("u", mfu)  # Displacement of the structure

# Truss problem
md.add_initialized_data("E", E)
md.add_initialized_data("A", A)
md.add_linear_term(
    mim,
    "((E*A)*[[+element_K(1)/Norm(element_K)*element_K(1)/Norm(element_K), -element_K(1)/Norm(element_K)*element_K(2)/Norm(element_K)],"
    + "[-element_K(2)/Norm(element_K)*element_K(1)/Norm(element_K), +element_K(2)/Norm(element_K)*element_K(2)/Norm(element_K)]]*Grad_u):Grad_Test_u",
)

md.add_explicit_rhs("u", np.array([0.0, 0.0, P, Q, 0.0, 0.0]))

md.add_Dirichlet_condition_with_simplification("u", BOUNDARY1)
md.add_Dirichlet_condition_with_simplification("u", BOUNDARY3)

#
# Model solve
#
md.solve()

#
# Solution export
#
U = md.variable("u")

#
# Exact solution
#
SM1 = E*A/L*np.array([[(8.0*np.sqrt(3.0)+9.0)/12.0, -np.sqrt(3)/4.0], [-np.sqrt(3)/4.0, 1.0/4.0]])
rhs1 = np.array([P, Q])
u2 = np.linalg.solve(SM1, rhs1)
Ue = np.array([0.0, 0.0, u2[0], u2[1], 0.0, 0.0])

assert np.allclose(U, Ue)
