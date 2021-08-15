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
""" 3D Truss problem test

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

#
# Numerical parameters
#
elements_degree = 1  # Degree of the finite element methods
b = 1.0
h = 2.0
x1 = 0.0
y1 = 0.0
z1 = 0.0
x2 = b
y2 = 0.0
z2 = 0.0
x3 = 0.0
y3 = b
z3 = 0.0
x4 = b
y4 = b
z4 = 0.0
x5 = 0.0
y5 = 0.0
z5 = h

#
# Mesh generation.
#
mesh = gf.Mesh("empty", 3)
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x1, x5], [y1, y5], [z1, z5]])
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x2, x5], [y2, y5], [z2, z5]])
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x3, x5], [y3, y5], [z3, z5]])
mesh.add_convex(gf.GeoTrans("GT_PK(1,1)"), [[x4, x5], [y4, y5], [z4, z5]])

#
# Boundary selection
#
pid1 = mesh.pid_from_coords([x1, y1, z1])
pid2 = mesh.pid_from_coords([x2, y2, z2])
pid3 = mesh.pid_from_coords([x3, y3, z3])
pid4 = mesh.pid_from_coords([x4, y4, z4])
pid5 = mesh.pid_from_coords([x5, y5, z5])
f1 = mesh.faces_from_pid(pid1)
f2 = mesh.faces_from_pid(pid2)
f3 = mesh.faces_from_pid(pid3)
f4 = mesh.faces_from_pid(pid4)
f5 = mesh.faces_from_pid(pid5)
BOUNDARY1 = 1
BOUNDARY2 = 2
BOUNDARY3 = 3
BOUNDARY4 = 4
BOUNDARY5 = 5
mesh.set_region(BOUNDARY1, f1)
mesh.set_region(BOUNDARY2, f2)
mesh.set_region(BOUNDARY3, f3)
mesh.set_region(BOUNDARY4, f4)
mesh.set_region(BOUNDARY5, f5)

#
# Definition of finite elements methods and integration method
#
mfu = gf.MeshFem(mesh, 3)  # Finite element for the elastic displacement
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
md.add_linear_term(mim, "((E*A)*" +
        "[[element_K(1)/Norm(element_K)*element_K(1)/Norm(element_K), element_K(1)/Norm(element_K)*element_K(2)/Norm(element_K), element_K(1)/Norm(element_K)*element_K(3)/Norm(element_K)]," +
        " [element_K(2)/Norm(element_K)*element_K(1)/Norm(element_K), element_K(2)/Norm(element_K)*element_K(2)/Norm(element_K), element_K(2)/Norm(element_K)*element_K(3)/Norm(element_K)]," +
        " [element_K(3)/Norm(element_K)*element_K(1)/Norm(element_K), element_K(3)/Norm(element_K)*element_K(2)/Norm(element_K), element_K(3)/Norm(element_K)*element_K(3)/Norm(element_K)]]" +
        "*Grad_u):Grad_Test_u")

F = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -P, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
md.add_explicit_rhs("u", F)

md.add_Dirichlet_condition_with_simplification("u", BOUNDARY1)
md.add_Dirichlet_condition_with_simplification("u", BOUNDARY2)
md.add_Dirichlet_condition_with_simplification("u", BOUNDARY3)
md.add_Dirichlet_condition_with_simplification("u", BOUNDARY4)

#
# Model solve
#
md.solve()

#
# Solution export
U = md.variable("u")

#
# Exact solution
#
SM = md.tangent_matrix()
print(SM.full())
rhs = md.rhs()
print(U)
