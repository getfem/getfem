#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2026 Yves Renard.
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
"""  Demo of getfem_plot : several plots 2D and 3D, mech plots and cut

  $Id$
"""
import os
import sys
import importlib.util

import numpy as np

import getfem as gf

import getfem.getfem_plot as gfp

try:
  import getfem.getfem_plot as gfp
except ImportError:
  import getfem_plot as gfp

# ---------------------------------------------------------------- 2D Laplacian
m = gf.Mesh("regular simplices", np.linspace(0, 2, 21), np.linspace(0, 1, 11))
m.set_region(1, m.outer_faces())
mf = gf.MeshFem(m, 1); mf.set_classical_fem(2)
mim = gf.MeshIm(m, 4)
md = gf.Model("real")
md.add_fem_variable("u", mf)
md.add_Laplacian_brick(mim, "u")
md.add_initialized_data("f", [1.0])
md.add_source_term_brick(mim, "u", "f")
md.add_Dirichlet_condition_with_multipliers(mim, "u", mf, 1)
md.solve()
U = md.variable("u")

gfp.plot_mesh(m, regions=[1], title="Mesh and Dirichlet boundary", show=True)
gfp.plot(mf, U, name="u", refine=4, title="Laplacian, P2", show=True)            # matplotlib
gfp.plot(mf, U, name="u", refine=4, contours=15, title="Iso-levels", show=True)

# ------------------------------------------------------- 3D linear elasticity
m3 = gf.Mesh("regular simplices", np.linspace(0, 4, 17), np.linspace(0, 1, 5),
             np.linspace(0, 1, 5))
m3.set_region(1, m3.outer_faces_with_direction([-1.0, 0.0, 0.0], 0.01))
mfu = gf.MeshFem(m3, 3); mfu.set_classical_fem(2)
mim3 = gf.MeshIm(m3, 4)
md = gf.Model("real")
md.add_fem_variable("u", mfu)
md.add_initialized_data("lambda", [1.0]); md.add_initialized_data("mu", [1.0])
md.add_isotropic_linearized_elasticity_brick(mim3, "u", "lambda", "mu")
md.add_initialized_data("f", [0.0, 0.0, -0.01])
md.add_source_term_brick(mim3, "u", "f")
md.add_Dirichlet_condition_with_multipliers(mim3, "u", mfu, 1)
md.solve()
Uel = md.variable("u")

# Surface only, coloured by |u|, on the deformed configuration, with axes and colorbar.
# The PyVista window can be rotated / zoomed / panned with the mouse.
gfp.plot(mfu, Uel, name="u", deformation=Uel, deformation_scale="20%", refine=2,
        title="Displacement norm (deformed)")

# Cut by the plane y = 0.5: skin of the remaining half + cut face, with arrows.
gfp.plot(mfu, Uel, name="u", slice_op=("boundary", ("planar", +1, [0, 0.5, 0], [0, 1, 0])),
        component=2, quiver=True, quiver_count=300, title="u_z, cut y = 0.5")

# Off-screen image (CI, papers)
gfp.plot(mfu, Uel, name="u", deformation=Uel, off_screen=True, screenshot="elasticity.png")

# ------------------------------------------------- mesh annotations (gf_plot_mesh)
mf1 = gf.MeshFem(m, 1); mf1.set_classical_fem(1)
gfp.plot_mesh(m, vertices=True, convexes=True, regions=[1], title="Vertices and convexes", show=True)
gfp.plot_mesh(mf1, dof=True, title="Dof numbers (P1)", show=True)
gfp.plot_mesh(m, quality=True, title="Mesh quality", show=True)
gfp.plot_mesh(m3, convexes=True, vertices=True, title="3D: numbers on the skin")     # pyvista
gfp.plot_mesh(mfu, deformation=Uel, deformation_scale="20%", vertices=True,
              title="Deformed mesh with vertex numbers")
