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
"""  Check of the visualisation module getfem_plot.

  The extraction of fields on slices is checked against exact values (P2
  fields reproducing quadratic functions), then all the plotting functions
  are run with the non interactive matplotlib backend (no display needed).
  PyVista, if installed, is only checked for the conversion of the data
  (no OpenGL rendering). The test is skipped if matplotlib is not installed.

  $Id$
"""
import os
import sys
import importlib.util

import numpy as np

import getfem as gf

try:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
except ImportError:
  print("matplotlib not installed: test skipped")
  sys.exit(77)

try:
  import getfem.getfem_plot as gfp
except ImportError:
  import getfem_plot as gfp

tol = 1e-9

def check(cond, msg):
  if not cond:
    print("FAILED:", msg)
    sys.exit(1)

def field_on_dofs(mf, funcs):
  """Values of the functions funcs (one per component) on the basic dofs."""
  P = np.asarray(mf.basic_dof_nodes())
  q = len(funcs)
  U = np.zeros(P.shape[1])
  for i in range(P.shape[1]):
    U[i] = funcs[i % q](*P[:, i])
  return U

def exact(points, funcs, sdim):
  X = [points[:, k] for k in range(sdim)]
  V = np.stack([f(*X) for f in funcs], axis=1)
  return V[:, 0] if len(funcs) == 1 else V


# ------------------------------------------------------------------ 2D mesh
m = gf.Mesh('regular simplices', np.linspace(0, 2, 9), np.linspace(0, 1, 5))
m.set_region(1, m.outer_faces())

mf = gf.MeshFem(m, 1); mf.set_classical_fem(2)
fs = [lambda x, y: x*x + x*y - y]
U = field_on_dofs(mf, fs)
d = gfp.slice_field(mf, U, refine=3, name="u", with_edges=True)
check(d.sdim == 2 and d.tdim == 2, "2D slice dimensions")
check(d.values.shape == (len(d.points),), "2D scalar field shape")
check(np.abs(d.values - exact(d.points, fs, 2)).max() < tol,
      "2D scalar P2 field not reproduced exactly")
check(d.edges is not None and len(d.edges[1]) > 0, "2D mesh edges")

mfv = gf.MeshFem(m, 2); mfv.set_classical_fem(2)
fv = [lambda x, y: x*y, lambda x, y: x - y*y]
V = field_on_dofs(mfv, fv)
d = gfp.slice_field(mfv, V, refine=2, name="v")
check(d.values.shape == (len(d.points), 2), "2D vector field shape")
check(np.abs(d.values - exact(d.points, fv, 2)).max() < tol,
      "2D vector P2 field not reproduced exactly (component ordering)")
check(np.abs(d.scalars() - np.linalg.norm(d.values, axis=1)).max() < tol,
      "norm of a vector field")

# deformation by the linear displacement (a x, b y) with an explicit scale
mfd = gf.MeshFem(m, 2); mfd.set_classical_fem(1)
D = field_on_dofs(mfd, [lambda x, y: 0.5*x, lambda x, y: -0.25*y])
d = gfp.slice_field(mf, U, deformation=D, deformation_mf=mfd,
                    deformation_scale=1.0)
check(abs(d.points[:, 0].max() - 3.0) < tol and
      abs(d.points[:, 1].max() - 0.75) < tol, "2D deformation")
check(np.abs(d.values - exact(d.points / [1.5, 0.75, 1.], fs, 2)).max() < tol,
      "field on the deformed configuration")
d = gfp.slice_field(mf, U, deformation=D, deformation_mf=mfd,
                    deformation_scale="10%")
diag = np.sqrt(2.0**2 + 1.0**2)
check(abs(d.meta["deformation_scale"] - 0.1*diag/np.abs(D).max()) < tol,
      "relative deformation scale")

# mesh annotations
X, lab = gfp._vertex_labels(m, None, False, None)
check(len(lab) == m.nbpts(), "vertex numbers")
md = gfp.slice_mesh(m)
X, lab = gfp._convex_labels(m, None, False, md)
check(len(lab) == m.nbcvs(), "convex numbers")
X, lab = gfp._dof_labels(mfv, None, False, None)
check(len(lab) == mfv.nb_basic_dof() // 2 and all("," in l for l in lab),
      "dof numbers grouped by node")
q = gfp._quality_values(m, md)
check(q.shape == (len(md.points),) and q.min() > 0. and q.max() <= 1. + tol,
      "mesh quality on the slice")

# rendering (matplotlib, no display)
fig, ax = gfp.plot_mesh(m, vertices=True, convexes=True, regions=[1],
                        title="mesh")
check(len(ax.get_legend().get_texts()) == 1, "region legend")
plt.close(fig)
fig, ax = gfp.plot_mesh(mfv, dof=True, quality=True)
plt.close(fig)
fig, ax = gfp.plot(mf, U, name="u", contours=10, title="u",
                   axis_labels=("x [m]", "y [m]", "z [m]"))
check(ax.get_xlabel() == "x [m]", "axis labels")
check(any(a.get_ylabel() == "u" for a in fig.axes), "colorbar label")
plt.close(fig)
fig, ax = gfp.plot(mfv, V, name="v", quiver=True, deformation=V,
                   deformation_mf=mfv)
check(any(a.get_ylabel() == "norm(v)" for a in fig.axes), "norm label")
plt.close(fig)
fig, ax = gfp.plot(mfv, V, name="v", component=1)
check(any(a.get_ylabel() == "v_y" for a in fig.axes), "component label")
plt.close(fig)

# several slices in the same axes: field, arrows only, boundary curves
check("matplotlib" in gfp.available_backends(), "available backends")
sl = gf.Slice(('none',), m, 2)
fig, ax = gfp.plot_slice(sl, mf, U, name="u")
n0 = len(ax.collections)
gfp.plot_slice(sl, mfv, V, name="v", faces=False, quiver=True, colorbar=False,
               ax=ax)
sb = gf.Slice(('none',), m, 1, m.region(1))
db = gfp.slice_to_data(sb)
check(db.tdim == 1 and len(db.simplexes[1]) > 0, "slice of boundary faces")
gfp.plot_slice(sb, line_color="red", ax=ax, colorbar=False)
check(len(ax.collections) == n0 + 2, "layers added to the same axes")
check(len(fig.axes) == 2, "only one colour bar")
plt.close(fig)


# ------------------------------------------------------------------ 3D mesh
m3 = gf.Mesh('regular simplices', np.linspace(0, 1, 4), np.linspace(0, 1, 4),
             np.linspace(0, 2, 5))
mf3 = gf.MeshFem(m3, 3); mf3.set_classical_fem(2)
f3 = [lambda x, y, z: z*x, lambda x, y, z: x - y, lambda x, y, z: y*y]
U3 = field_on_dofs(mf3, f3)
d3 = gfp.slice_field(mf3, U3, refine=2, name="u", with_edges=True)
check(d3.sdim == 3 and d3.tdim == 2, "3D: only the surface is extracted")
P = d3.points
on_skin = ((np.abs(P[:, 0]) < tol) | (np.abs(P[:, 0] - 1) < tol) |
           (np.abs(P[:, 1]) < tol) | (np.abs(P[:, 1] - 1) < tol) |
           (np.abs(P[:, 2]) < tol) | (np.abs(P[:, 2] - 2) < tol))
check(on_skin.all(), "3D surface points lie on the boundary")
check(np.abs(d3.values - exact(P, f3, 3)).max() < tol,
      "3D vector P2 field not reproduced exactly")

# cut by the plane z = 1: skin of one half and the cut face
dc = gfp.slice_field(mf3, U3, slice_op=('boundary', ('planar', +1, [0, 0, 1],
                                                    [0, 0, 1])), refine=1)
z = dc.points[:, 2]
check((z >= 1 - tol).all() or (z <= 1 + tol).all(), "planar cut keeps one half")
check((np.abs(z - 1) < tol).any(), "planar cut face")
check(np.abs(dc.values - exact(dc.points, f3, 3)).max() < tol, "field on a cut")

# annotations on the skin only
X, lab = gfp._vertex_labels(m3, None, True, None)
nskin = len(np.unique(np.asarray(m3.pid_in_faces(m3.outer_faces()))))
check(len(lab) == nskin, "3D vertex numbers on the skin")
md3 = gfp.slice_mesh(m3)
X, lab = gfp._convex_labels(m3, None, True, md3)
check(len(lab) == len(np.unique(np.asarray(m3.outer_faces())[0])),
      "3D convex numbers on the skin")

fig, ax = gfp.plot(mf3, U3, name="u", refine=1, backend="matplotlib")
plt.close(fig)
fig, ax = gfp.plot_mesh(m3, convexes=True, backend="matplotlib")
plt.close(fig)

# key saving the view of an interactive window (without PyVista: fake plotter)
import tempfile
class FakePlotter:
  iren = object()
  def __init__(self): self.keys = {}
  def add_key_event(self, key, f): self.keys[key] = f
  def screenshot(self, f, scale=1): open(f, "w").write("bitmap")
  def save_graphic(self, f): open(f, "w").write("vector")
with tempfile.TemporaryDirectory() as tmp:
  fp = FakePlotter()
  pattern = os.path.join(tmp, "view_{:02d}.png")
  open(pattern.format(1), "w").close()
  gfp._install_save_key(fp, "i", pattern)
  fp.keys["i"](); fp.keys["i"]()
  check(sorted(os.listdir(tmp)) == ["view_01.png", "view_02.png", "view_03.png"],
        "numbered saving of the view")
  check(open(pattern.format(1)).read() == "", "existing file not overwritten")
  gfp.save_view(fp, os.path.join(tmp, "view.pdf"))
  check(open(os.path.join(tmp, "view.pdf")).read() == "vector", "vector export")
fp = FakePlotter(); fp.iren = None
gfp._install_save_key(fp, "i", "x.png")
check(fp.keys == {}, "no key binding off-screen")

# PyVista: conversion only, no rendering
try:
  import pyvista
except Exception:
  pyvista = None
if pyvista is not None:
  grid = gfp.to_pyvista(d3)
  check(grid.n_points == len(d3.points), "pyvista points")
  check(grid.n_cells == len(d3.simplexes[2]), "pyvista cells")
  check("norm(u)" in grid.point_data and "u" in grid.point_data,
        "pyvista point data")
  check(grid.point_data["u"].shape == (len(d3.points), 3), "pyvista vectors")

print("check_plot: all tests passed")
