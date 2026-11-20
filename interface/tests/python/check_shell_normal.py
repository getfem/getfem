#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2026-2026 GetFEM contributors.
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
"""Normal and element_curvature inside the elements of a curve in 2D and
of a surface in 3D, the raytracing and projection transformations on such
elements (shells) and the sphere and circle meshes of
'structured_ball_shell'."""
import numpy as np

import getfem as gf

def check(name, val, ref, tol):
  err = abs(val - ref) / max(abs(ref), 1.)
  print("%-45s %12.6g (expected %12.6g)" % (name, val, ref))
  assert err < tol, "%s: %g instead of %g" % (name, val, ref)

R = 2.
for N in [2, 3]:
  m = gf.Mesh('import', 'structured_ball_shell',
              "GT='GT_QK(%i,2)';ORG=[%s];SIZES=[%g];NSUBDIV=[8]"
              % (N-1, ",".join(["0"]*N), R))
  mim = gf.MeshIm(m, 6)
  md = gf.Model('real')
  meas = 2*np.pi*R if N == 2 else 4*np.pi*R*R
  area = gf.asm_generic(mim, 0, "1", -1, md)
  check("%iD: measure" % N, area, meas, 1e-4)
  # the parametrization gives the outward normal
  check("%iD: int Normal.X/R" % N,
        gf.asm_generic(mim, 0, "Normal.X/%g" % R, -1, md), meas, 1e-4)
  check("%iD: int |Normal|" % N,
        gf.asm_generic(mim, 0, "Norm(Normal)", -1, md), area, 1e-12)
  check("%iD: int |Normal-X/R|" % N, gf.asm_generic(
        mim, 0, "Norm(Normal-X/%g)" % R, -1, md)/area, 0., 1e-3)
  check("%iD: int Normal" % N, np.linalg.norm(
        gf.asm_generic(mim, 0, "Normal", -1, md)), 0., 1e-8)
  # surface gradient of the normal: (I - n n^T)/R
  check("%iD: int Trace(element_curvature)" % N,
        gf.asm_generic(mim, 0, "Trace(element_curvature)", -1, md),
        (N-1)*meas/R, 2e-3)
  check("%iD: int |element_curvature - (I-X@X/R^2)/R|" % N,
        gf.asm_generic(mim, 0, "Norm(element_curvature-(Id(%i)-X@X/%g)/%g)"
                       % (N, R*R, R), -1, md)/area, 0., 1e-2)
  check("%iD: int element_curvature*Normal" % N, np.linalg.norm(
        gf.asm_generic(mim, 0, "element_curvature*Normal", -1, md)), 0., 1e-8)
  # no curvature on elements with a linear transformation (flat)
  ml = gf.Mesh('import', 'structured_ball_shell',
               "GT='GT_PK(%i,1)';ORG=[%s];SIZES=[%g];NSUBDIV=[4]"
               % (N-1, ",".join(["0"]*N), R))
  check("%iD: element_curvature on linear elements" % N,
        gf.asm_generic(gf.MeshIm(ml, 2), 0, "Norm(element_curvature)", -1,
                       md), 0., 1e-12)

# Raytracing and projection in 3D: from a sphere surface of radius 1 (slave shell)
# towards the inner face of a solid ball shell of inner radius 1.5 (master),
# and from a solid ball of radius 1 (slave faces) towards a sphere surface
# of radius 1.5 (master shell, two-sided); distances 0.5.
R1, R2 = 1., 1.5
SLAVE, MASTER = 1, 2
shell1 = gf.Mesh('import', 'structured_ball_shell',
                 "GT='GT_QK(2,2)';ORG=[0,0,0];SIZES=[%g];NSUBDIV=[4]" % R1)
shell1.set_region(SLAVE, np.array([shell1.cvid()]))
solid2 = gf.Mesh('import', 'structured_ball_shell',
                 "GT='GT_QK(3,2)';ORG=[0,0,0];SIZES=[%g,0.5];NSUBDIV=[4,1]"
                 % (R2 + 0.5))
solid2.set_region(MASTER, solid2.outer_faces_in_ball([0., 0., 0.], R2 + 0.25))
solid1 = gf.Mesh('import', 'structured_ball',
                 "GT='GT_QK(3,2)';ORG=[0,0,0];SIZES=[%g];NSUBDIV=[2,3]" % R1)
solid1.set_region(SLAVE, solid1.outer_faces())
shell2 = gf.Mesh('import', 'structured_ball_shell',
                 "GT='GT_QK(2,2)';ORG=[0,0,0];SIZES=[%g];NSUBDIV=[4]" % R2)
shell2.set_region(MASTER, np.array([shell2.cvid()]))

for kind in ["raytracing", "projection"]:
  for ms, mm, what in [(shell1, solid2, "shell slave"),
                       (solid1, shell2, "shell master")]:
    md = gf.Model('real')
    mfs = gf.MeshFem(ms, 3); mfs.set_classical_fem(2)
    mfm = gf.MeshFem(mm, 3); mfm.set_classical_fem(2)
    md.add_fem_variable("us", mfs)
    md.add_fem_variable("um", mfm)
    getattr(md, "add_%s_transformation" % kind)("t", 1.)
    getattr(md, "add_slave_contact_boundary_to_%s_transformation"
            % kind)("t", ms, "us", SLAVE)
    getattr(md, "add_master_contact_boundary_to_%s_transformation"
            % kind)("t", mm, "um", MASTER)
    mim = gf.MeshIm(ms, 6)
    meas = 4*np.pi*R1*R1
    # (the projection misses a few points near the edges of the elements)
    tol = 1e-3 if kind == "raytracing" else 5e-3
    check("%s, %s: int 1" % (kind, what), gf.asm_generic(
          mim, 0, "Interpolate_filter(t, 1, 1)", SLAVE, md), meas, tol)
    check("%s, %s: int distance" % (kind, what), gf.asm_generic(
          mim, 0, "Interpolate_filter(t, Norm(Interpolate(X,t)-X), 1)",
          SLAVE, md), (R2-R1)*meas, tol)
    # the target normal faces the slave point
    check("%s, %s: int Interpolate(Normal).Normal" % (kind, what),
          gf.asm_generic(mim, 0, "Interpolate_filter(t, "
                         "Interpolate(Normal,t).Normal, 1)", SLAVE, md),
          -meas if what == "shell slave" else meas, tol)
    # derivatives of the transformation: gradient of
    # E = int |y - x|^2/2 against finite differences, with a smooth
    # displacement of both bodies
    for v, mf in [("us", mfs), ("um", mfm)]:
      P = mf.basic_dof_nodes()
      f = [lambda p: np.sin(2*p[1]+p[2]), lambda p: np.cos(p[0]-p[2]),
           lambda p: np.sin(p[0]*p[1])]
      md.set_variable(v, np.array([0.02*f[i % 3](P[:, i])
                                   for i in range(mf.nbdof())]))
    E = ("Interpolate_filter(t, Norm_sqr(Interpolate(X,t)+Interpolate(um,t)"
         "-X-us)/2, 1)")
    X0 = md.from_variables()
    G = gf.asm_generic(mim, 1, E, SLAVE, md)
    dX = np.random.default_rng(0).standard_normal(len(X0))*1e-2
    h = 1e-5
    md.to_variables(X0 + h*dX); Ep = gf.asm_generic(mim, 0, E, SLAVE, md)
    md.to_variables(X0 - h*dX); Em = gf.asm_generic(mim, 0, E, SLAVE, md)
    check("%s, %s: gradient, relative error" % (kind, what),
          abs((Ep - Em)/(2*h) - G @ dX)/abs(G @ dX), 0., 1e-4)
    md.to_variables(X0)
    # tangent terms can be assembled
    K = gf.asm_generic(mim, 2, E, SLAVE, md)
    assert K.size()[0] == md.nbdof()

print("check_shell_normal: ok")
