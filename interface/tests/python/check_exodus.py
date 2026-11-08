#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2026 GetFEM contributors.
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
"""  Test the Exodus II import/export of the python-getfem interface.

  Writes a mesh and a nodal field to an Exodus file, reads them back, and
  checks the geometry, field values and a boundary region are recovered.
  Also checks the incremental 'append' time-series path. Requires GetFEM
  built with --enable-exodus.

  $Id$
"""
import numpy as np
import getfem as gf

print("Exodus python import/export test")


def field(P):
    # P is a 1D coordinate array
    v = 1.0 + 2.0 * P[0] + 0.3 * P[0] ** 2
    if P.size > 1:
        v += 3.0 * P[1] - 0.2 * P[1] ** 2
    if P.size > 2:
        v += 4.0 * P[2]
    return v


x = np.arange(0.0, 1.0 + 1e-9, 1.0 / 4)
m = gf.Mesh("cartesian", x, x, x)
m.set_region(7, m.outer_faces())
mf = gf.MeshFem(m, 1)
mf.set_classical_fem(2)
coords = mf.basic_dof_nodes()
U = np.array([field(coords[:, i]) for i in range(coords.shape[1])])

fname = "check_exodus.exo"
mf.export_to_exodus(fname, U, "u")

m2 = gf.Mesh("import", "exodus", fname)
assert m2.nbcvs() == m.nbcvs(), "number of elements"
assert m2.dim() == m.dim(), "dimension"
assert 7 in list(m2.regions()), "boundary region present"
n_orig = np.asarray(m.region(7)).shape[1]
n_imp = np.asarray(m2.region(7)).shape[1]
assert n_imp == n_orig, "boundary face count (%d vs %d)" % (n_imp, n_orig)

vals = np.asarray(m2.exodus_nodal_data(fname, "u", 0)).ravel()
pts = m2.pts()                                  # (dim, nbpt), node i == point i
assert vals.size == pts.shape[1], "variable length"
err = max(abs(vals[i] - field(pts[:, i])) for i in range(vals.size))
assert err < 1e-10, "field values (err=%g)" % err
print("  roundtrip + region: %d elements, %d/%d boundary faces, value error %g"
      % (m2.nbcvs(), n_imp, n_orig, err))

# incremental transient export: the first call creates the file, each later
# call appends one step (no rewrite); re-read every step and check the values
fname2 = "check_exodus_stream.exo"
mfq = gf.MeshFem(gf.Mesh("cartesian", x, x), 1)
mfq.set_classical_fem(1)
P = mfq.basic_dof_nodes()
times = [0.0, 0.5, 1.0]
for k, t in enumerate(times):
    Uk = np.array([(1.0 + t) * field(P[:, i]) for i in range(P.shape[1])])
    if k == 0:
        mfq.export_to_exodus(fname2, "time", t, Uk, "u")
    else:
        mfq.export_to_exodus(fname2, "append", "time", t, Uk, "u")

m3 = gf.Mesh("import", "exodus", fname2)
pts2 = m3.pts()
err2 = 0.0
for k, t in enumerate(times):
    v = np.asarray(m3.exodus_nodal_data(fname2, "u", k)).ravel()
    expected = np.array([(1.0 + t) * field(pts2[:, i]) for i in range(pts2.shape[1])])
    err2 = max(err2, np.max(np.abs(v - expected)))
assert err2 < 1e-10, "incremental append: %d steps (err=%g)" % (len(times), err2)
print("  incremental append: %d steps, value error %g" % (len(times), err2))

print("all Exodus python tests passed")
