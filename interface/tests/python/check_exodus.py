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
  checks the geometry and the field values are recovered. Requires GetFEM
  built with --enable-exodus.

  $Id$
"""
import numpy as np
import getfem as gf

NX = 4
errors = 0


def check(cond, msg):
    global errors
    if not cond:
        errors += 1
        print("  *** error has been detected: " + msg)


def field(P):
    # P is a 1D coordinate array
    v = 1.0 + 2.0 * P[0] + 0.3 * P[0] ** 2
    if P.size > 1:
        v += 3.0 * P[1] - 0.2 * P[1] ** 2
    if P.size > 2:
        v += 4.0 * P[2]
    return v


def make_mesh(kind):
    x = np.arange(0.0, 1.0 + 1e-9, 1.0 / NX)
    if kind == "tri":
        return gf.Mesh("regular_simplices", x, x)
    elif kind == "quad":
        return gf.Mesh("cartesian", x, x)
    elif kind == "tet":
        return gf.Mesh("regular_simplices", x, x, x)
    elif kind == "hex":
        return gf.Mesh("cartesian", x, x, x)
    raise ValueError(kind)


def roundtrip(kind, order):
    tag = "%s (order %d)" % (kind, order)
    fname = "check_exodus_%s_%d.exo" % (kind, order)
    m = make_mesh(kind)
    mf = gf.MeshFem(m, 1)
    mf.set_classical_fem(order)

    coords = mf.basic_dof_nodes()              # (dim, nbdof)
    U = np.array([field(coords[:, i]) for i in range(coords.shape[1])])

    mf.export_to_exodus(fname, U, "u")

    # read the mesh back through getfem
    m2 = gf.Mesh("import", "exodus", fname)
    check(m2.nbcvs() == m.nbcvs(), tag + ": number of elements")
    check(m2.dim() == m.dim(), tag + ": dimension")

    # read the field back through getfem and compare at every node
    vals = np.asarray(m2.exodus_nodal_data(fname, "u", 0)).ravel()
    pts = m2.pts()                             # (dim, nbpt), node i == point i
    check(vals.size == pts.shape[1], tag + ": variable length")
    err = 0.0
    for i in range(vals.size):
        err = max(err, abs(vals[i] - field(pts[:, i])))
    check(err < 1e-10, tag + ": field values (err=%g)" % err)
    print("  %s: %d elements, %d nodes, value error %g"
          % (tag, m2.nbcvs(), vals.size, err))


def region_roundtrip(kind):
    # a boundary face region must survive export (as an Exodus side set) and
    # import (back to a getfem face region)
    fname = "check_exodus_reg_%s.exo" % kind
    m = make_mesh(kind)
    m.set_region(7, m.outer_faces())
    n_orig = np.asarray(m.region(7)).shape[1]
    m.export_to_exodus(fname)
    m2 = gf.Mesh("import", "exodus", fname)
    present = 7 in list(m2.regions())
    check(present, "%s region: present" % kind)
    n_imp = np.asarray(m2.region(7)).shape[1] if present else 0
    check(n_imp == n_orig, "%s region: boundary face count (%d vs %d)"
          % (kind, n_imp, n_orig))
    print("  %s region: %d/%d boundary faces round-tripped"
          % (kind, n_imp, n_orig))


def volume_region_roundtrip(kind):
    # a convex (volume) region must survive export (Exodus element set) and
    # import (back to a getfem convex region)
    fname = "check_exodus_vol_%s.exo" % kind
    m = make_mesh(kind)
    m.set_region(5, m.cvid())          # all convexes -> region 5
    n_orig = m.nbcvs()
    m.export_to_exodus(fname)
    m2 = gf.Mesh("import", "exodus", fname)
    present = 5 in list(m2.regions())
    check(present, "%s volume region: present" % kind)
    n_imp = np.asarray(m2.region(5)).shape[1] if present else 0
    check(n_imp == n_orig, "%s volume region: convex count (%d vs %d)"
          % (kind, n_imp, n_orig))
    print("  %s volume region: %d/%d convexes round-tripped"
          % (kind, n_imp, n_orig))


def streaming_roundtrip():
    # incremental transient export: the first call creates the file, each later
    # call appends one step (no rewrite); re-read every step and check the values
    fname = "check_exodus_stream.exo"
    m = make_mesh("quad")
    mf = gf.MeshFem(m, 1); mf.set_classical_fem(1)
    P = mf.basic_dof_nodes()
    times = [0.0, 0.5, 1.0, 1.5]
    for k, t in enumerate(times):
        U = np.array([(1.0 + t) * field(P[:, i]) for i in range(P.shape[1])])
        if k == 0:
            mf.export_to_exodus(fname, "time", t, U, "u")
        else:
            mf.export_to_exodus(fname, "append", "time", t, U, "u")

    m2 = gf.Mesh("import", "exodus", fname)
    pts = m2.pts()
    err = 0.0
    for k, t in enumerate(times):
        vals = np.asarray(m2.exodus_nodal_data(fname, "u", k)).ravel()
        expv = np.array([(1.0 + t) * field(pts[:, i]) for i in range(pts.shape[1])])
        err = max(err, np.max(np.abs(vals - expv)))
    check(err < 1e-10, "incremental append: %d steps (err=%g)" % (len(times), err))
    print("  incremental append (export_to_exodus 'append'): %d steps, value error %g"
          % (len(times), err))


def independent_reader_check():
    # confirm the written file is genuine Exodus, read by an independent
    # NetCDF reader (skipped if scipy is unavailable)
    try:
        from scipy.io import netcdf_file
    except Exception:
        print("  (scipy not available, skipping independent-reader check)")
        return
    m = make_mesh("quad")
    m.set_region(3, m.outer_faces())   # boundary -> side set + node set
    mf = gf.MeshFem(m, 1)
    mf.set_classical_fem(2)
    coords = mf.basic_dof_nodes()
    U = np.array([field(coords[:, i]) for i in range(coords.shape[1])])
    # 'uncompressed' so the classic-only scipy reader can open it even when the
    # default output is compressed NetCDF4.
    mf.export_to_exodus("check_exodus_indep.exo", "uncompressed", "region field",
                        U, "u")
    f = netcdf_file("check_exodus_indep.exo", "r", mmap=False)
    et = f.variables["connect1"].elem_type
    et = et.decode() if hasattr(et, "decode") else et
    check(et == "QUAD9", "independent reader: elem_type (%s)" % et)
    check(f.dimensions["num_dim"] == 2, "independent reader: num_dim")
    check("coordx" in f.variables and "vals_nod_var1" in f.variables,
          "independent reader: required Exodus variables present")
    check("elem_ss1" in f.variables and "side_ss1" in f.variables,
          "independent reader: side set present")
    check("node_ns1" in f.variables, "independent reader: node set present")
    # identity id maps (help ParaView's GlobalElementId / GlobalNodeId)
    check("elem_num_map" in f.variables and "node_num_map" in f.variables,
          "independent reader: id maps present")
    ne = f.dimensions["num_elem"]
    check(np.array_equal(np.asarray(f.variables["elem_num_map"][:]),
                         np.arange(1, ne + 1)),
          "independent reader: elem_num_map identity")
    # the "region" element (cell) variable carrying the block id
    check(f.dimensions["num_elem_var"] == 1,
          "independent reader: one element variable")
    rname = np.asarray(f.variables["name_elem_var"][0]).tobytes().split(b"\x00", 1)[0]
    check(rname.decode("ascii", "ignore") == "region",
          "independent reader: element var named 'region'")
    reg = np.asarray(f.variables["vals_elem_var1eb1"][0])
    check(np.all(reg == f.variables["eb_prop1"][0]),
          "independent reader: region value == block id")
    f.close()
    print("  independent reader: valid Exodus II with side set, node set, "
          "id maps and region cell var")


def _exo_name(f, var, i):
    return np.asarray(f.variables[var][i]).tobytes().split(b"\x00", 1)[0] \
             .decode("ascii", "ignore")


def region_names_check():
    # the 'region names' option names the Exodus blocks and element sets
    try:
        from scipy.io import netcdf_file
    except Exception:
        return
    m = make_mesh("hex")
    cv = np.asarray(m.cvid()).ravel()
    m.set_region(1, cv[0::2])               # two disjoint volume regions
    m.set_region(2, cv[1::2])
    mf = gf.MeshFem(m, 1); mf.set_classical_fem(1)
    P = mf.basic_dof_nodes()
    U = np.array([field(P[:, i]) for i in range(P.shape[1])])
    mf.export_to_exodus("check_exodus_named.exo", "uncompressed",
                        "region names", [1, 2], ["inner", "outer"], U, "u")
    f = netcdf_file("check_exodus_named.exo", "r", mmap=False)
    ebp = [int(x) for x in f.variables["eb_prop1"][:]]
    want = {1: "inner", 2: "outer"}
    ok = all(_exo_name(f, "eb_names", i) == want[ebp[i]] for i in range(len(ebp)))
    check(ok, "region names: blocks named inner/outer")
    check("els_names" in f.variables, "region names: els_names present")
    f.close()
    print("  region names: blocks/sets named via the 'region names' option")


def streaming_region_check():
    # the region cell var must be written at every (appended) time step
    try:
        from scipy.io import netcdf_file
    except Exception:
        return
    fname = "check_exodus_stream_region.exo"
    m = make_mesh("quad")
    m.set_region(1, m.cvid())               # all convexes -> one volume region
    mf = gf.MeshFem(m, 1); mf.set_classical_fem(1)
    P = mf.basic_dof_nodes()
    times = [0.0, 0.5, 1.0]
    for k, t in enumerate(times):
        U = np.array([(1.0 + t) * field(P[:, i]) for i in range(P.shape[1])])
        if k == 0:
            mf.export_to_exodus(fname, "uncompressed", "region field",
                                "time", t, U, "u")
        else:
            mf.export_to_exodus(fname, "append", "time", t, U, "u")
    f = netcdf_file(fname, "r", mmap=False)
    reg = np.asarray(f.variables["vals_elem_var1eb1"][:])   # (nsteps, nelem)
    bid = int(f.variables["eb_prop1"][0])
    check(reg.shape[0] == len(times) and np.all(reg == bid),
          "streaming: region cell var at every step")
    f.close()
    print("  streaming: region cell var written at each appended step")


def compression_check():
    # Default output is compressed NetCDF4 only when the build proved usable
    # NetCDF4/HDF5 deflate support; otherwise it is classic 64-bit offset.
    import os
    x = np.arange(0.0, 1.0 + 1e-9, 1.0 / 12)
    m = gf.Mesh("cartesian", x, x, x)
    mf = gf.MeshFem(m, 1); mf.set_classical_fem(1)
    P = mf.basic_dof_nodes()
    U = np.array([field(P[:, i]) for i in range(P.shape[1])])
    mf.export_to_exodus("check_exodus_comp.exo", U, "u")                 # default
    mf.export_to_exodus("check_exodus_unc.exo", "uncompressed", U, "u")  # opt out
    sc = os.path.getsize("check_exodus_comp.exo")
    su = os.path.getsize("check_exodus_unc.exo")
    with open("check_exodus_comp.exo", "rb") as f:
        default_is_hdf5 = f.read(8) == b"\x89HDF\r\n\x1a\n"
    if default_is_hdf5:
        check(sc < su, "compression: default file smaller than uncompressed "
              "(%d vs %d)" % (sc, su))
    else:
        check(sc > 0 and su > 0, "compression: default classic file written")
    # default output must round-trip through GetFEM in either format
    m2 = gf.Mesh("import", "exodus", "check_exodus_comp.exo")
    pts = m2.pts()
    vals = np.asarray(m2.exodus_nodal_data("check_exodus_comp.exo", "u", 0)).ravel()
    err = max(abs(vals[i] - field(pts[:, i])) for i in range(vals.size))
    check(err < 1e-10, "compression: compressed file round-trips (err=%g)" % err)
    try:
        from scipy.io import netcdf_file
        netcdf_file("check_exodus_unc.exo", "r", mmap=False).close()
        ok_unc = True
    except Exception:
        ok_unc = False
    check(ok_unc, "compression: scipy reads the uncompressed file")
    mode = "NetCDF4 compressed" if default_is_hdf5 else "classic 64-bit offset"
    print("  compression: default %s %d B vs uncompressed %d B, round-trip err %g"
          % (mode, sc, su, err))


print("Exodus python import/export test")
for kind in ("tri", "quad", "tet", "hex"):
    for order in (1, 2):
        roundtrip(kind, order)
for kind in ("tri", "quad", "tet", "hex"):
    region_roundtrip(kind)
for kind in ("tri", "quad", "tet", "hex"):
    volume_region_roundtrip(kind)
streaming_roundtrip()
independent_reader_check()
region_names_check()
streaming_region_check()
compression_check()

if errors:
    print("FAILED (%d errors)" % errors)
else:
    print("all Exodus python tests passed")
assert errors == 0
