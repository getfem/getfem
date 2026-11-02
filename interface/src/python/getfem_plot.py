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
"""Visualisation of meshes and finite element fields.

Python counterpart of the Matlab/Octave functions gf_plot_mesh, gf_plot and
gf_plot_slice::

    import getfem as gf
    import getfem_plot as gfp

    gfp.plot_mesh(m, vertices=True, regions=[1])
    gfp.plot(mf, U, name="u", deformation=U, deformation_scale="10%")

Principle::

    (mesh_fem, U) --gf.Slice(+refinement)--> (points, simplexes, values) --> renderer

* The extraction layer only uses the Python API of GetFEM (gf.Slice,
  Slice.pts/splxs/edges/cvs, gf.compute_interpolate_on) and returns plain
  NumPy arrays in a ``SliceData`` object (see ``slice_field``, ``slice_mesh``).
* Two renderers, both optional dependencies imported only when used:

  - Matplotlib: 1D and 2D (publication quality figures), small 3D surfaces;
  - PyVista/VTK: 3D, interactive (rotation, zoom, pan), Jupyter (trame),
    off-screen rendering.

* In 3D only the surface is drawn: the default slice operation is
  ('boundary',), which only produces the 2D simplexes of the boundary skin.

Main functions: ``plot_mesh``, ``plot``, ``plot_slice``, ``slice_field``,
``slice_mesh``, ``to_pyvista``.  Indices are 0-based, as everywhere in the
Python interface.
"""

from __future__ import annotations

__all__ = ["plot", "plot_mesh", "plot_slice", "slice_field", "slice_mesh",
           "slice_to_data", "to_pyvista", "available_backends", "save_view",
           "SliceData"]

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

try:
    import getfem as gf
except ImportError:  # allows importing the module for documentation / tests
    gf = None


# ---------------------------------------------------------------------------
#  Extraction layer (no rendering dependency)
# ---------------------------------------------------------------------------

@dataclass
class SliceData:
    """Data extracted from a gf.Slice, independent of any renderer.

    points    : (n, 3) coordinates (padded with zeros when sdim < 3)
    sdim      : dimension of the ambient space of the slice (1, 2 or 3)
    simplexes : {d: (m, d+1) int array} simplexes of dimension d
    values    : None, (n,) scalar values or (n, q) vector/tensor values
    name      : name of the field (used for colorbar / legend)
    edges     : None or (pts (k, 3), segments (s, 2)) mesh edges to overlay
    """
    points: np.ndarray
    sdim: int
    simplexes: dict
    values: Optional[np.ndarray] = None
    name: str = ""
    edges: Optional[tuple] = None
    meta: dict = field(default_factory=dict)

    @property
    def tdim(self):
        """Highest dimension of the simplexes present in the slice."""
        ds = [d for d, s in self.simplexes.items() if len(s) > 0]
        return max(ds) if ds else 0

    @property
    def qdim(self):
        if self.values is None:
            return 0
        return 1 if self.values.ndim == 1 else self.values.shape[1]

    def scalars(self, component=None):
        """Scalar values to colour with: component index, 'norm' or None (auto)."""
        if self.values is None:
            return None
        if self.values.ndim == 1:
            return self.values
        if component is None or component == "norm":
            return np.linalg.norm(self.values, axis=1)
        return self.values[:, int(component)]

    def vectors3(self):
        """Vector values padded to 3 components (for arrows / glyphs)."""
        if self.values is None or self.values.ndim == 1:
            return None
        v = self.values[:, : min(3, self.values.shape[1])]
        out = np.zeros((v.shape[0], 3))
        out[:, : v.shape[1]] = v
        return out


def _pad3(P):
    """GetFEM points (sdim, n) -> (n, 3)."""
    P = np.asarray(P, dtype=float)
    if P.ndim == 1:
        P = P[None, :]
    d, n = P.shape
    X = np.zeros((n, 3))
    X[:, :d] = P.T
    return X, d


def _get_splxs(sl, d, npts):
    """Simplexes of dimension d of a slice, as an (m, d+1) int array."""
    try:
        res = sl.splxs(d)
    except Exception:
        return np.empty((0, d + 1), dtype=int)
    S = res[0] if isinstance(res, (tuple, list)) else res
    S = np.asarray(S)
    if S.size == 0:
        return np.empty((0, d + 1), dtype=int)
    S = S.astype(int).reshape(d + 1, -1).T
    if S.max() >= npts:  # safety net if the indices were 1-based
        S = S - 1
    return S


def _get_edges(sl):
    """Edges of the original mesh contained in the slice (Slice.edges)."""
    try:
        P, E1, E2 = sl.edges()
    except Exception:
        return None
    X, _ = _pad3(P)
    E = np.asarray(E1)
    if E.size == 0:
        return None
    E = E.astype(int).reshape(2, -1).T
    if E.max() >= len(X):
        E = E - 1
    return X, E


def _interpolate(mf, U, sl, npts):
    """Values of the field (mf, U) at the points of the slice -> (n,) or (n, q)."""
    # The interface returns an array of shape (qdim, npts) (or (npts,) for a
    # scalar field), the point index being the LAST one.  Flattening must keep
    # the Fortran ordering used by GetFEM, i.e. components fastest.
    V = np.asarray(gf.compute_interpolate_on(mf, U, sl), dtype=float)
    if V.size % npts != 0 or (V.ndim > 1 and V.shape[-1] != npts):
        raise ValueError("unexpected shape of the interpolated field "
                         f"{V.shape} for {npts} points")
    q = V.size // npts
    if q == 1:
        return V.reshape(npts)
    if V.ndim == 1:  # flat vector, components interleaved point by point
        return V.reshape(npts, q)
    # (qdim, npts) or (N, N, npts) for tensors: component indices first
    return V.reshape(q, npts, order="F").T


def slice_to_data(sl, mf=None, U=None, name="", with_edges=False):
    """Convert an existing gf.Slice (and optionally a field) into SliceData."""
    X, sdim = _pad3(sl.pts())
    n = X.shape[0]
    simplexes = {d: _get_splxs(sl, d, n) for d in range(0, min(sdim, 3) + 1)}
    values = _interpolate(mf, U, sl, n) if mf is not None else None
    edges = _get_edges(sl) if with_edges else None
    return SliceData(X, sdim, simplexes, values, name, edges)


def _default_op(mesh):
    return ("boundary",) if mesh.dim() == 3 else ("none",)


def _mesh_size(mesh):
    P = np.asarray(mesh.pts(), dtype=float)
    return float(np.linalg.norm(P.max(axis=1) - P.min(axis=1)))


def _make_slice(op, mesh, refine, cvlst=None, mf_def=None, D=None):
    """gf.Slice on the mesh, or on the mesh deformed by (mf_def, D)."""
    args = [op]
    if mf_def is not None:
        # the interface checks the *last* dimension of U against nb_dof:
        # a 1D numpy array is seen as a column, so pass it as a (1, nb_dof) row.
        args += [mf_def, np.asarray(D, dtype=float).reshape(1, -1)]
    else:
        args += [mesh]
    args.append(int(refine))
    if cvlst is not None:
        args.append(cvlst)
    return gf.Slice(*args)


def _deformation_scale(scale, D, mesh):
    """Matlab-like scale: a number, or 'p%' (max displacement = p% of mesh size)."""
    if isinstance(scale, str) and scale.strip().endswith("%"):
        pct = float(scale.strip()[:-1]) / 100.0
        dmax = float(np.max(np.abs(D))) if np.size(D) else 0.0
        return 0.0 if dmax == 0.0 else pct * _mesh_size(mesh) / dmax
    return float(scale)


def slice_field(mf, U, slice_op=None, refine=3, cvlst=None,
                deformation=None, deformation_mf=None, deformation_scale="10%",
                name="", with_edges=False):
    """Extract a field (mf, U) on a slice of its mesh.

    slice_op      : GetFEM slice operation, e.g. ('none',), ('boundary',),
                    ('planar', +1, x0, n), ('boundary', ('planar', +1, x0, n)),
                    ('isovalues', 0, mf2, U2, val) ... Default: ('boundary',)
                    in 3D (surface only), ('none',) otherwise.
    refine        : number of subdivisions of each convex (high order elements
                    and curved geometric transformations need refine > 1).
    cvlst         : convex list / region (as returned by mesh.region(r)).
    deformation   : displacement vector D defined on deformation_mf (default mf).
    deformation_scale : number or 'p%' (default '10%' as in Matlab gf_plot).
    """
    mesh = mf.linked_mesh()
    op = slice_op if slice_op is not None else _default_op(mesh)
    if deformation is not None:
        mf_def = deformation_mf if deformation_mf is not None else mf
        D = np.asarray(deformation, dtype=float)
        s = _deformation_scale(deformation_scale, D, mesh)
        sl = _make_slice(op, mesh, refine, cvlst, mf_def, D * s)
    else:
        s = 0.0
        sl = _make_slice(op, mesh, refine, cvlst)
    data = slice_to_data(sl, mf, U, name=name, with_edges=with_edges)
    data.meta["deformation_scale"] = s
    data.meta["slice"] = sl
    return data


def slice_mesh(mesh, slice_op=None, refine=1, cvlst=None,
               deformation=None, deformation_mf=None, deformation_scale=1.0):
    """Extract a mesh (or the surface of a 3D mesh) as SliceData with its edges."""
    op = slice_op if slice_op is not None else _default_op(mesh)
    if deformation is not None:
        if deformation_mf is None:
            raise ValueError("deformation requires deformation_mf (the mesh_fem of D)")
        D = np.asarray(deformation, dtype=float)
        s = _deformation_scale(deformation_scale, D, mesh)
        sl = _make_slice(op, mesh, refine, cvlst, deformation_mf, D * s)
        deform = (deformation_mf, D * s)
    else:
        sl = _make_slice(op, mesh, refine, cvlst)
        deform = None
    data = slice_to_data(sl, with_edges=True)
    data.meta.update(slice=sl, deform=deform, op=op)
    return data


# ---------------------------------------------------------------------------
#  Mesh annotations (vertex, convex and dof numbers, quality)
# ---------------------------------------------------------------------------

def _as_cvfids(cvlst):
    """cvlst -> (convex ids, face array or None)."""
    C = np.asarray(cvlst)
    if C.ndim == 2 and C.shape[0] == 2:
        return np.unique(C[0].astype(int)), C
    return np.unique(C.reshape(-1).astype(int)), None


def _visible_entities(mesh, cvlst, surface):
    """Convex ids and convex faces whose numbers should be displayed."""
    if cvlst is not None:
        return _as_cvfids(cvlst)
    if surface:  # 3D: only what lies on the skin
        F = np.asarray(mesh.outer_faces())
        return np.unique(F[0].astype(int)), F
    return np.asarray(mesh.cvid(), dtype=int).reshape(-1), None


def _displaced(points, deform):
    """Apply the (already scaled) deformation (mf, D) to points (k, 3)."""
    if deform is None or len(points) == 0:
        return points
    mf, D = deform
    sdim = mf.linked_mesh().dim()
    # points given as an (sdim, k) array (one point per column, as mesh.pts())
    V = np.asarray(gf.compute_interpolate_on(mf, D, points[:, :sdim].T), dtype=float)
    V = V.reshape(-1, len(points), order="F").T if V.ndim > 1 else V.reshape(len(points), -1)
    out = points.copy()
    out[:, : V.shape[1]] += V
    return out


def _vertex_labels(mesh, cvlst, surface, deform):
    cv, F = _visible_entities(mesh, cvlst, surface)
    pids = mesh.pid_in_faces(F) if F is not None else mesh.pid_in_cvids(cv)
    pids = np.unique(np.asarray(pids, dtype=int).reshape(-1))
    X, _ = _pad3(np.asarray(mesh.pts(pids)))
    return _displaced(X, deform), [str(p) for p in pids]


def _convex_labels(mesh, cvlst, surface, data):
    """Convex numbers at the convex centres (2D), or on the skin (3D surface)."""
    if surface:
        # place the number on the visible surface: for each convex of the slice,
        # the centroid of its simplex closest to the mean of its slice points.
        sl = data.meta["slice"]
        cvs = np.asarray(sl.cvs(), dtype=int).reshape(-1)
        S, cv2s = sl.splxs(2)
        S = np.asarray(S, dtype=int).reshape(3, -1).T
        cv2s = np.asarray(cv2s, dtype=int).reshape(-1)
        X = data.points
        pos, lab = [], []
        keep = None if cvlst is None else set(_as_cvfids(cvlst)[0].tolist())
        for i, cv in enumerate(cvs):
            Si = S[cv2s[i]:cv2s[i + 1]]
            if len(Si) == 0 or (keep is not None and cv not in keep):
                continue
            G = X[Si].mean(axis=1)                   # centroids of the simplexes
            j = np.argmin(np.linalg.norm(G - G.mean(axis=0), axis=1))
            pos.append(G[j]); lab.append(str(cv))
        return np.array(pos).reshape(-1, 3), lab
    cv, _ = _visible_entities(mesh, cvlst, False)
    pid, idx = mesh.pid_from_cvid(cv)
    pid = np.asarray(pid, dtype=int).reshape(-1)
    idx = np.asarray(idx, dtype=int).reshape(-1)
    P, _ = _pad3(np.asarray(mesh.pts()))
    centers = np.add.reduceat(P[pid], idx[:-1], axis=0) / np.diff(idx)[:, None]
    return _displaced(centers, data.meta.get("deform")), [str(c) for c in cv]


def _dof_labels(mf, cvlst, surface, deform):
    """Basic dof numbers at their nodes; dofs sharing a node are grouped."""
    mesh = mf.linked_mesh()
    if cvlst is None and not surface:
        dofs = np.arange(int(mf.nb_basic_dof()))
    else:
        cv, F = _visible_entities(mesh, cvlst, surface)
        dofs = mf.basic_dof_from_cv(F if F is not None else cv)
        dofs = np.unique(np.asarray(dofs, dtype=int).reshape(-1))
    if len(dofs) == 0:
        return np.empty((0, 3)), []
    X, _ = _pad3(np.asarray(mf.basic_dof_nodes(dofs)))
    _, first, inv = np.unique(np.round(X, 10), axis=0, return_index=True,
                              return_inverse=True)
    inv = inv.reshape(-1)
    labels = []
    for g in range(len(first)):
        ids = dofs[inv == g]
        labels.append(",".join(map(str, ids)) if len(ids) <= 3
                      else f"{ids[0]}..{ids[-1]}")
    return _displaced(X[first], deform), labels


def _quality_values(mesh, data):
    """Quality of each convex, as values on the slice points (0 <= Q <= 1)."""
    cvids = np.asarray(mesh.cvid(), dtype=int).reshape(-1)
    Q = np.zeros(cvids.max() + 1)
    Q[cvids] = np.asarray(mesh.quality(cvids), dtype=float).reshape(-1)
    # np.real: the interface currently returns a complex array for real input
    return np.real(np.asarray(data.meta["slice"].interpolate_convex_data(Q))).reshape(-1)


def _draw_labels(out, backend, X, labels, color, size, sdim, shift=False):
    """Draw text labels at points X; shift=True offsets them slightly (2D, as
    Matlab gf_plot_mesh) so that vertex / dof numbers do not hide the edges."""
    if len(labels) == 0:
        return
    if shift and backend != "pyvista" and sdim <= 2:
        ext = np.ptp(X[:, :2], axis=0)
        X = X.copy()
        X[:, :2] += np.where(ext > 0, ext, 1.0) / 60.0
    if backend == "pyvista":
        out.add_point_labels(X, labels, font_size=size, text_color=color,
                             shape_opacity=0.0, show_points=False,
                             always_visible=False)
        return
    ax = out[1]
    for x, t in zip(X, labels):
        if sdim == 3:
            ax.text(x[0], x[1], x[2], t, color=color, fontsize=size,
                    ha="center", va="center")
        else:
            ax.text(x[0], x[1], t, color=color, fontsize=size, ha="center",
                    va="center", clip_on=False)


# ---------------------------------------------------------------------------
#  Backend selection
# ---------------------------------------------------------------------------

_PYVISTA_ERROR = None


def _has_pyvista():
    """True if pyvista can be imported; otherwise keep the reason of the failure."""
    global _PYVISTA_ERROR
    try:
        import pyvista  # noqa: F401
        return True
    except Exception as e:  # ImportError, but also broken vtk / numpy ABI, OpenGL...
        import sys
        _PYVISTA_ERROR = f"{type(e).__name__}: {e} (python: {sys.executable})"
        return False


def available_backends():
    """List of the renderers that can be used: 'pyvista' and/or 'matplotlib'."""
    out = ["pyvista"] if _has_pyvista() else []
    try:
        import matplotlib  # noqa: F401
        out.append("matplotlib")
    except Exception:
        pass
    return out


def _choose_backend(backend, data):
    if backend in ("matplotlib", "mpl"):
        return "mpl"
    if backend in ("pyvista", "vtk"):
        if not _has_pyvista():
            raise ImportError("backend 'pyvista' requested but pyvista cannot be "
                              f"imported: {_PYVISTA_ERROR}")
        return "pyvista"
    if data.sdim == 3:
        if _has_pyvista():
            return "pyvista"
        warnings.warn("pyvista cannot be imported, falling back to matplotlib for 3D "
                      f"(slow, approximate depth sorting).\n  Reason: {_PYVISTA_ERROR}",
                      stacklevel=3)
    return "mpl"


AXIS_LABELS = ("x", "y", "z")


# ---------------------------------------------------------------------------
#  Matplotlib backend
# ---------------------------------------------------------------------------

def _mpl_axes(ax, sdim):
    import matplotlib.pyplot as plt
    if ax is not None:
        return ax.figure, ax
    if sdim == 3:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
    else:
        fig, ax = plt.subplots()
    return fig, ax


def _mpl_decorate(ax, data, title, axis_labels, show_axes):
    if title:
        ax.set_title(title)
    if not show_axes:
        ax.set_axis_off()
        return
    labels = axis_labels or AXIS_LABELS
    ax.set_xlabel(labels[0])
    if data.sdim >= 2 or data.tdim == 1:
        ax.set_ylabel(labels[1] if data.sdim >= 2 else (data.name or "u"))
    if data.sdim == 3:
        ax.set_zlabel(labels[2])
        X = data.points
        ax.set_box_aspect(np.ptp(X, axis=0) + 1e-12)
    elif data.sdim == 2:
        ax.set_aspect("equal")


def _mpl_edges(ax, data, color, lw):
    from matplotlib.collections import LineCollection
    if data.edges is None:
        return
    X, E = data.edges
    segs = X[E]
    if data.sdim == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        ax.add_collection3d(Line3DCollection(segs, colors=color, linewidths=lw))
    else:
        ax.add_collection(LineCollection(segs[:, :, :2], colors=color, linewidths=lw))


def _render_mpl(data, ax=None, component=None, cmap="viridis", clim=None,
                contours=0, quiver=False, quiver_count=400, quiver_scale=None,
                edge_color="k", edge_width=0.4, face_color="lightgray",
                line_width=2.0, line_color="k", quiver_color="k",
                title=None, colorbar=True, axis_labels=None, show_axes=True,
                show=False, show_faces=True, **_):
    import matplotlib.pyplot as plt
    previous = None  # limits of what is already drawn in the given axes
    if ax is not None and ax.has_data() and data.sdim <= 2:
        previous = (ax.get_xlim(), ax.get_ylim())
    fig, ax = _mpl_axes(ax, data.sdim)
    X = data.points
    c = data.scalars(component)
    label = _scalar_label(data, component)
    mappable = None

    if data.sdim == 3:  # surface in 3D (fallback, prefer pyvista)
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        T = data.simplexes.get(2, np.empty((0, 3), int))
        poly = Poly3DCollection(X[T], linewidths=0)
        if c is not None and show_faces:
            import matplotlib as mpl
            vals = c[T].mean(axis=1)
            norm = mpl.colors.Normalize(*(clim or (c.min(), c.max())))
            poly.set_facecolor(plt.get_cmap(cmap)(norm(vals)))
            mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        elif c is None and show_faces:
            poly.set_facecolor(face_color)
        else:
            poly.set_facecolor((0, 0, 0, 0))
        ax.add_collection3d(poly)
        L = data.simplexes.get(1, np.empty((0, 2), int))
        if len(L) and data.tdim == 1:  # curves in 3D (e.g. streamlines)
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            ax.add_collection3d(Line3DCollection(X[L], colors=line_color,
                                                 linewidths=line_width))
        _mpl_edges(ax, data, edge_color, edge_width)
        lo, hi = X.min(axis=0), X.max(axis=0)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])

    elif data.tdim == 2:  # 2D field / mesh
        import matplotlib.tri as mtri
        tri = mtri.Triangulation(X[:, 0], X[:, 1], data.simplexes[2])
        if c is not None:
            kw = {} if clim is None else {"vmin": clim[0], "vmax": clim[1]}
            if not show_faces:
                pass
            elif contours:
                mappable = ax.tricontourf(tri, c, levels=contours, cmap=cmap, **kw)
            else:
                mappable = ax.tripcolor(tri, c, shading="gouraud", cmap=cmap, **kw)
        elif data.edges is None:
            ax.triplot(tri, color=edge_color, lw=edge_width)
        elif show_faces:
            ax.tripcolor(tri, np.zeros(len(X)), cmap="Greys", vmin=-1, vmax=4,
                         shading="flat")
        _mpl_edges(ax, data, edge_color, edge_width)
        if quiver and data.qdim >= 2:
            idx = _subsample(X, quiver_count)
            V = data.values
            ax.quiver(X[idx, 0], X[idx, 1], V[idx, 0], V[idx, 1],
                      color=quiver_color, scale=quiver_scale, angles="xy")
        lo, hi = X[:, :2].min(axis=0), X[:, :2].max(axis=0)
        pad = 0.03 * (hi - lo).max()
        lo, hi = lo - pad, hi + pad
        if previous is not None:
            lo = np.minimum(lo, [previous[0][0], previous[1][0]])
            hi = np.maximum(hi, [previous[0][1], previous[1][1]])
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])

    elif data.tdim == 1:  # 1D: curve u(x), or mesh edges / boundary in 2D
        from matplotlib.collections import LineCollection
        S = data.simplexes[1]
        if data.sdim == 1 and c is not None:
            segs = np.stack([X[S, 0], c[S]], axis=-1)
            ax.add_collection(LineCollection(segs, colors="C0"))
        else:
            segs = X[S][:, :, :2]
            ax.add_collection(LineCollection(segs, colors=line_color,
                                             linewidths=line_width))
        ax.autoscale_view()

    if mappable is not None and colorbar:
        cb = fig.colorbar(mappable, ax=ax, shrink=0.8)
        cb.set_label(label)
    _mpl_decorate(ax, data, title, axis_labels, show_axes)
    if show:
        plt.show()
    return fig, ax


def _subsample(X, count):
    """About `count` well spread distinct points (for arrows)."""
    _, uniq = np.unique(np.round(X, 12), axis=0, return_index=True)
    uniq = np.sort(uniq)
    if len(uniq) <= count:
        return uniq
    return uniq[np.linspace(0, len(uniq) - 1, count).astype(int)]


def _scalar_label(data, component):
    name = data.name or "u"
    if data.qdim <= 1:
        return name
    if component is None or component == "norm":
        return f"norm({name})"
    c = int(component)
    return f"{name}_{'xyz'[c]}" if data.qdim <= 3 else f"{name}[{c}]"


# ---------------------------------------------------------------------------
#  PyVista backend
# ---------------------------------------------------------------------------

_VTK_CELL = {0: 1, 1: 3, 2: 5, 3: 10}  # vertex, line, triangle, tetra


def to_pyvista(data, component=None):
    """SliceData -> pyvista.UnstructuredGrid (point data: scalars and vectors)."""
    import pyvista as pv
    cells, types = [], []
    for d in sorted(data.simplexes):
        S = data.simplexes[d]
        if len(S) == 0 or (d == 0 and data.tdim > 0):
            continue
        cells.append(np.hstack([np.full((len(S), 1), d + 1), S]).ravel())
        types.append(np.full(len(S), _VTK_CELL[d], dtype=np.uint8))
    grid = pv.UnstructuredGrid(np.concatenate(cells) if cells else np.empty(0, int),
                               np.concatenate(types) if types else np.empty(0, np.uint8),
                               data.points)
    if data.values is not None:
        grid.point_data[_scalar_label(data, component)] = data.scalars(component)
        if data.qdim >= 2:
            grid.point_data[data.name or "u"] = data.vectors3()
    return grid


_VECTOR_FORMATS = (".svg", ".eps", ".ps", ".pdf", ".tex")


def _save_view(pl, filename, scale=1):
    """Save the current view of a PyVista plotter: bitmap (png, jpg, ...) or
    vector graphics (svg, eps, ps, pdf, tex, through VTK/gl2ps)."""
    import os
    ext = os.path.splitext(filename)[1].lower()
    if ext in _VECTOR_FORMATS:
        pl.save_graphic(filename)
    else:
        pl.screenshot(filename, scale=scale)
    return os.path.abspath(filename)


def _install_save_key(pl, key, pattern, scale=1):
    """Bind a key of the interactive window to the saving of the current view.

    pattern is a file name, possibly with a counter field, e.g.
    'getfem_plot_{:03d}.png': each key press then writes a new file without
    overwriting the existing ones.
    """
    if not key or getattr(pl, "iren", None) is None:
        return  # off-screen or notebook plotter: no keyboard
    import os
    counter = [0]

    def save():
        fname = pattern
        if "{" in pattern:
            while True:
                counter[0] += 1
                fname = pattern.format(counter[0])
                if not os.path.exists(fname):
                    break
        try:
            print("getfem_plot: view saved in", _save_view(pl, fname, scale))
        except Exception as e:
            print("getfem_plot: could not save %s (%s)" % (fname, e))

    try:
        pl.add_key_event(key, save)
    except Exception as e:
        warnings.warn("cannot bind key %r to save the view: %s" % (key, e))


def _render_pyvista(data, plotter=None, component=None, cmap="viridis", clim=None,
                    contours=0, quiver=False, quiver_count=400, quiver_scale=None,
                    edge_color="black", edge_width=1.0, face_color="lightgray",
                    line_width=3.0, line_color="black", quiver_color="black",
                    title=None, colorbar=True, axis_labels=None, show_axes=True,
                    axes_style="grid", show=True, notebook=None, off_screen=False,
                    screenshot=None, show_faces=True, save_key="i",
                    save_filename="getfem_plot_{:03d}.png", save_scale=1, **_):
    import pyvista as pv
    if plotter is not None:
        pl = plotter
    else:
        pl = pv.Plotter(notebook=notebook, off_screen=off_screen)
        if not off_screen and not notebook:
            _install_save_key(pl, save_key, save_filename, save_scale)
    grid = to_pyvista(data, component)
    label = _scalar_label(data, component)
    if data.values is not None and show_faces:
        pl.add_mesh(grid, scalars=label, cmap=cmap, clim=clim,
                    n_colors=contours if contours else 256,
                    show_scalar_bar=colorbar,
                    scalar_bar_args={"title": label, "vertical": True})
    elif data.values is None and show_faces:
        if data.tdim == 1:  # curves (e.g. streamlines), drawn as tubes
            pl.add_mesh(grid, color=line_color, line_width=line_width,
                        render_lines_as_tubes=True)
        else:
            pl.add_mesh(grid, color=face_color)
    if data.edges is not None:
        X, E = data.edges
        lines = np.hstack([np.full((len(E), 1), 2), E]).ravel()
        pl.add_mesh(pv.PolyData(X, lines=lines), color=edge_color, line_width=edge_width)
    if quiver and data.qdim >= 2:
        idx = _subsample(data.points, quiver_count)
        pts = pv.PolyData(data.points[idx])
        V = data.vectors3()[idx]
        pts["vec"] = V
        vmax = np.linalg.norm(V, axis=1).max() or 1.0
        factor = quiver_scale or 0.05 * float(np.ptp(data.points, axis=0).max()) / vmax
        pl.add_mesh(pts.glyph(orient="vec", scale="vec", factor=factor),
                    color=quiver_color)
    if show_axes:
        lx, ly, lz = axis_labels or AXIS_LABELS
        pl.add_axes(xlabel=lx, ylabel=ly, zlabel=lz)  # orientation triad
        if axes_style in ("grid", "bounds"):
            pl.show_bounds(grid="back" if axes_style == "grid" else False,
                           location="outer", xtitle=lx, ytitle=ly, ztitle=lz)
    if title:
        pl.add_title(title, font_size=12)
    if data.sdim <= 2:
        pl.view_xy()
    if screenshot:
        pl.show(screenshot=screenshot, auto_close=False)
    elif show and plotter is None:
        pl.show()
    return pl


def save_view(plotter, filename, scale=1):
    """Save the current view of a PyVista plotter returned by plot, plot_mesh
    or plot_slice (e.g. after pl.show(auto_close=False)).  The format is given
    by the extension: png, jpg, ... or svg, eps, ps, pdf, tex (vector).
    Returns the absolute path of the file."""
    return _save_view(plotter, filename, scale)


def _render(data, backend="auto", **kw):
    b = _choose_backend(backend, data)
    if b == "pyvista":
        kw.pop("ax", None)
        return _render_pyvista(data, **kw)
    kw.pop("plotter", None)
    kw.setdefault("show", False)
    for k in ("notebook", "off_screen", "screenshot", "axes_style",
              "save_key", "save_filename", "save_scale"):
        kw.pop(k, None)
    return _render_mpl(data, **kw)


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------

def plot(mf, U, slice_op=None, refine=3, cvlst=None, name="u",
         component=None, mesh_edges=True, faces=True,
         deformation=None, deformation_mf=None, deformation_scale="10%",
         backend="auto", **kw):
    """Plot a finite element field (mf, U).

    component : None / 'norm' (norm of a vector field), or a component index.
    quiver=True draws arrows for vector fields; contours=N draws N filled levels.
    mesh_edges : overlay the edges of the (possibly deformed) mesh.
    faces=False : do not draw the coloured surface (e.g. to draw arrows only).
    title, axis_labels=('x','y','z'), show_axes, colorbar, cmap, clim,
    edge_color, face_color, quiver_color, line_color, line_width.
    Matplotlib: ax=..., returns (fig, ax).  PyVista: plotter=..., show=...,
    notebook=..., off_screen=..., screenshot='file.png', axes_style='grid'|'bounds'|'triad';
    returns the Plotter (rotation / zoom with the mouse).  In the PyVista
    window, the key save_key ('i' by default) saves the current view in
    save_filename ('getfem_plot_{:03d}.png' by default, numbered so as not to
    overwrite; .svg/.eps/.pdf give vector graphics), magnified by save_scale.
    """
    data = slice_field(mf, U, slice_op, refine, cvlst, deformation, deformation_mf,
                       deformation_scale, name=name, with_edges=mesh_edges)
    return _render(data, backend, component=component, show_faces=faces, **kw)


def plot_slice(sl, mf=None, U=None, name="u", mesh_edges=False, faces=True,
               backend="auto", **kw):
    """Plot an existing gf.Slice, optionally coloured by the field (mf, U).

    Without field, the slice is drawn with face_color; a slice made of curves
    (e.g. gf.Slice('streamlines', ...)) is drawn with line_color and line_width
    (as tubes with PyVista).
    Several slices can be drawn in the same scene by passing plotter= (PyVista)
    or ax= (Matplotlib). Other keywords: as for plot().
    """
    data = slice_to_data(sl, mf, U, name=name, with_edges=mesh_edges)
    return _render(data, backend, show_faces=faces, **kw)


def plot_mesh(mesh, slice_op=None, refine=1, cvlst=None, regions=None,
              region_colors=None, vertices=False, convexes=False, dof=False,
              quality=False, edges=True, faces=True,
              deformation=None, deformation_mf=None, deformation_scale=1.0,
              label_size=10, backend="auto", **kw):
    """Plot a mesh (its surface in 3D), as Matlab gf_plot_mesh.

    mesh       : a Mesh, or a MeshFem (required for dof=True; its mesh is used).
    regions    : region numbers to highlight, each with a colour and a legend
                 entry, or a dict {region number: label for the legend}.
    vertices   : display the vertex numbers.
    convexes   : display the convex numbers.
    dof        : display the basic dof numbers at their nodes (MeshFem only);
                 dofs sharing a node (qdim > 1) are shown as "first..last".
    quality    : colour each convex by its quality estimate (0 bad, 1 good).
    edges, faces : draw the mesh edges / fill the faces.
    cvlst      : restrict to these convexes (or convex faces, 2-row array).
    In 3D only the skin is drawn and only the numbers lying on it are shown.
    Other keywords: as for plot() (title, axis_labels, cmap, ax, plotter, ...).
    """
    mf = mesh if hasattr(mesh, "linked_mesh") else None
    if mf is not None:
        mesh = mf.linked_mesh()
        if deformation is not None and deformation_mf is None:
            deformation_mf = mf
    if dof and mf is None:
        raise ValueError("dof=True requires a MeshFem instead of a Mesh")
    data = slice_mesh(mesh, slice_op, refine, cvlst, deformation, deformation_mf,
                      deformation_scale)
    if not edges:
        data.edges = None
    if quality:
        data.values = _quality_values(mesh, data)
        data.name = "quality"
        kw.setdefault("clim", (0.0, 1.0))
        kw.setdefault("cmap", "RdYlGn")
    b = _choose_backend(backend, data)
    show = kw.pop("show", b == "pyvista")
    screenshot = kw.pop("screenshot", None)
    out = _render(data, b, show=False, show_faces=faces, **kw)
    if regions:
        _plot_regions(mesh, regions, region_colors, refine, b, out,
                      deformation, deformation_mf, deformation_scale)

    surface = mesh.dim() == 3 and slice_op is None
    deform = data.meta.get("deform")
    if vertices:
        X, lab = _vertex_labels(mesh, cvlst, surface, deform)
        _draw_labels(out, b, X, lab, "black", label_size, data.sdim, shift=True)
    if convexes:
        X, lab = _convex_labels(mesh, cvlst, surface, data)
        _draw_labels(out, b, X, lab, "darkred", label_size, data.sdim)
    if dof:
        X, lab = _dof_labels(mf, cvlst, surface, deform)
        _draw_labels(out, b, X, lab, "darkblue", label_size, data.sdim, shift=True)

    if b == "pyvista":
        if screenshot:
            out.show(screenshot=screenshot, auto_close=False)
        elif show and kw.get("plotter") is None:
            out.show()
    elif show:
        import matplotlib.pyplot as plt
        plt.show()
    return out


def _plot_regions(mesh, regions, colors, refine, backend, out,
                  deformation, deformation_mf, deformation_scale):
    import matplotlib.pyplot as plt
    labels = dict(regions) if isinstance(regions, dict) else {r: f"region {r}"
                                                              for r in regions}
    # default colour given by the region number, so that the regions of
    # several meshes drawn in the same figure have distinct colours
    cycle = colors or [plt.get_cmap("tab10")(int(r) % 10) for r in labels]
    for r, col in zip(labels, cycle):
        rd = slice_mesh(mesh, ("none",), refine, mesh.region(r), deformation,
                        deformation_mf, deformation_scale)
        label = labels[r]
        d = rd.tdim
        S = rd.simplexes.get(d)
        if S is None or len(S) == 0:
            continue
        if backend == "pyvista":
            import pyvista as pv
            out.add_mesh(to_pyvista(SliceData(rd.points, rd.sdim, {d: S})),
                         color=col[:3] if isinstance(col, tuple) else col,
                         line_width=4, label=label)
        else:
            fig, ax = out
            from matplotlib.collections import LineCollection, PolyCollection
            P = rd.points[S]
            if rd.sdim == 3:
                from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
                coll = (Poly3DCollection(P, facecolor=col, alpha=0.9) if d == 2
                        else Line3DCollection(P, colors=col, linewidths=3))
                ax.add_collection3d(coll)
                coll.set_label(label)
            elif d == 1:
                ax.add_collection(LineCollection(P[:, :, :2], colors=[col],
                                                 linewidths=3, label=label))
            else:
                ax.add_collection(PolyCollection(P[:, :, :2], facecolors=[col],
                                                 alpha=0.6, label=label))
    if backend == "pyvista":
        out.add_legend()
    else:
        out[1].legend(loc="best")
