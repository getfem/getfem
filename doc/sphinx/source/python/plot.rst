.. $Id$

.. include:: ../replaces.txt

.. highlight:: python

.. _py-plot:

Visualisation of meshes and fields
==================================

The module ``getfem_plot`` displays meshes and finite element fields directly
from |py|, in the same spirit as the functions ``gf_plot_mesh``, ``gf_plot``
and ``gf_plot_slice`` of the Matlab/Octave interface. It relies on two
optional libraries, which are only imported when they are used:

* `Matplotlib <https://matplotlib.org>`_, for 1D and 2D plots (well suited to
  figures for papers) and small 3D surfaces;
* `PyVista <https://pyvista.org>`_ (based on VTK), for 3D plots that can be
  rotated, zoomed and panned with the mouse, also inside Jupyter notebooks, or
  rendered off-screen to image files.

They can be installed with::

  $ pip install matplotlib pyvista

In 3D, PyVista is used when it is available; otherwise the module falls back
on Matplotlib with a warning (Matplotlib is slow for large surfaces and only
sorts the faces approximately in depth). The renderer can be chosen
explicitly with ``backend="matplotlib"`` or ``backend="pyvista"``.


Quick start
-----------

::

  import numpy as np
  import getfem as gf
  import getfem_plot as gfp

  m = gf.Mesh('regular simplices', np.linspace(0, 2, 21), np.linspace(0, 1, 11))
  m.set_region(1, m.outer_faces())
  mf = gf.MeshFem(m, 1); mf.set_classical_fem(2)
  # ... compute a solution U on mf ...

  gfp.plot_mesh(m, regions=[1], vertices=True, show=True)
  gfp.plot(mf, U, name="u", refine=4, title="Solution", show=True)

For a 3D elasticity problem with a displacement ``U`` on a vector
``MeshFem`` ``mfu``::

  # norm of the displacement on the deformed surface (interactive window)
  gfp.plot(mfu, U, name="u", deformation=U, deformation_scale="20%")

  # vertical component on the half y > 0.5, with arrows
  gfp.plot(mfu, U, name="u", component=2, quiver=True,
           slice_op=('boundary', ('planar', +1, [0, 0.5, 0], [0, 1, 0])))

  # image file, without opening a window
  gfp.plot(mfu, U, name="u", off_screen=True, screenshot="u.png")

A complete example is given by the test ``check_plot.py`` in the directory
``interface/tests/python/`` of the |gf| distribution.


How it works
------------

A field is not drawn directly: it is first evaluated on a *slice* of the mesh
(see the ``Slice`` object), then the points, simplexes and values of the slice
are given to the renderer as |np| arrays. This has several consequences:

* In 3D, only the surface is drawn: the default slice operation is
  ``('boundary',)``, which keeps only the triangles of the skin of the
  domain. Any other slice operation can be passed with ``slice_op`` (planar
  cuts, iso-values, balls, intersections, ...).
* High order elements and curved geometric transformations are represented by
  refining each convex ``refine`` times (default 3 for fields).
* The points of the slice are not shared between convexes, so discontinuous
  fields keep their jumps.
* The deformation is applied by building the slice on the deformed mesh, so
  that the colour, the arrows and the mesh edges all follow it.


Plotting a field: ``plot``
--------------------------

``gfp.plot(mf, U, **options)`` draws the field ``U`` defined on ``mf``:

=========================== ==================================================
option                      meaning
=========================== ==================================================
``name``                    name of the field, used in the colour bar
``component``               ``None`` or ``'norm'`` (norm of a vector field,
                            the default), or the index of a component
``slice_op``, ``refine``    slice operation and refinement (see above)
``cvlst``                   restrict to a list of convexes or a region
``deformation``             displacement vector, defined on
                            ``deformation_mf`` (default ``mf``)
``deformation_scale``       a number, or ``'p%'`` for a maximal displacement
                            of p% of the size of the mesh (default ``'10%'``)
``mesh_edges``              draw the edges of the mesh (default ``True``)
``faces``                   draw the coloured surface (``False`` to draw only
                            the arrows or the edges)
``quiver``                  draw arrows for a vector field;
                            ``quiver_count`` sets their number
``contours``                draw this number of iso-value levels
``cmap``, ``clim``          colour map and range of the colour bar
``title``, ``axis_labels``  title and names of the axes
``edge_color``              colour of the mesh edges
``face_color``              colour of the faces when there is no data
``quiver_color``            colour of the arrows
``line_color``              colour of the curves (e.g. streamlines)
``line_width``              width of the curves
``show_axes``, ``colorbar`` display the axes / the colour bar
=========================== ==================================================

With Matplotlib, ``ax=`` draws in existing axes and ``(fig, ax)`` is returned;
use ``show=True`` to open the window. With PyVista, ``plotter=`` adds to an
existing ``pyvista.Plotter``, which is returned; ``axes_style`` is
``'grid'`` (default, graduated axes), ``'bounds'`` or ``'triad'``, and
``notebook``, ``off_screen`` and ``screenshot`` are passed to PyVista.

``gfp.plot_slice(sl, mf, U)`` does the same for any ``Slice`` built by the
user, or draws the slice itself when no field is given (a slice made of
curves, such as streamlines, is drawn as tubes with PyVista).

Several slices can be combined in the same scene by passing the PyVista
plotter (or the Matplotlib axes) returned by the first call::

  pl = gfp.plot_slice(sl, mfp, P, name="p", backend="pyvista", show=False)
  gfp.plot_slice(sl, mfu, U, faces=False, quiver=True, colorbar=False, plotter=pl)
  gfp.plot_slice(gf.Slice('streamlines', mfu, U, H), line_color="white", plotter=pl)
  pl.show()

See ``demo_stokes_3D_tank_draw.py`` in ``interface/tests/python/`` for a
complete example. ``gfp.available_backends()`` returns the renderers which are
installed.


Saving figures
--------------

In the interactive PyVista window, the key ``i`` saves the current view
(after any rotation or zoom) in the files ``getfem_plot_001.png``,
``getfem_plot_002.png``, ... of the current directory, without overwriting
existing files. This is set by the options ``save_key``, ``save_filename``
(possibly with a counter field ``{:03d}``) and ``save_scale`` (magnification
of the image with respect to the window)::

  gfp.plot(mfu, U, name="u", save_key="s", save_filename="tank_{:02d}.pdf")

The format is given by the extension: bitmap (``png``, ``jpg``, ...) or vector
graphics (``svg``, ``eps``, ``ps``, ``pdf``, ``tex``). In vector formats,
the colour gradients of the surfaces are less faithful than in a bitmap of high
resolution (``save_scale=3`` for instance).

From a script, ``gfp.save_view(pl, "fig.png")`` saves the view of a plotter
``pl``. To keep the view chosen interactively, show the window with
``pl.show(auto_close=False)``: after closing it, ``gfp.save_view`` saves the
last view as a bitmap image, and ``pl.camera_position`` can be stored to reproduce the same view
later off-screen. Images can also be produced without any window with
``off_screen=True, screenshot="file.png"``.

The Matplotlib windows have their own button to save the figure, and
``fig.savefig("file.pdf")`` can be used in scripts.


Plotting a mesh: ``plot_mesh``
------------------------------

``gfp.plot_mesh(m, **options)`` draws a mesh (its skin in 3D). ``m`` can also
be a ``MeshFem``, which is required for the dof numbers.

=========================== ==================================================
option                      meaning
=========================== ==================================================
``regions``                 regions to highlight, each with its own colour and
                            an entry in the legend
``vertices``                display the vertex numbers
``convexes``                display the convex numbers
``dof``                     display the basic dof numbers at their nodes (the
                            dofs sharing a node are displayed together)
``quality``                 colour each convex with its quality (0 to 1)
``edges``, ``faces``        draw the edges / fill the faces
``cvlst``                   restrict to these convexes, or to these faces with
                            a two-row array
``deformation``             as for ``plot``; the numbers follow the deformed
                            mesh
=========================== ==================================================

In 3D, only the numbers lying on the skin are displayed and PyVista hides the
ones which are behind the surface.


Using the extracted data
------------------------

The functions ``gfp.slice_field`` and ``gfp.slice_mesh`` only extract the data
and return a ``SliceData`` object, with the attributes ``points`` (n x 3
array), ``simplexes`` (dictionary dimension -> array of point indices),
``values`` (n or n x q array), ``edges`` and ``name``. It can be used with any
other plotting tool, and ``gfp.to_pyvista(data)`` converts it to a
``pyvista.UnstructuredGrid``.
