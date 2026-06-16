.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-export:

Export and view a solution
==========================

There are essentially four ways to view the result of getfem computations:

* Scilab, Octave or Matlab, with the interface.
* The open-source Paraview, Mayavi2, PyVista or any other VTK/VTU file viewer.
* The open-source OpenDX program.
* The open-source Gmsh program.

The objects that can be exported are, |m|, |mf| objects, and |smsl|.

Saving mesh and mesh_fem objects for the Matlab interface
---------------------------------------------------------

If you have installed the Scilab, Octave or Matlab interface, you can simply use
``mesh_fem::write_to_file`` and save the solution as a plain text file, and then,
load them with the interface. For example, supposing you have a solution ``U`` on a |mf|
``mf``,::

  std::fstream f("solution.U",std::ios::out);
  for (unsigned i=0; i < gmm::vect_size(U); ++i)
    f << U[i] << "\verb+\+n";

  // when the 2nd arg is true, the mesh is saved with the |mf|
  mf.write_to_file("solution.mf", true);

and then, under Scilab, Octave or Matlab:

.. code-block:: matlab

   >> U=load('solution.U');
   >> mf=gfMeshFem('load','solution.mf');
   >> gf_plot(mf,U,'mesh','on');

See the getfem-matlab interface documentation for more details.

Four file formats are supported for export: the `VTK`_ and `VTU`_ file
formats, the`OpenDX`_ file format and the `Gmsh`_ post-processing file
format. All four can be used for exporting either a |gf_m| or |gf_mf|, and
all except `VTU`_ can be used for exporting the more versatile |gf_smsl|.
The corresponding four classes: |gf_vtk_export|, |gf_vtu_export|,
|gf_dx_export| and |gf_pos_export| are contained in the file
:file:`getfem/getfem_export.h`.

Examples of use can be found in the examples of the tests directory.

In addition, if |gf| is configured with ``--enable-exodus``, the `Exodus II
<https://sandialabs.github.io/seacas-docs/>`_ format can be used to export and
re-import a |gf_m| or |gf_mf| (including transient results); see
:ref:`ud-export-exodus` below.

.. _ud-export_slices:

Producing mesh slices
---------------------

|gf| provides "slicers" objects which are dedicated to generating post-treatment
data from meshes and solutions. These slicers, defined in the file
:file:`getfem/getfem_mesh_slicers.h` take a |m| (and sometimes a |mf| with a
solution field) on input, and produce a set of simplices after applying some
operations such as *intersection with a plane*, *extraction of the mesh
boundary*, *refinement of each convex*, *extraction of isosurfaces*, etc. The
output of these slicers can be stored in a |gf_smsl| object (see the file
:file:`getfem/getfem_mesh_slice.h`). A |smsl| object may be considered as a P1
discontinuous FEM on a non-conformal mesh with fast interpolation ability. Slices
are made of segments, triangles and tetrahedrons, so the convexes of the original
mesh are always simplexified.

All slicer operation inherit from |gf_sl_a|, it is very easy to create a new
slicer. Example of slicers are (some of them use a |gf_sl_ddb| which is just a
reference to a |mf| ``mf`` and a field ``U`` on this |mf|).

.. cpp:function:: getfem::slicer_none()

   empty slicer.

.. cpp:function:: getfem::slicer_boundary(const mesh &m, ldots)

   extract the boundary of a mesh.

.. cpp:function:: getfem::slicer_apply_deformation(mesh_slice_cv_dof_data_base &)

   apply a deformation to the mesh , the deformation field is defined on a |mf|.

.. cpp:function:: getfem::slicer_half_space(base_node x0, base_node n, int orient)

   cut the mesh with a half space (if ``orient`` = -1 or +1), or a plane (if
   ``orient`` = 0), ``x0`` being a node of the plane, and ``n`` being a normal
   of the plane.

.. cpp:function::  getfem::slicer_sphere(base_node x0, scalar_type R, int orient)

   cut with the interior (``orient``=-1), boundary (``orient``=0) or exterior
   (``orient``=+1) or a sphere of center ``x0`` and radius ``R``.

.. cpp:function:: getfem::slicer_cylinder(base_node x0, base_node x1, scalar_type R, int orient)

   slice with the interior/boundary/exterior of a cylinder of axis ``(x0,x1)``
   and radius ``R``.

.. cpp:function:: getfem::slicer_isovalues(const mesh_slice_cv_dof_data_base& mfU, scalar_type val, int orient)

   cut with the isosurface defined by the scalar field ``mfU`` and ``val``.
   Keep only simplices where ::math:`u(x)<val` (``orient``=-1), :math:`u(x)=val`
   (``orient=0`` or :math:`u(x)>val`.

.. cpp:function:: getfem::slicer_mesh_with_mesh(const mesh& m2)

   cut the convexes with the convexes of the mesh ``m2``.

.. cpp:function:: getfem::slicer_union(const slicer_action &sA, const slicer_action &sB)

   merges the output of two slicer operations.

.. cpp:function:: getfem::slicer_intersect(slicer_action &sA, slicer_action &sB)

   intersect the output of two slicer operations.

.. cpp:function:: getfem::slicer_complementary(slicer_action &s)

   return the complementary of a slicer operation.

.. cpp:function:: getfem::slicer_build_edges_mesh(mesh& edges_m)

   slicer whose side-effect is to build the mesh ``edges_m`` with the edges of
   the sliced mesh.

.. cpp:function:: getfem::slicer_build_mesh(mesh &m)

   in some (rare) occasions , it might be useful to build a mesh from a slice.
   Note however that there is absolutely no guaranty that the mesh will be
   conformal (although it is often the case).

.. cpp:function:: getfem::slicer_build_stored_mesh_slice(stored_mesh_slice& sl)

   record the output of the slicing operation into a |smsl| object. Note that it
   is often more convenient to use the ``stored_mesh_slice::build(...)`` method to
   achieve the same result.

.. cpp:function:: getfem::slicer_explode(c)

   shrink or expand each convex with respect to its gravity center.

In order to apply these slicers, a ``getfem::mesh_slicer(mesh&)`` object should be
created, and the |gf_sl_a| are then stacked with
``mesh_slicer::push_back_action(slicer_action&)`` and
``mesh_slicer::push_front_action(slicer_action&)``. The slicing operation is
finally executed with ``mesh_slicer::exec(int nrefine)`` (or
``mesh_slicer::exec(int nrefine, const mesh_region &cvlst)`` to apply the operation
to a subset of the mesh, or its boundary etc.).

The ``nrefine`` parameter is very important, as the "precision" of the final result
will depend on it: if the data that is represented on the final slice is just P1
data on convexes with a linear geometric transformation, ``nrefine = 1`` is the
right choice, but for P2, P3, non linear transformation etc, it is better to refine
each convex of the original mesh during the slicing operation. This allows an
accurate representation of any finite element field onto a very simple structure
(linear segment/triangles/tetrahedrons with P1 discontinuous data on them) which is
what most visualization programs (gmsh, mayavi, opendx, scilab, octave, matlab, etc.) expect.

Example of use (cut the boundary of a mesh ``m`` with a half-space, and save the
result into a |smsl|)::

  getfem::slicer_boundary a0(m);
  getfem::slicer_half_space a1(base_node(0,0), base_node(1, 0), -1);
  getfem::stored_mesh_slice sl;
  getfem::slicer_build_stored_mesh_slice a2(sl);
  getfem::mesh_slicer slicer(m);
  slicer.push_back_action(a1);
  slicer.push_back_action(a2);
  int nrefine = 3;
  slicer.exec(nrefine);

In order to build a |gf_smsl| object during the slicing operation, the ``stored_mesh_slice::build()`` method is often more convenient than using explicitly the ``slicer_build_stored_mesh_slice`` slicer::

  getfem::stored_mesh_slice sl;
  sl.build(m, getfem::slicer_boundary(m),
           getfem::slicer_half_space(base_node(0,0), base_node(1, 0), -1),
           nrefine);

The simplest way to use these slices is to export them to |vtk|,
|opendx|, or |gmsh|.


Exporting |m|, |mf| or slices to VTK/VTU
-----------------------------------------

VTK/VTU files can handle data on segment, triangles, quadrangles,
tetrahedrons and hexahedrons of first or second degree.

For example, supposing that a |smsl| ``sl`` has already been built::

  // an optional the 2nd argument can be set to true to produce
  // a text file instead of a binary file
  vtk_export exp("output.vtk");
  exp.exporting(sl); // will save the geometrical structure of the slice
  exp.write_point_data(mfp, P, "pressure"); // write a scalar field
  exp.write_point_data(mfu, U, "displacement"); // write a vector field

In this example, the fields ``P`` and ``U`` are interpolated on the slice
nodes and then written into the VTK field.

It is also possible to export a |mf| ``mfu`` without having to build a slice::

  // an optional the 2nd argument can be set to true to produce
  // a text file instead of a binary file
  vtk_export exp("output.vtk");
  exp.exporting(mfu);
  exp.write_point_data(mfp, P, "pressure"); // write a scalar field
  exp.write_point_data(mfu, U, "displacement"); // write a vector field

An |mf| ``mfu`` can also be exported in the VTU format with::

  vtu_export exp("output.vtu);
  exp.exporting(mfu); // will save the geometrical structure of the mesh_fem
  exp.write_point_data(mfp, P, "pressure"); // write a scalar field
  exp.write_point_data(mfu, U, "displacement"); // write a vector field

Note however that when exporing a |mf| with ``vtk_export`` or ``vtu_export``
each convex/fem of ``mfu`` will be mapped to a VTK/VTU element type. As
VTK/VTU does not handle elements of degree greater than 2, there will be a
loss of precision for higher degree FEMs.

Exporting |m|, |mf| or slices to OpenDX
---------------------------------------

The OpenDX data file is more versatile than the VTK one. It is able to store more
that one mesh, any number of fields on these meshes etc. However, it does only
handle elements of degree 1 and 0 (segments, triangles, tetrahedrons, quadrangles
etc.). And each mesh can only be made of one type of element, it cannot mix
triangles and quadrangles in a same object. For that reason, it is generally
preferable to export |gf_smsl| objects (in which non simplex elements are
simplexified, and which allows refinement of elements) than |gf_mf| and |gf_m|
objects.

The basic usage is very similar to |gf_vtk_export|::

  getfem::dx_export exp("output.dx");
  exp.exporting(sl);
  exp.write_point_data(mfu, U, "displacement");

Moreover, |gf_dx_export| is able to reopen a '.dx' file and append new data into
it. Hence it is possible, if many time-steps are to be saved, to view intermediate
results in OpenDX during the computations. The prototype of the constructor is::

  dx_export(const std::string& filename, bool ascii = false, bool append = false);
  dx_export(std::ostream &os_, bool ascii = false);

An example of use, with multiple time steps (taken from
:file:`tests/dynamic_friction.cc`)::

  getfem::stored_mesh_slice sl;
  getfem::dx_export exp("output.dx", false);
  if (N <= 2) sl.build(mesh, getfem::slicer_none(),4);
  else        sl.build(mesh, getfem::slicer_boundary(mesh),4);
  exp.exporting(sl,true);

  // for each mesh object, a corresponding ``mesh'' object will be
  // created in the data file for the edges of the original mesh
  exp.exporting_mesh_edges();

  while (t <= T) {
    ...
    exp.write_point_data(mf_u, U0);
    exp.serie_add_object("deformation");
    exp.write_point_data(mf_vm, VM);
    exp.serie_add_object("von_mises_stress");
  }

In this example, an OpenDX "time series" is created, for each time step, two data
fields are saved: a vector field called "deformation", and a scalar field called
"von_mises_stress".

Note also that the ``dx_export::exporting_mesh_edges()`` function has been called.
It implies that for each mesh exported, the edges of the original mesh are also
exported (into another OpenDX mesh). In this example, you have access in OpenDX to
4 data fields: "deformation", "deformation_edges", "von_mises_stress" and
"von_mises_stress_edges".

The ``tests/dynamic_friction.net`` is an example of OpenDX program for these data
(run it with ``cd tests; dx -edit dynamic_friction.net`` , menu
"Execute/sequencer").

.. _ud-export-exodus:

Exodus II export and import
---------------------------

The `Exodus II <https://sandialabs.github.io/seacas-docs/>`_ format (defined by
Sandia's SEACAS project) is a finite-element database built on top of `NetCDF
<https://www.unidata.ucar.edu/software/netcdf/>`_. It stores a mesh together
with optional nodal results and, unlike the formats above, a single file
natively holds a transient (time-dependent) series.

|gf| support for Exodus is **optional**: it is only compiled when the library is
configured with ``--enable-exodus``, which requires the NetCDF library and
headers (the SEACAS library itself is *not* needed — the documented Exodus
layout is written directly through NetCDF). The classes ``getfem::exodus_export``
and ``getfem::exodus_import`` are declared in :file:`getfem/getfem_exodus.h`.

Writing a mesh and some fields is similar to the other exporters::

  #include "getfem/getfem_exodus.h"
  ...
  getfem::exodus_export exp("output.exo");
  exp.exporting(mf);                  // a mesh_fem (or a mesh)
  exp.write_mesh();
  exp.write_point_data(mf, U, "displacement");

As with the VTK export, the finite element and geometric transformations are
mapped to order 1 or 2 isoparametric Pk/Qk elements. The supported element types
are segments, triangles, quadrilaterals, tetrahedra, hexahedra and prisms (linear
and quadratic), written as the matching Exodus element types (``BAR2``/``BAR3``,
``TRI3``/``TRI6``, ``QUAD4``/``QUAD8``/``QUAD9``, ``TETRA4``/``TETRA10``,
``HEX8``/``HEX20``/``HEX27``, ``WEDGE6``/``WEDGE15``).

A transient series is written to a single file by opening a new time step before
writing each field::

  getfem::exodus_export exp("output.exo");
  exp.exporting(mf);
  exp.write_mesh();
  while (t <= T) {
    ...
    exp.set_time(t);                  // open a new time step at value t
    exp.write_point_data(mf, U, "u");
    exp.sync();                       // optional: flush this complete step
  }
  // the file is finalised when exp goes out of scope (or on exp.close())

Calling ``sync()`` after all fields for a time step have been written flushes the
completed step to disk. This is useful for long-running simulations: ParaView can
reload/refresh the same file and see the synced time steps while the simulation
continues. Incomplete steps are rejected instead of being deliberately synced
with fill values.

For best write performance when the list of result fields is known up front,
predeclare them before ``write_mesh()``::

  getfem::exodus_export exp("output.exo");
  exp.exporting(mf);
  exp.declare_point_data("u", mf.get_qdim());
  exp.write_mesh();

This lets the exporter create the transient variables during the initial NetCDF
definition phase and avoids a later file redefine. The Python/Matlab interface
uses this path automatically for ``export_to_exodus`` calls that create a new
file.

Appending time steps is supported only for files written by GetFEM with the same
exported mesh. New files store a compact GetFEM mesh fingerprint (dimension,
node coordinates, block layout and connectivity); append mode checks it before
writing so that a result cannot be appended to a file with a different block or
element layout.

|gf| regions are exported as Exodus sets, keyed by the region number: the face
entries of a region become a *side set*, its whole-convex entries become an
*element set*, and the nodes it touches become a *node set*. In addition, each
volume (whole-convex) region is written as its own Exodus *element block* whose
block id equals the region number, so a viewer such as ParaView can colour the
regions by the block id (``ObjectId``) or with ``vtkBlockColors``.

``enable_region_field()`` in C++ (or the ``'region field'`` keyword of the Python
``export_to_exodus``) additionally writes a ``region`` element (cell) variable
holding that block id per element, which gives a discrete cell field to colour
by. It is **off by default**: being constant in time, it would otherwise be
stored at every transient step. If it is explicitly enabled for a mesh-only
export, GetFEM writes one static Exodus time step containing only this cell
variable.

When GetFEM is built with usable NetCDF4/HDF5 deflate support, the Exodus writer
creates **compressed** NetCDF4 classic-model files by default (deflated numeric
arrays), which substantially reduces the size of large meshes and transient
series. Such files are still valid Exodus II and are read transparently by any
reader built with NetCDF4/HDF5 support, including modern ParaView. If that
support is not available at configure time, GetFEM defaults to classic
64-bit-offset files. To force classic output — for the widest reader
compatibility, or for very small meshes where the HDF5 container overhead can
make a compressed file *larger* — pass ``enable_compression(0)`` in C++ or the
``'uncompressed'`` keyword in the interface. The deflate level (1..9) can be set
with ``enable_compression(level)`` on builds with NetCDF4/HDF5 support.
Compression is fixed when the file is created and cannot be changed on
``append``.

By default, region blocks are named ``region_<id>`` and set names are left
empty. Names can be set explicitly with ``exp.set_region_name(id, "name")`` in
C++, or with the ``'region names'`` option of the Python ``export_to_exodus`` (an
id vector and a list of names of the same length).

When *importing* externally produced files, both the split (``coordx``/``coordy``/
``coordz``) and the older packed (``coord``) coordinate layouts are read, and
connectivity node ids are range-checked. Shell elements (``SHELL*``) are read as
2D surface elements. Side sets defined on shell elements are rejected with a clear
error because their top/bottom-face numbering differs from a 2D quad's edge
numbering.

The mesh is read back with ``getfem::import_mesh`` (side sets are restored as
face regions and element sets as convex regions)::

  getfem::mesh m;
  getfem::import_mesh("output.exo", "exodus", m);

The transient nodal variables can also be read back, which makes an
export/import round-trip verifiable::

  getfem::exodus_import imp("output.exo");
  getfem::mesh m;  imp.read_mesh(m);
  std::vector<double> U;
  imp.read_nodal_var("u", step, U);   // values at time step ``step``

The Python interface exposes the same operations:

.. code-block:: python

   mf.export_to_exodus('output.exo', U, 'u')        # compressed when supported
   mf.export_to_exodus('plain.exo', 'uncompressed', U, 'u')   # force classic
   # optionally name the blocks/sets matching regions 1 and 2:
   mf.export_to_exodus('named.exo', 'region names', [1, 2], ['left', 'right'], U, 'u')
   # request NetCDF4 compression when creating a new file:
   mf.export_to_exodus('small.exo', 'compress', U, 'u')
   m2 = gf.Mesh('import', 'exodus', 'output.exo')
   vals = m2.exodus_nodal_data('output.exo', 'u', 0)
