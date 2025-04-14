.. $Id$

.. include:: ../replaces.txt

.. highlight:: matlab

.. _mlab-mlabgf:

|gfm| organization
==================

The |gfm| toolbox is just a convenient interface to the |gf| library: you must
have a working |gf| installed on your computer. This toolbox provides a big
:envvar:`mex-file` (c++ binary callable from |octv| or |mlab|) and some additional
``m-files`` (documentation and extra-functionalities). All the functions of |gfm|
are prefixed by ``gf_`` (hence typing ``gf_`` at the |octv| or |mlab| prompt and then
pressing the ``<tab>`` key is a quick way to obtain the list of getfem
functions).


Functions
---------

* ``gf_workspace`` : workspace management.
* ``gf_util`` : miscellanous utility functions.
* ``gf_delete`` : destroy a |gf| object (|mlab_m| , |mlab_mf| , |mlab_mim| etc.).
* ``gf_cvstruct_get`` : retrieve informations from a |mlab_cs| object.
* ``gf_geotrans`` : define a geometric transformation.
* ``gf_geotrans_get`` : retrieve informations from a |mlab_gt| object.
* ``gf_mesh`` : creates a new |mlab_m| object.
* ``gf_mesh_get`` : retrieve informations from a |mlab_m| object.
* ``gf_mesh_set`` : modify a |mlab_m| object.
* ``gf_eltm`` : define an elementary matrix.
* ``gf_fem`` : define a |mlab_fem|.
* ``gf_fem_get`` : retrieve informations from a |mlab_fem| object.
* ``gf_integ`` : define a integration method.
* ``gf_integ_get`` : retrieve informations from an |mlab_int| object.
* ``gf_mesh_fem`` : creates a new |mlab_mf| object.
* ``gf_mesh_fem_get`` : retrieve informations from a |mlab_mf| object.
* ``gf_mesh_fem_set`` : modify a |mlab_mf| object.
* ``gf_mesh_im`` : creates a new |mlab_mim| object.
* ``gf_mesh_im_get`` : retrieve informations from a |mlab_mim| object.
* ``gf_mesh_im_set`` : modify a |mlab_mim| object.
* ``gf_slice`` : create a new |mlab_sl| object.
* ``gf_slice_get`` : retrieve informations from a |mlab_sl| object.
* ``gf_slice_set`` : modify a |mlab_sl| object.
* ``gf_spmat`` : create a |mlab_sm| object.
* ``gf_spmat_get`` : perform computations with the |mlab_sm|.
* ``gf_spmat_set`` : modify the |mlab_sm|.
* ``gf_precond`` : create a |mlab_pc| object.
* ``gf_precond_get`` : perform computations with the |mlab_pc|.
* ``gf_linsolve`` : interface to various linear solvers provided by getfem
  (|sLU|, conjugated gradient, etc.).
* ``gf_asm`` : assembly routines.
* ``gf_solve`` : various solvers for usual PDEs (obsoleted by the |mlab_mbr|
  objects).
* ``gf_compute`` : computations involving the solution of a PDE (norm,
  derivative, etc.).
* ``gf_mdbrick`` : create a ("model brick") |mlab_mbr| object.
* ``gf_mdbrick_get`` : retrieve information from a |mlab_mbr| object.
* ``gf_mdbrick_set`` : modify a |mlab_mbr| object.
* ``gf_mdstate`` : create a ("model state") |mlab_ms| object.
* ``gf_mdstate_get`` : retrieve information from a |mlab_ms| object.
* ``gf_mdstate_set`` : modify a |mlab_ms| object.
* ``gf_model`` : create a |mlab_md| object.
* ``gf_model_get`` : retrieve information from a |mlab_md| object.
* ``gf_model_set`` : modify a |mlab_md| object.
* ``gf_mumps_context`` : create a |mlab_mumps| object.
* ``gf_mumps_context_get`` : retrieve information from a |mlab_mumps| object.
* ``gf_mumps_context_set`` : modify a |mlab_mumps| object.
* ``gf_global_function`` : create a gfGlobalFunction object.
* ``gf_model_get`` : retrieve information from a gfGlobalFunction object.
* ``gf_model_set`` : modify a GlobalFunction object.
* ``gf_plot_mesh`` : plotting of mesh.
* ``gf_plot`` : plotting of 2D and 3D fields.
* ``gf_plot_1D`` : plotting of 1D fields.
* ``gf_plot_slice`` : plotting of a mesh slice.


Objects
-------

Various "objects" can be manipulated by the |gfm| toolbox, see fig.
:ref:`malb-fig-hierarchy`. The MESH and MESHFEM objects are the two most
important objects.

.. _malb-fig-hierarchy:
.. figure:: images/hierarchy.png
   :align: center

   |gfm| objects hierarchy.

* :envvar:`gfGeoTrans`: geometric transformations (defines the shape/position of
  the convexes), created with ``gf_geotrans``
* :envvar:`gfGlobalFunction`: represent a global function for the enrichment of finite element methods.
* :envvar:`gfMesh` : mesh structure (nodes, convexes, geometric transformations for
  each convex), created with ``gf_mesh``
* :envvar:`gfInteg` : integration method (exact, quadrature formula...).  Although
  not linked directly to GEOTRANS, an integration method is usually specific to a
  given convex structure. Created with ``gf_integ``
* :envvar:`gfFem` : the finite element method (one per convex, can be PK, QK,
  HERMITE, etc.). Created with ``gf_fem``
* :envvar:`gfCvStruct` : stores formal information convex structures (nb. of points,
  nb. of faces which are themselves convex structures).
* :envvar:`gfMeshFem` : object linked to a mesh, where each convex has been assigned
  an FEM. Created with ``gf_mesh_fem``.
* :envvar:`gfMeshImM` : object linked to a mesh, where each convex has been assigned
  an integration method. Created with ``gf_mesh_im``.
* :envvar:`gfMeshSlice` : object linked to a mesh, very similar to a
  P1-discontinuous |mlab_mf|. Used for fast interpolation and plotting.
* :envvar:`gfMdBrick` : |mlab_mbr| , an abstraction of a part of solver (for
  example, the part which build the tangent matrix, the part which handles the
  dirichlet conditions, etc.). These objects are stacked to build a complete
  solver for a wide variety of problems. They typically use a number of
  |mlab_mf|, |mlab_mim| etc. Deprecated object, replaced now by gfModel.
* :envvar:`gfMdState` : "model state", holds the global data for a stack of mdbricks
  (global tangent matrix, right hand side etc.). Deprecated object, replaced now by gfModel.
* :envvar:`gfModel` : "model", holds the global data, variables and description of a
  model. Evolution of "model state" object for 4.0 version of |gf|.

The |gfm| toolbox uses its own :envvar:`memory management`. Hence |gf| objects
are not cleared when a::

  >> clear all

is issued at the |octv| or |mlab| prompt, but instead the function::

  >> gf_workspace('clear all')

should be used. The various |gfm| object can be accessed via *handles* (or
*descriptors*), which are just |octv| / |mlab| structures containing 32-bits integer
identifiers to the real objects. Hence the |octv| or |mlab| command::

  >> whos

does not report the memory consumption of |gf| objects (except the marginal space
used by the handle). Instead, you should use::

  >> gf_workspace('stats')

There are two kinds of |gfm| objects:

* static ones, which can not be deleted: ELTM, FEM, INTEG, GEOTRANS and CVSTRUCT.
  Hopefully their memory consumption is very low.
* dynamic ones, which can be destroyed, and are handled by the ``gf_workspace``
  function: MESH, MESHFEM, MESHIM, SLICE, SPMAT, PRECOND.

The objects MESH and MESHFEM are not independent: a MESHFEM object is always
linked to a MESH object, and a MESH object can be used by several MESHFEM
objects. Hence when you request the destruction of a MESH object, its destruction
might be delayed until it is not used anymore by any MESHFEM (these objects
waiting for deletion are listed in the *anonymous workspace* section of
``gf_workspace('stats')``).
