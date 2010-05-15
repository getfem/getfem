.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-xfem:

Level-sets, Xfem, fictitious domains
====================================

|gf| offers (since v2.0) a certain number of functionalities concerning
level-sets, support for Xfem and fictitious domain methods and discontinuous field
across a level-set.

.. important::

   All the tools listed below needs the package `qhull <http://www.qhull.org>`_
   installed on your system. This package is widely available. It computes convex
   hull and delaunay triangulations in arbitrary dimension. Everything here is
   considered "work in progress", it is still subject to major changes if needed.

The program :file:`tests/crack.cc` is a good example of use of these tools.


Representation of level-sets
----------------------------

|gf| deals with level-set defined by piecewise polynomial function on a mesh. It
will be defined as the zero of this function. In the file
:file:`getfem/getfem_levelset.h` a level-set is represented by a function defined
on a lagrange fem of a certain degree on a mesh. The constructor to define a new
|gf_ls| is the following::

  getfem::level_set ls(mesh, degree = 1, with_secondary = false);

where ``mesh`` is a valid mesh of type |gf_m|, ``degree`` is the degree of the
polynomials (1 is the default value), and ``with_secondary`` is a boolean whose
default value is false. The secondary level-set is used to represent fractures (if
:math:`p(x)` is the primary levelset function and :math:`s(x)` is the secondary
levelset function, the crack is defined by :math:`p(x) = 0` and :math:`s(x) \leq
0`: the role of the secondary is to stop the crack).

Each level-set function is defined by a |mf| ``mf`` and the dof values over this
|mf|, in a vector. The object |gf_ls| contains a |mf| and the vectors of dof for
the corresponding function(s). The method ``ls.value(0)`` returns the vector of
dof for the primary level-set function, so that these values can be set. The
method ``ls.value(1)`` returns the dof vector for the secondary level-set function
if any. The method ``ls.get_mesh_fem()`` returns a reference on the |gf_mf|
object.


Mesh cut by level-sets
----------------------

In order to compute adapted integration methods and finite element methods to
represent a field which is discontinuous across a level-set, a certain number of
pre-computations have to be done at the mesh level. The file
:file:`getfem/getfem_mesh_level_set.h` defines the object |gf_mls| which handles
these pre-computations. The constructor of this object is the following::

  getfem::mesh_level_set mls(mesh);

where ``mesh`` is a valid mesh of type |gf_m|. In order to indicate that the mesh
is cut by a level-set, one has to call the method ``mls.add_level_set(ls)``, where
``ls`` is an object of type |gf_ls|. An arbitrary number of level-sets can be
added. To initialize the object or to actualize it when the value of the level-set
function is modified, one has to call the method ``mls.adapt()``.

In particular a subdivision of each element cut by the level-set is made with
simplices.


Adapted integration methods
---------------------------

For fields which are discontinuous across a level-set, integration methods have
to be adapted. The object |gf_mimls| defined in the file
:file:`getfem/getfem_mesh_im_level_set.h` defines a composite integration method
for the elements cut by the level-set. The constructor of this object is the
following::

  getfem::mesh_im_level_set mim(mls, where, regular_im = 0, singular_im = 0);

where ``mls`` is an object of type |gf_mls|, ``where`` is an enum for which
possible values are

* ``getfem::mesh_im_level_set::INTEGRATE_INSIDE`` (integrate over :math:`p(x)<0`),

* ``getfem::mesh_im_level_set::INTEGRATE_OUTSIDE`` (integrate over :math:`p(x)>0`),

* ``getfem::mesh_im_level_set::INTEGRATE_ALL``,

* ``getfem::mesh_im_level_set::INTEGRATE_BOUNDARY`` (integrate over :math:`p(x)=0`
  and :math:`s(x)\leq 0`)

The argument ``regular_im`` should be of type ``pintegration_method``, and will be
the integration method applied on each sub-simplex of the composite integration
for convexes cut by the levelset. The optional ``singular_im`` should be also of
type ``pintegration_method`` and is used for crack singular functions: it is
applied to sub-simplices which share a vertex with the crack tip (the specific
integration method ``IM_QUASI_POLAR(..)`` is well suited for this purpose).

The object |gf_mimls| can be used as a classical |gf_mim| object (for instance the
method ``mim.set_integration_method(...)`` allows to set the integration methods
for the elements which are not cut by the level-set).

To initialize the object or to actualize it when the value of the level-set
function is modified, one has to call the method ``mim.adapt()``.


Discontinuous field across some level-sets
------------------------------------------

The object |gf_mfls| is defined in the file
:file:`getfem/getfem_mesh_fem_level_set.h`. It is derived from |gf_mf| object
and can be used in the same way. It defines a finite element method with
discontinuity across the level-sets (it can deal with an arbitrary number of
level-sets). The constructor is the following::

  getfem::mesh_fem_level_set mfls(m, mf);

where ``m`` is a valid mesh of type |gf_m| and ``mf`` is the an object of type
|gf_mf| which defines the finite element method used for elements which are not
cut by the level-sets.

To initialize the object or to actualize it when the value of the level-set
function is modified, one has to call the method ``mfls.adapt()``.

To represent discontinuous fields, the finite element method is enriched with
discontinuous functions which are the product of a Heaviside function by the base
functions of the finite element method represented by ``mf`` (see [Xfem]_ for
more details).


Fictitious domain approach with Xfem
------------------------------------

An example of a Poisson problem with a Dirichlet condition posed on a boundary
independant of the mesh is present on the ``tests`` directory of the distribution.

See :file:`contrib/xfem_contact/xfem_dirichlet.cc` file.
