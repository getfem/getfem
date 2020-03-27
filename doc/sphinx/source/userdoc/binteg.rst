.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-binteg:

Selecting integration methods
=============================

The description of an integration method on a whole mesh is done thanks to the
structure |gf_mim|, defined in the file :file:`getfem/getfem_mesh_im.h`.
Basically, this structure describes the integration method on each element of the
mesh. One can instantiate a |gf_mim| object as follows::

  getfem::mesh_im mim(mymesh);

where ``mymesh`` is an already existing mesh. The structure will be linked to this
mesh and will react when modifications will be done on it (for example when the
mesh is refined, the integration method will be also refined).

It is possible to specify element by element the integration method, so that
element of mixed types can be treated, even if the dimensions are different.

To select a particular integration method on a given element, one can use::

  mim.set_integration_method(i, ppi);

where ``i`` is the index of the element and ``ppi`` is the descriptor of the
integration method. Alternative forms of this member function are::

  void mesh_im::set_integration_method(const dal::bit_vector &cvs,
                                        getfem::pintegration_method ppi);
  void mesh_im::set_integration_method(getfem::pintegration_method ppi);

which set the integration method for either the convexes listed in the |bv| cvs,
or all the convexes of the mesh.

The list of all available descriptors of integration methods is in the file
:file:`getfem/getfem_integration.h`. Descriptors for integration methods are
available thanks to the following function::

  getfem::pintegration_method ppi = getfem::int_method_descriptor("name of method");

where ``"name of method"`` is to be chosen among the existing methods. A name of a
method can be retrieved with::

  std::string im_name = getfem::name_of_int_method(ppi);

A non exhaustive list (see :ref:`ud-appendixb` or
:file:`getfem/getfem_integration.h` for exhaustive lists) of integration methods
is given below.

Examples of exact integration methods:

* ``"IM_NONE()"``:
  Dummy integration method (new in getfem++-1.7).

* ``"IM_EXACT_SIMPLEX(n)"``:
  Description of the exact integration of polynomials on the simplex of reference
  of dimension ``n``.

* ``"IM_PRODUCT(a, b)"``:
  Description of the exact integration on the convex which is the direct product
  of the convex in ``a`` and in ``b``.

* ``"IM_EXACT_PARALLELEPIPED(n)"``:
  Description of the exact integration of polynomials on the parallelepiped of
  reference of dimension ``n``

* ``"IM_EXACT_PRISM(n)"``:
  Description of the exact integration of polynomials on the prism of reference of
  dimension ``n``

Examples of approximated integration methods:

* ``"IM_GAUSS1D(k)"``:
  Description of the Gauss integration on a segment of order ``k``. Available for
  all odd values of ``k <= 99``.

* ``"IM_NC(n,k)"``:
  Description of the integration on a simplex of reference of dimension ``n`` for
  polynomials of degree ``k`` with the Newton Cotes method (based on Lagrange
  interpolation).

* ``"IM_PRODUCT(a,b)"``:
  Build a method doing the direct product of methods ``a`` and ``b``.

* ``"IM_TRIANGLE(2)"``:
  Integration on a triangle of order 2 with 3 points.

* ``"IM_TRIANGLE(7)"``:
  Integration on a triangle of order 7 with 13 points.

* ``"IM_TRIANGLE(19)"``:
  Integration on a triangle of order 19 with 73 points.

* ``"IM_QUAD(2)"``:
  Integration on quadrilaterals of order 2 with 3 points.

* ``"IM_GAUSS_PARALLELEPIPED(2,3)"``:
  Integration on quadrilaterals of order 3 with 4 points (shortcut for
  ``"IM_PRODUCT(IM_GAUSS1D(3),IM_GAUSS1D(3))"``).

* ``"IM_TETRAHEDRON(5)"``:
  Integration on a tetrahedron of order 5 with 15 points.

.. note::

    Note that ``"IM_QUAD(3)"`` is not able to integrate exactly the base functions
    of the ``"FEM_QK(2,3)"`` finite element! Since its base function are tensorial
    product of 1D polynomials of degree 3, one would need to use ``"IM_QUAD(7)"``
    (6 is not available). Hence ``"IM_GAUSS_PARALLELEPIPED(2,k)"`` should always
    be preferred over ``"IM_QUAD(2*k)"`` since it has less integration points.

An alternative way to obtain integration methods::

  getfem::pintegration_method ppi =
    getfem::classical_exact_im(bgeot::pgeometric_trans pgt);

  getfem::pintegration_method ppi =
    getfem::classical_approx_im(bgeot::pgeometric_trans pgt, dim_type d);

These functions return an exact (i.e. analytical) integration method, or select an
approximate integration method which is able to integrate exactly polynomials of
degree <= ``d`` (at least) for convexes defined with the specified geometric
transformation.


Methods of the |mim| object
---------------------------

Once an integration method is defined on a mesh, it is possible to obtain
information on it with the following methods (the list is not exhaustive).

.. function:: mim.convex_index()

   Set of indexes (a |dal_bv|) on which an integration method is defined.

.. function:: mim.linked_mesh()

   Gives a reference to the linked mesh.

.. function:: mim.int_method_of_element(i)

   Gives a descriptor on the integration method defined on element of index ``i``.

.. function:: mim.clear()

   Clear the structure. There are no further integration method defined on the
   mesh.
