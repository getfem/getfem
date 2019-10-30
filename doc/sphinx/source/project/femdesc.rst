.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-femdesc:

The FEM description in |gf|
===========================

The aim of this section is to briefly introduce the FEM description in |gf|
mainly in order to fix the notation used in the rest of the document (definition
of element, reference element, geometric transformation, gradient of the
geometric transformation ...).


Convex structures
-----------------

Finite element methods are defined on small convex domains called elements. The
simplest element on which a finite element method can be defined is a segment
(simplex of dimension 1), other possibilities are triangles, tetrahedrons
(simplices of dimension 2 and 3), prisms, parallelepiped, etc. In |gf|, a type of
element (for us, a convex) is described by the object |bg_cs| defined in the file
:file:`bgeot_convex_structure.h`.

It describes only the structure of the convex not the coordinates of the
vertices. This structure is not to be manipulated by itself, because it is not
necessary that more than one structure of this type describe the same type of
convex. What will be manipulated is a pointer on such a descriptor which has to
be declared with the type |bg_pcs|

The following functions give a pointer onto the descriptor of the usual type of
elements:

.. c:function:: bgeot::simplex_structure(dim_type d)

   description of a simplex of dimension ``d``.

.. c:function:: bgeot::parallelepiped_structure(dim_type d)

   description of a parallelepiped of dimension ``d``.

.. c:function:: bgeot::convex_product_structure(bgeot::pconvex_structure p1, bgeot::pconv$

   description of the direct product of ``p1`` and ``p2``.

.. c:function:: bgeot::prism_P1_structure(dim_type d)

   description of a prism of dimension ``d``.

For instance if one needs the description of a square, one can call
equivalently::

  p = bgeot::parallelepiped_structure(2);

or::

 p = bgeot::convex_product_structure(bgeot::simplex_structure(1),
                                     bgeot::simplex_structure(1));

The descriptor contains in particular the number of faces (``p->nb_faces()``),
the dimension of the convex (``p->dim()``), for the number of vertices
(``p->nb_points()``). Other information is the number of vertices of each face,
the description of a face and the eventual reference to a more basic description
(used for the description of geometric transformations).

.. _dp-fig-elem:
.. figure:: images/getfemelemelem.png
   :align: center
   :scale: 60

   usual elements


Convexes of reference
---------------------

A convex of reference is a particular real element, i.e. a structure of convex
with a list of vertices. It describes the particular element from which a finite
element method is defined. In the file :file:`bgeot_convex_ref.h` the object
|bg_cr| makes this description. The library keeps only one description for each
type of convex. So what will be manipulated is a pointer of type |bg_pcr| on the
descriptor.

The following functions build the descriptions:

.. c:function:: bgeot::simplex_of_reference(dim_type d)

   description of the simplex of reference of dimension ``d``.

.. c:function:: bgeot::simplex_of_reference(dim_type d, short_type k)

   description of the simplex of reference of dimension ``d`` with degree ``k``
   Lagrange grid.

.. c:function:: bgeot::convex_ref_product(pconvex_ref a, pconvex_ref b)

   description of the direct product of two convexes of reference.

.. c:function:: bgeot::parallelepiped_of_reference(dim_type d)

   description of the parallelepiped of reference of dimension ``d``.

The vertices correspond to the classical vertices for such reference element. For
instance the vertices for the triangle are :math:`(0, 0)`, :math:`(1, 0)` and
:math:`(0, 1)`. It corresponds to the configuration shown in Figure
:ref:`dp-fig-elem`

If ``p`` is of type |bg_pcr| then ``p->structure()`` is the corresponding convex
structure. Thus for instance ``p->structure()->nb_points()`` gives the number of
vertices. The function ``p->points()`` give the array of vertices and
``p->points()[0]`` is the first vertex. The function ``p->is_in(const base_node
&pt)`` return a real which is negative or null if the point ``pt`` is in the
element. The function ``p->is_in_face(short_type f, const base_node &pt)`` return
a real which is null if the point ``pt`` is in the face ``f`` of the element.
Other functions can be found in :file:`bgeot_convex_ref.h` and
:file:`bgeot_convex.h`.


Shape function type
-------------------

Most of the time the shape functions of finite element methods are polynomials,
at least on the convex of reference. But, the possibility is given to have other
types of elements. It is possible to define other kind of base functions such as
piecewise polynomials, interpolant wavelets, etc.

To be used by the finite element description, a shape function type must be able
to be evaluated on a point (``a = F.eval(pt)``, where ``pt`` is a ``base_node``)
and must have a method to compute the derivative with respect to the ith variable
(``F.derivative(i)``).

For the moment, only polynomials and piecewise polynomials are defined in the
files :file:`bgeot_poly.h` and :file:`bgeot_poly_composite.h`.

.. _dp-transgeo:

Geometric transformations
-------------------------

.. _dp-fig-transgeo:
.. figure:: images/getfemtransgeo.png
   :align: center
   :scale: 60

   geometric transformation

A geometric transformation is a polynomial application:

.. math::

   \tau : \widehat{T} \subset \Reel^P \longrightarrow T \subset \Reel^N,

which maps the reference element :math:`\widehat{T}` to the real element :math:`T`. The
geometric nodes are denoted:

.. math::

   g^i, i = 0, \ldots, n_g - 1.

The geometric transformation is described thanks to a :math:`n_g` components
polynomial vector (In fact, as an extention, non polynomial geometric
transformation can also be supported by |gf|, but this is very rarely used)

.. math::

   {\cal N}(\widehat{x}),

such that

.. math::

  \tau(\widehat{x}) = \sum_{i = 0}^{n_g - 1}{\cal N}_i(\widehat{x}) g^i.

Denoting

.. math::

   G = (g^0; g^1; ...; g^{n_g - 1}),

the :math:`N\times n_g` matrix containing of all the geometric nodes, one has

.. math::

   \fbox{$\tau(\widehat{x}) = G\cdot{\cal N}(\widehat{x})$.}

The derivative of :math:`\tau` is then

.. math::

   \fbox{$K(\widehat{x}) := \nabla\tau(\widehat{x}) = G\cdot\nabla {\cal N}(\widehat{x})$,}

where :math:`K(\widehat{x}) = \nabla\tau(\widehat{x})` is a :math:`N\times P` matrix and
:math:`\nabla {\cal N}(\widehat{x})` a :math:`n_g\times P` matrix. The (transposed)
pseudo-inverse of :math:`\nabla\tau(\widehat{x})` is a :math:`N\times P` matrix denoted
:math:`B(\widehat{x})`:

.. math::

   \fbox{$B(\widehat{x}) := K(\widehat{x})(K(\widehat{x})^T K(\widehat{x}))^{-1}$,}

Of course, when :math:`P=N`, one has :math:`B(\widehat{x})=K(\widehat{x})^{-T}`.

Pointers on a descriptor of a geometric transformation can be obtained by the
following function defined in the file :file:`bgeot_geometric_trans.h`::

  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor("name of trans");

where ``"name of trans"`` can be chosen among the following list.

* ``"GT_PK(n,k)"``

  Description of the simplex transformation of dimension ``n`` and degree ``k``
  (Most of the time, the degree 1 is used).

* ``"GT_QK(n,k)"``

  Description of the parallelepiped transformation of dimension ``n`` and degree
  ``k``.

* ``"GT_PRISM(n,k)"``

  Description of the prism transformation of dimension ``n`` and degree ``k``.

* ``"GT_PRODUCT(a,b)"``

  Description of the direct product of the two transformations ``a`` and ``b``.

* ``"GT_LINEAR_PRODUCT(a,b)"``

  Description of the direct product of the two transformations ``a`` and ``b``
  keeping a linear transformation (this is a restriction of the previous
  function). This allows, for instance, to use exact integrations on regular
  meshes with parallelograms.


Finite element methods description
----------------------------------

A finite element method is defined on a reference element
:math:`\widehat{T}\subset\Reel^P` by a set of :math:`n_d` nodes :math:`a^i` and
corresponding base functions

.. math::

   (\widehat{\varphi})^i : \widehat{T}\subset\Reel^P \longrightarrow \Reel^Q

Denoting

.. math::

   \psi^i(x) = (\widehat{\varphi})^i(\widehat{x}) = (\widehat{\varphi})^i(\tau^{-1}(x)),

a supplementary linear transformation is allowed for the real base function

.. math::

   \varphi^i(x) = \sum_{j = 0}^{n_d - 1} M_{ij} \psi^j(x),

where :math:`M` is a :math:`n_d \times n_d` matrix possibly depending on the
geometric transformation (i.e. on the real element). For basic elements as
Lagrange elements this matrix is the identity matrix (it is simply ignored). In
this case, we will say that the element is :math:`\tau`-equivalent.

This approach allows to define hermite elements (Argyris for instance) in a
generic way, even with non linear transformations (i.e. mainly for curved
boundaries). We denote :math:`[\widehat{\varphi}(\widehat{x})]` the :math:`n_d \times Q` matrix
whose ith line is :math:`(\widehat{\varphi})^i(\widehat{x})`. Whis this notation, for a function is
defined by

.. math::

   f(x) = \sum_{i = 0}^{n_d - 1} \alpha_i \varphi^i(x),

one has

.. math::

   \fbox{$f(\tau(\widehat{x})) = \alpha^T M [\widehat{\varphi}(\widehat{x})]$,}

where :math:`\alpha` is the vector whose ith component is :math:`\alpha_i`.

A certain number of description of classical finite element method are defined in
the file :file:`getfem_fem.h`. See :ref:`ud-appendixa` for an exhaustive list of
available finite element methods.

A pointer to the finite element descriptor of a method is obtained using the
function::

  getfem::pfem pfe = getfem::fem_descriptor("name of method");

We refer to the file :file:`getfem_fem.cc` for how to define a new finite element
method.
