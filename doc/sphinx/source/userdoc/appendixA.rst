.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. |nbsp| unicode:: U+00A0 .. non-breaking space

.. _ud-appendixa:

Appendix A. Finite element method list
======================================

  .. list-table:: Symbols representing degree of freedom types
     :widths: 30 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlistsymbols00.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols01.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols02.png
            :align: center
            :scale: 50
     * - Value of the function at the node.
       - Value of the gradient along of the first coordinate.
       - Value of the gradient along of the second coordinate.
     * - .. image:: images/getfemlistsymbols03.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols04.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols05.png
            :align: center
            :scale: 50
     * - Value of the gradient along of the third coordinate for 3D elements.
       - Value of the whole gradient at the node.
       - Value of the normal derivative to a face.
     * - .. image:: images/getfemlistsymbols06.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols07.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols08.png
            :align: center
            :scale: 50
     * - Value of the second derivative along the first coordinate (twice).
       - Value of the second derivative along the second coordinate (twice).
       - Value of the second cross derivative in 2D or second derivative
         along the third coordinate (twice) in 3D.
     * - .. image:: images/getfemlistsymbols09.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols10.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols11.png
            :align: center
            :scale: 50
     * - Value of the whole second derivative (hessian) at the node.
       - Scalar product with a certain vector (for instance an edge) for a
         vector elements.
       - Scalar product with the normal to a face for a vector elements.
     * - .. image:: images/getfemlistsymbols12.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistsymbols13.png
            :align: center
            :scale: 50
       -
     * - Bubble function on an element or a face, to be specified.
       - Lagrange hierarchical d.o.f. value at the node in a space of details.
       -

Let us recall that all finite element methods defined in |gf| are declared in the
file ``getfem_fem.h`` and that a descriptor on a finite element method is obtained
thanks to the function::

  getfem::pfem pf = getfem::fem_descriptor("name of method");

where ``"name of method"`` is a string to be choosen among the existing methods.


Classical :math:`P_K` Lagrange elements on simplices
----------------------------------------------------

.. _ud-fig-segmentpk:
.. figure:: images/getfemlistsegmentPk.png
   :align: center
   :scale: 60

   Examples of classical :math:`P_K` Lagrange elements on a segment

It is possible to define a classical :math:`P_K` Lagrange element of arbitrary
dimension and arbitrary degree. Each degree of freedom of such an element
corresponds to the value of the function on a corresponding node. The grid of
node is the so-called Lagrange grid. Figures :ref:`ud-fig-segmentpk`.

  .. _ud-fig-trianglepk:
  .. list-table:: Examples of classical :math:`P_K` Lagrange elements on a triangle.
     :widths: 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlisttriangleP1.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttriangleP2.png
            :align: center
            :scale: 50
     * - :math:`P_1`, 3 d.o.f., :math:`C^0`
       - :math:`P_2` element, 6 d.o.f., :math:`C^0`
     * - .. image:: images/getfemlisttriangleP3.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttriangleP6.png
            :align: center
            :scale: 50
     * - :math:`P_3`, 10 d.o.f., :math:`C^0`
       - :math:`P_6` element, 28 d.o.f., :math:`C^0`

The number of degrees of freedom for a classical :math:`P_K` Lagrange element of
dimension :math:`P` and degree :math:`K` is :math:`\dfrac{(P+K)!}{P!K!}`. For
instance, in dimension 2 :math:`(P = 2)`, this value is :math:`\dfrac{(K+1)
(K+2)}{2}` and in dimension 3 :math:`(P = 3)`, it is :math:`\dfrac{(K+1) (K+2)
(K+3)}{6}`.

  .. _ud-fig-tetrahedronpk:
  .. list-table:: Examples of classical :math:`P_K` Lagrange elements on a tetrahedron.
     :widths: 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlisttetrahedronP1.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttetrahedronP2.png
            :align: center
            :scale: 50
     * - :math:`P_1` element, 4 d.o.f., :math:`C^0`
       - :math:`P_2` element, 10 d.o.f., :math:`C^0`
     * - .. image:: images/getfemlisttetrahedronP4.png
            :align: center
            :scale: 50
       -
     * - :math:`P_4` element, 35 d.o.f., :math:`C^0`
       -

The particular way used in |gf| to numerate the nodes are also shown in figures
:ref:`segment<ud-fig-segmentpk>`, :ref:`triangle<ud-fig-trianglepk>` and
:ref:`tetrahedron<ud-fig-tetrahedronpk>`. Using another numeration, let

.. math::

  i_0, i_1, ... i_P,

be some indices such that

.. math::

  0 \leq i_0, i_1, ... i_P \leq K, \ \mbox{ and } \ \sum_{n = 0}^{P} i_n = K.

Then, the coordinate of a node can be computed as

.. math::

   a_{i_0, i_1, ... i_P} = \sum_{n = 0}^{P} \dfrac{i_n}{K}S_n, \ \ \mbox{ for } K \neq 0,

where :math:`S_0, S_1, ... S_N` are the vertices of the simplex (for :math:`K = 0`
the particular choice :math:`a_{0, 0, ... 0} =  \sum_{n = 0}^{P}
\dfrac{1}{P+1}S_n` has been chosen). Then each base function, corresponding of each
node :math:`a_{i_0, i_1, ... i_P}` is defined by

.. math::

  \phi_{i_0, i_1, ... i_P} = \prod_{n = 0}^{P} \prod_{j=0}^{i_n-1} \left(\dfrac{K \lambda_n - j}{j+1}\right).

where :math:`\lambda_n` are the barycentric coordinates, i.e. the polynomials of
degree 1 whose value is :math:`1` on the vertex :math:`S_n` and whose value is
:math:`0` on other vertices. On the reference element, one has

.. math::

  \lambda_n = x_n, \ \ 0 \leq n < P,


.. math::

  \lambda_P = 1 - x_0 - x_1 - ... - x_{P-1}.

When between two elements of the same degrees (even with different dimensions),
the d.o.f. of a common face are linked, the element is of class :math:`C^0`. This
means that the global polynomial is continuous. If you try to link elements of
different degrees, you will get some trouble with the unlinked d.o.f. This is not
automatically supported by |gf|, so you will have to support it (add constraints
on these d.o.f.).

For some applications (computation of a gradient for instance) one may not want
the d.o.f. of a common face to be linked. This is why there are two versions of
the classical :math:`P_K` Lagrange element.

  .. list-table:: Classical :math:`P_K` Lagrange element ``"FEM_PK(P, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K \leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`\dfrac{(K+P)!}{K! P!}`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

  .. list-table:: Discontinuous :math:`P_K` Lagrange element ``"FEM_PK_DISCONTINUOUS(P, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K \leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`\dfrac{(K+P)!}{K! P!}`
       - discontinuous
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`         

  .. list-table:: Discontinuous :math:`P_K` Lagrange element with internal dofs ``"FEM_PK_DISCONTINUOUS(P, K, alpha)"``. The method ``"FEM_PK_DISCONTINUOUS(P, K, 0)"`` is identical to ``"FEM_PK_DISCONTINUOUS(P, K)"``. For alpha > 0, ``"FEM_PK_DISCONTINUOUS(P, K, alpha)"`` corresponds to a Lagrange method with all finite element nodes in the interior of the domain located at the position :math:`(\mbox{alpha})g + (1-\mbox{alpha})a_i` for :math:`g` the centroid of the element and :math:`a_i` the node of the standard :math:`P_K` method.
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K \leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`\dfrac{(K+P)!}{K! P!}`
       - discontinuous
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes         

Even though Lagrange elements are defined for arbitrary degrees, choosing a
high degree can be problematic for a large number of applications due to
the "noisy" characteristic of the lagrange basis. These elements are
recommended for the basic interpolation but for p.d.e. applications elements
with hierarchical basis are preferable (see the corresponding section).

Classical Lagrange elements on other geometries
-----------------------------------------------

Classical Lagrange elements on parallelepipeds or prisms are obtained as tensor
product of Lagrange elements on simplices. When two elements are defined, one on a
dimension :math:`P^1` and the other in dimension :math:`P^2`, one obtains the base
functions of the tensorial product (on the reference element) as

.. math::

  \widehat{\varphi}_{ij}(x,y) = \widehat{\varphi}^1_i(x) \widehat{\varphi}^2_j(y), ~~ x \in \rm I\hspace{-0.15em}R^{P^1}, y \in  \rm I\hspace{-0.15em}R^{P^2},

where :math:`\widehat{\varphi}^1_i` and :math:`\widehat{\varphi}^2_i` are respectively the base functions
of the first and second element.

  .. _ud-fig-prodpkdeux:
  .. list-table:: Examples of classical :math:`Q_K` Lagrange elements in dimension 2.
     :widths: 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlistquadQ1.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistquadQ3.png
            :align: center
            :scale: 50
     * - :math:`Q_1` element, 4 d.o.f., :math:`C^0`
       - :math:`Q_3` element, 16 d.o.f., :math:`C^0`

The :math:`Q_K` element on a parallelepiped of dimension :math:`P` is obtained as
the tensorial product of :math:`P` classical :math:`P_K` elements on the segment.
Examples in dimension 2 are shown in figure :ref:`dimension 2<ud-fig-prodpkdeux>`
and in dimension 3 in figure :ref:`dimension 3<ud-fig-prodpktrois>`.

A prism in dimension :math:`P > 1` is the direct product of a simplex of dimension
:math:`P-1` with a segment. The :math:`P_K \otimes P_K` element on this prism is
the tensorial product of the classical :math:`P_K` element on a simplex of
dimension :math:`P-1` with the classical :math:`P_K` element on a segment. For
:math:`P=2` this coincide with a parallelepiped. Examples in dimension :math:`3`
are shown in figure :ref:`dimension 3<ud-fig-prodpktrois>`. This is also possible
not to have the same degree on each dimension. An example is shown on figure
:ref:`dimension 3, prism<ud-fig-prism_P2_p1>`.

  .. _ud-fig-prodpktrois:
  .. list-table:: Examples of classical Lagrange elements in dimension 3.
     :widths: 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlistcubeQ1.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistcubeQ3.png
            :align: center
            :scale: 50
     * - :math:`Q_1` element, 8 d.o.f., :math:`C^0`
       - :math:`Q_3` element, 64 d.o.f., :math:`C^0`
     * - .. image:: images/getfemlistprismP1.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlistprismP3.png
            :align: center
            :scale: 50
     * - :math:`P_1 \otimes P_1` element, 6 d.o.f., :math:`C^0`
       - :math:`P_3 \otimes P_3` element, 40 d.o.f., :math:`C^0`

:math:`.\\`

.. _ud-fig-prism_P2_p1:
.. figure:: images/getfemlistprismP2P1.png
   :align: center
   :scale: 60

   :math:`P_2 \otimes P_1` Lagrange element on a prism, 12 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: . :math:`Q_K` Lagrange element on parallelepipeds ``"FEM_QK(P, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`KP`, :math:`0 \leq K \leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`(K+1)^P`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

  .. list-table:: . :math:`P_K \otimes P_K` Lagrange element on prisms ``"FEM_PK_PRISM(P, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`2K`, :math:`0 \leq K \leq 255`
       - :math:`P`, :math:`~ 2 \leq P \leq 255`
       - :math:`(K+1)` :math:`\times~\dfrac{(K+P-1)!}{K! (P-1)!}`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

  .. list-table:: . :math:`P_{K_1} \otimes P_{K_2}` Lagrange element on prisms ``"FEM_PRODUCT(FEM_PK(P-1, K1), FEM_PK(1, K2))"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K_1+K_2`, :math:`0 \leq K_1,K_2 \leq 255`
       - :math:`P`, :math:`~ 2 \leq P \leq 255`
       - :math:`(K_2+1)` :math:`\times~\dfrac{(K_1+P-1)!}{K_1! (P-1)!}`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

.. figure:: images/getfemlistincomplete.png
   :align: center
   :scale: 60

   Incomplete :math:`Q_2` elements in dimension two and three, 8 or 20 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: Incomplete :math:`Q_2` Lagrange element on parallelepipeds (Quad 8 and Hexa 20 serendipity elements) ``"FEM_Q2_INCOMPLETE(P)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - 3
       - :math:`P`, :math:`~ 2 \leq P \leq 3`
       - :math:`8\ \text{for}\ P = 2~~~~~` :math:`20\ \text{for}\ P = 3`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes


Elements with hierarchical basis
--------------------------------

The idea behind hierarchical basis is the description of the solution at different
level: a rough level, a more refined level ... In the same discretization some
degrees of freedom represent the rough description, some other the more rafined
and so on. This corresponds to imbricated spaces of discretization. The
hierarchical basis contains a basis of each of these spaces (this is not the case
in classical Lagrange elements when the mesh is refined).

Among the advantages, the condition number of rigidity matrices can be greatly
improved, it allows local raffinement and a resolution with a multigrid approach.

Hierarchical elements with respect to the degree
+++++++++++++++++++++++++++++++++++++++++++++++++

.. _ud-fig-seg_hier:
.. figure:: images/getfemlistsegmenthier.png
   :align: center
   :scale: 60

   :math:`P_K` Hierarchical element on a segment, :math:`C^0`

:math:`.\\`

  .. list-table:: . :math:`P_{K}` Classical Lagrange element on simplices but with a hierarchical basis with respect to the degree ``"FEM_PK_HIERARCHICAL(P,K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K\leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`\dfrac{(K+P)!}{K! P!}`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

  .. list-table:: . :math:`Q_{K}` Classical Lagrange element on parallelepipeds but with a hierarchical basis with respect to the degree ``"FEM_QK_HIERARCHICAL(P,K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K\leq 255`
       - :math:`P`, :math:`~ 1 \leq P \leq 255`
       - :math:`(K+1)^P`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

:math:`.\\`

  .. list-table:: . :math:`P_{K}` Classical Lagrange element on prisms but with a hierarchical basis with respect to the degree ``"FEM_PK_PRISM_HIERARCHICAL(P,K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`, :math:`0 \leq K\leq 255`
       - :math:`P`, :math:`~ 2 \leq P \leq 255`
       - :math:`(K+1)` :math:`\times~\dfrac{(K+P-1)!}{K! (P-1)!}`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - Yes

some particular choices: :math:`P_4` will be built with the basis of the
:math:`P_1`, the additional basis of the :math:`P_2` then the additional basis of the :math:`P_4`.

:math:`P_6` will be built with the basis of the :math:`P_1`, the additional basis
:of the :math:`P_2` then the additional basis of the :math:`P_6` (not with the
:basis of the :math:`P_1`, the additional basis of the :math:`P_3` then the
:additional basis of the :math:`P_6`, it is possible to build the latter with
:``"FEM_GEN_HIERARCHICAL(a,b)"``)

Composite elements
++++++++++++++++++

The principal interest of the composite elements is to build hierarchical
elements. But this tool can also be used to build piecewise polynomial elements.

.. _ud-fig-triangle_comp:
.. figure:: images/getfemlisttriangleP1comp.png
   :align: center
   :scale: 60

   composite element ``"FEM_STRUCTURED_COMPOSITE(FEM_PK(2,1), 3)"``

:math:`.\\`

  .. list-table:: Composition of a finite element method on an element with ``S`` subdivisions ``"FEM_STRUCTURED_COMPOSITE(FEM1, S)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - degree of FEM1
       - dimension of FEM1
       - variable
       - variable
       - No :math:`(Q = 1)`
       - If ``FEM1`` is
       - piecewise

It is important to use a corresponding composite integration method.


Hierarchical composite elements
+++++++++++++++++++++++++++++++

.. _ud-fig-triangle_compdeux:
.. figure:: images/getfemlisttriangleP1comphier.png
   :align: center
   :scale: 60

   hierarchical composite element ``"FEM_PK_HIERARCHICAL_COMPOSITE(2,1,3)"``

:math:`.\\`

  .. list-table:: Hierarchical composition of a :math:`P_K` finite element method on a simplex with ``S`` subdivisions ``"FEM_PK_HIERARCHICAL_COMPOSITE(P,K,S)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`
       - :math:`P`
       - :math:`\dfrac{(SK+P)!}{(SK)! P!}`
       - variable
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - piecewise

:math:`.\\`

  .. list-table:: Hierarchical composition of a hierarchical :math:`P_K` finite element method on a simplex with ``S`` subdivisions ``"FEM_PK_FULL_HIERARCHICAL_COMPOSITE(P,K,S)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`
       - :math:`P`
       - :math:`\dfrac{(SK+P)!}{(SK)! P!}`
       - variable
       - No :math:`(Q = 1)`
       - Yes :math:`(M = Id)`
       - piecewise

Other constructions are possible thanks to ``"FEM_GEN_HIERARCHICAL(FEM1, FEM2)"``
and ``"FEM_STRUCTURED_COMPOSITE(FEM1, S)"``.

It is important to use a corresponding composite integration method.


Classical vector elements
----------------------------

Raviart-Thomas of lowest order elements
+++++++++++++++++++++++++++++++++++++++

.. _ud-fig-triangle_comptrois:
.. figure:: images/getfemlistRT0.png
   :align: center
   :scale: 60

   RT0 elements in dimension two and three. (P+1 dof, H(div))

:math:`.\\`

  .. list-table:: Raviart-Thomas of lowest order element on simplices ``"FEM_RT0(P)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`1`
       - :math:`P`
       - :math:`P+1`
       - H(div)
       - Yes :math:`(Q = P)`
       - No
       - Yes

:math:`.\\`

  .. list-table:: Raviart-Thomas of lowest order element on parallelepipeds (quadrilaterals, hexahedrals) ``"FEM_RT0Q(P)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`1`
       - :math:`P`
       - :math:`2P`
       - H(div)
       - Yes :math:`(Q = P)`
       - No
       - Yes


Nedelec (or Whitney) edge elements
++++++++++++++++++++++++++++++++++

.. _ud-fig-triangle_compquatre:
.. figure:: images/getfemlistnedelec.png
   :align: center
   :scale: 60

   Nedelec edge elements in dimension two and three. (P(P+1)/2 dof, H(rot))

:math:`.\\`

  .. list-table:: Nedelec (or Whitney) edge element `"FEM_NEDELEC(P)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`1`
       - :math:`P`
       - :math:`P(P+1)/2`
       - H(rot)
       - Yes :math:`(Q = P)`
       - No
       - Yes


Specific elements in dimension 1
--------------------------------


GaussLobatto element
++++++++++++++++++++

The 1D GaussLobatto :math:`P_K` element is similar to the classical :math:`P_K`
fem on the segment, but the nodes are given by the Gauss-Lobatto-Legendre
quadrature rule of order :math:`2K-1`. This FEM is known to lead to better
conditioned linear systems, and can be used with the corresponding quadrature to
perform mass-lumping (on segments or parallelepipeds).

The polynomials coefficients have been pre-computed with Maple (they require the
inversion of an ill-conditioned system), hence they are only available for the
following values of :math:`K`: :math:`1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
14, 16, 24, 32`. Note that for :math:`K=1` and :math:`K=2`, this is the classical
:math:`P1` and :math:`P2` fem.

  .. list-table:: GaussLobatto :math:`P_K` element on the segment ``"FEM_PK_GAUSSLOBATTO1D(K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`K`
       - :math:`1`
       - :math:`K+1`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes


Hermite element
+++++++++++++++

.. _ud-fig-segment_hermite:
.. figure:: images/getfemlistsegmenthermite.png
   :align: center
   :scale: 60

   :math:`P_3` Hermite element on a segment, 4 d.o.f., :math:`C^1`

Base functions on the reference element

.. math::

  \begin{array}{ll}
    \widehat{\varphi}_0 = (2x+1)(x-1)^2,&\ \ \ \widehat{\varphi}_1 = x(x-1)^2, \\
    \widehat{\varphi}_2 = x^2(3-2x),& \ \ \ \widehat{\varphi}_3 = x^2(x - 1).
  \end{array}

This element is close to be :math:`\tau`-equivalent but it is not. On the real
element the value of the gradient on vertices will be multiplied by the gradient
of the geometric transformation. The matrix :math:`M` is not equal to identity but
is still diagonal.

  .. list-table:: Hermite element on the segment ``"FEM_HERMITE(1)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`1`
       - :math:`4`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - Yes


Lagrange element with an additional bubble function
+++++++++++++++++++++++++++++++++++++++++++++++++++

.. _ud-fig-segment_bubble:
.. figure:: images/getfemlistsegmentbubble.png
   :align: center
   :scale: 60

   :math:`P_1` Lagrange element on a segment with additional internal bubble function, 3 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: Lagrange :math:`P_1` element with an additional internal bubble function ``"FEM_PK_WITH_CUBIC_BUBBLE(1, 1)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`2`
       - :math:`1`
       - :math:`3`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes


Specific elements in dimension 2
--------------------------------


Elements with additional bubble functions
+++++++++++++++++++++++++++++++++++++++++

  .. _ud-fig-triangle_p1_bubble:
  .. list-table:: Lagrange element on a triangle with additional internal bubble function
     :widths: 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlisttriangleP1bubble.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttriangleP2bubble.png
            :align: center
            :scale: 50
     * - :math:`P_1` with additional bubble function, 4 d.o.f., :math:`C^0`
       - :math:`P_2` with additional bubble function, 7 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: Lagrange :math:`P_1` or :math:`P_2` element with an additional internal bubble function ``"FEM_PK_WITH_CUBIC_BUBBLE(2, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`4` or :math:`7`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes

:math:`.\\`

.. _ud-fig-triangle_p1_bubblepie:
.. figure:: images/getfemlisttriangleP1linbubble.png
   :align: center
   :scale: 60

   :math:`P_1` Lagrange element on a triangle with additional internal piecewise linear bubble function

:math:`.\\`

  .. list-table:: Lagrange :math:`P_1` with an additional internal piecewise linear bubble function ``"FEM_P1_PIECEWISE_LINEAR_BUBBLE"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`1`
       - :math:`2`
       - :math:`4` or :math:`7`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Piecewise

:math:`.\\`

.. _ud-fig-triangle_p1_bubble_face:
.. figure:: images/getfemlisttriangleP1bubbleface.png
   :align: center
   :scale: 60

   :math:`P_1` Lagrange element on a triangle with additional bubble function on face 0, 4 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: Lagrange :math:`P_1` element with an additional bubble function on face 0 ``"FEM_P1_BUBBLE_FACE(2)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`2`
       - :math:`2`
       - :math:`4`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes

:math:`.\\`

.. _ud-fig-triangle_p1_p2_face:
.. figure:: images/getfemlisttriangleP1withP2face.png
   :align: center
   :scale: 60

   :math:`P_1` Lagrange element on a triangle with additional d.o.f on face 0, 4 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: . :math:`P_1` Lagrange element on a triangle with additional d.o.f on face 0 ``"FEM_P1_BUBBLE_FACE_LAG"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`2`
       - :math:`2`
       - :math:`4`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes


Non-conforming :math:`P_1` element
++++++++++++++++++++++++++++++++++

.. _ud-fig-triangle_non_conforming:
.. figure:: images/getfemlisttriangleP1nonconforming.png
   :align: center
   :scale: 60

   :math:`P_1` non-conforming element on a triangle, 3 d.o.f., discontinuous

:math:`.\\`

  .. list-table:: . :math:`P_1` non-conforming element on a triangle ``"FEM_P1_NONCONFORMING"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`1`
       - :math:`2`
       - :math:`3`
       - :math:`discontinuous`
       - No :math:`(Q = 1)`
       - Yes
       - Yes


Hermite element
+++++++++++++++

.. _ud-fig-triangle_hermite:
.. figure:: images/getfemlisttrianglehermite.png
   :align: center
   :scale: 60

   Hermite element on a triangle, :math:`P_3`, 10 d.o.f., :math:`C^0`

Base functions on the reference element:

.. math::

  \begin{array}{ll}
  \widehat{\varphi}_0 = (1-x-y)(1+x+y-2x^2-2y^2-11xy),~~ & (\widehat{\varphi}_0(0,0) = 1), \\
  \widehat{\varphi}_1 = x(1-x-y)(1-x-2y), & (\partial_x\widehat{\varphi}_1(0,0) = 1), \\
  \widehat{\varphi}_2 = y(1-x-y)(1-2x-y), & (\partial_y\widehat{\varphi}_2(0,0) = 1), \\
  \widehat{\varphi}_3 = -2x^3 + 7 x^2y + 7xy^2 + 3x^2 - 7xy, & (\widehat{\varphi}_3(1,0) = 1), \\
  \widehat{\varphi}_4 = x^3-2x^2y-2xy^2-x^2+2xy, & (\partial_x\widehat{\varphi}_4(1,0) = 1), \\
  \widehat{\varphi}_5 = xy(y+2x-1), & (\partial_y\widehat{\varphi}_5(1,0) = 1), \\
  \widehat{\varphi}_6 = 7x^2y + 7xy^2 - 2y^3+3y^2-7xy, & (\widehat{\varphi}_6(0,1) = 1), \\
  \widehat{\varphi}_7 = xy(x+2y-1), & (\partial_x\widehat{\varphi}_7(0,1) = 1), \\
  \widehat{\varphi}_8 = y^3-2x^2y-2xy^2-y^2+2xy, & (\partial_y\widehat{\varphi}_8(0,1) = 1), \\
  \widehat{\varphi}_9 = 27xy(1-x-y), & (\widehat{\varphi}_9(1/3,1/3) = 1), \\
  \end{array}

This element is not :math:`\tau`-equivalent (The matrix :math:`M` is not equal to
identity). On the real element linear combinations of :math:`\widehat{\varphi}_4` and
:math:`\widehat{\varphi}_7` are used to match the gradient on the corresponding vertex.
Idem for the two couples :math:`(\widehat{\varphi}_5`, :math:`\widehat{\varphi}_8)` and
:math:`(\widehat{\varphi}_6`, :math:`\widehat{\varphi}_9)` for the two other vertices.

  .. list-table:: Hermite element on a triangle ``"FEM_HERMITE(2)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`10`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - No
       - Yes


Morley element
++++++++++++++

.. _ud-fig-triangle_morley:
.. figure:: images/getfemlistmorley.png
   :align: center
   :scale: 60

   triangle Morley element, :math:`P_2`, 6 d.o.f., :math:`C^0`


This element is not :math:`\tau`-equivalent (The matrix :math:`M` is not equal to
identity). In particular, it can be used for non-conforming discretization of
fourth order problems, despite the fact that it is not :math:`{\cal C}^1`.

  .. list-table:: Morley element on a triangle ``"FEM_MORLEY"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`2`
       - :math:`2`
       - :math:`6`
       - discontinuous
       - No :math:`(Q = 1)`
       - No
       - Yes


Argyris element
+++++++++++++++

.. _ud-fig-argyris:
.. figure:: images/getfemlistargyris.png
   :align: center
   :scale: 60

   Argyris element, :math:`P_5`, 21 d.o.f., :math:`C^1`

The base functions on the reference element are:

.. math::

  \begin{array}{ll}
  \widehat{\varphi}_{0}(x,y) = 1 - 10x^3 - 10y^3 + 15x^4 - 30x^2y^2 + 15y^4 - 6x^5 + 30x^3y^2 + 30x^2y^3 - 6y^5, & (\widehat{\varphi}_0(0,0) = 1), \\
  \widehat{\varphi}_{1}(x,y) = x - 6x^3 - 11xy^2 + 8x^4 + 10x^2y^2 + 18xy^3 - 3x^5 + x^3y^2 - 10x^2y^3 - 8xy^4, & (\partial_x\widehat{\varphi}_1(0,0) = 1),\\
  \widehat{\varphi}_{2}(x,y) = y - 11x^2y - 6y^3 + 18x^3y + 10x^2y^2 + 8y^4 - 8x^4y - 10x^3y^2 + x^2y^3 - 3y^5, & (\partial_y\widehat{\varphi}_2(0,0) = 1),\\
  \widehat{\varphi}_{3}(x,y) = 0.5x^2 - 1.5x^3 + 1.5x^4 - 1.5x^2y^2 - 0.5x^5 + 1.5x^3y^2 + x^2y^3, & (\partial^2_{xx}\widehat{\varphi}_3(0,0) = 1),\\
  \widehat{\varphi}_{4}(x,y) = xy - 4x^2y - 4xy^2 + 5x^3y + 10x^2y^2 + 5xy^3 - 2x^4y - 6x^3y^2 - 6x^2y^3 - 2xy^4, & (\partial^2_{xy}\widehat{\varphi}_{4}(0,0) = 1),\\
  \widehat{\varphi}_{5}(x,y) = 0.5y^2 - 1.5y^3 - 1.5x^2y^2 + 1.5y^4 + x^3y^2 + 1.5x^2y^3 - 0.5y^5, & (\partial^2_{yy}\widehat{\varphi}_{5}(0,0) = 1),\\
  \widehat{\varphi}_{6}(x,y) = 10x^3 - 15x^4 + 15x^2y^2 + 6x^5 - 15x^3y^2 - 15x^2y^3, & (\widehat{\varphi}_6(1,0) = 1),\\
  \widehat{\varphi}_{7}(x,y) = -4x^3 + 7x^4 - 3.5x^2y^2 - 3x^5 + 3.5x^3y^2 + 3.5x^2y^3, & (\partial_x\widehat{\varphi}_7(1,0) = 1),\\
  \widehat{\varphi}_{8}(x,y) = -5x^2y + 14x^3y + 18.5x^2y^2 - 8x^4y - 18.5x^3y^2 - 13.5x^2y^3, & (\partial_y\widehat{\varphi}_8(1,0) = 1),\\
  \widehat{\varphi}_{9}(x,y) = 0.5x^3 - x^4 + 0.25x^2y^2 + 0.5x^5 - 0.25x^3y^2 - 0.25x^2y^3, & (\partial^2_{xx}\widehat{\varphi}_{9}(1,0) = 1),\\
  \widehat{\varphi}_{10}(x,y) = x^2y - 3x^3y - 3.5x^2y^2 + 2x^4y + 3.5x^3y^2 + 2.5x^2y^3, & (\partial^2_{xy}\widehat{\varphi}_{10}(1,0) = 1),\\
  \widehat{\varphi}_{11}(x,y) = 1.25x^2y^2 - 0.75x^3y^2 - 1.25x^2y^3, & (\partial^2_{yy}\widehat{\varphi}_{11}(1,0) = 1),\\
  \widehat{\varphi}_{12}(x,y) = 10y^3 + 15x^2y^2 - 15y^4 - 15x^3y^2 - 15x^2y^3 + 6y^5, & (\widehat{\varphi}_{12}(0,1) = 1),\\
  \widehat{\varphi}_{13}(x,y) = -5xy^2 + 18.5x^2y^2 + 14xy^3 - 13.5x^3y^2 - 18.5x^2y^3 - 8xy^4, & (\partial_x\widehat{\varphi}_{13}(0,1) = 1),\\
  \widehat{\varphi}_{14}(x,y) = -4y^3 - 3.5x^2y^2 + 7y^4 + 3.5x^3y^2 + 3.5x^2y^3 - 3y^5, & (\partial_y\widehat{\varphi}_{14}(0,0) = 1),\\
  \widehat{\varphi}_{15}(x,y) = 1.25x^2y^2 - 1.25x^3y^2 - 0.75x^2y^3, & (\partial^2_{xx}\widehat{\varphi}_{15}(0,1) = 1),\\
  \widehat{\varphi}_{16}(x,y) = xy^2 - 3.5x^2y^2 - 3xy^3 + 2.5x^3y^2 + 3.5x^2y^3 + 2xy^4, & (\partial^2_{xy}\widehat{\varphi}_{16}(0,1) = 1),\\
  \widehat{\varphi}_{17}(x,y) = 0.5y^3 + 0.25x^2y^2 - y^4 - 0.25x^3y^2 - 0.25x^2y^3 + 0.5y^5, & (\partial^2_{yy}\widehat{\varphi}_{17}(0,1) = 1),\\
  \widehat{\varphi}_{18}(x,y) = \sqrt{2}(-8x^2y^2 + 8x^3y^2 + 8x^2y^3), & ~\hspace{-10.5em}(\sqrt{0.5}(\partial_{x}\widehat{\varphi}_{18}(0.5,0.5) + \partial_{y}\widehat{\varphi}_{18}(0.5,0.5)) = 1),\\
  \widehat{\varphi}_{19}(x,y) = -16xy^2 + 32x^2y^2 + 32xy^3 - 16x^3y^2 - 32x^2y^3 - 16xy^4, & (-\partial_{x}\widehat{\varphi}_{19}(0,0.5) = 1),\\
  \widehat{\varphi}_{20}(x,y) = -16x^2y + 32x^3y + 32x^2y^2 - 16x^4y - 32x^3y^2 - 16x^2y^3, & (-\partial_{y}\widehat{\varphi}_{20}(0.5,0) = 1),\\
  \end{array}

This element is not :math:`\tau`-equivalent (The matrix :math:`M` is not equal to
identity). On the real element linear combinations of the transformed base
functions :math:`\widehat{\varphi}_i` are used to match the gradient, the second
derivatives and the normal derivatives on the faces. Note that the use of the
matrix :math:`M` allows to define Argyris element even with nonlinear geometric
transformations (for instance to treat curved boundaries).

  .. list-table:: Argyris element on a triangle ``"FEM_ARGYRIS"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`5`
       - :math:`2`
       - :math:`21`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - Yes


Hsieh-Clough-Tocher element
+++++++++++++++++++++++++++

.. _ud-fig-HCT_tr:
.. figure:: images/getfemlistHCT.png
   :align: center
   :scale: 60

   Hsieh-Clough-Tocher (HCT) element, :math:`P_3`, 12 d.o.f., :math:`C^1`


This element is not :math:`\tau`-equivalent. This is a composite element.
Polynomial of degree 3 on each of the three sub-triangles (see figure
:ref:`ud-fig-HCT_tr` and [ciarlet1978]_). It is strongly advised to use a
``"IM_HCT_COMPOSITE"`` integration method with this finite element. The numeration
of the dof is the following: 0, 3 and 6 for the lagrange dof on the first second
and third vertex respectively; 1, 4, 7 for the derivative with respects to the
first variable; 2, 5, 8 for the derivative with respects to the second variable
and 9, 10, 11 for the normal derivatives on face 0, 1, 2 respectively.

  .. list-table:: HCT element on a triangle ``"FEM_HCT_TRIANGLE"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`12`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - piecewise

:math:`.\\`

.. _ud-fig-reduced_HCT_tr:
.. figure:: images/getfemlistreducedHCT.png
   :align: center
   :scale: 60

   Reduced Hsieh-Clough-Tocher (reduced HCT) element, :math:`P_3`, 9 d.o.f., :math:`C^1`

This element exists also in its reduced form, where the normal derivatives are
assumed to be polynomial of degree one on each edge (see figure
:ref:`ud-fig-reduced_HCT_tr`)

  .. list-table:: Reduced HCT element on a triangle ``"FEM_REDUCED_HCT_TRIANGLE"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`9`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - piecewise


A composite :math:`C^1` element on quadrilaterals
+++++++++++++++++++++++++++++++++++++++++++++++++

.. _ud-fig-QC1_tr:
.. figure:: images/getfemlistquadc1composite.png
   :align: center
   :scale: 60

   Composite element on quadrilaterals, piecewise :math:`P_3`, 16 d.o.f., :math:`C^1`


This element is not :math:`\tau`-equivalent. This is a composite element.
Polynomial of degree 3 on each of the four sub-triangles (see figure
:ref:`ud-fig-QC1_tr`). At least on the reference element it corresponds to the
Fraeijs de Veubeke-Sander element (see  [ciarlet1978]_). It is strongly advised
to use a ``"IM_QUADC1_COMPOSITE"`` integration method with this finite element.

  .. list-table:: . :math:`C^1` composite element on a quadrilateral (FVS) ``"FEM_QUADC1_COMPOSITE"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`16`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - piecewise

:math:`.\\`

.. _ud-fig-reduced_QC1_tr:
.. figure:: images/getfemlistreducedquadc1composite.png
   :align: center
   :scale: 60

   Reduced composite element on quadrilaterals, piecewise :math:`P_3`, 12 d.o.f., :math:`C^1`


This element exists also in its reduced form, where the normal derivatives are
assumed to be polynomial of degree one on each edge (see figure
:ref:`ud-fig-reduced_QC1_tr`)

  .. list-table:: Reduced :math:`C^1` composite element on a quadrilateral (reduced FVS) ``"FEM_REDUCED_QUADC1_COMPOSITE"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`2`
       - :math:`12`
       - :math:`C^1`
       - No :math:`(Q = 1)`
       - No
       - piecewise


Specific elements in dimension 3
--------------------------------

Lagrange elements on 3D pyramid
+++++++++++++++++++++++++++++++

|gf| proposes some Lagrange pyramidal elements of degree 0, 1 and two based on [GR-GH1999]_ and [BE-CO-DU2010]_. See these references for more details. The proposed element can be raccorded to standard :math:`P_1` or :math:`P_2` Lagrange fem on the triangular faces and to a standard :math:`Q_1` or :math:`Q_2` Lagrange fem on the quatrilateral face.

.. _ud-fig-pyramid-lagrange:
.. list-table:: Lagrange element on a pyramidal element of order 0, 1 and 2
   :widths: 20 20 20
   :header-rows: 0
   :class: figure

   * - .. image:: images/getfemlistpyramidP0.png
          :align: center
          :scale: 50
     - .. image:: images/getfemlistpyramidP1.png
          :align: center
          :scale: 50
     - .. image:: images/getfemlistpyramidP2.png
          :align: center
          :scale: 50

   * - Degree 0 pyramidal element with 1 dof
     - Degree 1 pyramidal element with 5 dof
     - Degree 2 pyramidal element with 14 dof


The associated geometric transformations are ``"GT_PYRAMID(K)"`` for K = 1 or 2. The associated integration methods ``"IM_PYRAMID(im)"`` where ``im`` is an integration method on a hexahedron (or alternatively ``"IM_PYRAMID_COMPOSITE(im)"`` where ``im`` is an integration method on a tetrahedron, but it is theoretically less accurate)
The shape functions are not polynomial ones but rational fractions. For the first degree the shape functions read:

.. math::

   \begin{array}{l}
   \widehat{\varphi}_{0}(x,y,z) =  \frac{1}{4}\left(1-x-y-z+\dfrac{xy}{1-z}\right), \\
   \widehat{\varphi}_{1}(x,y,z) =  \frac{1}{4}\left(1+x-y-z-\dfrac{xy}{1-z}\right), \\
   \widehat{\varphi}_{2}(x,y,z) =  \frac{1}{4}\left(1-x+y-z-\dfrac{xy}{1-z}\right), \\
   \widehat{\varphi}_{3}(x,y,z) =  \frac{1}{4}\left(1+x+y-z+\dfrac{xy}{1-z}\right), \\
   \widehat{\varphi}_{4}(x,y,z) =  z.\\
   \end{array}

For the second degree, setting

.. math::

   \xi_0 = \dfrac{1-z-x}{2}, ~~~\xi_1 = \dfrac{1-z-y}{2}, ~~~\xi_2 = \dfrac{1-z+x}{2}, ~~~\xi_3 = \dfrac{1-z+y}{2}, ~~~\xi_4 = z,

the shape functions read:

.. math::

   \begin{array}{l}
   \widehat{\varphi}_{0}(x,y,z) = \dfrac{\xi_0 \xi_1}{(1-\xi_4)^2}((1-\xi_4-2\xi_0)(1-\xi_4-2\xi_1) -\xi_4(1-\xi_4)), \\
   \widehat{\varphi}_{1}(x,y,z) = 4\dfrac{\xi_0\xi_1\xi_2}{(1-\xi_4)^2}(2\xi_1-(1-\xi_4)), \\
   \widehat{\varphi}_{2}(x,y,\xi_4) = \dfrac{\xi_1 \xi_2}{(1-\xi_4)^2}((1-\xi_4-2\xi_1)(1-\xi_4-2\xi_2) -\xi_4(1-\xi_4)), \\
   \widehat{\varphi}_{3}(x,y,z) = 4\dfrac{\xi_3\xi_0\xi_1}{(1-\xi_4)^2}(2\xi_0-(1-\xi_4)), \\
   \widehat{\varphi}_{4}(x,y,z) = 16\dfrac{\xi_0\xi_1\xi_2\xi_3}{(1-\xi_4)^2}, \\
   \widehat{\varphi}_{5}(x,y,z) = 4\dfrac{\xi_1\xi_2\xi_3}{(1-\xi_4)^2}(2\xi_2-(1-\xi_4)), \\
   \widehat{\varphi}_{6}(x,y,z) = \dfrac{\xi_3 \xi_0}{(1-\xi_4)^2}((1-\xi_4-2\xi_3)(1-\xi_4-2\xi_0) -\xi_4(1-\xi_4)), \\
   \widehat{\varphi}_{7}(x,y,z) = 4\dfrac{\xi_2\xi_3\xi_0}{(1-\xi_4)^2}(2\xi_3-(1-\xi_4)), \\
   \widehat{\varphi}_{8}(x,y,z) = \dfrac{\xi_2 \xi_3}{(1-\xi_4)^2}((1-\xi_4-2\xi_2)(1-\xi_4-2\xi_3) -\xi_4(1-\xi_4)), \\
   \widehat{\varphi}_{9}(x,y,z) = 4\dfrac{\xi_4}{1-\xi_4}\xi_0\xi_1, \\
   \widehat{\varphi}_{10}(x,y,z) = 4\dfrac{\xi_4}{1-\xi_4}\xi_1\xi_2,  \\
   \widehat{\varphi}_{11}(x,y,z) = 4\dfrac{\xi_4}{1-\xi_4}\xi_3\xi_0,  \\
   \widehat{\varphi}_{12}(x,y,z) = 4\dfrac{\xi_4}{1-\xi_4}\xi_2\xi_3,  \\
   \widehat{\varphi}_{13}(x,y,z) = \xi_4(2\xi_4-1). \\
   \end{array}

.. list-table:: Continuous Lagrange element of order 0, 1 or 2 ``"FEM_PYRAMID_LAGRANGE(K)"``
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - degree
     - dimension
     - d.o.f. number
     - class
     - vector
     - :math:`\tau`-equivalent
     - Polynomial

   * - :math:`0`
     - :math:`3`
     - :math:`1`
     - discontinuous
     - No :math:`(Q = 1)`
     - Yes
     - No

   * - :math:`1`
     - :math:`3`
     - :math:`5`
     - :math:`C^0`
     - No :math:`(Q = 1)`
     - Yes
     - No

   * - :math:`2`
     - :math:`3`
     - :math:`14`
     - :math:`C^0`
     - No :math:`(Q = 1)`
     - Yes
     - No



.. list-table:: Discontinuous Lagrange element of order 0, 1 or 2 ``"FEM_PYRAMID_DISCONTINUOUS_LAGRANGE(K)"``
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - degree
     - dimension
     - d.o.f. number
     - class
     - vector
     - :math:`\tau`-equivalent
     - Polynomial

   * - :math:`0`
     - :math:`3`
     - :math:`1`
     - discontinuous
     - No :math:`(Q = 1)`
     - Yes
     - No

   * - :math:`1`
     - :math:`3`
     - :math:`5`
     - discontinuous
     - No :math:`(Q = 1)`
     - Yes
     - No

   * - :math:`2`
     - :math:`3`
     - :math:`14`
     - discontinuous
     - No :math:`(Q = 1)`
     - Yes
     - No


	

Elements with additional bubble functions
+++++++++++++++++++++++++++++++++++++++++

  .. _ud-fig-tetrahedron_p1_bubble:
  .. list-table:: Lagrange element on a tetrahedron with additional internal bubble function
     :widths: 30 30 30
     :header-rows: 0
     :class: figure

     * - .. image:: images/getfemlisttetrahedronP1bubble.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttetrahedronP2bubble.png
            :align: center
            :scale: 50
       - .. image:: images/getfemlisttetrahedronP3bubble.png
            :align: center
            :scale: 50
     * - :math:`P_1` with additional bubble function, 5 d.o.f., :math:`C^0`
       - :math:`P_2` with additional bubble function, 11 d.o.f., :math:`C^0`
       - :math:`P_3` with additional bubble function, 21 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: :math:`P_K` Lagrange element with an additional internal bubble function ``"FEM_PK_WITH_CUBIC_BUBBLE(3, K)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`4`
       - :math:`3`
       - :math:`5`, :math:`11` or :math:`21`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes

:math:`.\\`

.. _ud-fig-tetrahedron_p1_bubble_face:
.. figure:: images/getfemlisttetrahedronP1bubbleface.png
   :align: center
   :scale: 60

   :math:`P_1` Lagrange element on a tetrahedron with additional bubble function on face 0, 5 d.o.f., :math:`C^0`

:math:`.\\`

  .. list-table:: Lagrange :math:`P_1` element with an additional bubble function on face 0 ``"FEM_P1_BUBBLE_FACE(3)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`3`
       - :math:`5`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - Yes
       - Yes


Hermite element
+++++++++++++++

.. _ud-fig-tetrahedron_hermite:
.. figure:: images/getfemlisttetrahedronhermite.png
   :align: center
   :scale: 60

   Hermite element on a tetrahedron, :math:`P_3`, 20 d.o.f., :math:`C^0`

Base functions on the reference element:

.. math::

  \begin{array}{ll}
  \widehat{\varphi}_{0}(x,y) = 1 - 3x^2 - 13xy - 13xz - 3y^2 - 13yz - 3z^2 + 2x^3 + 13x^2y + 13x^2z & \\
  ~~~~~~~~~~~~~~~ + 13xy^2 + 33xyz + 13xz^2 + 2y^3 + 13y^2z + 13yz^2 + 2z^3, & (\widehat{\varphi}_0(0,0,0) = 1),\\
  \widehat{\varphi}_{1}(x,y) = x - 2x^2 - 3xy - 3xz + x^3 + 3x^2y + 3x^2z + 2xy^2 + 4xyz + 2xz^2, & (\partial_x\widehat{\varphi}_1(0,0,0) = 1),\\
  \widehat{\varphi}_{2}(x,y) = y - 3xy - 2y^2 - 3yz + 2x^2y + 3xy^2 + 4xyz + y^3 + 3y^2z + 2yz^2, & (\partial_y\widehat{\varphi}_2(0,0,0) = 1),\\
  \widehat{\varphi}_{3}(x,y) = z - 3xz - 3yz - 2z^2 + 2x^2z + 4xyz + 3xz^2 + 2y^2z + 3yz^2 + z^3, & (\partial_z\widehat{\varphi}_3(0,0,0) = 1),\\
  \widehat{\varphi}_{4}(x,y) = 3x^2 - 7xy - 7xz - 2x^3 + 7x^2y + 7x^2z + 7xy^2 + 7xyz + 7xz^2, & (\widehat{\varphi}_4(1,0,0) = 1),\\
  \widehat{\varphi}_{5}(x,y) = -x^2 + 2xy + 2xz + x^3 - 2x^2y - 2x^2z - 2xy^2 - 2xyz - 2xz^2, & (\partial_x\widehat{\varphi}_5(1,0,0) = 1),\\
  \widehat{\varphi}_{6}(x,y) = -xy + 2x^2y + xy^2, & (\partial_y\widehat{\varphi}_6(1,0,0) = 1),\\
  \widehat{\varphi}_{7}(x,y) = -xz + 2x^2z + xz^2, & (\partial_z\widehat{\varphi}_7(1,0,0) = 1),\\
  \widehat{\varphi}_{8}(x,y) = -7xy + 3y^2 - 7yz + 7x^2y + 7xy^2 + 7xyz - 2y^3 + 7y^2z + 7yz^2, & (\widehat{\varphi}_8(0,1,0) = 1),\\
  \widehat{\varphi}_{9}(x,y) = -xy + x^2y + 2xy^2, & (\partial_x\widehat{\varphi}_9(0,1,0) = 1),\\
  \widehat{\varphi}_{10}(x,y) = 2xy - y^2 + 2yz - 2x^2y - 2xy^2 - 2xyz + y^3 - 2y^2z - 2yz^2, & (\partial_y\widehat{\varphi}_{10}(0,1,0) = 1),\\
  \widehat{\varphi}_{11}(x,y) = -yz + 2y^2z + yz^2, & (\partial_z\widehat{\varphi}_{11}(0,1,0) = 1),\\
  \widehat{\varphi}_{12}(x,y) = -7xz - 7yz + 3z^2 + 7x^2z + 7xyz + 7xz^2 + 7y^2z + 7yz^2 - 2z^3, & (\widehat{\varphi}_{12}(0,0,1) = 1),\\
  \widehat{\varphi}_{13}(x,y) = -xz + x^2z + 2xz^2, & (\partial_x\widehat{\varphi}_{13}(0,0,1) = 1),\\
  \widehat{\varphi}_{14}(x,y) = -yz + y^2z + 2yz^2, & (\partial_y\widehat{\varphi}_{14}(0,0,1) = 1),\\
  \widehat{\varphi}_{15}(x,y) = 2xz + 2yz - z^2 - 2x^2z - 2xyz - 2xz^2 - 2y^2z - 2yz^2 + z^3, & (\partial_z\widehat{\varphi}_{15}(0,0,1) = 1),\\
  \widehat{\varphi}_{16}(x,y) = 27xyz, & (\widehat{\varphi}_{16}(1/3,1/3,1/3) = 1),\\
  \widehat{\varphi}_{17}(x,y) = 27yz - 27xyz - 27y^2z - 27yz^2, & (\widehat{\varphi}_{17}(0,1/3,1/3) = 1),\\
  \widehat{\varphi}_{18}(x,y) = 27xz - 27x^2z - 27xyz - 27xz^2, & (\widehat{\varphi}_{18}(1/3,0,1/3) = 1),\\
  \widehat{\varphi}_{19}(x,y) = 27xy - 27x^2y - 27xy^2 - 27xyz, & (\widehat{\varphi}_{19}(1/3,1/3,0) = 1),\\
  \end{array}

This element is not :math:`\tau`-equivalent (The matrix :math:`M` is not equal to
identity). On the real element linear combinations of :math:`\widehat{\varphi}_8`,
:math:`\widehat{\varphi}_{12}` and :math:`\widehat{\varphi}_{16}` are used to match the gradient on
the corresponding vertex. Idem on the other vertices.

  .. list-table:: Hermite element on a tetrahedron ``"FEM_HERMITE(3)"``
     :widths: 10 10 10 10 10 10 10
     :header-rows: 1

     * - degree
       - dimension
       - d.o.f. number
       - class
       - vector
       - :math:`\tau`-equivalent
       - Polynomial

     * - :math:`3`
       - :math:`3`
       - :math:`20`
       - :math:`C^0`
       - No :math:`(Q = 1)`
       - No
       - Yes
