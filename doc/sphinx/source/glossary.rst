
.. include:: replaces.txt

********
Glossary
********

.. if you add new entries, keep the alphabetical sorting!

.. glossary::

   Convex
      See **element**

   Cubature method
      A cubature method on an **element** consists in a set of nodes
      (generally called gauss points) and corresponding loads which
      define a approximated integration method. In |Gf| it is defined
      on the **reference elements**.

   Degree of freedom
      The degrees of freedom for a finite element method is the coefficients
      which multiply the shape functions in order to describe a
      (scalar or vector) field. Generally, they are the unknowns of the
      problem in general.

   Element
      An element is a small piece of a domain with a special shape (a segment,
      a triangle, a quadrilateron, an tetrahedron, a hexahedron or a prism)
      for dimensions less or equal to three. A mesh is the union of
      non intersecting elements.

   Finite element method (fem)
      A finite element method is defined on a real element. It consist on a
      certain number of degrees of freedom linked to the corresponding shape
      functions and a manner to glue the degrees of freedom from a element
      to a neighbour element.

   Integration method
      See **cubature method**.

   Quadrature method
      See **cubature method**.

   Mesh
      The mesh is composed of **elements**. in |gf|, these elements are
      often called **convexes**. A mesh can be composed of elements of different
      dimensions (triangles, segments, quadrilaters, tetrahedra,
      hexahedra ...).

   Mesh_Fem
      The mesh_fem object is a mesh with a **finite element method** defined
      on each **element**. This
      represent a finite element space on which a unknown or a data on the
      considered domain will be described.

   Mesh_Im
      The mesh_im object is a mesh with a **cubature method** defined on
      each **element**. It is used in assembly procedures.

   Reference element
      A reference element or a convex of reference is a special **element**
      on which the elementary computations (integrals) are performed.
      For instance, the reference segment in |gf| is the segment [0,1].
      The reference triangle is the triangle (0,0), (0,1), (1,0). etc.


