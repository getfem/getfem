.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc_fem:

Fem module
----------

Description
^^^^^^^^^^^

The Fem module is the part of |gf| which describes the finite elements at the
element level and the degrees of freedom. Finite element methods can be of
different types. They could be scalar or vectorial, polynomial, piecewise
polynomial or non-polynomial, equivalent via the geometric transformation or not.
Moreover, the description of the degrees of freedom have to be such that it is
possible to gather the compatible degrees of freedom between two neighbour
elements in a generic way (for instance connecting a Lagrange 2D element to
another Lagrange 1D element).


Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`bgeot_poly.h` and :file:`bgeot_poly_composite.h` and :file:`bgeot_poly.cc` and :file:`bgeot_poly_composite.cc`, "Some classes to represent polynomials and piecewise polynomials in order to describe shape functions of a finite element method on the reference element."
   :file:`getfem_fem.h` and :file:`getfem_fem.cc` and :file:`getfem_fem_composite.cc`, "Descriptors for finite element and a degree of freedom. Polynomial finite elements are defined in :file:`getfem_fem.cc` and piecewise polynomial finite elements in :file:`getfem_fem_composite.cc`"
   :file:`getfem_fem_global_function.h` and :file:`getfem_fem_global_function.cc`, "Defines a fem with base functions defined as global functions given by the user. Useful for enrichment with singular functions and for implementation of meshless methods."
   :file:`getfem_projected_fem.h` and :file:`getfem_projected_fem.cc`, "Defines a fem which is the projection of a finite element space (represented by a mesh_fem) on a different mesh. Note that the high-generic assembly language offers also this functionality by means of the interpolated transformations."
   :file:`getfem_interpolated_fem.h` and :file:`getfem_interpolated_fem.cc`, "Dfines a fem which is the interpolation of a finite element space (represented by a mesh_fem) on a different mesh. Note that the high-generic assembly language offers also this functionality by means of the interpolated transformations."




State
^^^^^



The two files :file:`getfem_fem.cc` and :file:`getfem_fem_composite.cc` mainly
contains all the finite element description for basic elements. A exhaustive list
of the defined finite elements is given in :ref:`ud-appendixa`.

Some other files define some specific finite element such as
:file:`getfem_fem_level_set.h` which is a complex construction which allows to
"cut" a existing element by one or several level sets.

The manner to describe the degrees of freedom globally satisfies the needing
(connecting dof from an element to another in a generic way) but is a little bit
obscure and too much complicated.

Conversely, the way to represent non-equivalent elements with the supplementary
matrix ``M`` has proven its efficiency on several elements (Hermites elements,
Argyris, etc.).

Perspectives
^^^^^^^^^^^^

The principal dissatisfaction of this module is that description of the degrees
of freedom is not completely satisfactory. It is the principal reason why one
documentation on how to build an element from A to Z was not made for the moment
because description of the degrees of freedom was conceived to be temporary. An
effort of design is thus to be provided to completely stabilize this module
mainly thus with regard to the description of degrees of freedom but also perhaps
the description of finite elements which could be partially externalized in a
similar way to the cubature methods , at least for the simplest finite elements
(equivalent and polynomial finite elements).

