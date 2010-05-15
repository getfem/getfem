.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-ifem:

Incorporate new finite element methods in |gf|
==============================================

Basically, It is sufficient to describe an element on the reference element, i.e.
to describe each base function of each degree of freedom. Intrinsically vectorial
elements are supported (see for instance Nedelec and Raviart-Thomas elements).
Finite element methods that are not equivalent via the geometric transformation
(not :math:`\tau`-equivalent in |gf| jargon, such as vectorial elements, Hermite
elements ...) an additional linear transformation of the degrees of freedom
depending on the real element should be described (see the implementation of
Argyris element for instance).

Please read :ref:`dp` for more details and see the files
:file:`getfem/getfem_fem.h`, :file:`getfem_fem.cc` for practical implementation.
