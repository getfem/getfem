.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-appendixb:

Appendix B. Cubature method list
================================

The integration methods are of two kinds. Exact integrations of polynomials and
approximated integrations (cubature formulas) of any function. The exact
integration can only be used if all the elements are polynomial and if the
geometric transformation is linear.

A descriptor on an integration method is given by the function::

  ppi = getfem::int_method_descriptor("name of method");

where ``"name of method"`` is a string to be choosen among the existing methods.

The program ``integration`` located in the ``tests`` directory lists and checks
the degree of each integration method.
