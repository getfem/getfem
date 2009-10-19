.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-appendixa:

Appendix A. Finite element method list
======================================

Let us recall that all finite element methods defined in |gf| are declared in the
file ``getfem_fem.h`` and that a descriptor on a finite element method is obtained
thanks to the function::

  getfem::pfem pf = getfem::fem_descriptor("name of method");

where ``"name of method"`` is a string to be choosen among the existing methods.


Dof graphical codification
--------------------------

.. _fig-symbols:
.. figure:: images/getfemlistsymbols.png
   :align: center
   :width: 450pt

   Symbols representing degree of freedom types


