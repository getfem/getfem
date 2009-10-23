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


Classical :math:`P_K` Lagrange elements on simplices
----------------------------------------------------

It is possible to define a classical :math:`P_K` Lagrange element of arbitrary 
dimension and arbitrary degree. Each degree of freedom of such an element 
corresponds to the value of the function on a corresponding node. The grid of 
node is the so-called Lagrange grid. Figures :ref:`ud-fig-segmentpk`,
:ref:`ud-fig-trianglepk`

.. _ud-fig-segmentpk:
.. figure:: images/getfemlistsegmentPk.png
   :align: center
   :width: 300pt

   Examples of classical :math:`P_K` Lagrange elements on a segment

.. _ud-fig-trianglepk:
.. tabularcolumns:: cccc
.. csv-table:: Examples of classical :math:`P_K` Lagrange elements on a triangle.
   :class: figure
   :widths: 1, 2, 2, 1
   :delim: &

   & .. image:: images/getfemlisttriangleP1.png & .. image:: images/getfemlisttriangleP2.png &
   & :math:`P_1`, 3 d.o.f., :math:`C^0` & :math:`P_2` element, 6 d.o.f., :math:`C^0` &
   & .. image:: images/getfemlisttriangleP3.png & .. image:: images/getfemlisttriangleP6.png &
   & :math:`P_3`, 10 d.o.f., :math:`C^0` & :math:`P_6` element, 28 d.o.f., :math:`C^0` &

.. _ud-fig-symbols:
.. figure:: images/getfemlistsymbols.png
   :align: center
   :width: 450pt

   Symbols representing degree of freedom types
