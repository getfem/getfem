.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-contact-friction-large:

Large sliding/large deformation contact with friction bricks
------------------------------------------------------------

These bricks present some algorithms for contact and friction in the large sliding/large deformation framework. Of course, their computational cost is greatly higher than small sliding-small deformation bricks.

The multi-contact frame object
++++++++++++++++++++++++++++++

A \gf object is dedicated to the computation of effective contact surfaces which is shared by all the bricks. This object stores the different potential contact surfaces. On most of methods, potential contact surface are classified into two categories: master and slave surface (see  :ref:`figure<ud-fig-masterslave>`).

.. ud-fig-masterslave:
.. figure:: images/getfemusermasterslave.png
   :align: center
   :scale: 60

The slave surface is the "contactor" and the master one the "target". Rigid obstacle are also considered. They are always master surfaces.  The basic rule is that the contact is considered between a slave surface and a master one. However, the multi-contact frame object and the \gf bricks allow multi-contact situations, including contact between two master surfaces, self-contact of a master surface and an arbitrary number of slave and master surfaces. 

Basically, in order to detect the contact pairs, Gauss points or f.e.m. nodes of slave surfaces are projected on master surfaces (see  :ref:`figure<ud-fig-masterslave>`). If self-contact is considered, Gauss points or f.e.m. nodes of master surface are also projected on master surfaces.

The use of multi-contact frame object
*************************************



The contact pair detection algorithm
************************************

A contact pair is ...

scheme of contact pair detection

General heuristic

Avoid false contact detection ...

two basic algoritms: Delaunay and influence boxes ...