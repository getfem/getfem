.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-iinteg:

Incorporate new approximated integration methods in |gf|
========================================================

A perl script automatically incorporates new cubature methods from a description
file. You can see in the directory ``cubature`` such description files (with
extension ``.IM``) . For instance for ``IM_TETRAHEDRON(5)`` the following file
describes the method::

  NAME = IM_TETRAHEDRON(5)
  N = 3
  GEOTRANS = GT_PK(3,1)
  NBPT = 4
  0, 0.25, 0.25, 0.25, 0.008818342151675485
  1, 0.31979362782962991, 0.31979362782962991, 0.31979362782962991, 0.011511367871045398
  1, 0.091971078052723033, 0.091971078052723033, 0.091971078052723033, 0.01198951396316977
  1, 0.056350832689629156, 0.056350832689629156, 0.44364916731037084, 0.008818342151675485
  NBF = 4 IM_TRIANGLE(5)
  IM_TRIANGLE(5)
  IM_TRIANGLE(5)
  IM_TRIANGLE(5)

where ``NAME`` is the name of the method in |gf| (constant integer parameter are
allowed), ``N`` is the dimension, ``GEOTRANS`` describes a valid geometric
transformation of |gf|. This geometric transformation just defines the reference
element on which the integration method is described. ``NBPT`` is the number of
integration node definitions. Integration node definitions include a symmetry
definition such that the total number of integration nodes would be greater than
``NBPT``.

Composition of the integration node definition:

* an integer: 0 = no symmetry, 1 = full symmetric (x6 for a triangle, x4 for a
  quadrangle, x24 for a tetrahedron ...),

* the ``N`` coordinates of the integration node,

* the load.

``NBF`` is the number of faces of the reference element (should
correspond to ``GEOTRANS``). Then follows an already existing
integration method for each face (each on a line). This is necessary
to make integrations on boundaries.

The file format is inspired from [EncyclopCubature]_.
