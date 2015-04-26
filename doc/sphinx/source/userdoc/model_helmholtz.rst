.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-helmholtz:


Helmholtz brick
---------------

This brick represents the complex or real Helmholtz problem:

.. math::

   \Delta u + k^2 u = \ldots

where :math:`k` the wave number is a real or complex value. For a complex
version, a complex model has to be used (see :file:`tests/helmholtz.cc`).

The function adding a Helmholtz brick to a model is::

  getfem::add_Helmholtz_brick(md, mim, varname, dataexpr, region);

where ``varname`` is the variable on which the Helmholtz term is added and
``dataexpr`` is the wave number.
