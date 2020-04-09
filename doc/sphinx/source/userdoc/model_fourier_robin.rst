.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. index:: models, model bricks

.. _ud-model-fourier-robin:




Fourier-Robin brick
-------------------

This brick can be used to add boundary conditions of Fourier-Robin type like:

.. math::

   \frac{\partial u}{\partial \nu} = Qu

for scalar problems, or

.. math::

   \sigma\cdot \nu = Qu

for linearized elasticity problems. ``Q`` is a scalar field in the scalar case or
a matrix field in the vectorial case. This brick works for both real or complex
terms in scalar or vectorial problems.

The function adding this brick to a model is::

  add_Fourier_Robin_brick(md, mim, varname, dataexpr, region);

where ``dataexpr`` is the data of the model which represents the coefficient
:math:`Q`.  It can be an arbitrary valid expression of GWFL, the generic weak form language (except for the complex version for which it should be a data of the model)

Note that an additional right hand side can be added with a source term brick.

