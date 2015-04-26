.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-mass:


Mass brick
----------

This brick represents a weak term of the form

.. math::

   \int_{\Omega} \rho u\cdot v\ dx + \ldots

It mainly represents a mass term for transient problems but can also be used for
other applications (it can be used on a boundary). Basically, this brick adds a
mass matrix on the tangent linear system with respect to a certain variable.

The function which adds this brick to a model is::

  ind_brick = getfem::add_mass_brick
              (md, mim, varname, dataexpr_rho="", region = size_type(-1));

where ``dataexpr_rho`` is an optional expression representing the density
:math:`\rho`. If it is omitted, the density is assumed to be equal to one.
