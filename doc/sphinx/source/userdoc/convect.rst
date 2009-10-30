.. $Id: interMM.rst 3275 2009-10-29 19:14:14Z lsaavedr $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-intermm:

A pure convection method 
========================

A method to compute a pure convection is defined in the file :file:`getfem/getfem_convect.h`. The call of the function is

getfem::convect(mf, U, mf_v, V, dt, nt);

where ``mf`` is a variable of type |gf_mf|, ``U`` is a vector which represent the field to be convected, ``mf_v`` is a |gf_mf| for the velocity field, ``V`` is the dof vector for the velocity field, ``dt`` is the pseudo times of convection and ``nt`` the number of iterations for the computation of characteristics.

The method integrate the partial differential equation

.. math::

   \frac{\partial U}{\partial\t} + V\cdot\nabla U = 0,

on the time intervall :math:`[0, dt]`.

The method used is of Galerkin-Characteristic kind. It is a very simple version which is inconditionnally stable but rather dissipative. See the book of Zienkiewicz and Taylor <<The finite element method>> 5th edition volume 3 : Fluids Dynamics, section 2.6 and also the Freefem++ documentation on convect command.

The defined method works only if ``mf`` is a pure Lagrange finite element method for the moment. The principle is to convect bacward the finite element nodes. This convection is made with ``nt`` steps. Then the solution is interploated on the convected nodes. On the boundary where there is an input convection a simple extrapolation of the field is done on the nearest element.

In order to make the extrapolation not too expensive, the product :math:`dt\ V` should not be too large.

