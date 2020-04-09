.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-computed:

Compute derivatives
===================

The file :file:`getfem/getfem_derivatives.h` defines the following function to
compute the gradient of a solution::

  getfem::compute_gradient(mf1, mf2, U, V);

where ``mf1`` is a variable of type |mf| and describes the finite element method
on which the solution is defined, ``mf2`` describes the finite element method to
compute the gradient, ``U`` is a vector representing the solution and should be
of size ``mf1.nb_dof()``, ``V`` is the vector on which the gradient will be
computed and should be of size ``N * mf2.nb_dof()``, with ``N`` the dimension of
the domain.

.. important:

   This function only works when ``mf2`` is a Lagrange element. This element
   should be, most of the time, a discontinuous Lagrangian element, because for
   usual element (for instance ``getfem::FEM_PK_DISCONTINUOUS(n, k)``), the
   gradient is not continuous.
