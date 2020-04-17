.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-computel2h1:

Compute :math:`L^2` and :math:`H^1` norms
=========================================

The file :file:`getfem/getfem_assembling.h` defines the functions to compute
:math:`L^2` and :math:`H^1` norms of a solution. The following functions compute
the different norms::

  getfem::asm_L2_norm(mim, mf, U, region = mesh_region::all_convexes());
  getfem::asm_H1_semi_norm(mim, mf, U, region = mesh_region::all_convexes());
  getfem::asm_H1_norm(mim, mf, U, region = mesh_region::all_convexes());

where ``mim`` is a |gf_mim| used for the integration, ``mf`` is a |gf_mf| and
describes the finite element method on which the solution is defined, ``U`` is the
vector of values of the solution on each degree of freedom of ``mf`` and ``region`` is an optional parameter which specify the mesh region on which the norm is computed. The size of
``U`` should be ``mf.nb_dof()``.

In order to compare two solutions, it is often simpler and faster to use the
following function than to interpolate one |mf| on another::

  getfem::asm_L2_dist(mim, mf1, U1, mf2, U2, region = mesh_region::all_convexes());
  getfem::asm_H1_dist(mim, mf1, U1, mf2, U2, region = mesh_region::all_convexes());

These functions return the :math:`L^2` and :math:`H^1` norms of :math:`u_1-u_2`.
