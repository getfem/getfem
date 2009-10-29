.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-intermm:

Interpolation on different meshes
=================================

The file :file:`getfem/getfem_interpolation.h` defines the function 
``getfem:interpolation(...)`` to interpolate a solution from a given mesh/finite 
element method on another mesh and/or another Lagrange finite element method::

  getfem::interpolation(mf1, mf2, U, V, extrapolation = 0);

where ``mf1`` is a variable of type |gf_mf| and describes the finite element 
method on which the source field ``U`` is defined, ``mf2`` is the finite element 
method on which ``U`` will be interpolated. ``extrapolation`` is a optional 
parameter. The values are ``0`` not to allow the extrapolation, ``1`` for an 
extrapolation of the exterior points near the boundary and ``1`` for the 
extrapolation of all exterior points (could be expensive).


The dimension of ``U`` should be a multiple of ``mf1.nb_dof()``, and the 
interpolated data ``V`` should be correctly sized (multiple of ``mf2.nb_dof()``).

... important::

    ``mf2`` should be of Lagrange type for the interpolation to make sense but the
    meshes linked to ``mf1`` and ``mf2`` may be different (and this is the
    interest of this function). There is no restriction for the dimension of the
    domain (you can interpolate a 2D mesh on a line etc.).

If you need to perform more than one interpolation between the same finite element
methods, it might be more efficient to use the function::

  getfem::interpolation(mf1, mf2, M, extrapolation = 0);

where ``M`` is a row matrix which will be filled with the linear map representing
the interpolation (i.e. such that ``V = MU``). The matrix should have the correct
dimensions (i.e. ``mf2.nb_dof()``x``mf1.nb_dof()``). Once this matrix is built,
the interpolation is done with a simple matrix multiplication::

  gmm::mult(M, U, V);
