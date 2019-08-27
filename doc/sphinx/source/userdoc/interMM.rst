.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: none

.. _ud-intermm:

Interpolation of arbitrary quantities
=====================================

Once a solution has been computed, it is quite easy to extract any quantity of interest on it with the interpolation functions for instance for post-treatment.

Basic interpolation
*******************

The file :file:`getfem/getfem_interpolation.h` defines the function
``getfem::interpolation(...)`` to interpolate a solution from a given mesh/finite
element method on another mesh and/or another Lagrange finite element method::

  getfem::interpolation(mf1, mf2, U, V, extrapolation = 0);

where ``mf1`` is a variable of type |gf_mf| and describes the finite element
method on which the source field ``U`` is defined, ``mf2`` is the finite element
method on which ``U`` will be interpolated. ``extrapolation`` is an optional
parameter. The values are ``0`` not to allow the extrapolation, ``1`` for an
extrapolation of the exterior points near the boundary and ``2`` for the
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


Interpolation based on the high-level weak form language
********************************************************

It is possible to extract some arbitrary expressions on possibly several fields thanks to the weak form language and the interpolation functions.

This is specially dedicated to the model object (but it can also be used with a ga_workspace object). For instance if ``md`` is a valid object containing some defined variables ``u`` (vectorial) and ``p`` (scalar), one can interpolate on a Lagrange finite element method an expression such as ``p*Trace(Grad_u)``. The resulting expression can be scalar, vectorial or tensorial. The size of the resulting vector is automatically adapted.

The high-level generic interpolation functions are defined in the file :file:`getfem/getfem_generic_assembly.h`.

There is different interpolation functions corresponding to the interpolation on a Lagrange fem on the same mesh, the interpolation on a cloud on points or on a ``getfem::im_data`` object.

Interpolation on a Lagrange fem::

  void getfem::ga_interpolation_Lagrange_fem(workspace, mf, result);

where ``workspace`` is a ``getfem::ga_workspace`` object which aims to store the different variables and data (see  :ref:`ud-gasm-high`), ``mf`` is the ``getfem::mesh_fem`` object reresenting the Lagrange fem on which the interpolation is to be done and ``result`` is a ``beot::base_vector`` which store the interpolatin. Note that the workspace should contain the epression to be interpolated. ::

  void getfem::ga_interpolation_Lagrange_fem(md, expr, mf, result, rg=mesh_region::all_convexes());

where ``md`` is a ``getfem::model`` object (containing the variables and data), ``expr`` (std::string object) is the expression to be interpolated, ``mf`` is the ``getfem::mesh_fem`` object reresenting the Lagrange fem on which the interpolation is to be done, ``result`` is the vector in which the interpolation is stored and ``rg`` is the optional mesh region.

Interpolation on a cloud of points::

  void getfem::ga_interpolation_mti(md, expr, mti, result, extrapolation = 0, rg=mesh_region::all_convexes(), nbpoints = size_type(-1));

where ``md`` is a ``getfem::model`` object (containing the variables and data), ``expr`` (std::string object) is the expression to be interpolated, ``mti`` is a ``getfem::mesh_trans_inv`` object which stores the cloud of points (see :file:`getfem/getfem_interpolation.h`), ``result`` is the vector in which the interpolation is stored, ``extrapolation`` is an option for extrapolating the field outside the mesh for outside points, ``rg`` is the optional mesh region and ``nbpoints`` is the optional maximal number of points.

Interpolation on an im_data object (on the Gauss points of an integration method)::

  void getfem::ga_interpolation_im_data(md, expr, im_data &imd,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());

where ``md`` is a ``getfem::model`` object (containing the variables and data), ``expr`` (std::string object) is the expression to be interpolated, ``imd`` is a ``getfem::im_data`` object which refers to a integration method (see :file:`getfem/getfem_im_data.h`), ``result`` is the vector in which the interpolation is stored and ``rg`` is the optional mesh region.
