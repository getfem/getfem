.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-internmm:

Interpolation of a finite element method on non-matching meshes
===============================================================

A special finite element method is defined in
:file:`getfem/getfem_interpolated_fem.h` which is not a real finite element
method, but a pseudo-fem which interpolates a finite element method defined on
another mesh. If you need to assemble a matrix with finite element methods
defined on different meshes, you may use the "interpolated fem" for that
purpose::

  getfem::new_interpolated_fem(getfem::mesh_fem mf, getfem::mesh_im mim);

Because each base function of the finite element method has to be interpolated,
such a computation can be a heavy procedure. By default, the interpolated fem
object store the interpolation data.

The interpolation is made on each Gauss point of the integration methods of
``mim``, so that you have to use these integration methods in the assembling
procedures.

For instance if you need to compute the mass matrix between two different finite
element methods defined on two different meshes, this is an example of code which
interpolate the second FEM. on the mesh of the first FEM., assuming that ``mf``
describes the finite element method and ``mim`` is the chosen integration method::

  getfem::mesh_fem mf_interpole(mfu.linked_mesh());
  pfem ifem = getfem::new_interpolated_fem(mf, mim);
  dal::bit_vector nn = mfu.convex_index();
  mf_interpole.set_finite_element(nn, ifem);
  getfem::asm_mass_matrix(SM1, mim, mfu, mf_interpole);
  del_interpolated_fem(ifem);

The object pointed by ``ifem`` contains all the information concerning the
interpolation. It could use a lot of memory. As pfem is a smart pointer (a boost
`intrusive_ptr <http://www.boost.org/libs/smart_ptr/intrusive_ptr.html>`_), the
interpolated fem will be automatically destroyed when the last pointer on it is
destroyed. To obtain a better accuracy, it is better to refine the integration
method (with ``IM_STRUCTURED_COMPOSITE`` for instance) rather than increase its
order.


mixed methods with different meshes
-----------------------------------
  to be described ...


mortar methods
--------------
  to be described ...
