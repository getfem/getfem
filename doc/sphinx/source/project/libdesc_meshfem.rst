.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc-mesh-fem:


MeshFem module
--------------

Description
^^^^^^^^^^^

The MeshFem module aims to represent a finite element method (space) with respect to a given mesh. The mesh_fem object will be permanently linked to the given mesh and will be able to react to  changes in the mesh (addition or deletion of elements, in particular). A mesh_fem object may associate a different finite element method on each element of the mesh even though of course, the most common case it that all the element share the same finite element method.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15


   :file:`getfem_mesh_fem.h` and :file:`getfem_mesh_fem.cc`, "Defines the structure representing a finite element on a whole mesh. Each element of the mesh is associated with a finite element method. This is a quite complex structure which perform dof identification and numbering, allows a global linear reduction."
   :file:`getfem_mesh_fem_global_function.h` and :file:`getfem_mesh_fem_global_function.cc`, "Defines mesh_fem with fem defined as a fem_global_function. It provides convenience methods for updating the list of base functions in the linked fem_global_function."
   :file:`getfem_mesh_fem_product.h` and :file:`getfem_mesh_fem_product.cc`, "Produces a mesh_fem object which is a kind of direct product of two finite element method. Useful for Xfem enrichment."
   :file:`getfem_mesh_fem_sum.h` and :file:`getfem_mesh_fem_sum.cc`, "Produces a mesh_fem object which is a kind of direct sum of two finite element method. Useful for Xfem enrichment."
   :file:`getfem_partial_mesh_fem.h` and :file:`getfem_partial_mesh_fem.cc`, "Produces a mesh_fem with a reduced number of dofs"
   :file:`getfem_interpolation.h` and :file:`getfem_interpolation.cc`, "Interpolation between two finite element methods, possibly between different meshes. The interpolation facilities of the high-level generic assembly can be used instead."
   :file:`getfem_derivatives.h`, "Interpolation of some derivatives of a finite element field on a (discontinuous) Lagrange finite element. The interpolation facilities of the high-level generic assembly can be used instead."
   :file:`getfem_inter_element.h` and :file:`getfem_inter_element.cc`, "An attempt to make framework for inter-element computations (jump in normal derivative for instance). To be continuated and perhaps integrated into the generic assembly language."
   :file:`getfem_error_estimate.h` and :file:`getfem_error_estimate.cc`, "An attempt to make framework for computation of error estimates. To be continuated and perhaps integrated into the generic assembly language."
   :file:`getfem_crack_sif.h`, "Crack support functions for computation of SIF(stress intensity factors)."
   :file:`getfem_torus.h` and :file:`getfem_torus.cc`, "Adapt a mesh_fem object which extends a 2D dimensional structure with a radial dimension."

State
^^^^^

Stable. Not evolving so much.

Perspectives
^^^^^^^^^^^^

Parallelisation of dof numbering to be done. An optimal (an simple) algorithm
exists.


