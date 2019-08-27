.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc_levelset:

Level-set module
----------------

Description
^^^^^^^^^^^

Define level-set objects and cut meshes, integration method and finite element method with respect to one or several level-set.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_level_set.h` and :file:`getfem_level_set.cc`, "Define a level-set function (scalar field defined on a Lagrange fem) with an optional secondary level-set function."
   :file:`getfem_mesh_level_set.h` and :file:`getfem_mesh_level_set.cc`, "Cut a mesh with respect to one or several level-sets."
   :file:`getfem_fem_level_set.h` and :file:`getfem_fem_level_set.cc`, "Define a special finite element method which depends on the element and which is cut by one or several level-sets."
   :file:`getfem_mesh_fem_level_set.h` and :file:`getfem_mesh_fem_level_set.cc`, "Produces a mesh_fem object with shape functions cut by one or several level-sets."
   :file:`getfem_mesh_im_level_set.h` and :file:`getfem_mesh_im_level_set.cc`, "Produce a mesh_im representing an intergration method cut by the level set and being on on side of level-set, the oter side, both or only on the levelset itself."
   :file:`getfem_level_set_contact.h` and :file:`getfem_level_set_contact.cc`, "A level set based large sliding contact algorithm for an easy analysis of implant positioning."
   :file:`getfem_convect.h`, "Compute the convection of a quantity with respect to a vector field. Used to computate the evolution of a level-set function for instance. Galerkin characteristic method."

State
^^^^^

Stable.



Perspectives
^^^^^^^^^^^^

Clarify the algorithm computing the different zones.

