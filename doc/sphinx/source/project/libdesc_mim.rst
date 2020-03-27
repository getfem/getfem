.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc-mim:


MeshIm module
-------------

Description
^^^^^^^^^^^

Defines an integration method on a whole mesh.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_mesh_im.h` and :file:`getfem_mesh_im.cc`, "Object which defines an integration method on each element of the mesh. Reacts to the main mesh changes (add or deletion of elements)."
   :file:`getfem_im_data.h` and :file:`getfem_im_data.cc`, "Define an object representing a scalar, a vector or a tensor on each Gauss point of a mesh_im object. Used for instance in plasticity approximation. Interpolation of arbitrary expressions can be made thanks to the weak form language."


State
^^^^^

Stable, not evolving so much.

Perspectives
^^^^^^^^^^^^

