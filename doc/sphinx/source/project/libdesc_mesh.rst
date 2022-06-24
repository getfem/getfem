.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc_mesh:


Mesh module
-----------



Description
^^^^^^^^^^^

This part of the library has the role to store and manage the meshes, i.e. a
collection of elements (real elements) connected to each other by some of their
faces. For that, it develops concepts of elements, elements of reference,
structure of meshes, collection of nodes, geometric transformations, subpart of
the boundary or subzone of the mesh.

There is no really effective meshing capabilities available for the moment in
|gf|. The meshes of complex objects must be imported from existing meshers such
as `Gmsh`_ or `GiD`_. Some importing functions of meshes have been written and
can be easily extended for other formats.

The object which represents a mesh declared in the file :file:`getfem_mesh.h` and
which is used as a basis for handling of the meshes in |gf| manages also the
possibility for the structures depending on a mesh (see MESHFEM and MESHIM
modules) to react to the evolution of the mesh (addition or removal of elements,
etc.).


Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`bgeot_convex_structure.h` and :file:`bgeot_convex_structure.cc`, "Describes the structure of an element disregarding the coordinates of its vertices."
   :file:`bgeot_mesh_structure.h` and :file:`bgeot_mesh_structure.cc`, "Describes the structure of a mesh disregarding the coordinates of the nodes."
   :file:`bgeot_node_tab.h` and :file:`bgeot_node_tab.cc`, "A node container allowing the fast search of a node and store nodes identifying the too much close nodes."
   :file:`bgeot_convex.h`, "Describes an element with its vertices."
   :file:`bgeot_convex_ref.h` and :file:`bgeot_convex_ref.cc` and :file:`bgeot_convex_structure.cc`, "Describe reference elements."
   :file:`bgeot_mesh.h`, "Describes a mesh with the collection of node (but without the description of geometric transformations)."
   :file:`getfem_mesh_region.h` and :file:`getfem_mesh_region.cc`, "Object representing a mesh region (boundary or part of a mesh)."
   :file:`bgeot_geometric_trans.h` and :file:`bgeot_geometric_trans.cc`, "Describes geometric transformations."
   :file:`bgeot_geotrans_inv.h` and :file:`bgeot_geotrans_inv.cc`, "A tool to invert geometric transformations."
   :file:`getfem_mesh.h` and :file:`getfem_mesh.cc`, "Fully describes a mesh (with the geometric transformations, subparts of the mesh, support for parallelization). Includes the Bank algorithm to refine a mesh."
   :file:`getfem_deformable_mesh.h`, "defines an object capable to deform a mesh with respect to a displacement field and capable to restore it"
   :file:`getfem_mesher.h` and :file:`getfem_mesher.cc`, "An experimental mesher, in arbitrary dimension. To be used with care and  quite slow (because of node optimization). It meshes geometries defined by some level sets."
   :file:`getfem_import.h` and :file:`getfem_import.cc`, "Import mesh files in various formats"
   :file:`getfem_regular_meshes.h` and :file:`getfem_regular_meshes.cc`, "Produces structured meshes"
   :file:`getfem_mesh_slicers.h` and :file:`getfem_mesh_slicers.cc`, "A slice is built from a mesh, by applying some slicing operations (cut the mesh with a plane, intersect with a sphere, take the boundary faces, etc..). They are used for post-treatment (exportation of results to VTK or OpenDX,  etc.)."
   :file:`getfem_mesh_slice.h` and :file:`getfem_mesh_slice.cc`, "Store mesh slices."


State
^^^^^

Stable and not evolving so much.

Perspectives
^^^^^^^^^^^^

For the moment, the module is split into two parts which lie into two different
namespaces. Of course, It would be more coherent to gather the module in only one
namespace (``getfem``).

.. note::

   The file :file:`bgeot_mesh.h` could be renamed :file:`getfem_basic_mesh.h`.

A  bibliographical review on how to efficiently store a mesh and implement the main operations (add a node, an element, deal with faces, find the neighbour elements, the isolated faces ...) would be interesting to make the mesh structure evolve.

A senstive algorithm is the one (in bgeot_node_tab.cc) which identify the too much close nodes. More investigations (and documentation) are probably necessary.




