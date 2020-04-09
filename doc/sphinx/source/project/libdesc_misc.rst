.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc_misc:


Miscellaneous algorithms
------------------------

Description
^^^^^^^^^^^

A set of miscellaneous basic algorithms and definitions used in |gf|.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15


   :file:`bgeot_comma_init.h`, "Allow to init  container with a list of values, from boost init.hpp."
   :file:`bgeot_ftool.h` and :file:`bgeot_ftool.cc`, "Small language allowing to read a parameter file with a Matlab syntax like. Used also for structured meshes."
   :file:`bgeot_kdtree.h` and :file:`bgeot_kdtree.cc`, "Balanced N-dimensional tree. Store a list of points and allows a quick search of points lying in a given box."
   :file:`bgeot_rtree.h` and :file:`bgeot_rtree.cc`, "Rectangle tree. Store a list of N-dimensional rectangles and allows a quick search of rectangles containing a given point."
   :file:`permutations.h`, "Allows to iterate on permutations. Only used in :file:`getfem_integration.cc`."
   :file:`bgeot_small_vector.h` and :file:`bgeot_small_vector.cc`, "Defines a vector of low dimension mainly used to represent mesh nodes. Optimized operations."
   :file:`bgeot_tensor.h`, "Arbitrary order tensor. Used in assembly."
   :file:`bgeot_sparse_tensors.h` and :file:`bgeot_sparse_tensors.cc`, "Arbitrary order sparse tensor. Used in the low-level generic assembly."
   :file:`getfem_omp.h` and :file:`getfem_omp.cc`, "Tools for multithreaded, OpenMP and Boost based parallelization."
   :file:`getfem_export.h` and :file:`getfem_export.cc`, "Export in pos and vtk formats"
   :file:`getfem_superlu.h` and :file:`getfem_superlu.cc`, "Interface with Superlu (the included version or an external one)"



State
^^^^^



Perspectives
^^^^^^^^^^^^

