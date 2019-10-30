.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc_low_assemb:



The low-level generic assembly module in |gf|
---------------------------------------------

Description
^^^^^^^^^^^

First version of the generic assembly. Base on tensor reduction. Not very convenient for nonlinear terms. The high-level generic assembly have to be prefered now.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_mat_elem_type.h` and :file:` getfem_mat_elem_type.cc, "Defines base type for components of an elementary matrix."
   :file:`getfem_mat_elem.h` and :file:` getfem_mat_elem.cc, "Describes an compute elementary matrices."
   :file:`getfem_assembling_tensors.h` and :file:`getfem_assembling_tensors.cc`, "Performs the assembly."
   :file:`getfem_assembling.h`, "Various assembly terms (linear elasticity, generic elliptic term, Dirichlet condition ..."

State
^^^^^

Stable.

Perspectives
^^^^^^^^^^^^

Will not evolve since the efforts are now focused on the high-level generic assembly.
