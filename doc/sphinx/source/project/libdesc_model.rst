.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc-model:

Model module
------------


Description
^^^^^^^^^^^

Describe a model (variable, data and equation terms linking the variables).

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_models.h` and :file:`getfem_models.cc`, "Defines the object models, its internal and the standard bricks (linear elasticity, generic elliptic brick, Dirichlet boundary conditions ...)."
   :file:`getfem_model_solvers.h` and :file:`getfem_model_solvers.cc`, "Defines the the standard solvers for the model object."
   :file:`getfem_contact_and_friction_common.h` and :file:`getfem_contact_and_friction_common.cc`, "Common algorithms for contact/friction conditions on deformable bodies"
   :file:`getfem_contact_and_friction_integral.h` and :file:`getfem_contact_and_friction_integral.cc`, "Small sliding contact/friction bricks of integral type."
   :file:`getfem_contact_and_friction_large_sliding.h` and :file:`getfem_contact_and_friction_large_sliding.cc` , "Large sliding contact/friction bricks."
   :file:`getfem_contact_and_friction_nodal.h` and :file:`getfem_contact_and_friction_nodal.cc`, "Small sliding nodal contact/friction bricks."
   :file:`getfem_Navier_Stokes.h` , "An attempt for Navier-Stokes bricks. To be improved."
   :file:`getfem_fourth_order.h` and :file:`getfem_fourth_order.cc` , "Bilaplacian and Kirchhoff-Love plate bricks"
   :file:`getfem_linearized_plates.h` and :file:`getfem_linearized_plates.cc` , "Mindlin-Reissner plate brick"
   :file:`getfem_nonlinear_elasticity.h` and :file:`getfem_nonlinear_elasticity.cc` , "Large deformation elasticity bricks."
   :file:`getfem_plasticity.h` and :file:`getfem_plasticity.cc` , "Plasticity bricks."

State
^^^^^

Constant evolution to includes next models.

Perspectives
^^^^^^^^^^^^

More plate, load and shell bricks, plasticity in large deformation, ...
