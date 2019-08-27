.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model:

============================================
The model description and basic model bricks
============================================


The model description of |gf| allows
to quickly build some fem applications on complex linear or nonlinear PDE coupled
models. The principle is to propose predefined bricks which can be assembled to
describe a complex situation. A brick can describe either an equation (Poisson
equation, linear elasticity ...) or a boundary condition (Dirichlet, Neumann ...)
or any relation between two variables. Once a brick is written, it is possible to
use it in very different situations. This allows a reusability of the produced
code and the possibility of a growing library of bricks. An effort has been made in
order to facilitate as much as possible the definition of a new brick. A brick is
mainly defined by its contribution in the tangent linear system to be solved.

This model description is an evolution of the model bricks of previous versions of
|gf|. Compared to the old system, it is more flexible, more general, allows the
coupling of model (multiphysics) in a easier way and facilitates the writing of new
components. It also facilitate the write of time integration schemes for evolving
PDEs.

The kernel of the model description is contained in the file
:file:`getfem/getfem_models.h`. The two main objects are the |mo| and the |br|.



.. toctree::
   :maxdepth: 2

   model_object
   model_generic_assembly
   model_generic_elliptic
   model_dirichlet
   model_source_term
   model_solvers
   model_poisson
   model_Nitsche
   model_constraint
   model_explicit
   model_helmholtz
   model_fourier_robin
   model_linear_elasticity
   model_mass
   model_bilaplacian
   model_Mindlin_plate
   model_time_integration
   model_contact_friction
   model_contact_friction_large_sliding

