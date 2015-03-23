.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-linear-elasticity:


Isotropic linearized elasticity brick
-------------------------------------

This brick represents a term

.. math::

   -div(\sigma) = \ldots

with

.. math::

   \sigma &= \lambda\mbox{tr}(\varepsilon(u))I + 2\mu\varepsilon(u) \\
   \varepsilon(u) &= (\nabla u + \nabla u^T)/2

:math:`\varepsilon(u)` is the small strain tensor, :math:`\sigma` is the stress
tensor, :math:`\lambda` and :math:`\mu` are the Lamé coefficients. This represents
the system of linearized isotropic elasticity. It can also be used with
:math:`\lambda=0` together with the linear incompressible brick to build the
Stokes problem.

The function which adds this brick to a model is::

  ind_brick = getfem::add_isotropic_linearized_elasticity_brick
              (md, mim, varname, dataname_lambda, dataname_mu,
               region = size_type(-1));

where ``dataname_lambda`` and ``dataname_mu`` are the data of the model
representing the Lamé coefficients (constant or described on a finite element
method).

The function::

  getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
    (md, varname, dataname_lambda, dataname_mu, mf_vm, VM, tresca_flag = false);

compute the Von Mises criterion (or Tresca if ``tresca_flag`` is set to true) on
the displacement field stored in ``varname``. The stress is evaluated on the |mf|
``mf_vm`` and stored in the vector ``VM``.

The program :file:`tests/elastostatic.cc` can be taken as a model of use of this
brick.


Linear incompressibility (or nearly incompressibility) brick
------------------------------------------------------------

This brick adds a linear incompressibility condition (or a nearly incompressible
condition) in a problem of type:

.. math::

   \mbox{div}(u) = 0,\quad (\mbox{ or } \mbox{div}(u) = \varepsilon p)

This constraint is enforced with Lagrange multipliers representing the pressure,
introduced in a mixed formulation.

The function adding this incompressibility condition is::

  ind_brick = getfem::add_linear_incompressibility
              (md, mim, varname, multname_pressure, region = size_type(-1),
               dataname_penal_coeff = std::string());

where ``varname`` is the variable on which the incompressibility condition is
prescribed, ``multname_pressure`` is a variable which should be described on a
scalar fem representing the multiplier (the pressure) and ``dataname_penal_coeff``
is an optional penalization coefficient (constant or described on a finite element
method) for the nearly incompressible condition.

In nearly incompressible homogeneous linearized elasticity, one has
:math:`\varepsilon = 1 / \lambda` where :math:`\lambda` is one of the Lamé
coefficient and :math:`\varepsilon` the penalization coefficient.

For instance, the following program defines a Stokes problem with a source term
and an homogeneous Dirichlet condition on boundary 0. ``mf_u``, ``mf_data`` and
``mf_p`` have to be valid finite element description on the same mesh. ``mim``
should be a valid integration method on the same mesh::

  typedef std::vector<getfem::scalar_type> plain_vector;
  size_type N = mf_u.linked_mesh().dim();

  getfem::model Stokes_model;

  laplacian_model.add_fem_variable("u", mf_u);

  getfem::scalar_type mu = 1.0;
  Stokes_model.add_initialized_data("lambda", plain_vector(1, 0.0));
  Stokes_model.add_initialized_data("mu", plain_vector(1, mu));

  getfem::add_isotropic_linearized_elasticity_brick(Stokes_model, mim,
                                                    "u", "lambda", "mu");

  laplacian_model.add_fem_variable("p", mf_p);
  getfem::add_linear_incompressibility(Stokes_model, mim, "u", "p");

  plain_vector F(mf_data.nb_dof()*N);
  for (int i = 0; i < mf_data.nb_dof()*N; ++i) F(i) = ...;
  Stokes_model.add_initialized_fem_data("VolumicData", mf_data, F);
  getfem::add_source_term_brick(Stokes_model, mim, "u", "VolumicData");

  getfem::add_Dirichlet_condition_with_multipliers(Stokes_model, mim,
                                                   "u", mf_u, 1);

  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(Stokes_model, iter);

  plain_vector U(mf_u.nb_dof());
  gmm::copy(Stokes_model.real_variable("u"), U);

An example for a nearly incompressibility condition can be found in the program
:file:`tests/elastostatic.cc`.

