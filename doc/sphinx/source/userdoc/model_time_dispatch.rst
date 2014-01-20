.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-time-dispatchers:



The time dispatchers: integration of transient problems
-------------------------------------------------------

The role of time dispatchers is to allow the integration of transient problems 
with some pre-defined time integration schemes. The principle of the time 
dispatchers is to dispatch the terms of a brick on the different time steps of the 
considered time integration scheme. When time derivative terms are present in the 
model (this is generally the case except for quasistatic models), the time 
dispatcher will be associated to a specific brick representing this time 
derivative term (:math:`\partial u / \partial t` or :math:`\partial^2 u / \partial 
t^2` for instance). For this, a number of tools are available in |gf| to help the 
construction of a time dispatcher. Mainly they are the two following:

* The variables can be duplicated to take into account the different versions 
  corresponding to each time iteration. For instance, for simplest time 
  integration schemes, two versions :math:`U^n` and :math:`U^{n+1}` of a variable 
  :math:`U` are stored. The addition of a variable :math:`u` with two versions can 
  be done with the method of the model object::

    model.add_fem_variable("u", mf_u, 2);

  where :math:`2` is here the number of versions. The variable which is actually 
  computed have always the index 0 and will be accessed with 
  ``model.real_variable("u", 0)`` or simply with ``model.real_variable("u")``. It 
  will generally represent the version :math:`U^{n+1}`. The version :math:`U^{n}` 
  (corresponding to the previous time step) will be accessed with 
  ``model.real_variable("u", 1)``. Generally, it will be necessary to set this 
  version with ``model.set_real_variable("u", 1)`` to define the initial condition 
  of the model. At the end of each iteration, the different versions of a variable 
  are automatically shifted (version 0 :math:`\rightarrow` version 1 ...).

* The right hand side of a brick is dispatched into several right hand sides for 
  each time iteration which are stored. To avoid unnecessary computation, the time 
  dispatcher can shift these extra right hand sides at the end of each time 
  iteration.


Theta-method dispatcher
-----------------------

This is the simplest time dispatcher. The use of this dispatcher will be described 
in details. Since the use of the other dispatchers is similar, only their 
specificities will be described later on.

The principle of the :math:`\theta`-method is to dispatch the term :math:`F` into 
:math:`(\theta) F^{n+1} + (1-\theta) F^{n},`

For specific values of :math:`\theta` one obtains some classical schemes: backward 
Euler for :math:`\theta = 1`, forward Euler for :math:`\theta = 0` and 
Crank-Nicolson scheme for :math:`\theta = 1/2` (which is an order two scheme).

For instance, if the dispatcher is applied to a brick representing a linear 
elliptic term :math:`KU` where :math:`K` is the stiffness matrix and :math:`U` the 
unknown, it will be transformed into :math:`(\theta) KU^{n+1} + (1-\theta) 
KU^{n}`.

Since :math:`U^{n+1}` is the real unknown, the effect will be to multiply by 
:math:`\theta` the stiffness matrix and to add to the right hand side the term 
:math:`(1-\theta) KU^{n}`. This means also that :math:`U^{n}` have to be 
initialized (with something like ``gmm::copy(U0, model.real_variable("u",1))``). 
It represents an initial data for the problem. Remember this principle: each time 
you apply a time dispatcher to a brick, the corresponding variables have to have 
the right number of versions (see above) and should be initialized before the 
first time iteration.

You can apply the dispatcher to a brick having only a right hand side (a source 
term for instance). It is not necessary if the term is constant in time.

When a brick represents a constraint (Dirichlet condition, incompressibility ...) 
this is not mandatory to apply the dispatcher. Of course, the result will not be 
exactly the same if you apply or not the dispatcher. If you do not apply it, the 
constraint will be applied to the current variable (:math:`U^{n+1}` for the 
:math:`\theta`-method). If you apply it, the constraint will be in a sense applied 
to :math:`(\theta) U^{n+1} + (1-\theta) U^{n}`. If the constraint is applied 
thanks to a multiplier, this multiplier will need to have different versions and 
will need to have an initial value.

In order to apply the :math:`\theta`-method dispatcher to a set of brick you must 
execute::

  model.add_initialized_scalar_data("theta", theta);
  getfem::add_theta_method_dispatcher(model, transient_bricks, "theta");

where ``transient_bricks`` is a ``dal::bit_vector`` containing the indices of the 
corresponding bricks. The value of :math:`\theta` can be modified from an 
iteration to another.

The global structure of the loop solving the different time steps should be the 
following::

  gmm::iteration solver_iter(residual, 0, 40000);

  // Set here the initial values.

  model.first_iter(); // initialize the iterations.

  for (scalar_type t = 0; t < T; t += dt) {

    solver_iter.init();
    getfem::standard_solve(model, solver_iter); // solve an iteration.

    model.next_iter(); // shift the variables and additional right hand sides.
  }

where ``model.first_iter()`` should be called before the first iteration to 
initialize the right hand side of the time dispatchers. The initial data should be 
set before the call to ``model.first_iter()``. The method ``model.next_iter()`` is 
to be called at the end of each iteration. It calls the dispatcher to shift there 
additional right hand side and it shifts the version of the variables.


Basic first order time derivative brick
+++++++++++++++++++++++++++++++++++++++

A term like :math:`\rho \partial u / \partial t` will be represented in the model 
by :math:`(MU^{n+1} - MU^{n}) / dt`, where :math:`M` is the mass matrix and 
:math:`dt` is the time step. The :math:`\theta`-method is compatible with this. A 
brick is dedicated to represent this term. It can be added to the model by the 
function::

  getfem::add_basic_d_on_dt_brick(model, mim, varname, dataname_dt,
                                  dataname_rho = std::string(),
                                  region = size_type(-1));

where ``varname`` is the name of the variable on which the time derivative is 
applied (should have at least two versions), ``dataname_dt`` is the name of the 
data corresponding to the time step (added by 
``model.add_initialized_scalar_data("dt", dt)`` for instance) which could be 
modified from an iteration to another and ``dataname_rho`` is an optional 
parameter (whose default value is 1) corresponding to the term :math:`\rho` in 
:math:`\rho \partial u / \partial t`.

NOTE that the time dispatcher should not be applied to this brick !

A good model of the use of this brick and the :math:`\theta`-method time 
dispatcher can be found in the test program ``tests/heat_equation.cc``.


Basic second order time derivative brick
++++++++++++++++++++++++++++++++++++++++

This brick represents a second order time derivative like :math:`\rho \partial^2 u 
/ \partial t^2`. The problem with such a term is that the :math:`\theta`-method 
should be applied both on :math:`u` and :math:`\partial u / \partial t` which 
means that :math:`\partial u / \partial t` is a natural unknown of the problem. 
The easiest way is then to add the time derivative of the variable :math:`u` has an independent variables of the model (a drawback, of course, is that one has twice 
as much unknowns). This Basic second order time derivative brick does not apply 
this strategy. The time derivative :math:`\partial u / \partial t` is considered 
as a data which is updated at a post-treatment stage (in some cases, this strategy 
cannot be applied if the time derivative appears to be a required unknown of the 
model).

The term :math:`\rho \partial^2 u / \partial t^2` will be represented by 
:math:`(MU^{n+1} - MU^{n}) / (\alpha dt^2) - M V^n / (\alpha dt) ~~~~~~~~(*)`, 
where :math:`M` is the mass matrix, :math:`dt` is the time step, :math:`\alpha` is 
a parameter which is equal to :math:`\theta` for the :math:`\theta`-method and 
:math:`V^n` the time derivative at the previous time step. This means in 
particular that :math:`V` should be added as a data on the model with (at least) 
two versions.

The function adding the brick is::

  getfem::add_basic_d2_on_dt2_brick(model, mim, varname, dataname_V,
             dataname_dt, dataname_alpha, dataname_rho = std::string(),
             region = size_type(-1));

where ``varname`` is the name of the variable on which the second order time 
derivative is applied, ``dataname_V`` is the data representing the time 
derivative, ``dataname_dt`` is the name of the data corresponding to the time step 
(added by ``model.add_initialized_scalar_data("dt", dt)`` for instance) which 
could be modified from an iteration to another, ``dataname_alpha`` is the name of 
the data containing the parameter :math:`\alpha` in (*) and ``dataname_rho`` is an 
optional parameter (whose default value is 1) corresponding to the term 
:math:`\rho` in :math:`\rho \partial^2 u / \partial t^2`.

At the end of each iteration, the data ``dataname_V`` should be updated (before 
the call to ``model.next_iter()`` by the call to::

  getfem::velocity_update_for_order_two_theta_method
      (model, varname, dataname_V, dataname_dt, dataname_alpha);

A good model of the use of this brick and the :math:`\theta`-method time 
dispatcher can be found in the test program ``tests/wave_equation.cc``.


Midpoint dispatcher
-------------------

The principle of the midpoint scheme is to dispacth a term :math:`F(U)` into 
:math:`F((U^{n+1}+U^{n})/2),`

It is different from the Crank-Nicolson scheme (:math:`\theta`-method for 
:math:`\theta=1/2`) only for nonlinear terms.

The real unknown remains :math:`U^{n+1}`. the effect will be to multiply by 
:math:`1/2` the stiffness (or tangent) matrix and to add to a right hand side the 
term :math:`(KU^{n}/2` for a linear matrix term :math:`K`. As for the 
:math:`\theta`-method, the variables have to have two version and the second 
version have to be initialized.

You can apply the dispatcher to a brick having only a right hand side (a source 
term for instance). It is not necessary if the term is constant in time.

NOTE that if the brick depend on a data which is not constant in time, the data 
either have to have to versions (and the mean of the two versions are taken into 
account) or evaluated at the middle of the time step.

When a brick represents a constraint (Dirichlet condition, incompressibility ...) 
this is not mandatory to apply the dispatcher. Of course, the result will not be 
exactly the same if you apply or not the dispatcher. If you do not apply it, the 
constraint will be applied to the current variable :math:`U^{n+1}`. If you apply 
it, the constraint will be applied to :math:`(U^{n+1} + U^{n})/2`. If the 
constraint is applied thanks to a multiplier, this multiplier will need to have 
different versions and will need to have an initial value.

In order to apply the midpoint dispatcher to a set of brick you must execute::

  getfem::add_midpoint_dispatcher(model, transient_bricks);

where ``transient_bricks`` is a ``dal::bit_vector`` containing the indices of the 
corresponding bricks.


Basic first order time derivative brick
+++++++++++++++++++++++++++++++++++++++

The same brick as for the :math:`\theta`-method can be used to represent a first 
order time derivative.


Basic second order time derivative brick
++++++++++++++++++++++++++++++++++++++++

The same brick as for the :math:`\theta`-method can be used to represent a second 
order time derivative. The value of :math:`\alpha` should be :math:`1/2`.


Newmark scheme
--------------

For a system

.. math::

   M\ddot{U} + K(U) = F,

the Newmark scheme of parameter :math:`\beta` and :math:`\gamma` is defined by

.. math::

   M(U^{n+1} - U^{n}) = dt M V^n + dt^2/2( 2\beta(F^{n+1}-K(U^{n+1})) + (1-2\beta)(F^{n}-K(U^{n}))),\\
   M(V^{n+1} - V^{n}) = dt ( 2\gamma(F^{n+1}-K(U^{n+1})) + (1-2\gamma)(F^{n}-K(U^{n}))),

where :math:`V` represents the time derivative of :math:`U`.

The implementation of the Newmark scheme proposed is not optimal and should be 
adapted. It can be optained using the basic second order time derivative brick 
(see :math:`\theta`-method) and the :math:`\theta`-method time dispatcher used 
with :math:`\theta = 2\beta`. Additionally, one has to use the following function 
which computes the time derivative of the variable as a post-computation::

  getfem::velocity_update_for_Newmark_scheme
      (model, id2dt2, varname, dataname_V, dataname_dt, dataname_alpha);

where ``id2dt2`` is the index of the basic second order time derivative brick (see 
the section on the :math:`\theta`-method for more details and the implementation 
in the test program ``tests/wave_equation.cc``).

This implementation of the Newmark-scheme is not optimal since the latter function 
inverts the mass matrix to compute the time derivative using a conjugate gradient. 
This linear system solve could be avoided by keeping the multiplication of the 
mass matrix with the time derivative as a data, with an adaptation of the time 
derivative brick.

