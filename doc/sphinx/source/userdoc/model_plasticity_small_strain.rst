.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-plasticity-small-strain:



Small strain plasticity
-----------------------

A framework for the approximation of plasticity models in |gf|. See in :file:`src/getfem_plasticity.cc` and :file:`interface/src/gf_model_set.cc` for the brick implementation and to extend the implementation to new plasticity models.


Theoretical background
++++++++++++++++++++++

We present a short introduction to small strain plasticity. We refer mainly to [SI-HU1998]_ and [SO-PE-OW2008]_ for a more detailed presentation.

Additive decomposition of the small strain tensor
=================================================

Let :math:`\Omega \subset \R^3` be the reference configuration of a deformable body and :math:`u : \Omega \rightarrow \R^3` be the displacement field. Small strain plasticity is based on the additive decomposition of the small strain tensor :math:`\varepsilon(u) = \Frac{\nabla u + \nabla u^T}{2}` in

.. math::
   \varepsilon(u) = \varepsilon^e + \varepsilon^p

where :math:`\varepsilon^e` is the elastic part of the strain tensor and :math:`\varepsilon^p` the plastic one.

Internal variables, free energy potential and elastic law
=========================================================

We consider

.. math::
   \alpha : \Omega \rightarrow \R^{d_{\alpha}},

a vector field of :math:`d_{\alpha}` strain type internal variables (:math:`d_{\alpha} = 0` if no internal variables are considered). We consider also a free energy potential

.. math::
   \psi(\varepsilon^e, \alpha),

such that corresponding stress type variables are determined by

.. math::
   \sigma = \Frac{\partial \psi}{\partial \varepsilon^e}(\varepsilon^e, \alpha), ~~~~ A =  \Frac{\partial \psi}{\partial \alpha}(\varepsilon^e, \alpha),

where :math:`\sigma` is the Cauchy stress tensor and :math:`A` the stress type internal variables. The plastic dissipation is given by

.. math::

   \sigma:\dot{\varepsilon}^p - A.\dot{\alpha} \ge 0.

In the standard cases, :math:`\psi(\varepsilon^e, \alpha)` is decomposed into

.. math:: \psi(\varepsilon^e, \alpha) = \psi^e(\varepsilon^e) + \psi^p(\alpha).

In the case of linearized elasticity, one has :math:`\psi^e(\varepsilon^e) = \frac{1}{2} ({\cal A}\varepsilon^e) :\varepsilon^e` where :math:`{\cal A}` is the fourth order elasticity tensor. For isotropic linearized elasticity this expression reduces to :math:`\psi^e(\varepsilon^e) = \mu \mbox{dev}(\varepsilon^e) : \mbox{dev}(\varepsilon^e) + \frac{1}{2} K (\mbox{tr}(\varepsilon^e))^2` where :math:`\mu` is the shear modulus and :math:`K = \lambda + 2\mu/3` is the bulk modulus.



Plastic potential, yield function and plastic flow rule
=======================================================

Plastic yielding is supposed to occur when the stress attains a critical value. This limit is determined by a yield function :math:`f(\sigma, A)` and the condition

.. math:: f(\sigma, A) \le 0.

The surface :math:`f(\sigma, A) = 0` is the yield surface where the plastic deformation may occur.

Let us also consider the plastic potential :math:`\Psi(\sigma, A)`, (convex with respect to its both variables) which determines the plastic flow direction in the sense that the flow rule is defined as

.. math:: \dot{\varepsilon}^p = \gamma \Frac{\partial \Psi}{\partial \sigma}(\sigma, A), ~~~~~~ \dot{\alpha} = -\gamma \Frac{\partial \Psi}{\partial A}(\sigma, A),

with the additional complementarity condition

.. math:: f(\sigma, A) \le 0, ~~~ \gamma \ge 0, ~~~ f(\sigma, A) \gamma = 0.

The variable :math:`\gamma` is called the plastic multiplier. Note that when :math:`\psi(\varepsilon^e, \alpha), f(\sigma, A) \mbox{ or } \Psi(\sigma, A)` are not differentiable, subdifferentials have to be used. Associated plasticity corresponds to the choice :math:`\Psi(\sigma, A) = f(\sigma, A)`.

Initial boundary value problem
==============================

The weak formulation of a dynamic elastoplastic problem can be written, for an arbitrary kinematically admissible test function :math:`v`, as follows:

.. math::

   \left| \begin{array}{l}
   \ds \int_{\Omega} \rho \ddot{u}\cdot v + \sigma : \nabla v dx =  \int_{\Omega} f\dot v dx + \int_{\Gamma_N} g\dot v dx, \\
   u(0,x) = u_0(x), ~~~\dot{u}(0) = \mathrm{v}_0(x), \\
   \varepsilon^p(0,x) = \varepsilon^p_0, ~~~ \alpha(0,x) = \alpha_0,
   \end{array} \right.

for :math:`u_0, \mathrm{v}_0, \varepsilon^p_0, \alpha_0` being initial values and :math:`f` and :math:`g` being prescribed forces in the interior of domain :math:`\Omega` and on the part of the boundary :math:`\Gamma_N`.

Note that plasticity models are often applied on quasi-static problems which correspond to the term :math:`\rho \ddot{u}` being neglected.

Given a time step :math:`\Delta t = t_{n+1} -t_n`, from time :math:`t_n` to :math:`t_{n+1}`, we will denote in the sequel :math:`u_n, \varepsilon^p_n  \mbox{ and } \alpha_n` the approximations at time :math:`t_n` of :math:`u(t_n), \varepsilon^p_n(t_n) \mbox{ and } \alpha(t_n)` respectively. These approximations correspond to the chosen time integration scheme (for instance one of the proposed schemes in :ref:`ud-model-time-integration`) which can be different than the time integration scheme used for the integration of the flow rule (see below).


Flow rule integration
+++++++++++++++++++++

The plastic flow rule has to be integrated with its own time integration scheme. Among standards schemes, the backward Euler scheme, the :math:`\theta`-scheme (or generalized trapezoidal rule) and the generalized mid-point scheme are the most commonly used in that context. We make here the choice of the :math:`\theta`-scheme (:math:`\theta = 1` corresponds to the backward Euler scheme as a special case).


Let :math:`u_{n+1}` be the displacement at the considered time step and  :math:`u_{n}` at the previous one.

The :math:`\theta`-scheme for the integration of the plastic flow rules reads as

.. math:: \varepsilon^p_{n+1} - \varepsilon^p_{n} = (1-\theta)\Delta t \gamma_n \Frac{\partial \Psi}{\partial \sigma}(\sigma_{n}, A_{n}) + \theta \Delta t \gamma_{n+1} \Frac{\partial \Psi}{\partial \sigma}(\sigma_{n+1}, A_{n+1}),
  :label: thetascheme1

.. math:: \alpha_{n+1} - \alpha_n = -(1-\theta)\Delta t \gamma_n \Frac{\partial \Psi}{\partial A}(\sigma_{n}, A_{n}) - \theta\Delta t \gamma_{n+1} \Frac{\partial \Psi}{\partial A}(\sigma_{n+1}, A_{n+1}),
  :label: thetascheme2

with the complementary condition

.. math:: f(\sigma_{n+1}, A_{n+1}) \le 0, ~~~ \gamma_{n+1} \ge 0, ~~~ f(\sigma_{n+1}, A_{n+1}) \gamma_{n+1} = 0.

where :math:`0 < \theta \le 1` is the parameter of the :math:`\theta`-scheme. We exclude :math:`\theta = 0` because we will not consider explicit integration of plasticity. Let us recall that :math:`\theta = 1` corresponds to the backward Euler scheme and :math:`\theta = 1/2` to the Crank-Nicolson scheme (or trapezoidal rule) which is a second order consistent scheme. Note that the complementarity condition for the quantities at time step :math:`n` is prescribed at the previous time step (:math:`\sigma_{n}, \alpha_n, \mbox{and } \gamma_n` are supposed to be already determined).

A solution would be to solve the whole problem with all the unknows, that is :math:`u_{n+1},  \gamma_{n+1}, \varepsilon^p_{n+1} \mbox{ and } A_{n+1}`. This is of course possible but would be a rather expensive strategy because of the resulting high number of degrees of freedom. A classical strategy (the return mapping one for instance, see [SO-PE-OW2008]_ or the closest point projection one) consist in integrating locally the plastic flow on each Gauss point of the considered integration method separately, or more precisely to consider on each Gauss point the maps

.. math::
   {\mathscr E}^p : (u_{n+1}, \zeta_n, \eta_n) \mapsto \varepsilon^p_{n+1}

   {\mathscr A} : (u_{n+1}, \zeta_{n}, \eta_n) \mapsto \alpha_{n+1}



with :math:`\eta_n, \zeta_{n}` the right hand side of equations :eq:`thetascheme1`, :eq:`thetascheme2`, i.e.

.. math::
   \zeta_n = \varepsilon^p_{n} + (1-\theta)\Delta t \gamma_n \Frac{\partial \Psi}{\partial \sigma}(\sigma_{n}, A_{n}) ,

   \eta_n = \alpha_n - (1-\theta)\Delta t \gamma_n \Frac{\partial \Psi}{\partial A}(\sigma_{n}, A_{n})

This means in particular that :math:`(\varepsilon^p_{n+1}, \alpha_{n+1}) = ({\mathscr E}^p(u_{n+1},  \zeta_n, \eta_n), {\mathscr A}(u_{n+1}, \zeta_{n}, \eta_n))` is the solution to equations :eq:`thetascheme1` and :eq:`thetascheme2`. Both these maps and their tangent moduli (usually called consistent tangent moduli) are then used in the global solve of the problem with a Newton method and for :math:`u_{n+1}` the unique remaining variable. The advantage of the return mapping strategy is that the unique variable of the global solve is the displacement :math:`u_{n+1}`. A nonlinear solve on each Gauss point is often necessary which is usualy performed with a local Newton method.

In |gf| we propose both the return mapping strategy and also an alternative strategy developed below which is mainly inspired from  [PO-NI2016]_,  [SE-PO-WO2015]_ and [HA-WO2009]_ and allow more simple tangent moduli. It consists in keeping (a multiple of) :math:`\gamma_{n+1}` as an additional unknown with respect to :math:`u_{n+1}`. As we will see, this will allow a more generic treatment of the yield functions, the price for the simplicity being this additional unknown scalar field.

First, we consider an additional (and optional) given function :math:`\alpha(\sigma_{n+1}, A_{n+1}) > 0` whose interest will appear later on (it will allow simple local inverses) and the new unknown scalar field

.. math:: \xi_{n+1} = \Frac{\gamma_{n+1}}{\alpha(\sigma_{n+1}, A_{n+1})} ,

so that our two main unknows are now :math:`u_{n+1} \mbox{ and } \xi_{n+1}`. The discretized plastic flow rule integration now reads:

.. math:: \varepsilon^p_{n+1} - \varepsilon^p_{n} = (1-\theta)\alpha(\sigma_n,A_n)\Delta t \xi_n \Frac{\partial \Psi}{\partial \sigma}(\sigma_{n}, A_{n}) + \theta \alpha(\sigma_{n+1},A_{n+1}) \Delta t \xi_{n+1} \Frac{\partial \Psi}{\partial \sigma}(\sigma_{n+1}, A_{n+1}),
  :label: flowrule1

.. math:: \alpha_{n+1} - \alpha_n = (1-\theta) \alpha(\sigma_n,A_n)\Delta t \xi_n \Frac{\partial \Psi}{\partial A}(\sigma_{n}, A_{n}) + \theta \alpha(\sigma_{n+1},A_{n+1}) \Delta t \xi_{n+1} \Frac{\partial \Psi}{\partial A}(\sigma_{n+1}, A_{n+1}),
  :label: flowrule2

.. math:: f(\sigma_{n+1}, A_{n+1}) \le 0, ~~~ \xi_{n+1} \ge 0, ~~~ f(\sigma_{n+1}, A_{n+1}) \xi_{n+1} = 0.
  :label: flowrule3


For :math:`u_{n+1} \mbox{ and } \xi_{n+1}` be given, we define the two maps

.. math::
   \tilde{\mathscr E}^p : (u_{n+1}, \theta \Delta t \xi_{n+1}, \zeta_{n}, \eta_n) \mapsto \varepsilon^p_{n+1}

   \tilde{\mathscr A} : (u_{n+1}, \theta \Delta t \xi_{n+1}, \zeta_{n}, \eta_n) \mapsto \alpha_{n+1}


where the pair :math:`(\varepsilon^p_{n+1}, \alpha_{n+1}) = (\tilde{\mathscr E}^p(u_{n+1}, \theta \xi_{n+1}, \zeta_{n}, \eta_n), \tilde{\mathscr A}(u_{n+1}, \theta \xi_{n+1}, \zeta_{n}, \eta_n))` is the solution to equations :eq:`flowrule1`, :eq:`flowrule2` (without the consideration of  :eq:`flowrule3`). We will see later, that, at least for simple isotropic plastic flow rules, these maps have a simple expression, even sometimes a linear one with respect to :math:`u_{n+1}`.

Still :math:`u_{n+1} \mbox{ and } \xi_{n+1}` be given the stress :math:`\sigma_{n+1}` reads

.. math:: \sigma_{n+1} = \Frac{\partial \psi^e}{\partial \varepsilon^e}(\varepsilon(u_{n+1}) -\varepsilon^p_{n+1}).

.. math:: A_{n+1} = \Frac{\partial \psi^p}{\partial \alpha}(\alpha_{n+1}).

The complementarity equation :eq:`flowrule3` is then prescribed with the use of a well chosen complementarity function, as in [HA-WO2009]_ for :math:`r > 0` such as:

.. math:: \ds \int_{\Omega} (\xi_{n+1} - (\xi_{n+1} + r f(\sigma_{n+1}, A_{n+1}))_+) \lambda dx = 0,   \forall \lambda

or

.. math:: \ds \int_{\Omega} (f(\sigma_{n+1} + (-f(\sigma_{n+1}, A_{n+1}) - \xi_{n+1}/r)_+ , A_{n+1}) ) \lambda dx = 0,   \forall \lambda

NOTE : The notation :math:`\Delta \xi_{n+1} = \Delta t \xi_{n+1}` is often used in the litterature. The choice here is to preserve the distinction between the two quantities, mainly because ot the possible use of adaptative time step : when the time step is changing, the value :math:`\xi_n` has to be multiplied by the new time step, so that it is preferable to store :math:`\xi_n` instead of :math:`\Delta \xi_{n}` when using the :math:`\theta`-scheme.


Plane strain approximation
==========================

A plane strain approximation is a 2D problem which corresponds to the deformation of a long cylindrical object where the strain in the length direction (assumed to be along the :math:`z` axis) is considered small compared to the ones in the other directions and is neglected. It result in a plane strain tensor of the form

.. math:: \varepsilon(u) = \left(\hspace{-0.5em}\begin{array}{ccc} \varepsilon_{1,1} & \varepsilon_{1,2} & 0 \\ \varepsilon_{1,2} & \varepsilon_{2,2} & 0 \\ 0 & 0 & 0 \end{array}\hspace{-0.5em}\right).

We denote

.. math:: \bar{\varepsilon}(u) =  \left(\hspace{-0.5em}\begin{array}{cc} \varepsilon_{1,1} & \varepsilon_{1,2} \\ \varepsilon_{1,2} & \varepsilon_{2,2} \end{array}\hspace{-0.5em}\right)

the non neglected components of the strain tensor.
In the decomposition of plastic and elastic part of the strain tensor, we assume

.. math:: \varepsilon^p_{1,3} = \varepsilon^p_{2,3} = \varepsilon^e_{1,3} = \varepsilon^e_{2,3} = 0

and

.. math:: \varepsilon^e_{3,3} + \varepsilon^p_{3,3} = \varepsilon_{3,3} = 0.

The adaptation to the plane strain approximation to plastic model is most of the time an  easy task. An isotropic linearized elastic response reads

.. math:: \sigma = \lambda \mbox{tr}(\varepsilon(u)) I + 2\mu(\varepsilon(u) - \varepsilon^p),

and thus

.. math:: \bar{\sigma} = \lambda \mbox{tr}(\bar{\varepsilon}(u)) \bar{I} + 2\mu(\bar{\varepsilon}(u) -\bar{\varepsilon}^p),

The nonzero :math:`\sigma_{3,3}` component of the stress tensor is given by

.. math:: \sigma_{3,3} = \lambda \mbox{tr}(\bar{\varepsilon}(u)) - 2\mu \varepsilon^p_{3,3}

Note that in the common case where isochoric plastic strain is assumed, one has

.. math:: \mbox{ tr}(\varepsilon^p) = 0 ~~~~ \Rightarrow  ~~~ \varepsilon^p_{3,3} = - (\varepsilon^p_{1,1} + \varepsilon^p_{2,2}).



Plane stress approximation
==========================

The plane stress approximation describe generally the 2D membrane deformation of a thin plate. It consist in prescribing the stress tensor to have only in-plane nonzero components, i.e.

.. math:: \sigma = \left(\hspace{-0.5em}\begin{array}{ccc} \sigma_{1,1} & \sigma_{1,2} & 0 \\ \sigma_{1,2} & \sigma_{2,2} & 0 \\ 0 & 0 & 0 \end{array}\hspace{-0.5em}\right).

We will still denote

.. math:: \bar{\sigma} =  \left(\hspace{-0.5em}\begin{array}{cc} \sigma_{1,1} & \sigma_{1,2} \\ \sigma_{1,2} & \sigma_{2,2} \end{array}\hspace{-0.5em}\right)

the in-plane components of the stress tensor. For elastoplasticity, it consists generally to apply the 2D plastic flow rule, prescribing the out-plane components of the stress tensor to be zero with the additionnal variables :math:`\varepsilon^e_{1,3}`, :math:`\varepsilon^e_{2,3}`, :math:`\varepsilon^e_{3,3}` being unknown (see for instance [SO-PE-OW2008]_).

For an isotropic linearized elastic response, one has :math:`\sigma = \lambda \mbox{tr}(\varepsilon^e) + 2\mu\varepsilon^e` such that

.. math:: \varepsilon^e = \left(\hspace{-0.5em}\begin{array}{ccc} \varepsilon^e_{1,1} & \varepsilon^e_{1,2} & 0 \\ \varepsilon^e_{1,2} & \varepsilon^e_{2,2} & 0 \\ 0 & 0 & \varepsilon^e_{3,3} \end{array}\hspace{-0.5em}\right).

with

.. math:: \varepsilon^e_{3,3} = -\Frac{\lambda}{\lambda+2\mu}(\varepsilon^e_{1,1} + \varepsilon^e_{2,2})

so that

.. math:: \bar{\sigma} = \lambda^* \mbox{tr}(\bar{\varepsilon}^e) + 2\mu\bar{\varepsilon}^e ~~~\mbox{ with } \lambda^* = \Frac{2\mu\lambda}{\lambda+2\mu}
   :label: plane_stress_iso

Moreover

.. math:: \|\mbox{Dev}(\sigma)\| = \left(\|\bar{\sigma}\|^2 - \Frac{1}{3}(\mbox{tr}(\bar{\sigma}))^2\right)^{1/2}.
   :label: plane_stress_dev

Note that in the case where isochoric plastic strain is assumed, one still has

.. math:: \mbox{ tr}(\varepsilon^p) = 0 ~~~~ \Rightarrow  ~~~ \varepsilon^p_{3,3} = - (\varepsilon^p_{1,1} + \varepsilon^p_{2,2}).


Some classical laws
+++++++++++++++++++


Tresca : :math:`\rho(\sigma) \le \sigma_y` where :math:`\rho(\sigma)` spectral radius of the Cauchy stress tensor and :math:`\sigma_y` the uniaxial yield stress (which may depend on some hardening internal variables).

Von Mises :  :math:`\|\mbox{Dev}(\sigma)\| \le \sqrt{\frac{2}{3}}\sigma_y` where
:math:`\mbox{Dev}(\sigma) = \sigma - \frac{1}{3}\mbox{tr}(\sigma)I` the deviatoric part of :math:`\sigma` and :math:`\|\sigma\| = \sqrt{\sigma:\sigma}`.


Perfect isotropic associated elastoplasticity with Von-Mises criterion (Prandl-Reuss model)
===========================================================================================

There is no internal variables and we consider an isotropic elastic response. The flow rule reads

.. math:: \dot{\varepsilon}^p = \gamma \Frac{\mbox{Dev}(\sigma)}{\|\mbox{Dev}(\sigma)\|}

This corresponds to :math:`\Psi(\sigma) = f(\sigma) = \|\mbox{Dev}(\sigma)\| - \sqrt{\frac{2}{3}}\sigma_y`.


The :math:`\theta`-scheme for the integration of the plastic flow rule reads:

.. math:: \varepsilon^p_{n+1} - \varepsilon^p_{n} = (1-\theta)\alpha(\sigma_{n}) \Delta t \xi_n \Frac{\mbox{Dev}(\sigma_{n})}{\|\mbox{Dev}(\sigma_{n})\|} + \theta\alpha(\sigma_{n+1}) \Delta t \xi_{n+1} \Frac{\mbox{Dev}(\sigma_{n+1})}{\|\mbox{Dev}(\sigma_{n+1})\|}.

Choosing the factor :math:`\alpha(\sigma_{n}) = \|\mbox{Dev}(\sigma_{n})\|` and still with :math:`\xi_n = \Frac{\gamma_n}{\alpha(\sigma_{n})}` this gives the equation

.. math::  \varepsilon^p_{n+1} - \varepsilon^p_{n} = (1-\theta)\Delta t \xi_n \mbox{Dev}(\sigma_{n}) + \theta \Delta t \xi_{n+1} \mbox{Dev}(\sigma_{n+1}).

Since :math:`\mbox{Dev}(\sigma_{n+1}) = 2\mu\mbox{Dev}(\varepsilon(u_{n+1})) - 2\mu\varepsilon^p_{n+1}` this directly gives:

.. math:: \tilde{\mathscr E}^p(u_{n+1}, \theta \Delta t \xi_{n+1}, \zeta_{n}) = \zeta_n + \left(1-\Frac{1}{1+2\mu\theta\Delta t\xi_{n+1}}\right)(\mbox{Dev}(\varepsilon(u_{n+1})) - \zeta_n),

which is a linear expression with respect to :math:`u_{n+1}` (but not with respect to :math:`\xi_{n+1}`).

Moreover, :math:`\zeta_n` is defined by

.. math:: \zeta_n = \varepsilon^p_n+(1-\theta)\Delta t \xi_n (\mbox{Dev}(\sigma_n)) = \varepsilon^p_n+(1-\theta)\Delta t \xi_n 2\mu \left(\mbox{Dev}(\varepsilon(u_{n}))-\varepsilon^p_n\right).

**Elimination of the multiplier (for the return mapping approach)**

One has

.. math:: \|\mbox{Dev}(\sigma_{n+1})\| = 2\mu\|\mbox{Dev}(\varepsilon(u_{n+1})) -\varepsilon^p_{n+1}\| = \Frac{2\mu}{1+2\mu\theta\Delta t \xi_{n+1}}\|\mbox{Dev}(\varepsilon(u_{n+1})) - \zeta_n\|,

Thus, denoting :math:`B = \mbox{Dev}(\varepsilon(u_{n+1})) - \zeta_n`, either

.. math:: 2\mu\|B\| \le \sqrt{\frac{2}{3}}\sigma_y,

and :math:`\xi_{n+1} = 0`, i.e. we are in the elastic case, or  :math:`\|\mbox{Dev}(\sigma_{n+1})\| =  \sqrt{\frac{2}{3}}` and one obtains

.. math:: 1+2\mu\theta\Delta t \xi_{n+1} = \Frac{2\mu\|B\|}{\sqrt{\frac{2}{3}}\sigma_y},

and thus

.. math:: \varepsilon^p_{n+1} = \zeta_n + \left( 1 - \sqrt{\frac{2}{3}}\Frac{\sigma_y}{2\mu\|B\|}\right) B.

The two options can be summarized by

.. math:: \varepsilon^p_{n+1} = {\mathscr E}^p(u_{n+1}, \zeta_{n}) = \zeta_n + \left( 1 - \sqrt{\frac{2}{3}}\Frac{\sigma_y}{2\mu\|B\|}\right)_+ B.

The multiplier :math:`\xi_{n+1}` (needed for the :math:`\theta`-scheme for :math:`\theta \ne 1`) is given by

.. math:: \xi_{n+1} = \Frac{1}{\theta\Delta t}\left(\sqrt{\frac{3}{2}}\Frac{\|B\|}{\sigma_y} - \Frac{1}{2\mu}\right)_+.


**Plane strain approximation**

The plane strain approximation has the same expression replacing the 3D strain tensors by the in-plane ones :math:`\bar{\varepsilon}^p` and  :math:`\bar{\varepsilon}(u_{n+1})`.

.. math:: \bar{\tilde{\mathscr E}}^p(u_{n+1}, \theta \Delta t \xi_{n+1}, \bar{\zeta}_{n}) = \bar{\zeta}_n + \left(1-\Frac{1}{1+2\mu\theta\Delta t\xi_{n+1}}\right)(\overline{\mbox{Dev}}(\bar{\varepsilon}(u_{n+1})) - \bar{\zeta}_n),

where :math:`\overline{\mbox{Dev}}(\bar{\varepsilon}) = \bar{\varepsilon} - \Frac{\mbox{tr}(\bar{\varepsilon})}{3} \bar{I}` is the 2D restriction of the 3D deviator.

Moreover, for the yield condition,

.. math:: \|\mbox{Dev}(\sigma)\|^2 = 4\mu^2\left(\|\overline{\mbox{Dev}}\bar{\varepsilon}(u) - \bar{\varepsilon}^p\|^2 + \left(\Frac{\mbox{tr}(\bar{\varepsilon}(u))}{3} -\mbox{tr}(\bar{\varepsilon}^p) \right)^2\right).

And for the elimination of the multiplier,

.. math:: \bar{\mathscr E}^p(\bar{u}_{n+1}, \bar{\varepsilon}^p_{n}) = \bar{\zeta}^p_{n} + \left( 1 - \sqrt{\frac{2}{3}}\Frac{\sigma_y}{2\mu\|B\|}\right)_+ \bar{B}

with :math:`\bar{B} = \overline{\mbox{Dev}}(\bar{\varepsilon}(u_{n+1}))-\bar{\zeta}_{n}` and :math:`\|B\|^2 = \|\overline{\mbox{Dev}}(\bar{\varepsilon}(u_{n+1})) - \bar{\zeta}_n\|^2 + \left(\Frac{\mbox{tr}(\bar{\varepsilon}(u_{n+1}))}{3} -\mbox{tr}(\bar{\zeta}_n) \right)^2`.

**Plane stress approximation**

For plane stress approximation, using :eq:`plane_stress_iso` we deduce from the expression of the 3D case

.. math::  \bar{\varepsilon}^p_{n+1} = \Frac{1}{1+2\mu\theta\Delta \xi}\left(\bar{\zeta}_{n} +2\mu\theta\Delta \xi\left(\bar{\varepsilon}(u_{n+1}) - \Frac{2\mu}{3(\lambda+2\mu)}(\mbox{tr}(\bar{\varepsilon}(u_{n+1})) - \mbox{tr}(\bar{\varepsilon}_{n+1}^p))\bar{I}\right) \right),

since :math:`\mbox{Dev}(\varepsilon(u)) = \varepsilon(u) - \Frac{2\mu}{3(\lambda+2\mu)}(\mbox{tr}(\bar{\varepsilon}(u)) - \mbox{tr}(\bar{\varepsilon}^p))`. Of course, this relation still has to be inverted. Denoting :math:`\alpha = 1+2\mu\theta\Delta \xi`, :math:`\beta = \Frac{4\mu^2\theta\Delta \xi}{3\lambda+6\mu}` and :math:`C = \bar{\zeta}_{n} +2\mu\theta\Delta \xi\left(\bar{\varepsilon}(u_{n+1}) - \Frac{2\mu}{3(\lambda+2\mu)}(\mbox{tr}(\bar{\varepsilon}(u_{n+1}))))\bar{I}\right)` one obtains

.. math:: \bar{\varepsilon}^p_{n+1} = \Frac{\beta \mbox{tr}(C)}{\alpha(\alpha-2\beta)}\bar{I} + \Frac{1}{\alpha}C.

Moreover, for the yield condition, expression :eq:`plane_stress_dev` can be used.

Isotropic elastoplasticity with linear isotropic and kinematic hardening and Von-Mises criterion
================================================================================================

We consider an isotropic elastic reponse and the internal variable :math:`\alpha : \Omega \rightarrow \R` being the accumulated plastic strain which satisfies

.. math:: \dot{\alpha} = \sqrt{\Frac{2}{3}}\gamma

For :math:`H_i` the isotropic hardening modulus, the linear hardening consists in

.. math:: \psi^p(\alpha) = \frac{1}{2}H_i\alpha^2

i.e. :math:`A = H_i\alpha` and a uniaxial yield stress defined by


.. math:: \sigma_y(a) = \sigma_{y0} + A = \sigma_{y0} + H_i\alpha,

for :math:`\sigma_{y0}` the initial uniaxial yield stress. The yield function (and plastic potential since this is an associated plastic model) can be defined by

.. math:: \Psi(\sigma, A) = f(\sigma, A) = \|\mbox{Dev}(\sigma - \frac{2}{3}H_k\varepsilon^p)\| - \sqrt{\frac{2}{3}}(\sigma_{y0} + A),

where :math:`H_k` is the kinematic hardening modulus. The same computation as in the previous section leads to

.. math:: \tilde{\mathscr E}^p(u_{n+1}, \theta\Delta t \xi_{n+1}, \zeta_n) = \zeta_n + \Frac{1}{2(\mu+H_k/3)}\left(1 - \Frac{1}{1+2(\mu+H_k/3)\theta\Delta t\xi_{n+1}}\right)(2\mu\mbox{Dev}(\varepsilon(u_{n+1}))-2(\mu+H_k/3)\zeta_n)


.. math:: \begin{array}{rcl} \tilde{\mathscr A}(u_{n+1}, \theta \Delta t \xi_{n+1}, \zeta_{n}, \eta_n) &=& \eta_n + \sqrt{\Frac{2}{3}} \theta \Delta t \xi_{n+1}\|\mbox{Dev}(\sigma_{n+1} - \frac{2}{3}H_k\varepsilon^p_{n+1})\| \\ &=& \eta_n + \sqrt{\Frac{2}{3}} \theta \Delta t \xi_{n+1}\|2\mu\mbox{Dev}(\varepsilon(u_{n+1})) - 2(\mu+H_k/3)\varepsilon^p_{n+1}\| \\ &=&  \eta_n + \sqrt{\Frac{2}{3}} \Frac{\theta \Delta t \xi_{n+1}}{1+2(\mu+H_k/3)\theta\Delta t\xi_{n+1}}\|2\mu\mbox{Dev}(\varepsilon(u_{n+1})) - 2(\mu+H_k/3)\zeta_{n}\| \\ &=& \eta_n + \sqrt{\Frac{2}{3}}\Frac{1}{2(\mu+H_k/3)}\left(1 - \Frac{1}{1+2(\mu+H_k/3)\theta\Delta t\xi_{n+1}}\right) \|2\mu\mbox{Dev}(\varepsilon(u_{n+1})) - 2(\mu+H_k/3)\zeta_{n}\|\end{array}

where :math:`\zeta_n` and :math:`\eta_n` are defined by

.. math:: \zeta_n = \varepsilon^p_n+(1-\theta)\Delta t \xi_n (\mbox{Dev}(\sigma_n)-\frac{2}{3}H_k\varepsilon^n_p) = \varepsilon^p_n+(1-\theta)\Delta t \xi_n \left(2\mu\mbox{Dev}(\varepsilon(u_{n}))-2(\mu+H_k/3)\varepsilon^n_p\right),

.. math:: \eta_n  = \alpha_n+(1-\theta)\sqrt{\Frac{2}{3}}\Delta t \xi_n \|\mbox{Dev}(\sigma_n)-\frac{2}{3}H_k\varepsilon^n_p\| =  \alpha_n+(1-\theta)\sqrt{\Frac{2}{3}}\Delta t \xi_n \|2\mu\mbox{Dev}(\varepsilon(u_{n}))-2(\mu+H_k/3)\varepsilon^n_p\|.

Note that the isotropic hardening modulus do not intervene in :math:`\tilde{\mathscr E}^p(u_{n+1}, \theta \Delta \xi, \varepsilon^p_{n})` but only in :math:`f(\sigma, A)`.

**Elimination of the multiplier (for the return mapping approach)**

Denoting :math:`\delta = \Frac{1}{1+2(\mu+H_k/3)\theta\Delta t\xi_{n+1}}`, :math:`\beta = \Frac{1-\delta}{2(\mu+H_k/3)}` and :math:`B = 2\mu\mbox{Dev}(\varepsilon(u_{n+1}))-2(\mu+H_k/3)\zeta_n` the expression for :math:`\varepsilon^p_{n+1}` and :math:`\alpha_{n+1}` becomes

.. math:: \varepsilon^p_{n+1} = \zeta_n+\beta B, ~~~ \alpha_{n+1} = \eta_n + \sqrt{\Frac{2}{3}}\beta \|B\|,
  :label: hardeningepsalp

and the plastic constraint

.. math:: \delta \|B\| \le \sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \alpha_{n+1}).

Thus, either we are in the elastic case, i.e. :math:`\xi_{n+1} = 0, \delta = 1` and

.. math:: \|B\| \le \sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \eta_n),

or we are in the plastic case and :math:`\xi_{n+1} > 0, \delta < 1`, :math:`\delta \|B\| = \sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \alpha_{n+1})` and :math:`(1-\delta)` solves the equation

.. math:: \|B\| - (1-\delta)\|B\| = \sqrt{\Frac{2}{3}}\left(\sigma_{y0}+H_i \eta_n + \sqrt{\Frac{2}{3}} \Frac{H_i}{2(\mu+H_k/3)}(1-\delta)\|B\|\right),

which leads to

.. math:: 1-\delta = \Frac{2(\mu+H_k/3)}{\|B\|(2\mu+\frac{2}{3}(H_k+H_i))}\left(\|B\|-\sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \eta_n) \right)

The two cases can be summarized by

.. math:: \beta = \Frac{1}{\|B\|(2\mu+\frac{2}{3}(H_k+H_i))}\left(\|B\|-\sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \eta_n) \right)_+

which directly gives :math:`{\mathscr E}^p(u_{n+1}, \zeta_n, \eta_n)` and :math:`{\mathscr A}(u_{n+1}, \zeta_n, \eta_n)` thanks to :eq:`hardeningepsalp`. The multiplier :math:`\xi_{n+1}` being given by

.. math:: \xi_{n+1} = \Frac{1}{(2(\mu+H_k/3))\theta\Delta t}(\Frac{1}{\delta}-1) = \Frac{1}{\theta\Delta t}~\Frac{\beta}{1-2(\mu+H_k/3)\beta}.


**Plane strain approximation**

Still denoting  :math:`\delta = \Frac{1}{1+2(\mu+H_k/3)\theta\Delta t\xi_{n+1}}`, :math:`\beta = \Frac{1-\delta}{2(\mu+H_k/3)}`, :math:`B = 2\mu\mbox{Dev}(\varepsilon(u_{n+1}))-2(\mu+H_k/3)\zeta_n` and :math:`\overline{B} = 2\mu\overline{Dev}(\bar{\varepsilon}(u_{n+1}))-2(\mu+H_k/3)\bar{\zeta}_n` its in-plane part, one has

.. math:: \bar{\tilde{\mathscr E}}^p(u_{n+1}, \theta\Delta t \xi_{n+1}, \bar{\zeta}_n) = \bar{\zeta}_n + \beta \overline{B},


.. math:: \tilde{\mathscr A}(u_{n+1}, \theta \Delta t \xi_{n+1}, \zeta_{n}, \eta_n) = \eta_n + \sqrt{\Frac{2}{3}}\beta\|B\|,

with

.. math:: \|B\|^2 = \|2\mu\overline{\mbox{Dev}}(\bar{\varepsilon}(u_{n+1})) - 2(\mu+H_k/3)\bar{\zeta}_n\|^2 + \left(2\mu\Frac{\mbox{tr}(\bar{\varepsilon}(u_{n+1}))}{3} -2(\mu+H_k/3)\mbox{tr}(\bar{\zeta}_n) \right)^2.

The yield condition still reads

.. math:: \delta \|B\| \le \sqrt{\Frac{2}{3}}(\sigma_{y0}+H_i \alpha_{n+1}).

and for the elimination of the multiplier, :math:`\beta` has the same expression as in the previous section adapting the value of :math:`\|B\|`. The expressions of :math:`\bar{\zeta}_n` and :math:`\eta_n` have to be adapted accoringly.





Souza-Auricchio elastoplasticity law (for shape memory alloys)
==============================================================

See for instance [GR-ST2015]_ for the justification of the construction of this flow rule. A Von-Mises stress criterion together with an isotropic elastic response, no internal variables and a special type of kinematic hardening is considered with a constraint :math:`\|\varepsilon^p\| \le c_3`. The plastic potential and yield function have the form

.. math:: \Psi(\sigma) = f(\sigma)  = \left\|\mbox{Dev}\left(\sigma - c_1\Frac{\varepsilon^p}{\|\varepsilon^p\|} - c_2\varepsilon^p - \delta \Frac{\varepsilon^p}{\|\varepsilon^p\|}\right)\right\| - \sqrt{\frac{2}{3}}\sigma_{y},

with the complementarity condition

.. math::
   \delta \ge 0, \|\varepsilon^p\| \le c_3,  \delta (\|\varepsilon^p\|-c_3) = 0,

where :math:`c_1, c_2 \mbox{ and } c_3` are some physical parameters. Note that :math:`\Frac{\varepsilon^p}{\|\varepsilon^p\|}` has to be understood to be the whole unit ball for :math:`\varepsilon^p = 0`.


to be done ...






Elasto-plasticity bricks
+++++++++++++++++++++++++

See the test programs :file:`tests/plasticity.cc`, :file:`interface/tests/matlab/demo_plasticity.m`, :file:`interface/tests/matlab/demo_plasticity.py` and in :file:`contrib/test_plasticity`.

Generic brick
=============

There are two versions of the generic brick. A first one when the plastic multiplier is kept as a variable of the problem where the added term is of the form:

.. math:: \int_{\Omega} \sigma_{n+1} : \nabla \delta u dx +  \ds \int_{\Omega} (\xi_{n+1} - (\xi_{n+1} + r f(\sigma_{n+1}, A_{n+1}))_+) \delta\xi dx = 0,

with :math:`r > 0` having a specific value chosen by the brick (in terms of the elasticity coefficients), and when the return mapping strategy is selected (plastic multiplier is just a data), just the added term:

.. math:: \int_{\Omega} \sigma_{n+1} : \nabla v dx.

The function which adds the brick to a model `md` is ::

   getfem::add_small_strain_elastoplasticity_brick
     (md, mim, lawname, unknowns_type,
      const std::vector<std::string> &varnames,
      const std::vector<std::string> &params, region = size_type(-1));


where `lawname` is the name of an implemented plastic law, `unknowns_type`
indicates the choice between a discretization where the plastic multiplier
is an unknown of the problem or (return mapping approach) just a data of
the model stored for the next iteration. Remember that in both cases, a
multiplier is stored anyway. `varnames` is a set of variable and data
names with length which may depend on the plastic law (at least the
displacement, the plastic multiplier and the plastic strain).
`params` is a list of expressions for the parameters (at least elastic
coefficients and the yield stress). These expressions can be some data
names (or even variable names) of the model but can also be any scalar
valid expression of GWFL, the generic weak form language (such as "1/2",
"2+sin(X[0])", "1+Norm(v)" ...). The last two parameters optionally
provided in `params` are the `theta` parameter of the `theta`-scheme
(generalized trapezoidal rule) used for the plastic strain integration
and the time-step`dt`. The default value for `theta` if omitted is 1,
which corresponds to the classical Backward Euler scheme which is first
order consistent. `theta=1/2` corresponds to the Crank-Nicolson scheme
(trapezoidal rule) which is second order consistent. Any value
between 1/2 and 1 should be a valid value. The default value of `dt` is
'timestep' which simply indicates the time step defined in the model
(by md.set_time_step(dt)). Alternatively it can be any expression
(data name, constant value ...). The time step can be altered from one
iteration to the next one. `region` is a mesh region.

The available plasticity laws are:

- "Prandtl Reuss" (or "isotropic perfect plasticity").
  Isotropic elasto-plasticity with no hardening. The variables are the
  displacement, the plastic multiplier and the plastic strain.
  The displacement should be a variable and have a corresponding data
  having the same name preceded by "Previous\_" corresponding to the
  displacement at the previous time step (typically "u" and "Previous_u").
  The plastic multiplier should also have two versions (typically "xi"
  and "Previous_xi") the first one being defined as data if
  `unknowns_type = DISPLACEMENT_ONLY` or as a variable if
  `unknowns_type = DISPLACEMENT_AND_PLASTIC_MULTIPLIER`.
  The plastic strain should represent a n x n data tensor field stored
  on mesh_fem or (preferably) on an im_data (corresponding to `mim`).
  The data are the first Lame coefficient, the second one (shear modulus)
  and the uniaxial yield stress. IMPORTANT: Note that this law implements
  the 3D expressions. If it is used in 2D, the expressions are just
  transposed to the 2D. For the plane strain approximation, see below.
- "plane strain Prandtl Reuss"
  (or "plane strain isotropic perfect plasticity")
  The same law as the previous one but adapted to the plane strain
  approximation. Can only be used in 2D.
- "Prandtl Reuss linear hardening"
  (or "isotropic plasticity linear hardening").
  Isotropic elasto-plasticity with linear isotropic and kinematic
  hardening. An additional variable compared to "Prandtl Reuss" law:
  the accumulated plastic strain. Similarly to the plastic strain, it
  is only stored at the end of the time step, so a simple data is
  required (preferably on an im_data).
  Two additional parameters: the kinematic hardening modulus and the
  isotropic one. 3D expressions only.
- "plane strain Prandtl Reuss linear hardening"
  (or "plane strain isotropic plasticity linear hardening").
  The same law as the previous one but adapted to the plane strain
  approximation. Can only be used in 2D.

IMPORTANT : remember that `small_strain_elastoplasticity_next_iter` has
to be called at the end of each time step, before the next one
(and before any post-treatment : this sets the value of the plastic
strain and plastic multiplier).

Additionaly, the following function allow to pass from a time step to another for the small strain plastic brick: ::

   getfem::small_strain_elastoplasticity_next_iter
     (md, mim, lawname, unknowns_type,
      const std::vector<std::string> &varnames,
      const std::vector<std::string> &params, region = size_type(-1));

The parameters have to be exactly the same as the ones of the
`add_small_strain_elastoplasticity_brick`,  so see the documentation of
this function for any explanations.
Basically, this brick computes the plastic strain and the plastic
multiplier and stores them for the next step. Additionaly, it copies
the computed displacement to the data that stores the displacement
of the previous time step (typically "u" to "Previous\_u").
It has to be called before any use of
`compute_small_strain_elastoplasticity_Von_Mises`.

The function ::

  getfem::compute_small_strain_elastoplasticity_Von_Mises
    (md, mim, lawname, unknowns_type,
     const std::vector<std::string> &varnames,
     const std::vector<std::string> &params,
     const mesh_fem &mf_vm, model_real_plain_vector &VM,
     region = size_type(-1));

computes the Von Mises stress field with respect to
a small strain elastoplasticity term, approximated on `mf_vm`,
and stores the result into `VM`.  All other parameters have to be
exactly the same as for `add_small_strain_elastoplasticity_brick`.
Remember that `small_strain_elastoplasticity_next_iter` has to be called
before any call of this function.


A specific brick based on the low-level generic assembly for perfect plasticity
================================================================================

This is an previous version of a elastoplasticity brick which is restricted to  isotropic perfect plasticity and is based on the low-level generic assembly. Its specificity which could be interesting for testing is that the flow rule is integrated on  finite element nodes (not on Gauss points).

The function adding this brick to a model is: ::

      getfem::add_elastoplasticity_brick
          (md, mim, ACP, varname, previous_varname, datalambda, datamu, datathreshold, datasigma, region);

where:
      - ``varname`` represents the main displacement unknown on which the brick is added (u).
      - ``previous_varname`` is the displacement at the previous time step.
      - ``datalambda`` and ``datamu`` are the data corresponding to the Lame coefficients.
      - ``datathreshold`` represents the plastic threshold of the studied material.
      - ``datasigma`` represents the stress constraint values supported by the material. It should be composed of 2 iterates for the time scheme needed for the Newton algorithm used. Note that the finite element method on which ``datasigma`` is defined should be able to represent the derivative of ``varname``.
      - ``ACP`` corresponds to the type of projection to be used. It has an `abstract_constraints_projection` type and for the moment, only exists the `VM_projection` corresponding to the Von Mises one.


Be careful: ``datalambda``, ``datamu`` and ``datathreshold`` could be constants or described on the same finite element method.

This function assembles the tangent matrix and the right hand side vector which will be solved using a Newton algorithm.


Additionaly, The function: ::

      getfem::elastoplasticity_next_iter
          (md, mim, varname, previous_varname, ACP, datalambda, datamu, datathreshold, datasigma);

computes the new stress constraint values supported by the material after a load or an unload (once a solve has been done earlier) and upload the variables ``varname`` and ``datasigma`` as follows:

.. math::

   u^{n+1} \Rightarrow u^n \ \ \ \ \ \textrm{ and } \ \ \ \ \ \sigma^{n+1} \Rightarrow \sigma^n

Then, :math:`u^n` and :math:`\sigma^n` contains the new values computed and one can restart the process.


The function: ::

      getfem::compute_elastoplasticity_Von_Mises_or_Tresca
          (md, datasigma, mf_vm, VM, tresca=false);

computes the Von Mises (or Tresca if ``tresca`` = true) criterion on the stress tensor stored in ``datasigma`` . The stress is evaluated on the `mesh_fem` ``mf_vm`` and stored into the vector ``VM``.
Of course, this function can be used if and only if the previous function ``elastoplasticity_next_iter`` has been called earlier.

The function: ::

      getfem::compute_plastic_part
          (md, mim, mf_pl, varname, previous_varname, ACP, datalambda, datamu, datathreshold, datasigma, Plast);

computes on ``mf_pl`` the plastic part of the material, that could appear after a load and an unload, into the vector ``Plast``.

Note that ``datasigma`` should be the vector containing the new stress constraint values, i.e. after a load or an unload of the material.



