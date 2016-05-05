.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-plasticity-small-strain:



Small strain plasticity
-----------------------

Work in progress. Not available for the moment ...

A framework for the approximation of plasticity models in |gf|.


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

Plastic yielding is supposed to occur when the stress attains a critical value. This is determinated by a yield function :math:`f(\sigma, A)` and the condition

.. math:: f(\sigma, A) \le 0.

The surface :math:`f(\sigma, A) = 0` is the yield surface where the plastic deformation may occur.

Let us also consider the plastic potential :math:`\Psi(\sigma, A)`, (convex with respect to its both variables) which determines the plastic flow direction in the sense that the flow rule is defined as

.. math:: \dot{\varepsilon}^p = \dot{\gamma} \Frac{\partial \Psi}{\partial \sigma}(\sigma, A), ~~~~~~ \dot{\alpha} = \dot{\gamma} \Frac{\partial \Psi}{\partial A}(\sigma, A),

with the additional complementarity condition

.. math:: f(\sigma, A) \le 0, ~~~ \dot{\gamma} \ge 0, ~~~ f(\sigma, A) \dot{\gamma} = 0.

The variable :math:`\dot{\gamma}` is called the plastic multiplier. Note that when :math:`\psi(\varepsilon^e, \alpha), f(\sigma, A) \mbox{ or } \Psi(\sigma, A)` are not differentiable, subdifferentials have to be used. Associated plasticity corresponds to the choice :math:`\Psi(\sigma, A) = f(\sigma, A)`.

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

The plastic flow rule has to be integrated with its own time integration scheme. Among standards schemes, the backward Euler scheme, the :math:`\theta`-scheme and the generalized mid-point scheme are the most commonly used in that context. We make here the choice of the generalized mid-point scheme.


Let :math:`u_{n+1}` be the displacement at the considered time step and  :math:`u_{n}` at the previous one. For a quantity :math:`B` we denote :math:`B_{n+\theta} = \theta B_{n+1} + (1-\theta)B_n` the convex combination of the quantity at iterations :math:`n` and :math:`n+1`.

The mid-point scheme for the integration of the plastic flow rules reads as

.. math:: \varepsilon^p_{n+1} - \varepsilon^p_{n} = \Delta \gamma \Frac{\Psi}{\partial \sigma}(\sigma_{n+\theta}, A_{n+\theta}),

.. math:: \alpha_{n+1} - \alpha_n = \Delta \gamma \Frac{\Psi}{\partial A}(\sigma_{n+\theta}, A_{n+\theta}),

with the complementary condition

.. math:: f(\sigma_{n+\theta}, A_{n+\theta}) \le 0, ~~~ \Delta\gamma \ge 0, ~~~ f(\sigma_{n+\theta}, A_{n+\theta}) \Delta \gamma = 0.

where :math:`0 < \theta \le 1` is the parameter of the mid-point scheme. We exclude :math:`\theta = 0` because we will not consider explicit integration of plasticity. Let us recall that :math:`\theta = 1` corresponds to the backward Euler scheme and :math:`\theta = 1/2` to the mid-point scheme which is a second order consistent scheme.

A solution would be to solve the whole problem with all the unknows, that is :math:`u_{n+1}, \Delta \gamma, \varepsilon^p_{n+1} \mbox{ and } A_{n+1}`. This is of course possible but would be a rather expensive strategy because of the resulting high number of degrees of freedom. A classical strategy (the return mapping one for instance, see [SO-PE-OW2008]_ or the closes point projection one) consist in integrating locally the plastic flow on each Gauss point of the considered integration method separately, or more precisely to consider on each Gauss point the maps

.. math::
   {\mathscr E}^p : (u_{n+\theta}, \varepsilon^p_{n}, \alpha_n) \mapsto \varepsilon^p_{n+\theta}

   {\mathscr A} : (u_{n+\theta}, \varepsilon^p_{n}, \alpha_n) \mapsto \alpha_{n+\theta}

which results from the local flow rule integration  (the pair :math:`(\varepsilon^p_{n+\theta}, \alpha_{n+\theta}) = ({\mathscr E}^p(u_{n+\theta},  \varepsilon^p_{n}, \alpha_n), {\mathscr A}(u_{n+\theta}, \varepsilon^p_{n}, \alpha_n))` is the solution to equations :eq:`flowrule1`, :eq:`flowrule2` and  :eq:`flowrule3`). Both these maps and their tangent moduli (usually called consistent tangent moduli) are then used in the global solve of the problem with a Newton method and for :math:`u_{n+1}` the unique remaining variable. The advantage of the return mapping strategy is that the unique variable of the global solve is the displacement :math:`u_{n+1}`. A nonlinear solve on each Gauss point is often necessary which is usualy performed with a local Newton method.

In |gf| we propose both the return mapping trategy and also an alternative strategy developped below which is mainly inspired from  [PO-NI2016]_,  [SE-PO-WO2015]_ and [HA-WO2009]_ and allow more simple tangent moduli. It consists in keeping (a multiple of) :math:`\Delta \gamma` as an additional unknown with respect to :math:`u_{n+1}`. As we will see, this will allow a more generic treatment of the yield functions, the price for the simplicity being this additional unknown scalar field.

First, we consider an additional (and optional) given function :math:`\alpha(\sigma_{n+\theta}, A_{n+\theta}) > 0` whose interest will appear later on (it will allow simple local inverses) and the new unknown scalar field

.. math:: \Delta \xi = \Frac{\Delta \gamma}{\alpha(\sigma_{n+\theta}, A_{n+\theta})} ,

so that our two main unknows are now :math:`u_{n+1} \mbox{ and } \Delta \xi`. The plastic flow rule integration may now read:

.. math:: \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \alpha(\sigma_{n+\theta}, A_{n+\theta}) \theta \Delta \xi \Frac{\Psi}{\partial \sigma}(\sigma_{n+\theta}, A_{n+\theta}).
   :label: flowrule1

.. math::  \alpha_{n+\theta} - \alpha_n = \alpha(\sigma_{n+\theta}, A_{n+\theta}) \theta \Delta \xi \Frac{\Psi}{\partial A}(\sigma_{n+\theta}, A_{n+\theta}),
   :label: flowrule2

.. math:: f(\sigma_{n+\theta}, A_{n+\theta}) \le 0, ~~~ \Delta\xi \ge 0, ~~~ f(\sigma_{n+\theta}, A_{n+\theta}) \Delta \xi = 0.
   :label: flowrule3

For :math:`u_{n+1} \mbox{ and } \Delta \xi` be given, we define the two maps

.. math::
   \tilde{\mathscr E}^p : (u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) \mapsto \varepsilon^p_{n+\theta}

   \tilde{\mathscr A} : (u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) \mapsto \alpha_{n+\theta}


where the pair :math:`(\varepsilon^p_{n+\theta}, \alpha_{n+\theta}) = (\tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n), \tilde{\mathscr A}(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n))` is the solution to equations :eq:`flowrule1`, :eq:`flowrule2` (without the consideration of  :eq:`flowrule3`). We will see later, that, at least for simple isotropic plastic flow rules, these maps have a simple expression, even sometimes a linear one with respect to :math:`u_{n+\theta}`.

Still :math:`u_{n+1} \mbox{ and } \Delta \xi` be given the stress :math:`\sigma_{n+1}` reads

.. math:: \sigma_{n+1} = \Frac{\partial \psi^e}{\partial \varepsilon^e}(\varepsilon(u_{n+1}) -\varepsilon^p_{n+1}).

.. math:: A_{n+1} = \Frac{\partial \psi^p}{\partial \alpha}(\alpha_{n+1}).

The complementarity equation :eq:`flowrule3` is then prescribed with the use of a well chosen complementarity function, as in [HA-WO2009]_ for :math:`r > 0` such as:

.. math:: \ds \int_{\Omega} (\Delta \xi - (\Delta \xi + r f(\sigma_{n+\theta}, A_{n+\theta}))_+) \lambda dx = 0,   \forall \lambda

or

.. math:: \ds \int_{\Omega} (f(\sigma_{n+\theta} + (-f(\sigma_{n+\theta}, A_{n+\theta}) - \Delta \xi/r)_+ , A_{n+\theta}) ) \lambda dx = 0,   \forall \lambda


pb : need of :math:`A_{n+\theta}` 



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

Note that in the case where isochoric plastic strain is assumed, one still has

.. math:: \mbox{ tr}(\varepsilon^p) = 0 ~~~~ \Rightarrow  ~~~ \varepsilon^p_{3,3} = - (\varepsilon^p_{1,1} + \varepsilon^p_{2,2}).


Some classical laws
+++++++++++++++++++


Tresca : :math:`\rho(\sigma) \le \sigma_y` where :math:`\rho(\sigma)` spectral radius of the Cauchy stress tensor and :math:`\sigma_y` the uniaxial yield stress (which may depend on some hardening internal variables.

Von Mises :  :math:`\|\mbox{Dev}(\sigma)\| \le \sqrt{\frac{2}{3}}\sigma_y` where
:math:`\mbox{Dev}(\sigma) = \sigma - \frac{1}{3}\mbox{tr}(\sigma)I` the deviatoric part of :math:`\sigma` and :math:`\|\sigma\| = \sqrt{\sigma:\sigma}`.


Perfect isotropic associated elastoplasticity with Von-Mises criterion (Prandl-Reuss model)
===========================================================================================

There is no internal variables and we consider an isotropic elastic response. The flow rule reads

.. math:: \dot{\varepsilon}^p = \dot{\gamma} \Frac{\mbox{Dev}(\sigma)}{\|\mbox{Dev}(\sigma)\|}

This corresponds to :math:`\Psi(\sigma) = f(\sigma) = \|\mbox{Dev}(\sigma)\| - \sqrt{\frac{2}{3}}\sigma_y`.


The generalized mid-point scheme for the integration of the plastic flow rule reads:

.. math:: \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \alpha(\sigma_{n+\theta}, A_{n+\theta}) \Delta \xi \Frac{\mbox{Dev}(\sigma_{n+\theta})}{\|\mbox{Dev}(\sigma_{n+\theta})\|}.

Choosing the factor :math:`\alpha(\sigma_{n+\theta}) = \|\mbox{Dev}(\sigma_{n+\theta})\|` and still with :math:`\Delta \xi = \Frac{\Delta \gamma}{\alpha(\sigma_{n+\theta})}` this gives the equation

.. math::  \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \Delta \xi \mbox{Dev}(\sigma_{n+\theta}).

Since :math:`\mbox{Dev}(\sigma_{n+\theta}) = 2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - 2\mu\varepsilon^p_{n+\theta}` this directly gives:

.. math:: \tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = \Frac{1}{1+2\mu\theta\Delta \xi}(\varepsilon^p_{n} + 2\mu\theta\Delta \xi \mbox{Dev}(\varepsilon(u_{n+\theta}))),

which is a linear expression with respect to :math:`u_{n+1}` (but not with respect to :math:`\Delta \xi`).

**Closest point projection approach (elimination of the multiplier)**

The flow rule can be written in term of differential inclusion

.. math:: \dot{\varepsilon}^p \in  \partial I_K(\sigma),

where :math:`K = \left\{ \sigma : \|\mbox{Dev}(\sigma)\| \le \sqrt{\frac{2}{3}}\sigma_y\right\}` is the set of admissible stres tensors, :math:`I_K` is the indicator function of this set and :math:`\partial I_K` its sub-differential (normal cone to :math:`K`). This can be equivalently written

.. math:: \sigma\in K, ~~~(\tau - \sigma):\dot{\varepsilon}^p \le 0 ~~~ \forall \tau \in K.

or

.. math:: {\cal A}(\bar{\tau} - \varepsilon(u) + \varepsilon^p):\dot{\varepsilon}^p \le 0 ~~~ \forall \bar{\tau} \in {\cal A}^{-1}K,

where :math:`{\cal A}` is the fourth order elasticity tensor. Now, in term of projection it can be expressed as

.. math:: \varepsilon(u) - \varepsilon^p = P^{{\cal A}}_{{\cal A}^{-1}K}(\varepsilon(u) - \varepsilon^p + r\dot{\varepsilon}^p),

for any :math:`r > 0` and :math:`P^{{\cal A}}_{{\cal A}^{-1}K}` being the orthogonal projection on :math:`{\cal A}^{-1}K` with respect to the scalar product induced by :math:`{\cal A}`.

The generalize mid-point scheme reads

.. math:: \varepsilon(u_{n+\theta}) - \varepsilon^p_{n+\theta} = P^{{\cal A}}_{{\cal A}^{-1}K}\left(\varepsilon(u_{n+\theta}) - \varepsilon^p_{n+\theta} + \Frac{r}{\theta}(\varepsilon^p_{n+\theta}-\varepsilon^p_{n})\right).

With the choice :math:`\Frac{r}{\theta} = 1` this gives

.. math:: \varepsilon(u_{n+\theta}) = \varepsilon^p_{n} + (I - P^{{\cal A}}_{{\cal A}^{-1}K})(\varepsilon(u_{n+\theta})-\varepsilon^p_{n}).

Since :math:`P^{{\cal A}}_{{\cal A}^{-1}K}` can be expressed as

.. math:: P^{{\cal A}}_{{\cal A}^{-1}K}(\varepsilon) = \Frac{\mbox{tr}(\varepsilon)}{3}I + \min\left(\Frac{1}{2\mu}\sqrt{\frac{2}{3}}\sigma_y, \|\mbox{Dev}(\varepsilon)\|\right) \Frac{\mbox{Dev}(\varepsilon)}{\|\mbox{Dev}(\varepsilon)\|},

We finally find with :math:`B = \mbox{Dev}(\varepsilon(u_{n+\theta}))-\varepsilon^p_{n}`

.. math:: \varepsilon(u_{n+\theta}) = {\mathscr E}^p(u_{n+\theta}, \varepsilon^p_{n}) = \varepsilon^p_{n} + \left( 1 - \sqrt{\frac{2}{3}}\Frac{\sigma_y}{2\mu\|B\|}\right)_+ B

**Plane strain approximation**

The plane strain approximation has the same expression replacing the 3D strain tensors by the in-plane ones :math:`\bar{\varepsilon}^p` and  :math:`\bar{\varepsilon}(u_{n+\theta})`.

.. math:: \bar{{\mathscr E}}^p(\bar{u}_{n+\theta}, \theta \Delta \xi, \bar{\varepsilon}^p_{n}) = \Frac{1}{1+2\mu\theta\Delta \xi}(\bar{\varepsilon}^p_{n} + 2\mu\theta\Delta \xi \mbox{Dev}^*(\bar{\varepsilon}(\bar{u}_{n+\theta}))),

where :math:`\mbox{Dev}^*(\bar{\varepsilon}) = \bar{\varepsilon} - \Frac{\mbox{tr}(\bar{\varepsilon})}{3} \bar{I}` is still the 3D deviator.

Moreover, for the yield condition, 

.. math:: \mbox{Dev}(\sigma) = 2\mu\mbox{Dev}(\varepsilon(u) - \varepsilon^p) = 2\mu\left(\varepsilon(u) - \varepsilon^p - \Frac{\mbox{tr}(\bar{\varepsilon}(u)) - \mbox{tr}(\bar{\varepsilon}^p)}{3} I\right)

.. math:: \begin{array}{rcl} \|\mbox{Dev}(\sigma)\| &=& 2\mu\sqrt{\left\|\bar{\varepsilon}(u) - \bar{\varepsilon}^p - \Frac{\mbox{tr}(\bar{\varepsilon}(u)) - \mbox{tr}(\bar{\varepsilon}^p)}{3} \bar{I}\right\|^2 + \Frac{(\mbox{tr}(\bar{\varepsilon}(u)) - \mbox{tr}(\bar{\varepsilon}^p))^2}{9}} \\ &=& \sqrt{\left\|\bar{\sigma} - \Frac{3\lambda+2\mu}{6(\lambda+\mu)}\mbox{tr}(\bar{\sigma})\bar{I} \right\|^2 + \Frac{\mu^2}{9(\lambda+\mu)^2}\mbox{tr}(\bar{\sigma})^2 } \end{array}

**Plane stress approximation**

For plane stress approximation, we use :eq:`plane_stress_iso` which gives

.. math::  \bar{\varepsilon}^p_{n+\theta} - \bar{\varepsilon}^p_{n} = \theta \Delta \xi \mbox{Dev}^*(\bar{\sigma}_{n+\theta}) =  \theta \Delta \xi \mbox{Dev}^*(\lambda^*\mbox{tr}(\bar{\varepsilon}^e_{n+\theta})\bar{I} + 2\mu \bar{\varepsilon}^e_{n+\theta}) = \theta \Delta \xi\left(\Frac{\lambda^*-2\mu}{3}\mbox{tr}(\bar{\varepsilon}^e_{n+\theta})\bar{I} + 2\mu\bar{\varepsilon}^e_{n+\theta}\right)

thus with :math:`\beta = \Frac{\lambda^*-2\mu}{3}` one has

.. math::  (1+2\mu\theta \Delta \xi)\bar{\varepsilon}^p_{n+\theta} + \beta\theta \Delta \xi \mbox{tr}(\bar{\varepsilon}^p_{n+\theta})\bar{I} = \bar{\varepsilon}^p_{n} + \theta \Delta \xi\left(\beta\mbox{tr}(\bar{\varepsilon}(u_{n+\theta}))\bar{I} + 2\mu\bar{\varepsilon}(u_{n+\theta})\right)

By inverting this relation we find for :math:`A = \bar{\varepsilon}^p_{n} + \theta \Delta \xi\left(\beta\mbox{tr}(\bar{\varepsilon}(u_{n+\theta}))\bar{I} + 2\mu\bar{\varepsilon}(u_{n+\theta})\right)`

.. math::  \bar{{\mathscr E}}^p(\bar{u}_{n+\theta}, \theta \Delta \xi, \bar{\varepsilon}^p_{n}) = \Frac{1}{1+2\mu\theta\Delta \xi} A - \left( \Frac{\beta\theta\Delta \xi}{(1+2\mu\theta\Delta \xi)(1+(2\mu+2\beta)\theta\Delta \xi)} \right) \mbox{tr}(A)\bar{I}


Isotropic elastoplasticity with linear isotropic and kinematic hardening and Von-Mises criterion
================================================================================================

We consider an isotropic elastic reponse and the internal variable :math:`\alpha : \Omega \rightarrow \R` being the accumulated plastic strain which satisfies

.. math:: \dot{\alpha} = \sqrt{\Frac{2}{3}}\dot{\gamma}

For :math:`H_i` the isotropic hardening modulus, the linear hardening consists in

.. math:: \psi^p(\alpha) = \frac{1}{\sqrt{6}}H_i\alpha^2

i.e. :math:`A = \sqrt{\frac{2}{3}}H_i\alpha` and a uniaxial yield stress defined by


.. math:: \sigma_y(a) = \sigma_{y0} + \sqrt{\frac{3}{2}}A = \sigma_{y0} + H_i\alpha,

for :math:`\sigma_{y0}` the initial uniaxial yield stress. The yield function (and plastic potential since this is an associated plastic model) can be defined by

.. math:: \Psi(\sigma, A) = f(\sigma, A) = \|\mbox{Dev}(\sigma - H_k\varepsilon^p)\| - \sqrt{\frac{2}{3}}\sigma_{y0} - A,

where :math:`H_k` is the kinematic hardening modulus. The same computation as in the previous section leads to

.. math:: \tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = \Frac{1}{1+(2\mu+H_k)\theta\Delta \xi}(\varepsilon^p_{n} + 2\mu\theta\Delta \xi \mbox{Dev}(\varepsilon(u_{n+\theta}))),

.. math:: \begin{array}{rcl} \tilde{\mathscr A}(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) &=& \alpha_n + \sqrt{\Frac{2}{3}}\theta \Delta \xi\|\mbox{Dev}(\sigma_{n+\theta} - H_k\varepsilon^p_{n+\theta})\| = \alpha_n + \sqrt{\Frac{2}{3}}\theta \Delta \xi\|2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - (2\mu+H_k)\varepsilon^p_{n+\theta}\| \\ &=&  \alpha_n + \sqrt{\Frac{2}{3}}\Frac{\theta \Delta \xi}{1+(2\mu+H_k)\theta\Delta \xi}\|2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - (2\mu+H_k)\varepsilon^p_{n}\| \\ &=& \alpha_n + \sqrt{\Frac{2}{3}}\|\varepsilon^p_{n+\theta}- \varepsilon^p_{n}\|.\end{array}

Note that the isotropic hardening modulus do not intervene in :math:`\tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n})` but only in :math:`f(\sigma, A)`.

**Closest point projection approach**

Using the same approach as for the perfect isotropic plasticity, the flow rule can be written

.. math:: \dot{\varepsilon}^p \in  \partial I_K(\sigma-H_k \varepsilon^p),

where :math:`K = \left\{ \sigma : \|\mbox{Dev}(\sigma)\| \le \sqrt{\frac{2}{3}}(\sigma_{y0}+H_i\alpha)\right\}`. Writting

.. math:: ({\cal A} + H_k I)(\bar{\tau} - ({\cal A} + H_k I)^{-1}{\cal A}\varepsilon(u) + \varepsilon^p) : \dot{\varepsilon}^p \le 0 ~~ \forall \bar{\tau} \in ({\cal A} + H_k I)^{-1}K

we conclude applying the generalized mid-point scheme by

.. math:: \varepsilon^p_{n+\theta} =  \varepsilon^p_{n} + \left(I-P^{{\cal A} + H_k I}_{({\cal A} + H_k I)^{-1}K}\right)\left(({\cal A} + H_k I)^{-1}{\cal A}\varepsilon(u_{n+\theta}) - \varepsilon^p_n\right).

Now, since

.. math:: \left(I-P^{{\cal A} + H_k I}_{({\cal A} + H_k I)^{-1}K}\right)(\varepsilon) =  \left(1 - \sqrt{\frac{2}{3}}\Frac{(\sigma_{y0}+H_i\alpha_{n+\theta})}{(2\mu+H_k)\|\mbox{Dev}(\varepsilon)\|}\right)_+ \mbox{Dev}(\varepsilon),

and :math:`\mbox{Dev}(({\cal A} + H_k I)^{-1}{\cal A}\varepsilon(u_{n+\theta})) = \Frac{2\mu}{2\mu+H_k}\mbox{Dev}(\varepsilon(u_{n+\theta}))` and with :math:`B = 2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - (2\mu+H_k)\varepsilon^p_n` one has

.. math:: {\mathscr E}^p(u_{n+\theta}, \varepsilon^p_{n}) = \varepsilon^p_{n+\theta} =  \varepsilon^p_{n}  + \Frac{1}{2\mu+H_k}\left(1 - \sqrt{\frac{2}{3}}\Frac{(\sigma_{y0}+H_i\alpha_{n+\theta})}{\|B\|}\right)_+ B.

The problem is not completely solved since :math:`\alpha_{n+\theta}` is still undetermined. However

.. math:: \alpha_{n+\theta} = \alpha_{n} + \sqrt{\Frac{2}{3}}\|\varepsilon^p_{n+\theta} - \varepsilon^p_{n}\| = \alpha_{n} + \sqrt{\Frac{2}{3}}\Frac{1}{2\mu+H_k}\left(\|B\| - \sqrt{\frac{2}{3}}(\sigma_{y0}+H_i\alpha_{n+\theta})\right)_+.

Thus

.. math:: {\mathscr A}(u_{n+\theta}, \varepsilon^p_{n}) = \alpha_{n+\theta} = \max\left(\alpha_{n}, \Frac{\sqrt{\Frac{3}{2}}(2\mu+H_k)\alpha_n+\|B\| - \sqrt{\frac{2}{3}}\sigma_{y0}}{\sqrt{\Frac{3}{2}}(2\mu+H_k)+\sqrt{\frac{2}{3}}H_i}\right),

which complete the expression.

**Plane strain approximation**


The plane strain approximation has the same expression replacing the 3D strain tensors by the in-plane ones :math:`\bar{\varepsilon}^p` and  :math:`\bar{\varepsilon}(u_{n+\theta})`.

Souza-Auricchio elastoplasticity law (for shape memory alloys)
==============================================================

See for instance [GR-ST2015]_ for the justification of the construction of this flow rule. A Von-Mises stress criterion together with an isotropic elastic response, no internal variables and a special type of kinematic hardening is considered with a constraint :math:`\|\varepsilon^p\| \le c_3`. The plastic potential and yield function have the form

.. math:: \Psi(\sigma) = f(\sigma)  = \left\|\mbox{Dev}\left(\sigma - c_1\Frac{\varepsilon^p}{\|\varepsilon^p\|} - c_2\varepsilon^p - \delta \Frac{\varepsilon^p}{\|\varepsilon^p\|}\right)\right\| - \sqrt{\frac{2}{3}}\sigma_{y},

with the complementarity condition

.. math::
   \delta \ge 0, \|\varepsilon^p\| \le c_3,  \delta (\|\varepsilon^p\|-c_3) = 0,

where :math:`c_1, c_2 \mbox{ and } c_3` are some physical parameters. Note that :math:`\Frac{\varepsilon^p}{\|\varepsilon^p\|}` has to be understood to be the whole unit ball for :math:`\varepsilon^p = 0`.


The integration of the flow rule reads

.. math::
   \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \Delta \xi \mbox{Dev}\left(\sigma_{n+\theta} - (c_1 + \delta)\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} - c_2\varepsilon^p_{n+\theta}\right).
   :label: souza_auri_comp

which can be transformed in

.. math::
   (1+(2\mu+c_2)\theta\Delta \xi)\varepsilon^p_{n+\theta} + \theta\Delta \xi(c_1+\delta)\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} = \varepsilon^p_{n} + \theta\Delta \xi 2\mu \mbox{Dev}(\varepsilon(u_{n+\theta})).

With

.. math:: B = \varepsilon^p_{n} + \theta\Delta \xi 2\mu \mbox{Dev}(\varepsilon(u_{n+\theta})),

we conclude that :math:`\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} = \Frac{B}{\|B\|}` and then :math:`\varepsilon^p_{n+\theta} = 0` for :math:`\|B\| \le \theta\Delta \xi c_1` and

.. math:: (1+(2\mu+c_2)\theta\Delta \xi)\varepsilon^p_{n+\theta} = \Frac{B}{\|B\|} (\|B\| - \theta\Delta \xi(c_1+\delta)).

Since :math:`\|\varepsilon^p_{n+\theta}\| = c_3` for :math:`\delta > 0` (complementarity condition), we can deduce the following expression for :math:`\varepsilon^p_{n+\theta}`: 

.. math:: \tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = B \min\left(\Frac{c_3}{\|B\|}, \Frac{\left(1 - \Frac{\theta\Delta \xi c_1}{\|B\|}\right)_+}{1+(2\mu+c_2)\theta\Delta \xi}\right).

The yield condition reads then

.. math:: \begin{array}{l} f(\sigma_{n+\theta}) = \left\|2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - (2\mu+c_2)\varepsilon^p_{n+\theta} - \max\left(c_1, \Frac{\|B\|-c3}{\theta\Delta \xi} - (2\mu+c2)c_3\right)\varepsilon^p_{n+\theta} \right\| \\ ~~~ - \sqrt{\frac{2}{3}}\sigma_{y} \le 0,\end{array}

.. or using :eq:`souza_auri_comp`
.. .. math:: \|\tilde{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) - \varepsilon^p_{n}\| \le \theta\Delta \xi\sqrt{\frac{2}{3}}\sigma_{y}.
.. stupid ? Yes a priori ! 

**Plane strain approximation**


The plane strain approximation has the same expression replacing the 3D strain tensors by the in-plane ones :math:`\bar{\varepsilon}^p` and  :math:`\bar{\varepsilon}(u_{n+\theta})`.



Some classical modelizations
++++++++++++++++++++++++++++






Elasto-plasticity bricks 
+++++++++++++++++++++++++



Generic brick
=============

The generic brick add the following terms:

.. math:: \int_{\Omega} \sigma_{n+1} : \nabla v dx +  \ds \int_{\Omega} (\Delta \xi - (\Delta \xi + r f(\sigma_{n+\theta}, A_{n+\theta}))_+) \lambda dx = 0,   \forall \lambda


Other bricks
============

to be done: ::

      getfem::add_elastoplasticity_brick
          (md, mim, ACP, varname, datalambda, datamu, datathreshold, datasigma, region);





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


Other useful functions
**********************

The function: ::

      getfem::elastoplasticity_next_iter
          (md, mim, varname, previous_varname, ACP, datalambda, datamu, datathreshold, datasigma);

computes the new stress constraint values supported by the material after a load or an unload (once a solve has been done earlier) and upload the variables ``varname`` and ``datasigma`` as follows:

.. math::
   
   u^{n+1} \Rightarrow u^n \ \ \ \ \ \textrm{ and } \ \ \ \ \ \sigma^{n+1} \Rightarrow \sigma^n

Then, :math:`u^n` and :math:`\sigma^n` contains the new values computed and one can restart the process.



########################


The function: ::

      getfem::compute_elastoplasticity_Von_Mises_or_Tresca
          (md, datasigma, mf_vm, VM, tresca=false);

computes the Von Mises (or Tresca if ``tresca`` = true) criterion on the stress tensor stored in ``datasigma`` . The stress is evaluated on the `mesh_fem` ``mf_vm`` and stored into the vector ``VM``.
Of course, this function can be used if and only if the previous function ``elastoplasticity_next_iter`` has been called earlier.



##########################


The function: ::

      getfem::compute_plastic_part
          (md, mim, mf_pl, varname, previous_varname, ACP, datalambda, datamu, datathreshold, datasigma, Plast);

computes on ``mf_pl`` the plastic part of the material, that could appear after a load and an unload, into the vector ``Plast``. 

Note that ``datasigma`` should be the vector containing the new stress constraint values, i.e. after a load or an unload of the material.



