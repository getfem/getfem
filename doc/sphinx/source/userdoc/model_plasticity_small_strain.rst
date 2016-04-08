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

We present a short introduction to small strain plasticity (mainly elastoplasticity). We refer mainly to [SI-HU1998]_ and [SO-PE-OW2008]_ for a more detailed presentation.

Additive decomposition of the small strain tensor
=================================================

Let :math:`\Omega \subset \R^3` be the reference configuration of a deformable body and :math:`u : \Omega \rightarrow \R^3` be the displacement. Small strain plasticity is based on the additive decomposition of the small deformation tensor :math:`\varepsilon(u) = \Frac{\nabla u + \nabla u^T}{2}` in

.. math::
   \varepsilon(u) = \varepsilon^e + \varepsilon^p

where :math:`\varepsilon^e` is the elastic part of the deformation and :math:`\varepsilon^p` the plastic one.

Internal variables, free energy potential and elastic law
=========================================================

We consider

.. math::
   \alpha : \Omega \rightarrow \R^{d_{\alpha}},

a vector field of :math:`d_{\alpha}` strain type internal variables (:math:`d_{\alpha} = 0` if no internal variables are considered). We consider also a free energy potential

.. math::
   \psi(\varepsilon^e, \alpha),

such that the stress type variables are determined by

.. math::
   \sigma = \Frac{\partial \psi}{\partial \varepsilon^e}(\varepsilon^e, \alpha), ~~~~ A =  \Frac{\partial \psi}{\partial \alpha}(\varepsilon^e, \alpha),

where :math:`\sigma` is the Cauchy stress tensor and :math:`A` the stress type internal variables. The plastic dissipation being given by

.. math::

   \sigma:\dot{\varepsilon}^p - A.\dot{\alpha} \ge 0.

In the standard cases, :math:`\psi(\varepsilon^e, \alpha)` is decomposed into

.. math:: \psi(\varepsilon^e, \alpha) = \psi^e(\varepsilon^e) + \psi^p(\alpha).

In the case of linearized elasticity, one has :math:`\psi^e(\varepsilon^e) = \frac{1}{2} ({\cal A}\varepsilon^e) :\varepsilon^e` for :math:`{\cal A}` the fourth order elasticity tensor en more precisely :math:`\psi^e(\varepsilon^e) = \mu \mbox{dev}(\varepsilon^e) : \mbox{dev}(\varepsilon^e) + \frac{1}{2} K (tr(\varepsilon^e))^2` for isotropic linearized elasticity for :math:`\mu, K = \lambda + 2\mu/3` the shear and bulk modulus, respectively.



Plastic potential, yield function and plastic flow rule
=======================================================

The plastic deformation is supposed to occurs when the stress attains a critical value. This is determinated by a yield function :math:`f(\sigma, A)` and the condition

.. math:: f(\sigma, A) \le 0.

The surface :math:`f(\sigma, A) = 0` being the yield surface where the plastic deformation may occur.

Let us consider also the plastic potential :math:`\Psi(\sigma, A)`, (convex with respect to both its two variables) which determine the plastic flow direction in the sense that the flow rule reads as

.. math:: \dot{\varepsilon}^p = \dot{\gamma} \Frac{\Psi}{\partial \sigma}(\sigma, A), ~~~~~~ \dot{\alpha} = -\dot{\gamma} \Frac{\Psi}{\partial A}(\sigma, A),

with the additional complementary condition

.. math:: f(\sigma, A) \le 0, ~~~ \dot{\gamma} \ge 0, ~~~ f(\sigma, A) \dot{\gamma} = 0.

The variable :math:`\dot{\gamma}` is called the plastic multiplier. Note that when :math:`\psi(\varepsilon^e, \alpha), f(\sigma, A) \mbox{ or } \Psi(\sigma, A)` are not differentiable, subdifferentials have to be used. Associated plasticity corresponds to the choice :math:`\Psi(\sigma, A) = f(\sigma, A)`.

Initial boundary value problem
++++++++++++++++++++++++++++++

The weak formulation of a dynamic elastoplastic problem can be written as follows for an arbitrary cinematically admissible test function :math:`v`:

.. math::

   \left| \begin{array}{l}
   \ds \int_{\Omega} \rho \ddot{u}\cdot v + \sigma : \nabla v dx =  \int_{\Omega} f\dot v dx + \int_{\Gamma_N} g\dot v dx, \\
   u(0,x) = u_0(x), ~~~\dot{u}(0) = v_0(x), \\
   \varepsilon^p(0,x) = \varepsilon^p_0, ~~~ \alpha(0,x) = \alpha_0,
   \end{array} \right.

for :math:`u_0, v_0, \varepsilon^p_0, \alpha_0` the initial values and :math:`g` the force prescribed on the part of the boundary :math:`\Gamma_N`.

Note that plasticity models are often applied on quasitistic problem which corresponds to neglect the term :math:`\rho \ddot{u}`.

Given a time step :math:`\Delta t` we will denote in the sequel :math:`u_n, \varepsilon^p_n  \mbox{ and } \alpha_n` the approximation at time :math:`n\Delta t` of :math:`u(t), \varepsilon^p_n \mbox{ and } \alpha(t)` respectively. This approximation is given by a chose time integration scheme (for instance one of the proposed schemes in :ref:`ud-model-time-integration`) which can be different than the time integration scheme used for the integration of the flow rule (see below).


Flow rule integration
+++++++++++++++++++++

The plastic flow rule have to be integrated with its own time integration scheme.  Among standards schemes, backward Euler scheme, :math:`\theta`-scheme and generalized mid-point scheme are the most commonly used in that context. We make here the choice of the generalized mid-point scheme.


Let :math:`u_{n+1}` be the displacement at the considered time step and  :math:`u_{n}` at the previous one. For a quantity :math:`B` we denote :math:`B_{n+\theta} = \theta B_{n+1} + (1-\theta)B_n` the convex combination of the quantity at iterations :math:`n` and :math:`n+1`.

The mid-point scheme for the integration of the plastic flow rules reads as

.. math:: \varepsilon^p_{n+1} - \varepsilon^p_{n} = \Delta \gamma \Frac{\Psi}{\partial \sigma}(\sigma_{n+\theta}, A_{n+\theta}),

.. math:: \alpha_{n+1} - \alpha_n = -\Delta \gamma \Frac{\Psi}{\partial A}(\sigma_{n+\theta}, A_{n+\theta}),

with the complementary condition

.. math:: f(\sigma_{n+\theta}, A_{n+\theta}) \le 0, ~~~ \Delta\gamma \ge 0, ~~~ f(\sigma_{n+\theta}, A_{n+\theta}) \Delta \gamma = 0.

where :math:`0 < \theta \le 1` is the parameter of the mid-point scheme. We exclude :math:`\theta = 0` because we will not consider explicit integration of plasticity. Let us recall that :math:`\theta = 1` corresponds to the backward Euler scheme and :math:`\theta = 1/2` to the mid-point scheme which is a second order consistent scheme.

A solution would be to solve the whole problem with all the unknows, that is :math:`u_{n+1}, \Delta \gamma, \varepsilon^p_{n+1} \mbox{ and } A_{n+1}`. This is of course possible but would be a rather expensive strategy because of the resulting high number of degrees of freedom. A classical strategy (the return mapping one for instance, see [SO-PE-OW2008]_) consist in integrating locally the plastic flow on each Gauss point of the considered integration method separately, or more precisely to consider on each Gauss point the map :math:`F : u_{n+1} \mapsto \sigma_{n+1}` which results from the local flow rule integration. Both this map and its tangent modulus (usually called consistent tangent modulus) is then used in the global solve of the problem with a Newton method and for :math:`u_{n+1}` the unique remaining variable. The advantage of the return mapping strategy is that the unique variable of the global solve is the displacement :math:`u_{n+1}`. However, a nonlinear solve on each Gauss point is necessary which is often performed with a local Newton method, and the expression of the consistent tangent modulus is often quite complex.

The approach developped below is mainly inspired from  [PO-NI2016]_,  [SE-PO-WO2015]_ and [HA-WO2009]_ and allow more simple tangent moduli. It consists in keeping (a multiple of) :math:`\Delta \gamma` as an additional unknown with respect to :math:`u_{n+1}`. As we will see, this will allow a more generic treatment of the yield functions, the price for the simplicity being this additional unknown scalar field.

First, we consider an additional (and optional) given function :math:`\alpha(\sigma_{n+\theta}, A_{n+\theta}) > 0` whose interest will appear later on (it will allow simple local inverses) and the new unknown scalar field

.. math:: \Delta \xi = \Frac{\Delta \gamma}{\alpha(\sigma_{n+\theta}, A_{n+\theta})} ,

so that our two main unknows are now :math:`u_{n+1} \mbox{ and } \Delta \xi`. The plastic flow rule integration may now read:

.. math:: \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \alpha(\sigma_{n+\theta}, A_{n+\theta}) \theta \Delta \xi \Frac{\Psi}{\partial \sigma}(\sigma_{n+\theta}, A_{n+\theta}).
   :label: flowrule1

.. math::  \alpha_{n+\theta} - \alpha_n = -\alpha(\sigma_{n+\theta}, A_{n+\theta}) \theta \Delta \xi \Frac{\Psi}{\partial A}(\sigma_{n+\theta}, A_{n+\theta}),
   :label: flowrule2

.. math:: f(\sigma_{n+\theta}, A_{n+\theta}) \le 0, ~~~ \Delta\xi \ge 0, ~~~ f(\sigma_{n+\theta}, A_{n+\theta}) \Delta \xi = 0.
   :label: flowrule3

For :math:`u_{n+1} \mbox{ and } \Delta \xi` be given, we define the two maps

.. math::
   {\mathscr E}^p : (u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) \mapsto \varepsilon^p_{n+\theta}

   {\mathscr A} : (u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) \mapsto \alpha_{n+\theta}


where the pair :math:`(\varepsilon^p_{n+\theta}, \alpha_{n+\theta}) = ({\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n), {\mathscr A}(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n))` is the solution to equations :eq:`flowrule1`, :eq:`flowrule2` (without the consideration of  :eq:`flowrule3`). We will see later, that, at least for simple isotropic plastic flow rules, these maps have a simple expression, even sometimes a linear one with respect to :math:`u_{n+\theta}`. Most of the time, the map :math:`{\mathscr A}` can be optained as a simple post-treatment of the map :math:`{\mathscr E}^p`.

Still :math:`u_{n+1} \mbox{ and } \Delta \xi` be given the stress :math:`\sigma_{n+1}` reads

.. math:: \sigma_{n+1} = \Frac{\partial \psi^e}{\partial \varepsilon^e}(\varepsilon(u_{n+1}) -\varepsilon^p_{n+1}).

.. math:: A_{n+1} = \Frac{\partial \psi^p}{\partial \alpha}(\alpha_{n+1}).

The complementarity equation :eq:`flowrule3` is then prescribed with the use of a well chosen complementarity function, as in [HA-WO2009]_ for :math:`r > 0` such as:

.. math:: \ds \int_{\Omega} ((-f(\sigma_{n+\theta}, A_{n+\theta}) - r\Delta \xi)_+ + f(\sigma_{n+\theta}, A_{n+\theta}) ) \lambda dx = 0,   \forall \lambda


pb : need of :math:`A_{n+\theta}` ... Would be more simple in :math:`f(\sigma_{n+1}, A_{n+1})`

and ... the following should be more simple :

.. math:: \ds \int_{\Omega} (\Delta \xi - (\Delta \xi + r f(\sigma_{n+\theta}, A_{n+\theta}))_+) \lambda dx = 0,   \forall \lambda

Plane strain approximation
++++++++++++++++++++++++++

To be described

Plane stress approximation
++++++++++++++++++++++++++

To be described

Some classical laws
+++++++++++++++++++

Von Mises vs Tresca threshold ...


Tresca : :math:`\rho(\sigma) \le \sigma_y` where :math:`\rho(\sigma)` spectral radius of the Cauchy stress tensor and :math:`\sigma_y` the uniaxial yield stress (which may depend on some hardening internal variables.

Von Mises :  :math:`\|\mbox{Dev}(\sigma)\| \le \sqrt{\frac{2}{3}}\sigma_y` where
:math:`\mbox{Dev}(\sigma) = \sigma - \frac{1}{3}\mbox{Tr}(\sigma)I` the deviatoric part of :math:`\sigma` and :math:`\|\sigma\| = \sqrt{\sigma:\sigma}`.


Perfect isotropic associated elastoplasticity with Von-Mises criterion (Prandl-Reuss model)
===========================================================================================

There is no internal variables and we consider an isotropic elastic response. The flow rule reads

.. math:: \dot{\varepsilon}^p = \dot{\gamma} \Frac{\mbox{Dev}(\sigma)}{\|\mbox{Dev}(\sigma)\|}

This corresponds to :math:`\Psi(\sigma) = f(\sigma) = \|\mbox{Dev}(\sigma)\| - \sqrt{\frac{2}{3}}\sigma_y`.


The mid-point scheme for the integration of the plastic flow rule reads:

.. math:: \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \alpha(\sigma_{n+\theta}, A_{n+\theta}) \Delta \xi \sqrt{\frac{3}{2}}\Frac{\mbox{Dev}(\sigma_{n+\theta})}{\|\mbox{Dev}(\sigma_{n+\theta})\|}.

Choosing the factor :math:`\alpha(\sigma_{n+\theta}) = \|\mbox{Dev}(\sigma_{n+\theta})\|` and still with :math:`\Delta \xi = \Frac{\Delta \gamma}{\alpha(\sigma_{n+\theta})}` this gives the equation

.. math::  \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \Delta \xi \mbox{Dev}(\sigma_{n+\theta}).

Since :math:`\mbox{Dev}(\sigma_{n+\theta}) = 2\mu\mbox{Dev}(\varepsilon(u_{n+\theta})) - 2\mu\varepsilon^p_{n+\theta}` this directly gives:

.. math:: {\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = \Frac{1}{1+2\mu\theta\Delta \xi}(\varepsilon^p_{n} + 2\mu\theta\Delta \xi \mbox{Dev}(\varepsilon(u_{n+\theta}))),

which is a linear expression with respect to :math:`u_{n+1}` (but not with respect to :math:`\Delta \xi`).






Isotropic elastoplasticity with linear isotropic and kinematic hardening and Von-Mises criterion
================================================================================================

We consider an isotropic elastic reponse and the internal variable :math:`\alpha : \Omega \rightarrow \R` being the accumulated plastic strain which satisfies

.. math:: \dot{\alpha} = \dot{\gamma}

For :math:`H_i` the isotropic hardening modulus, the linear hardening consists in

.. math:: \psi^p(\alpha) = \frac{1}{\sqrt{6}}H_i\alpha^2

i.e. :math:`A = \sqrt{\frac{2}{3}}H_i\alpha` and a uniaxial yield stress defined by


.. math:: \sigma_y(a) = \sigma_{y0} + \sqrt{\frac{3}{2}}A = \sigma_{y0} + H_i\alpha,

for :math:`\sigma_{y0}` the initial uniaxial yield stress. The yield function (and plastic potential since this is a associated plastic model) can be defined by

.. math:: \Psi(\sigma, A) = f(\sigma, A) = \|\mbox{Dev}(\sigma - H_k\varepsilon^p)\| - \sqrt{\frac{2}{3}}\sigma_{y0} + A,

where :math:`H_k` is the kinematic hardening modulus. The same computation as in the previous section leads to

.. math:: {\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = \Frac{1}{1+(2\mu+H_k)\theta\Delta \xi}(\varepsilon^p_{n} + (2\mu)\theta\Delta \xi \mbox{Dev}(\varepsilon(u_{n+\theta}))),

.. math:: {\mathscr A}(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}, \alpha_n) = \alpha_n + \|\varepsilon^p_{n+\theta}-\varepsilon^p_{n}\| = \alpha_n + \| {\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n})-\varepsilon^p_{n}\|

Note that the isotropic Hardening modulus do not intervene in :math:`{\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n})` but only in :math:`f(\sigma, A)`.

Souza-Auricchio elastoplasticity law (for shape memory alloys)
==============================================================

See for instance [GR-ST2015]_ for the justification of the construction of this flow rule. A Von-Mises stress criterion together with an isotropic elastic response, no internal variables and a special type of kinematic hardening is considered with a constraint :math:`\|\varepsilon^p\| \le c_3`. The yield function has the form

.. math:: \Psi(\sigma, A) = f(\sigma, A) = \left\|\mbox{Dev}\left(\sigma - c_1\Frac{\varepsilon^p}{\|\varepsilon^p\|} - c_2\varepsilon^p - \delta \Frac{\varepsilon^p}{\|\varepsilon^p\|}\right)\right\| - \sqrt{\frac{2}{3}}\sigma_{y},

with the complementarity condition

.. math::
   \delta \ge 0, \|\varepsilon^p\| \le c_3,  \delta \|\varepsilon^p\| = 0,

where :math:`c_1, c_2 \mbox{ and } c_3` are some physical parameters. Note that :math:`\Frac{\varepsilon^p}{\|\varepsilon^p\|}` has to be understood as the wall unit ball for :math:`\varepsilon^p = 0`.


The integration of the flow rule reads

.. math::
   \varepsilon^p_{n+\theta} - \varepsilon^p_{n} = \theta \Delta \xi \mbox{Dev}\left(\sigma_{n+\theta} - (c_1 + \delta)\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} - c_2\varepsilon^p_{n+\theta}\right).

which can be transformed in

.. math::
   (1+(2\mu+c_2)\theta\Delta \xi)\varepsilon^p_{n+\theta} + \theta\Delta \xi(c_1+\delta)\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} = \varepsilon^p_{n} + \theta\Delta \xi 2\mu \mbox{Dev}(\varepsilon(u_{n+\theta})).

With

.. math:: B = \varepsilon^p_{n} + \theta\Delta \xi 2\mu \mbox{Dev}(\varepsilon(u_{n+\theta})),

we conclude that :math:`\Frac{\varepsilon^p_{n+\theta}}{\|\varepsilon^p_{n+\theta}\|} = \Frac{B}{\|B\|}` and then :math:`\varepsilon^p_{n+\theta} = 0` for :math:`\|B\| \le c_1` and

.. math:: (1+(2\mu+c_2)\theta\Delta \xi)\varepsilon^p_{n+\theta} = \Frac{B}{\|B\|} (\|B\| - \theta\Delta \xi(c_1+\delta))

Since :math:`\varepsilon^p_{n+\theta} = c_3` for :math:`\delta > 0` (complementarity condition), we can deduce the follwoing expression for :math:`\varepsilon^p_{n+\theta}`: 

.. math:: {\mathscr E}^p(u_{n+\theta}, \theta \Delta \xi, \varepsilon^p_{n}) = \Frac{B}{\|B\|} \max\left(0, \min\left(c_3, \Frac{\|B\| - \theta\Delta \xi(c_1+\delta)}{1+(2\mu+c_2)\theta\Delta \xi}\right)\right).




Some classical modelizations
++++++++++++++++++++++++++++






Elasto-plasticity bricks 
+++++++++++++++++++++++++

to be done: ::

      getfem::add_elastoplasticity_brick
          (md, mim, ACP, varname, datalambda, datamu, datathreshold, datasigma, region);




