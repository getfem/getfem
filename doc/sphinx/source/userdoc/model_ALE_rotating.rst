.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. index:: models, model bricks

.. _ud-model-ALE_rotating:


ALE Support for object having a large rigid body motion
*******************************************************

ALE terms for rotating objects
------------------------------

This section present a set of bricks facilitating the use of an ALE formulation for rotating bodies having a rotational symmetry (typically a train wheel).

Theoretical background
++++++++++++++++++++++

This strategy consists in adopting an intermediary description between an Eulerian and a Lagrangian ones for a rotating body having a rotational symmetry. This intermediary description consist in a rotating axes with respect to the reference configuration. See for instance [Dr-La-Ek2014]_ and [Nackenhorst2004]_.

It is supposed that the considered body is submitted approximately to a rigid body motion

.. math::
  \tau(X) = R(t)X + Z(t)

and may have additonal deformation (exptected smaller) with respect to this rigid motion, where :math:`R(t)` is a rotation matrix

.. math::
  R(t) = \left(\begin{array}{ccc}
  \cos(\theta(t)) & \sin(\theta(t)) & 0 \\
  -\sin(\theta(t)) & \cos(\theta(t)) & 0 \\
  0 & 0 & 1
  \end{array} \right),

and :math:`Z(t)` is a translation. Note that, in order to be consistent with a positive translation for a positive angle for a rolling contact, the rotation is **clockwise**. This illustrated in the following figure:

.. _ud-fig-rotating_cylinder:

.. figure:: images/ALE_rotating_body.png
    :width: 40%
    :align: center
    :alt: alternate text
    :figclass: align-center

Note that the description is given for a three-dimensional body. For two-dimensional bodies, the third axes is neglected so that :math:`R(t)` is a :math:`2\times 2` rotation matrix. Let us denote :math:`r(t)` the rotation:

.. math::
   r(t,X) = R(t)X, ~~~~~~~~~ \mbox{ and }
   A = \left(\begin{array}{ccc}
   0 & 1 & 0 \\
   -1 & 0 & 0 \\
   0 & 0 & 0
   \end{array} \right).

We have then

.. math::
 \dot{r}(t,X) = \dot{\theta}AR(t)X

If :math:`\varphi(t, X)` is the deformation of the body which maps the reference configuration :math:`\Omega^0` to the deformed configuration :math:`\Omega_t` at time :math:`t`, the ALE description consists in the decomposition of the deformation of the cylinder in

.. math::
   \varphi(t, X) = (\tau(t) \circ \bar{\varphi}(t) \circ r(t))(X) = \bar{\varphi}(t, r(t, X)) + Z(t)

With :math:`\bar{X} = R(t)X` the new considered deformation is

.. math::
  \bar{\varphi}(t,\bar{X}) = \varphi(X) - Z(t)


Thanks to the rotation symmetry of the reference configuration :math:`\Omega^0:`, we note that :math:`\bar{\Omega}^0 = r(t, \Omega^0)` is independant of :math:`t` and will serve as the new reference configuration. This is illustrated in the following figure:

.. _ud-fig-rotating_cylinder_conf:

.. figure:: images/ALE_rotating_conf.png
    :width: 80%
    :align: center
    :alt: alternate text
    :figclass: align-center

The denomination ALE of the method is justified by the fact that :math:`\bar{\Omega}^0` is an intermediate configuration which is of Euler type for the rigid motion and a Lagrangian one for the additional deformation of the solid. If we denote

.. math::
  \bar{u}(t,\bar{X}) = \bar{\varphi}(t, \bar{X}) - \bar{X}

the displacement with respect to this intermediate configuration, the advantage is that if this additional displacement with respect to the rigid body motion is small, it is possible to use a small deformation model (for instance linearized elasticity).

Due to the objectivity properties of standard constitutive laws, the expression of these laws in the intermediate configuration is most of the time identical to the expression in a standard reference configuration except for the expression of the time derivative which are modified because the change of coordinate is  nonconstant in time :

.. math::

  \dfrac{\partial \varphi}{\partial t} = \dfrac{\partial \bar{\varphi}}{\partial t} + \dot{\theta} \nabla \bar{\varphi} A \bar{X} + \dot{Z}(t),

  \dfrac{\partial^2 \varphi}{\partial t^2} = \dfrac{\partial^2 \bar{\varphi}}{\partial t^2} + 2\dot{\theta} \nabla\dfrac{\partial \bar{\varphi}}{\partial t}A \bar{X} + \dot{\theta}^2\mbox{div}((\nabla\bar{\varphi}A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta}\nabla\bar{\varphi}A \bar{X} + \ddot{Z}(t).

Note that the term :math:`\dot{\theta} A \bar{X} = \left(\hspace*{-0.5em}\begin{array}{c} \dot{\theta}\bar{X}_2 \\ -\dot{\theta}\bar{X}_1 \\ 0 \end{array}\hspace*{-0.5em}\right)` is the rigid motion velocity vector. Now, If :math:`\Theta(t,X)` is a quantity attached to the material points (for instance the temperature), then, with :math:`\bar{\Theta}(t,\bar{X}) = \Theta(t,X)` , one simply has

.. math::

  \dfrac{\partial \Theta}{\partial t} = \dfrac{\partial \bar{\Theta}}{\partial t} + \dot{\theta} \nabla \bar{\Theta} A \bar{X}

This should not be forgotten that a correction has to be provided for each evolving variable for which the time derivative intervene in the considered model (think for instance to platic flow for plasticity). So that certain model bricks canot be used directly (plastic bricks for instance).

|gf| bricks for structural mechanics are mainly considering the displacement as the amin unknown. The expression for the displacement is the following:

.. math::
  \dfrac{\partial u}{\partial t} = \dfrac{\partial \bar{u}}{\partial t} + \dot{\theta} (I_d + \nabla \bar{u}) A \bar{X} + \dot{Z}(t),

  \dfrac{\partial^2 u}{\partial t^2} = \dfrac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\dfrac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{Z}(t).





Weak formulation of the transient terms
#######################################

Assuming :math:`\rho^0` the density in the reference configuration having a rotation symmetry, the term corresponding to acceleration in the weak formulation reads (with :math:`v(X) = \bar{v}(\bar{X})` a test function):

.. math::
   \int_{\Omega^0} \rho^0 \dfrac{\partial^2 u}{\partial t^2}\cdot vdX =

   \int_{\bar{\Omega}^0} \rho^0 \left[\dfrac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\dfrac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{Z}(t) \right] \cdot \bar{v} d\bar{X}.

The third term in the right hand side can be integrated by part as follows:

.. math::
   \begin{array}{rcl}
    \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} &=& -  \int_{\bar{\Omega}^0} (\dot{\theta}^2 (I_d + \nabla\bar{u})A \bar{X})) \cdot (\nabla (\rho^0 \bar{v}) A \bar{X}) d\bar{X} \\
   && + \int_{\partial \bar{\Omega}^0} \rho^0 \dot{\theta}^2 (((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \bar{N} \cdot \bar{v} d\bar{\Gamma}.
  \end{array}

Since :math:`\bar{N}` the outward unit normal vector on :math:`\partial \bar{\Omega}^0` is orthogonal to :math:`A \bar{X}` the boundary term is zero and :math:`\nabla (\rho^0 \bar{v}) = \bar{v} \otimes \nabla \rho^0   + \rho^0 \nabla \bar{v}` and since :math:`\nabla \rho^0.(A\bar{X}) = 0` because of the assumption on :math:`\rho^0` to have a rotation symmetry, we have

.. math::
   \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} = - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2 (A^2 \bar{X})\cdot \bar{v} d\bar{X}.

Thus globally

.. math::
   \begin{array}{rcl}
    \int_{\Omega^0} \rho^0 \dfrac{\partial^2 u}{\partial t^2}\cdot vdX &=&
    \int_{\bar{\Omega}^0} \rho^0 \left[\dfrac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\dfrac{\partial \bar{u}}{\partial t}A \bar{X} + \ddot{\theta} \nabla\bar{u} A \bar{X}   \right] \cdot \bar{v} d\bar{X}\\
   && - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 (\dot{\theta}^2 A^2 \bar{X} + \ddot{\theta} A\bar{X} + \ddot{Z}(t))\cdot \bar{v} d\bar{X}.
   \end{array}

Note that two terms can deteriorate the coercivity of the problem and thus its well posedness and the stability of time integration schemes: the second one (convection term) and the fifth one. This may oblige to use additional stabilization techniques for large values of the angular velocity :math:`\dot{\theta}`.


The available bricks ...
++++++++++++++++++++++++

To be adapted  ::

  ind = getfem::brick_name(parmeters);

where ``parameters`` are the parameters ...





ALE terms for a uniformly translated part of an object
------------------------------------------------------

This section present a set of bricks facilitating the use of an ALE formulation for an object being potentially infinite in one direction and which whose part of interests (on which the computation is considered) is translated uniformly in that direction (typically a bar).

Theoretical background
++++++++++++++++++++++

Let us consider an object whose reference configuration :math:`\Omega^0 \in \rm I\hspace{-0.15em}R^{d}` is infinite in the direction :math:`E_1`, i.e. :math:`\Omega^0 = \rm I\hspace{-0.15em}R \times \omega^0` where :math:`\omega^0 \in \rm I\hspace{-0.15em}R^{d-1}`. At a time :math:`t`, only a "windows" of this object is considered

.. math::
   \Omega^{0t} = (\alpha + z(t), \beta + z(t)) \times \omega^0

where :math:`z(t)` represents the translation.

If :math:`\varphi(t, X)` is the deformation of the body which maps the reference configuration :math:`\Omega^0` to the deformed configuration :math:`\Omega_t` at time :math:`t`, the ALE description consists in considering the intermediary reference configuration

.. math::
   \bar{\Omega}^{0} = (\alpha, \beta) \times \omega^0

and :math:`\bar{\varphi}(t, X) : \rm I\hspace{-0.15em}R_+ \times \bar{\Omega}^{0} \rightarrow \rm I\hspace{-0.15em}R^d` defined by

.. math::
  \bar{\varphi}(t,\bar{X}) = \varphi(t,X), ~~~\mbox{ with } \bar{X} = X - Z(t),

where :math:`Z(t) = z(t)E_1`. The interest of :math:`\bar{\Omega}^{0}` is of course to be time independant. Of course, some special boundary conditions have to be defined on :math:`\{\alpha\} \times \omega^0` and :math:`\{\beta\} \times \omega^0` (absorbing or periodic boundary conditions) in order to approximate the fact that the body is infinite.

.. _ud-fig-translating_bar:

.. figure:: images/ALE_translation_body.png
    :width: 40%
    :align: center
    :alt: alternate text
    :figclass: align-center


If we denote

.. math::
  \bar{u}(t,\bar{X}) = \bar{\varphi}(t, \bar{X}) - X = u(t, X),

the displacement on the intermediary configuration, then it is easy to check that

.. math::
   \dfrac{\partial \varphi}{\partial t} = \dfrac{\partial \bar{u}}{\partial t} - \nabla \bar{u} \dot{Z}

   \dfrac{\partial^2 \varphi}{\partial t^2} = \dfrac{\partial^2 \bar{u}}{\partial t^2} - \nabla\dfrac{\partial \bar{u}}{\partial t}\dot{Z} + \dfrac{\partial^2 \bar{u}}{\partial \dot{Z}^2} - \nabla\bar{u}\ddot{Z}.




Weak formulation of the transient terms
#######################################

Assuming :math:`\rho^0` the density in the reference being invariant with the considered translation, the term corresponding to acceleration in the weak formulation reads (with :math:`v(X) = \bar{v}(\bar{X})` a test function and after integration by part):

.. math::
   \int_{\Omega^0} \rho^0 \dfrac{\partial^2 u}{\partial t^2}\cdot vdX =

   \int_{\bar{\Omega}^{0}} \rho^0 \left[\dfrac{\partial^2 \bar{u}}{\partial t^2} - 2\nabla\dfrac{\partial \bar{u}}{\partial t}\dot{Z} - \nabla\bar{u}\ddot{Z}\right]\cdot \bar{v}  - \rho^0 (\nabla\bar{u}\dot{Z}).(\nabla\bar{v}\dot{Z}) d\bar{X} + \int_{\partial \bar{\Omega}^0} \rho^0 (\nabla\bar{u}\dot{Z}).\bar{v}(\dot{Z}.\bar{N}) d\bar{\Gamma},

where :math:`\bar{N}` is the outward unit normal vector on :math:`\partial \bar{\Omega}^0`. Note that the last term vanishes on :math:`(\alpha, \beta) \times \partial \omega^0` but not necessarily on :math:`\{\alpha\} \times \omega^0` and :math:`\{\beta\} \times \omega^0`.






