.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-ALE_rotating:

ALE terms for rotating objects
------------------------------

This section present a set of bricks facilitating the use of an ALE formulation for rotating bodies having a rotational symmetry (typically a train wheel).

Theoretical background
++++++++++++++++++++++

This strategy consists in adopting an intermediary description between an Eulerian and a Lagrangian ones for a rotating body having a rotational symmetry. This intermediary description consist in a rotating axes with respect to the reference configuration. See for instance [Dr-La-Ek2014]_ and [Nackenhorst2004]_.

It is supposed that the considered body is submitted approximately to a rigid body motion

.. math::
  \tau(X) = R(t)X + z(t)

and may have additonal deformation (exptected smaller) with respect to this rigid motion, where :math:`R(t)` is a rotation matrix

.. math::
  R(t) = \left(\begin{array}{ccc}
  \cos(\theta(t)) & \sin(\theta(t)) & 0 \\
  -\sin(\theta(t)) & \cos(\theta(t)) & 0 \\
  0 & 0 & 1
  \end{array} \right),

and :math:`z(t)` is a translation. Note that, in order to be consistent with a positive translation for a positive angle, the rotation is **clockwise**. This illustrated in the following figure:

.. _ud-fig-rotating_cylinder:

.. figure:: images/ALE_rotating_body.png
    :height: 200pt
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
   \varphi(t, X) = (\tau(t) \circ \bar{\varphi}(t) \circ r(t))(X) = \bar{\varphi}(t, r(t, X)) + z(t)

With :math:`\bar{X} = R(t)X` the new considered deformation is

.. math::
  \bar{\varphi}(t,\bar{X}) = \varphi(X) - z(t)


Thanks to the rotation symmetry of the reference configuration :math:`\Omega^0:`, we note that :math:`\bar{\Omega}^0 = r(t, \Omega^0)` is independant of :math:`t` and will serve as the new reference configuration. This is illustrated in the following figure: 

.. _ud-fig-rotating_cylinder_conf:

.. figure:: images/ALE_rotating_conf.png
    :height: 200pt
    :align: center
    :alt: alternate text
    :figclass: align-center

The denomination ALE of the method is justified by the fact that :math:`\bar{\Omega}^0` is an intermediate configuration which is of Euler type for the rigid motion and a Lagrangian one for the additional deformation of the solid. If we denote

.. math::
  \bar{u}(t,\bar{X}) = \bar{\varphi}(t, \bar{X}) - \bar{X}

the displacement with respect to this intermediate configuration, the advantage is that if this additional displacement with respect to the rigid body motion is small, it is possible to use a small deformation model (for instance linearized elasticity).

Due to the objectivity properties of standard consistutive laws, the expression od these laws in the intermediate configuration is most of the time identical to the expression in a standard reference configuration except for the expression of the time derivative which are modified because the change of coordinate is  nonconstant in time :

.. math::

  \Frac{\partial \varphi}{\partial t} = \Frac{\partial \bar{\varphi}}{\partial t} + \dot{\theta} \nabla \bar{\varphi} A \bar{X} + \dot{z}(t),

  \Frac{\partial^2 \varphi}{\partial t^2} = \Frac{\partial^2 \bar{\varphi}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{\varphi}}{\partial t}A \bar{X} + \dot{\theta}^2\mbox{div}((\nabla\bar{\varphi}A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta}\nabla\bar{\varphi}A \bar{X} + \ddot{z}(t).

Note that the term :math:`\dot{\theta} A \bar{X} = \left(\hspace*{-0.5em}\begin{array}{c} \dot{\theta}\bar{X}_2 \\ -\dot{\theta}\bar{X}_1 \\ 0 \end{array}\hspace*{-0.5em}\right)` is the rigid motion velocity vector. Now, If :math:`\Theta(t,X)` is a quantity attached to the material points (for instance the temperature), then, with :math:`\bar{\Theta}(t,\bar{X}) = \Theta(t,X)` , one simply has

.. math::

  \Frac{\partial \Theta}{\partial t} = \Frac{\partial \bar{\Theta}}{\partial t} + \dot{\theta} \nabla \bar{\Theta} A \bar{X}

This should not be forgotten that a correction has to be provided for each evolving variable for which the time derivative intervene in the considered model (think for instance to platic flow for plasticity). So that certain model bricks canot be used directly (plastic bricks for instance).

|gf| bricks for structural mecanics are mainly considering the displacement as the amin unknown. The expression for the displacement is the following:

.. math::
  \Frac{\partial u}{\partial t} = \Frac{\partial \bar{u}}{\partial t} + \dot{\theta} (I_d + \nabla \bar{u}) A \bar{X} + \dot{z}(t),

  \Frac{\partial^2 u}{\partial t^2} = \Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{z}(t).
  
  



Weak formulation of the transient terms
#######################################

Assuming :math:`\rho^0` the density in the reference configuration having a rotation symmetry, the term corresponding to acceleration in the weak formulation reads (with :math:`v(X) = \bar{v}(\bar{X})` a test function):

.. math::
   \int_{\Omega^0} \rho^0 \Frac{\partial^2 u}{\partial t^2}\cdot vdX =

   \int_{\bar{\Omega}^0} \rho^0 \left[\Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{z}(t) \right] \cdot \bar{v} d\bar{X}.
   
The third term in the right hand side can be integrated by part as follows:

.. math::
   \begin{array}{rcl}
   \ds \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} &=& - \ds \int_{\bar{\Omega}^0} (\dot{\theta}^2 (I_d + \nabla\bar{u})A \bar{X})) \cdot (\nabla (\rho^0 \bar{v}) A \bar{X}) d\bar{X} \\
   &&\ds + \int_{\partial \bar{\Omega}^0} \rho^0 \dot{\theta}^2 (((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \bar{N} \cdot \bar{v} d\bar{\Gamma}.
  \end{array}

Since :math:`\bar{N}` the outward unit normal vector on :math:`\partial \bar{\Omega}^0` is orthogonal to :math:`A \bar{X}` the boundary term is zero and :math:`\nabla (\rho^0 \bar{v}) = \bar{v} \otimes \nabla \rho^0   + \rho^0 \nabla \bar{v}` and since :math:`\nabla \rho^0.(A\bar{X}) = 0` because of the assumption on :math:`\rho^0` to have a rotation symmetry, we have

.. math::
   \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} = - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2 (A^2 \bar{X})\cdot \bar{v} d\bar{X}.

Thus globally

.. math::
   \begin{array}{rcl}
   \ds \int_{\Omega^0} \rho^0 \Frac{\partial^2 u}{\partial t^2}\cdot vdX &=&
   \ds \int_{\bar{\Omega}^0} \rho^0 \left[\Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} + \ddot{\theta} \nabla\bar{u} A \bar{X}   \right] \cdot \bar{v} d\bar{X}\\
   &&\ds - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 (\dot{\theta}^2 A^2 \bar{X} + \ddot{\theta} A\bar{X} + \ddot{z}(t))\cdot \bar{v} d\bar{X}.
   \end{array}

Note that two terms can deteriorate the coercivity of the problem and thus its well posedness and the stability of time integration schemes: the second one (convection term) and the fifth one. This may oblige to use additional stabilization techniques for large values of the angular velocity :math:`\dot{\theta}`.


The available bricks ...
++++++++++++++++++++++++

To be adapted  ::

  ind = getfem::brick_name(parmeters);

where ``parameters`` are the parameters ... 

=======
.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-ALE_rotating:

ALE terms for rotating objects
------------------------------

This section present a set of bricks facilitating the use of an ALE formulation for rotating bodies having a rotational symmetry (typically a train wheel).

Theoretical background
++++++++++++++++++++++

This strategy consists in adopting an intermediary description between an Eulerian and a Lagrangian ones for a rotating body having a rotational symmetry. This intermediary description consist in a rotating axes with respect to the reference configuration. See for instance [Dr-La-Ek2014]_ and [Nackenhorst2004]_.

It is supposed that the considered body is submitted approximately to a rigid body motion

.. math::
  \tau(X) = R(t)X + z(t)

and may have additonal deformation (exptected smaller) with respect to this rigid motion, where :math:`R(t)` is a rotation matrix

.. math::
  R(t) = \left(\begin{array}{ccc}
  \cos(\theta(t)) & \sin(\theta(t)) & 0 \\
  -\sin(\theta(t)) & \cos(\theta(t)) & 0 \\
  0 & 0 & 1
  \end{array} \right),

and :math:`z(t)` is a translation. This illustrated in the following figure:

.. _ud-fig-rotating_cylinder:

.. figure:: images/ALE_rotating_body.png
    :height: 200pt
    :align: center
    :alt: alternate text
    :figclass: align-center

Note that the description is given for a three-dimensional body. For two-dimensional bodies, the third axes is neglected so that :math:`R(t)` is a :math:`2\times 2` rotation matrix. Let us denote :math:`r(t)` the rotation:

.. math::
   r(t,X) = R(t)X, ~~~~~~~~~ \mbox{ and }
   A = \left(\begin{array}{ccc}
   0 & -1 & 0 \\
   1 & 0 & 0 \\
   0 & 0 & 0
   \end{array} \right).

We have then

.. math::
 \dot{r}(t,X) = \dot{\theta}AR(t)X

If :math:`\varphi(t, X)` is the deformation of the body which maps the reference configuration :math:`\Omega^0` to the deformed configuration :math:`\Omega_t` at time :math:`t`, the ALE description consists in the decomposition of the deformation of the cylinder in 

.. math::
   \varphi(t, X) = (\tau(t) \circ \bar{\varphi}(t) \circ r(t))(X) = \bar{\varphi}(t, r(t, X)) + z(t)

With :math:`\bar{X} = R(t)X` the new considered deformation is

.. math::
  \bar{\varphi}(t,\bar{X}) = \varphi(X) - z(t)


Thanks to the rotation symmetry of the reference configuration :math:`\Omega^0:`, we note that :math:`\bar{\Omega}^0 = r(t, \Omega^0)` is independant of :math:`t` and will serve as the new reference configuration. This is illustrated in the following figure: 

.. _ud-fig-rotating_cylinder_conf:

.. figure:: images/ALE_rotating_conf.png
    :height: 200pt
    :align: center
    :alt: alternate text
    :figclass: align-center

The denomination ALE of the method is justified by the fact that :math:`\bar{\Omega}^0` is an intermediate configuration which is of Euler type for the rigid motion and a Lagrangian one for the additional deformation of the solid. If we denote

.. math::
  \bar{u}(t,\bar{X}) = \bar{\varphi}(t, \bar{X}) - \bar{X}

the displacement with respect to this intermediate configuration, the advantage is that if this additional displacement with respect to the rigid body motion is small, it is possible to use a small deformation model (for instance linearized elasticity).

Due to the objectivity properties of standard consistutive laws, the expression od these laws in the intermediate configuration is most of the time identical to the expression in a standard reference configuration except for the expression of the time derivative which are modified because the change of coordinate is  nonconstant in time :

.. math::

  \Frac{\partial \varphi}{\partial t} = \Frac{\partial \bar{\varphi}}{\partial t} + \dot{\theta} \nabla \bar{\varphi} A \bar{X} + \dot{z}(t),

  \Frac{\partial^2 \varphi}{\partial t^2} = \Frac{\partial^2 \bar{\varphi}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{\varphi}}{\partial t}A \bar{X} + \dot{\theta}^2\mbox{div}((\nabla\bar{\varphi}A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta}\nabla\bar{\varphi}A \bar{X} + \ddot{z}(t).

Note that the term :math:`\dot{\theta} A \bar{X} = \left(\begin{array}{c} -\dot{\theta}\bar{X}_2 \\ \dot{\theta}\bar{X}_1 \\ 0 \end{array}\right)` is the rigid motion velocity vector. Now, If :math:`\Theta(t,X)` is a quantity attached to the material points (for instance the temperature), then, with :math:`\bar{\Theta}(t,\bar{X}) = \Theta(t,X)` , one simply has

.. math::

  \Frac{\partial \Theta}{\partial t} = \Frac{\partial \bar{\Theta}}{\partial t} + \dot{\theta} \nabla \bar{\Theta} A \bar{X}

This should not be forgotten that a correction has to be provided for each evolving variable for which the time derivative intervene in the considered model (think for instance to platic flow for plasticity). So that certain model bricks canot be used directly (plastic bricks for instance).

|gf| bricks for structural mecanics are mainly considering the displacement as the amin unknown. The expression for the displacement is the following:

.. math::
  \Frac{\partial u}{\partial t} = \Frac{\partial \bar{u}}{\partial t} + \dot{\theta} (I_d + \nabla \bar{u}) A \bar{X} + \dot{z}(t),

  \Frac{\partial^2 u}{\partial t^2} = \Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{z}(t).
  
  



Weak formulation of the transient terms
#######################################

Asuming :math:`\rho^0` the density in the reference configuration having a rotation symmetry, the term corresponding to acceleration in the weak formulation reads (with :math:`v(X) = \bar{v}(\bar{X})` a test function):

.. math::
   \int_{\Omega^0} \rho^0 \Frac{\partial^2 u}{\partial t^2}\cdot vdX =

   \int_{\bar{\Omega}^0} \rho^0 \left[\Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} +  \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) )  + \ddot{\theta} (I_d + \nabla\bar{u}) A \bar{X}  + \ddot{z}(t) \right] \cdot \bar{v} d\bar{X}.
   
The third term in the right hand side can be integrated by part as follows:

.. math::
   \begin{array}{rcl}
   \ds \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} &=& - \ds \int_{\bar{\Omega}^0} (\dot{\theta}^2 (I_d + \nabla\bar{u})A \bar{X})) \cdot (\nabla (\rho^0 \bar{v}) A \bar{X}) d\bar{X} \\
   &&\ds + \int_{\partial \bar{\Omega}^0} \rho^0 \dot{\theta}^2 (((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \bar{N} \cdot \bar{v} d\bar{\Gamma}.
  \end{array}

Since :math:`\bar{N}` the outward unit normal vector on :math:`\partial \bar{\Omega}^0` is orthogonal to :math:`A \bar{X}` the boundary term is zero and :math:`\nabla (\rho^0 \bar{v}) = \bar{v} \otimes \nabla \rho^0   + \rho^0 \nabla \bar{v}` and since :math:`\nabla \rho^0.(A\bar{X}) = 0` because of the asumption on :math:`\rho^0` to have a rotation symmetry, we have

.. math::
   \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2\mbox{div}(((I_d + \nabla\bar{u})A \bar{X}) \otimes (A \bar{X}) ) \cdot \bar{v} d\bar{X} = - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2 (A^2 \bar{X})\cdot \bar{v} d\bar{X}.

Thus globally

.. math::
   \begin{array}{rcl}
   \ds \int_{\Omega^0} \rho^0 \Frac{\partial^2 u}{\partial t^2}\cdot vdX &=&
   \ds \int_{\bar{\Omega}^0} \rho^0 \left[\Frac{\partial^2 \bar{u}}{\partial t^2} + 2\dot{\theta} \nabla\Frac{\partial \bar{u}}{\partial t}A \bar{X} + \ddot{\theta} \nabla\bar{u} A \bar{X}   \right] \cdot \bar{v} d\bar{X}\\
   &&\ds - \int_{\bar{\Omega}^0} \rho^0 \dot{\theta}^2(\nabla\bar{u}A \bar{X}) \cdot (\nabla \bar{v} A \bar{X}) d\bar{X} - \int_{\bar{\Omega}^0} \rho^0 (\dot{\theta}^2 A^2 \bar{X} + \ddot{\theta} A\bar{X} + \ddot{z}(t))\cdot \bar{v} d\bar{X}.
   \end{array}

Note that two terms can deteriorate the coercivity of the problem and thus its well posedness and the stability of time integration schemes: the second one (convection term) and the fifth one. This may oblige to use additional stabilization techniques for large values of the angular velocity :math:`\dot{\theta}`.


The available bricks ...
++++++++++++++++++++++++

To be adapted  ::

  ind = getfem::brick_name(parmeters);

where ``parameters`` are the parameters ... 

