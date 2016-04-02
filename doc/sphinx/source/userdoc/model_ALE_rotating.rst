.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-ALE_rotating:

ALE terms for rotating objects
------------------------------

This section present a set of bricks facilitating the use of an ALE formulation for rotating bodies having a rotational symmetry (typically a train wheel).

Work in progress ...

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
   r(t,X) = R(t)X, ~~~~~~~~~ mbox{ and }
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
  
  



Main invariants and derivatives
###############################




the available bricks ...
++++++++++++++++++++++++

To be adapted ..

This brick represents a large strain elasticity problem. It is defined in the files :file:`getfem/getfem_nonlinear_elasticity.h` and :file:`getfem/getfem_nonlinear_elasticity.cc`. The function adding this brick to a model is ::

  ind = getfem::add_nonlinear_elasticity_brick
    (md, mim, varname, AHL, dataname, region = -1);

where ``AHL`` is an object of type ``getfem::abstract_hyperelastic_law`` which represents the considered hyperelastic law. It has to be chosen between: ::

  getfem::SaintVenant_Kirchhoff_hyperelastic_law AHL;
  getfem::Ciarlet_Geymonat_hyperelastic_law AHL;
  getfem::Mooney_Rivlin_hyperelastic_law AHL(compressible, neohookean);
  getfem::plane_strain_hyperelastic_law AHL(pAHL);
  getfem::generalized_Blatz_Ko_hyperelastic_law AHL;

The Saint-Venant Kirchhoff law is a linearized law defined with the two Lame coefficients, Ciarlet Geymonat law is defined with the two Lame coefficients and an additional coefficient (:math:`\lambda, \mu, a`).



Here is the list of nonlinear operators in the language which can be useful for nonlinear elasticity::

  Det(M)                                % determinant of the matrix M
  Trace(M)                              % trace of the matrix M
  Matrix_i2(M)                          % second invariant of M (in 3D): (sqr(Trace(m)) - Trace(m*m))/2
  Matrix_j1(M)                          % modified first invariant of M: Trace(m)pow(Det(m),-1/3).
  Matrix_j2(M)                          % modified second invariant of M: Matrix_I2(m)*pow(Det(m),-2/3).
  Right_Cauchy_Green(F)                 % F' * F
  Left_Cauchy_Green(F)                  % F * F'
  Green_Lagrangian(F)                   % (F'F - Id(meshdim))/2
  Cauchy_stress_from_PK2(sigma, Grad_u) % (Id+Grad_u)*sigma*(I+Grad_u')/det(I+Grad_u)

