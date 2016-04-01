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
   :align: center
   :scale: 80

Note that the description is given for a three-dimensional body. For two-dimensional bodies, the third axes is neglected so that :math:`R(t)` is a :math:`2\times 2` rotation matrix.

Denoting :math:`r(t)` the rotation

.. math::
   r(t,X) = R(t)X, ~~~~~~~~~
   A(t) = \left(\begin{array}{ccc}
   0 & -1 & 0 \\
   1 & 0 & 0 \\
   0 & 0 & 0
   \end{array} \right), ~~~~~~~~~
   \mbox{ and } B(t) = A^2(t) = \left(\begin{array}{ccc}
   -1 & 0 & 0 \\
   0 & -1 & 0 \\
   0 & 0 & 0
   \end{array} \right),

such that

.. math::
 \dot{r}(t,X) = \dot{\theta}A(t)R(t)X



The ALE description consists in the decomposition of the motion of the cylinder :math:`\varphi(t, X)` in 

.. math::
   \varphi(t, X) = (\tau(t) \circ \bar{\varphi}(t) \circ r(t))(X) = \bar{\varphi}(t, r(t, X)) + z(t)

With :math:`\bar{X} = R(t)X` the new considered deformation is

.. math::
  \bar{\varphi}(\bar{X}) = \varphi(X) - z(t)


.. _ud-fig-rotating_cylinder_conf:

.. figure:: images/ALE_rotating_conf.png
   :align: center
   :scale: 80


+ analyse ... et graphique eventuel


Main invariants and derivatives
###############################




the available bricks ...
++++++++++++++++++++++++

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

