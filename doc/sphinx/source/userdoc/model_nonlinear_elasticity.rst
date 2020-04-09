.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlight:: none

.. index:: models, model bricks

.. _ud-model-nonlinear-elasticity:

Finite strain Elasticity bricks
-------------------------------

This brick implements some classical hyperelastic constitutive law for large deformation elasticity.

Some recalls on finite strain elasticity
++++++++++++++++++++++++++++++++++++++++

Let :math:`\Omega` be the reference configuration and :math:`\Omega_t` the deformed configuration of an elastic media. Then for :math:`X \in \Omega` we will denote by :math:`\Phi(x) = u(X) + X` the deformation. the vector field :math:`u` is the displacement with respect to the initial position.

The Cauchy-Green tensor is defined by

.. math::

  C = \nabla\Phi^T\nabla\Phi

The deformation tensor (Green-Lagrange)

.. math::

  E = \frac{1}{2}\left(\nabla\Phi^T\nabla\Phi - I)\right)
    = \frac{1}{2}\left({\nabla u^T}{\nabla u} + {\nabla u^T} + {\nabla u}\right)


(In the case of linear elasticity, :math:`{\nabla u^T}{\nabla u}` is neglected).

One has

.. math::

  C = \nabla\Phi^T\nabla\Phi = 2 E + I.

Both tensors :math:`E` and :math:`C` are used to describe finite strain elasticity constitutive laws.

Main invariants and derivatives
###############################

The description of finite strain elasticity constitutive laws often requires the principal invariants of the deformation tensors:

:math:`i_1,i_2,i_3` are the invariants of orders :math:`1,2` and :math:`3`:

.. math::

  i_1( E) = \mbox{tr } E \hspace{5cm} &i_1( C) = 2\mbox{tr } E + 3\\
  i_2( E) = \frac{(\mbox{tr } E)^2 - \mbox{tr } E^2}{2}\quad\hspace{3cm}& i_2( C)=4i_2( E)+4i_1( E)+3\\
  i_3( E) = \det E \hspace{5cm} &i_3( C) = 8i_3( E) + 4i_2( E) + 2i_1( E) + 1

The derivatives of the invariants with respect to the tensor :math:`E` in the direction :math:`H` are:

.. math::

  &\frac{\partial i_1}{\partial E}(E;H) = I:H = \mbox{tr } H\\
  &\frac{\partial i_2}{\partial E}(E;H) = (i_1( E)I -  E^T):H = (\mbox{tr }  E)(\mbox{tr } H) -  E^T:H\\
  &\frac{\partial i_3}{\partial E}(E;H) = i_3( E)(E^{-T}):H  = (i_2( E)I - i_1( E) E +  E^2):H \mbox{ in 3D}.

We will write

.. math::

  &\frac{\partial i_1}{\partial E}(E) = I\\
  &\frac{\partial i_2}{\partial E}(E) = i_1( E)I -  E^T\\
  &\frac{\partial i_3}{\partial E}(E) = i_3( E)E^{-T}.

Let us also recall that

.. math::

  \frac{\partial (M^{-1})}{\partial M}(M;H) = -M^{-1}HM^{-1}


The second derivatives of the invariants are fourth order tensors defined by

.. math::

  &\frac{\partial^2 i_1}{\partial E^2}(E) = 0\\
  &\frac{\partial^2 i_2}{\partial E^2}(E)_{ijkl} = \delta_{ij}\delta_{kl} - \delta_{il}\delta_{jk} \\
  &\frac{\partial^2 i_3}{\partial E^2}(E)_{ijkl} = i_3(E) (E^{-1}_{ji}E^{-1}_{lk} - E^{-1}_{jk}E^{-1}_{li}).


The notation :math:`A:B` denotes the Frobenius product :math:`A:B = \displaystyle\sum_{ij}A_{ij}B_{ij}`. This product has the following properties:

.. math::

  A:B &= \mbox{tr }(A^TB) = \mbox{tr }(AB^T) = \mbox{tr }(BA^T) = \mbox{tr }(B^TA),\\
  A:BC &= B^TA:C,\\
  A:BC &= AC^T:B,\\
  \mbox{tr }(ABC) &= \mbox{tr }(B^TA^TC^T)


Note also that

.. math::

  \frac{\partial i_j}{\partial E}(C;H) = 2 \frac{\partial i_j}{\partial C}(C;H).

This property enables us to write the constitutive laws as a function of the Cauchy-Green tensor invariants, especially for the case of the generalized Blatz-Ko strain energy.


Potential elastic energy and its derivative
###########################################

The stress in the reference configuration can be describe by the second Piola-Kirchhoff stress tensor :math:`{\hat{\hat{\sigma}}} = \nabla\Phi^{-1}\sigma\nabla\Phi^{-t}~\det \nabla\Phi` where :math:`\sigma` is the Cauchy stress tensor in the deformed configuration :math:`\Omega_t`. An hyper-elastic constitutive law is given by

.. math::

  {\hat{\hat{\sigma}}} &= \frac{\partial}{\partial E} {W}(E) = 2\frac{\partial}{\partial C} {W}(C)

where :math:`{W}` is the density of strain energy of the material. The total strain energy is given by

.. math::

  \mathcal{I}(u) = \int_{\Omega} W( E(u)) dX

and the derivative of the energy in a direction :math:`v` can be writen

.. math::

  D\mathcal{I}(u;v) = \int_{\Omega} \frac{\partial W}{\partial E}( E(u)):(I+{\nabla u^T}){\nabla v}  dX

because in particular

.. math::

  D E(u;v) &= \frac{1}{2}({\nabla u^T}{\nabla v} + {\nabla v^T}{\nabla u} + {\nabla v^T} + {\nabla v})\\
  &= \frac{1}{2}({\nabla v^T}(I+{\nabla u}) + (I+{\nabla u^T}){\nabla v})

and :math:`A:B = A:(B+B^T)/2` when A is symmetric which is the case for :math:`{\hat{\hat{\sigma}}}`.

Another way is to consider the static equilibrium which can be written as follows in the reference configuration:

.. math::
  -\mbox{div } \left((I+{\nabla u}){\hat{\hat{\sigma}}}\right) = f.


Integrating by parts, one obtains:

.. math::

  \int_{\Omega}(I + {\nabla u}){\hat{\hat{\sigma}}} : {\nabla v}  dX = l(v).


Tangent matrix
##############

The displacement :math:`u` is fixed. In order to obtain the tangent matrix, one subsitutes :math:`u` with :math:`u+h`

.. math::

  \int_\Omega(I + {\nabla u} + {\nabla h}){\hat{\hat{\sigma}}}( E(u)+ E(h) + \frac{1}{2}({\nabla h^T}{\nabla u}+{\nabla u^T}{\nabla h})) : {\nabla v}  dX = l(v)

and considers the linear part w.r.t. :math:`h`, which is

.. math::

  \int_\Omega{\nabla h}~{\hat{\hat{\sigma}}}( E(u)) : {\nabla v}  dX +\\
  \int_\Omega \frac{\partial^2 W}{\partial E^2}\left(\frac{{\nabla h}+{\nabla h^T}+{\nabla h^T}{\nabla u}+{\nabla u^T}{\nabla h}}{2}\right) : (I+{\nabla u}^T){\nabla v}  dX


which is symmetric w.r.t. :math:`v` and :math:`h`. It can be rewritten as

.. math::

  \int_\Omega {\nabla h}~{\hat{\hat{\sigma}}}( E(u)) : {\nabla v}  + \mathcal{A}((I+{\nabla u^T}){\nabla h}):(I+{\nabla u}^T){\nabla v}~ dX

where :math:`\mathcal{A}` is the symmetric :math:`3\times3\times3\times3` tensor given by :math:`\mathcal{A}_{ijkl} = ((\frac{\partial^2 W}{\partial E^2})_{ijkl} + (\frac{\partial^2 W}{\partial E^2})_{ijlk})/2`.

Some classical constitutive laws
################################


``Linearized: Saint-Venant Kirchhoff  law (small deformations)``
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

.. math::

  {W} &= \frac{\lambda}{2}i_1( E)^2 + \mu i_1( E^2)\\
  {\hat{\hat{\sigma}}}   &= \lambda i_1( E)I + 2\mu E\\
  \mathcal{A} &= \lambda i_1(H)I + \mu (H + H^T)

``Three parameters Mooney-Rivlin law``
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Compressible material.

.. math::

  {W} = c_1(j_1( C) - 3) + c_2(j_2( C)-3) + d_1(i_3( C)^{1/2}-1)^2

where :math:`c_1`, :math:`c_2` and :math:`d_1` are given coefficients and

.. math::

  j_1(C) &= i_1(C) i_3(C)^{-1/3}\\
  j_2(C) &= i_2(C) i_3(C)^{-2/3}\\
  \frac{\partial j_1}{\partial C}(C) &= i_3(C)^{-1/3}\left(\frac{\partial i_1}{\partial C}(C) - \frac{i_1(C)}{3i_3(C)} \frac{\partial i_3}{\partial C}(C)\right)\\
  \frac{\partial j_2}{\partial C}(C) &= i_3(C)^{-2/3}\left(\frac{\partial i_2}{\partial C}(C) - \frac{2i_2(C)}{3i_3(C)} \frac{\partial i_3}{\partial C}(C)\right)\\
  \frac{\partial^2 j_1}{\partial C^2}(C) &= i_3(C)^{-1/3}\left(\frac{4i_1(C)}{9i_3(C)^2} \frac{\partial i_3}{\partial C}(C) \otimes \frac{\partial i_3}{\partial C}(C) - \frac{1}{3i_3(C)}\left(\frac{\partial i_3}{\partial C}(C) \otimes \frac{\partial i_1}{\partial C}(C)\right.\right. \\
  & ~~~~~~~~~~~~~~~~\left.\left. + \frac{\partial i_1}{\partial C}(C) \otimes \frac{\partial i_3}{\partial C}(C)\right) - \frac{i_1(C)}{3i_3(C)} \frac{\partial^2 i_3}{\partial C^2}(C)\right)\\
  \frac{\partial^2 j_2}{\partial C^2}(C) &= i_3(C)^{-2/3}\left(\frac{\partial^2 i_2}{\partial C^2}(C) + \frac{10i_2(C)}{9i_3(C)^2} \frac{\partial i_3}{\partial C}(C) \otimes \frac{\partial i_3}{\partial C}(C) \right. \\
  & ~~~~~~~~~~~~~~~~\left. - \frac{2}{3i_3(C)}(\frac{\partial i_3}{\partial C}(C) \otimes \frac{\partial i_2}{\partial C}(C) + \frac{\partial i_2}{\partial C}(C) \otimes \frac{\partial i_3}{\partial C}(C)) - \frac{2i_2(C)}{3i_3(C)} \frac{\partial^2 i_3}{\partial C^2}(C)\right)

and then

.. math::

  {\hat{\hat{\sigma}}}   &= 2c_1 \frac{\partial j_1}{\partial C}(C) + 2c_2 \frac{\partial j_2}{\partial C}(C)  + 2d_1\left(1-i_3(C)^{-1/2}\right)\frac{\partial i_3}{\partial C}(C) \\
  \mathcal{B} &= 4 c_1 \frac{\partial^2 j_1}{\partial C^2}(C) + 4c_2 \frac{\partial^2 j_2}{\partial C^2}(C) + 4d_1\left(\left(1-i_3(C)^{-1/2}\right)\frac{\partial^2 i_3}{\partial C^2}(C) + \frac{1}{2}i_3(C)^{-3/2} \frac{\partial i_3}{\partial C}(C) \otimes \frac{\partial i_3}{\partial C}(C)\right) \\
  \mathcal{A}_{ijkl} &= (\mathcal{B}_{ijkl} + \mathcal{B}_{jikl})/2

Incompressible material.

.. math::

  {d_1} = 0
  \intertext{with the additional constraint:}
  i_3( C) = 1

The incompressibility constraint :math:`i_3( C) = 1` is handled with a Lagrange multiplier :math:`p` (the pressure)

constraint: :math:`\sigma = -pI \rm I\hspace{-0.15em}Rightarrow {\hat{\hat{\sigma}}} = -p\nabla\Phi\nabla\Phi^{-T}\det\nabla\Phi`

.. math::

  1 - i_3(\nabla\Phi) &= 0 \\
  -\int_{\Omega_0} (\det\nabla\Phi  -1) q  dX &= 0 ~~~ \forall q


.. math::

  B &= -\int_{\Omega_0} p(\nabla\Phi)^{-T} \det \nabla\Phi : \nabla v  dX \\
  K &= \int_{\Omega_0} \left( p(\nabla\Phi)^{-T}(\nabla h)^{T}(\nabla\Phi)^{-T}\det\nabla\Phi : \nabla v  dX -
  p(\nabla\Phi)^{-T}(\det \nabla\Phi(\nabla\Phi)^{-T}:\nabla h) : \nabla v \right)  dX\\
  &= \int_{\Omega_0} p(\nabla h^T\nabla\Phi^{-T}):(\nabla\Phi^{-1}\nabla v)\det\nabla\Phi dX - \int_{\Omega_0} p(\nabla\Phi^{-T}:\nabla h)(\nabla\Phi^{-T}:\nabla v)\det\nabla\Phi dX


``Ciarlet-Geymonat law``
<<<<<<<<<<<<<<<<<<<<<<<<

.. math::

  {W} &= a\; i_1(C) + (\frac{\mu}{2} - a)i_2(C) + (\frac{\lambda}{4} - \frac{\mu}{2} + a)i_3(C) - (\frac{\mu}{2}+\frac{\lambda}{4})\log \det(C)

with  :math:`\lambda, \mu` the Lame coefficients and :math:`\max(0,\frac{\mu}{2}-\frac{\lambda}{4})<a<\frac{\mu}{2}` (see [ciarlet1988]_).


``Generalized Blatz-Ko law``
<<<<<<<<<<<<<<<<<<<<<<<<<<<<

.. math::

 {W} &= (ai_1(C) + bi_3(C)^{1/2} + c\frac{\i_2(C)}{\i_3(C)} + d)^n

Since :math:`\frac{\partial}{\partial C} {W}(C) = \displaystyle\sum_{j}\frac{\partial W}{\partial i_j(C)} \frac{\partial i_j(C)}{\partial C}`, and :math:`\frac{\partial^2}{\partial C^2} {W}(C) = \displaystyle\sum_{j} \displaystyle\sum_{k} \frac{\partial^2 W}{\partial i_j(C) \partial i_k(C)} \frac{\partial i_k(C)}{\partial C} \otimes \frac{\partial i_j(C)}{\partial C} + \displaystyle\sum_{j} \frac{\partial W}{\partial i_j(C)} \frac{\partial^2 i_j(C)}{\partial C^2}` we must compute the derivatives of the strain energy function with respect to the Cauchy-Green tensor invariants (we don't need to compute the invariants derivatives with respect to :math:`E` since :math:`\frac{\partial i_j}{\partial E}(C;H) = 2 \frac{\partial i_j}{\partial C}(C;H)`) :

.. math::
  \begin{array}{l}
  \frac{\partial W}{\partial i_1(C)} = naZ^{n-1}
  ~~~~\mbox{with } Z = (ai_1(C) + bi_3(C)^{1/2} + c\frac{\i_2(C)}{\i_3(C)} + d)\\
  \frac{\partial W}{\partial i_2(C)} = n\frac{c}{i_3(C)}Z^{n-1}\\
  \frac{\partial W}{\partial i_3(C)} = n(\frac{b}{2i_3(C)^{1/2}}-\frac{ci_2(C)}{i_3(C)^2})Z^{n-1}\\
  \frac{\partial W^2}{\partial^2 i_1(C)} = n(n-1)A^2Z^{n-2}\\
  \frac{\partial W^2}{\partial i_1(C) \partial i_2(C)} = n(n-1)A\frac{c}{i_3(C)}Z^{n-2}\\
  \frac{\partial W^2}{\partial i_1(C) \partial i_3(C)} = n(n-1)A(\frac{b}{2i_3(C)^{1/2}}-\frac{ci_2(C)}{i_3(C)^2})Z^{n-2}\\
  \frac{\partial W^2}{\partial^2 i_2(C)} = n(n-1)\frac{c^2}{i_3(C)^2}Z^{n-2}\\
  \frac{\partial W^2}{\partial i_2(C) \partial i_3(C)} = n(n-1)(\frac{b}{2i_3(C)^{1/2}}-\frac{ci_2(C)}{i_3(C)^2})Z^{n-2} - n\frac{c^2}{i_3(C)^2}Z^{n-1}\\
  \frac{\partial W^2}{\partial i_3(C)^2} = n(n-1)(\frac{b}{2i_3(C)^{1/2}}-\frac{ci_2(C)}{i_3(C)^2})^2Z^{n-2} + n(-\frac{b}{4i_3(C)^{3/2}}+2\frac{ci_2(C)}{i_3(C)^4})Z^{n-1}
  \end{array}

``Plane strain hyper-elasticity``
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

All previous models are valid in volumic domains. Corresponding plane strain 2D models can be obtained by restricting the stress tensor and the fourth order tensor :math:`\mathcal{A}` to their plane components.



Add an nonlinear elasticity brick to a model
++++++++++++++++++++++++++++++++++++++++++++

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

The Mooney-Rivlin law accepts two optional flags, the first one determines if the material will be compressible (:math:`d_1 \neq 0`) and the second one determines if the material is neo Hookean (:math:`c_2 = 0`). Depending on these flags one to three coefficients may be necessary. By default it is defined as incompressible and non neo Hookean, thus it needs two material coefficients (:math:`c_1`, :math:`c_2`). In this case, it is to be used with the large strain incompressibility condition.

The plane strain hyperelastic law takes a pointer on a hyperelastic law as a parameter and performs a 2D plane strain approximation.

``md`` is the model variable, ``mim`` the integration method, ``varname`` the string being the name of the variable on which the term is added, ``dataname`` the string being the name of the data in the model representing the coefficients of the law (can be constant or describe on a finite element method) and ``region`` is the region on which the term is considered (by default, all the mesh).


The program :file:`nonlinear_elastostatic.cc` in :file:`tests` directory and :file:`demo_nonlinear_elasticity.m` in :file:`interface/tests/matlab` directory are some examples of use of this brick with or without an incompressibility condition.


Note that the addition of a new hyperelastic constitutive law consists in furnishing the expression of the strain energy, the stress tensor and the derivative of the stress tensor. See the file  :file:`getfem/getfem_nonlinear_elasticity.cc` for more details. In particular, expression of the invariants and their derivatives are available.


A function which computes the Von Mises or Tresca stresses is also available: ::

  VM = compute_Von_Mises_or_Tresca
    (md, varname, AHL, dataname, mf_vm, VM, tresca)

It returns a vector of the degrees of freedom of the Von Mises or Tresca stress on the finite element method mf_vm. ``tresca`` is a boolean whose value should be ``true`` for Tresca stress and ``false`` for Von Mises stress.



Add a large strain incompressibility brick to a model
+++++++++++++++++++++++++++++++++++++++++++++++++++++


This brick adds an incompressibility condition in a large strain problem of type

.. math::

 \mbox{det}(I+\nabla u) = 1,

A Lagrange multiplier representing the pressure is introduced in a mixed formulation. The function adding this brick to a model is ::

  ind = add_nonlinear_incompressibility_brick
    (md, mim, varname, multname, region = -1)


where ``md`` is the model, ``mim`` the integration method, ``varname`` the variable of the model on which the incompressibility condition is added, ``multanme`` the multiplier variable corresponding to the pressure (be aware that at least a linear Ladyzhenskaja-Babuska-Brezzi inf-sup condition is satisfied between the f.e.m. of the variable and the one of the multiplier). ``region`` is an optional parameter correponding to the mesh region on which the term is considered (by default, all the mesh).


High-level generic assembly versions
++++++++++++++++++++++++++++++++++++

The generic weak form language (GWFL) gives access to the hyperelastic potential and constitutive laws implemented in |gf|. This allows to directly use them in the language, for instance using a generic assembly brick in a model or for interpolation of certain quantities (the stress for instance).

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

The potentials::

  Saint_Venant_Kirchhoff_potential(Grad_u, [lambda; mu])
  Plane_Strain_Saint_Venant_Kirchhoff_potential(Grad_u, [lambda; mu])
  Generalized_Blatz_Ko_potential(Grad_u, [a;b;c;d;n])
  Plane_Strain_Generalized_Blatz_Ko_potential(Grad_u, [a;b;c;d;n])
  Ciarlet_Geymonat_potential(Grad_u, [lambda;mu;a])
  Plane_Strain_Ciarlet_Geymonat_potential(Grad_u, [lambda;mu;a])
  Incompressible_Mooney_Rivlin_potential(Grad_u, [c1;c2])
  Plane_Strain_Incompressible_Mooney_Rivlin_potential(Grad_u, [c1;c2])
  Compressible_Mooney_Rivlin_potential(Grad_u, [c1;c2;d1])
  Plane_Strain_Compressible_Mooney_Rivlin_potential(Grad_u, [c1;c2;d1])
  Incompressible_Neo_Hookean_potential(Grad_u, [c1])
  Plane_Strain_Incompressible_Neo_Hookean_potential(Grad_u, [c1])
  Compressible_Neo_Hookean_potential(Grad_u, [c1;d1])
  Plane_Strain_Compressible_Neo_Hookean_potential(Grad_u, [c1;d1])
  Compressible_Neo_Hookean_Bonet_potential(Grad_u, [lambda;mu])
  Plane_Strain_Compressible_Neo_Hookean_Bonet_potential(Grad_u, [lambda;mu])
  Compressible_Neo_Hookean_Ciarlet_potential(Grad_u, [lambda;mu])
  Plane_Strain_Compressible_Neo_Hookean_Ciarlet_potential(Grad_u, [lambda;mu])


The second Piola-Kirchhoff stress tensors::

  Saint_Venant_Kirchhoff_PK2(Grad_u, [lambda; mu])
  Plane_Strain_Saint_Venant_Kirchhoff_PK2(Grad_u, [lambda; mu])
  Generalized_Blatz_Ko_PK2(Grad_u, [a;b;c;d;n])
  Plane_Strain_Generalized_Blatz_Ko_PK2(Grad_u, [a;b;c;d;n])
  Ciarlet_Geymonat_PK2(Grad_u, [lambda;mu;a])
  Plane_Strain_Ciarlet_Geymonat_PK2(Grad_u, [lambda;mu;a])
  Incompressible_Mooney_Rivlin_PK2(Grad_u, [c1;c2])
  Plane_Strain_Incompressible_Mooney_Rivlin_PK2(Grad_u, [c1;c2])
  Compressible_Mooney_Rivlin_PK2(Grad_u, [c1;c2;d1])
  Plane_Strain_Compressible_Mooney_Rivlin_PK2(Grad_u, [c1;c2;d1])
  Incompressible_Neo_Hookean_PK2(Grad_u, [c1])
  Plane_Strain_Incompressible_Neo_Hookean_PK2(Grad_u, [c1])
  Compressible_Neo_Hookean_PK2(Grad_u, [c1;d1])
  Plane_Strain_Compressible_Neo_Hookean_PK2(Grad_u, [c1;d1])
  Compressible_Neo_Hookean_Bonet_PK2(Grad_u, [lambda;mu])
  Plane_Strain_Compressible_Neo_Hookean_Bonet_PK2(Grad_u, [lambda;mu])
  Compressible_Neo_Hookean_Ciarlet_PK2(Grad_u, [lambda;mu])
  Plane_Strain_Compressible_Neo_Hookean_Ciarlet_PK2(Grad_u, [lambda;mu])


Note that the derivatives with respect to the material parameters have not been implemented apart for the Saint Venant Kirchhoff hyperelastic law. Therefore, it is not possible to make the parameter depend on other variables of a model (derivatives are not necessary complicated to implement but for the moment, only a wrapper with old implementations has been written).

Note that the coupling of models is to be done at the weak formulation level. In a general way, it is recommended not to use the potential to define a problem. Main couplings cannot be obtained at the potential level. Thus the use of potential should be restricted to the actual computation of the potential.

An example of use to add a Saint Venant-Kirchhoff hyperelastic term to a variable ``u`` in a model or a ga_workspace is given by the addition of the following assembly string::

  "((Id(meshdim)+Grad_u)*(Saint_Venant_Kirchhoff_PK2(Grad_u,[lambda;mu]))):Grad_Test_u"

Note that in that case, ``lambda`` and ``mu`` have to be declared data of the model/ga_workspace. It is of course possible to replace them by explicit constants or expressions depending on several data.

Concerning the incompressible Mooney-Rivlin law, it has to be completed by an incompressibility term. For instance by adding the following incompressibility brick::

 ind = add_finite_strain_incompressibility_brick(md, mim, varname, multname, region = -1);

This brick just adds the term ``p*(1-Det(Id(meshdim)+Grad_u))`` if ``p`` is the multiplier and ``u`` the variable which represents the displacement.

The addition of an hyperelastic term to a model can also be done thanks to the following function::

  ind = add_finite_strain_elasticity_brick(md, mim, lawname, varname, params,
                                           region = size_type(-1));

where ``md`` is the model, ``mim`` the integration method, ``varname`` the variable of the model representing the large strain displacement, ``lawname`` is the constitutive law name which could be ``Saint_Venant_Kirchhoff``, ``Generalized_Blatz_Ko``, ``Ciarlet_Geymonat``, ``Incompressible_Mooney_Rivlin``, ``Compressible_Mooney_Rivlin``, ``Incompressible_Neo_Hookean``, ``Compressible_Neo_Hookean``, ``Compressible_Neo_Hookean_Bonet`` or ``Compressible_Neo_Hookean_Ciarlet``. ``params`` is a string representing the parameters of the law defined as a small vector or a vector field.

The Von Mises stress can be interpolated with the following function::

  void compute_finite_strain_elasticity_Von_Mises(md, varname, lawname, params, mf_vm, VM,
                                                  rg=mesh_region::all_convexes());

where ``md`` is the model, ``varname`` the variable of the model representing the large strain displacement, ``lawname`` is the constitutive law name (see previou brick), ``params`` is a string representing the parameters of the law, ``mf_vm`` a (preferably discontinuous) Lagrange  finite element method on which the interpolation will be done and ``VM`` a vector of type ``model_real_plain_vector`` in which the interpolation will be stored.
