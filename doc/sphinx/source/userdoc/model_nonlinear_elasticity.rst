.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-nonlinear-elasticity:

Nonlinear Elasticity brick
--------------------------

This brick implements some classical hyperelastic constitutive law for large deformation elasticity.

Some recalls on nonlinear elasticity
++++++++++++++++++++++++++++++++++++

Let :math:`\Omega` be the reference configuration and :math:`\Omega_t` the deformed configuration of an elastic media. Then for :math:`x \in \Omega` we will denote by :math:`\Phi(x) = u(x) + x` the deformation. the vector field :math:`u` is the displacement with respect to the initial position.

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

Both tensors :math:`E` and :math:`C` are used to describe nonlinear elasticity constitutive laws.

Main invariants and derivatives
###############################

The description of nonlinear elasticity constitutive laws often requires the principal invariants of the deformation tensors:

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

This property enables us to write the constitutive laws as a function of the Cauchy-Green tensor invariants, espacially for the case of the generalized Blatz-Ko strain energy.


Potential elastic energy and its derivative
###########################################

The stress in the reference configuration can be describe by the second Piola-Kirchhoff stress tensor :math:`{\hat{\hat{\sigma}}} = \nabla\Phi^{-1}\sigma\nabla\Phi^{-t}~\det \nabla\Phi` where :math:`\sigma` is the Cauchy stress tensor in the deformed configuration :math:`\Omega_t`. An hyper-elastic constitutive law is given by

.. math::

  {\hat{\hat{\sigma}}} &= \frac{\partial}{\partial E} {W}(E) = 2\frac{\partial}{\partial C} {W}(C)

where :math:`{W}` is the density of strain energy of the material. The total strain energy is given by

.. math::

  \mathcal{I}(u) = \int_\Omega W( E(u)) dx

and the derivative of the energy in a direction :math:`v` can be writen

.. math::

  D\mathcal{I}(u;v) = \int_{\Omega} \frac{\partial W}{\partial E}( E(u)):(I+{\nabla u^T}){\nabla v}  dx

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

  \int_\Omega(I + {\nabla u}){\hat{\hat{\sigma}}} : {\nabla v}  dx = l(v).


Tangent matrix
##############

The displacement :math:`u` is fixed. In order to obtain the tangent matrix, one subsitutes :math:`u` with :math:`u+h`

.. math::

  \int_\Omega(I + {\nabla u} + {\nabla h}){\hat{\hat{\sigma}}}( E(u)+ E(h) + \frac{1}{2}({\nabla h^T}{\nabla u}+{\nabla u^T}{\nabla h})) : {\nabla v}  dx = l(v)

and considers the linear part w.r.t. :math:`h`, which is

.. math::

  \int_\Omega{\nabla h}~{\hat{\hat{\sigma}}}( E(u)) : {\nabla v}  dx +\\
  \int_\Omega \frac{\partial^2 W}{\partial E^2}\left(\frac{{\nabla h}+{\nabla h^T}+{\nabla h^T}{\nabla u}+{\nabla u^T}{\nabla h}}{2}\right) : (I+{\nabla u}^T){\nabla v}  dx


which is symmetric w.r.t. :math:`v` and :math:`h`. It can be rewritten as

.. math::

  \int_\Omega {\nabla h}~{\hat{\hat{\sigma}}}( E(u)) : {\nabla v}  + \mathcal{A}((I+{\nabla u^T}){\nabla h}):(I+{\nabla u}^T){\nabla v}~ dx

where :math:`\mathcal{A}` is the symmetric :math:`3\times3\times3\times3` tensor given by :math:`\mathcal{A}_{ijkl} = ((\frac{\partial^2 W}{\partial E^2})_{ijkl} + (\frac{\partial^2 W}{\partial E^2})_{jikl})/2`.

Some classical constitutive laws
################################


``Linearized: Saint-Venant Kirchhoff  law (small deformations)``


.. math::

  {W} &= \frac{\lambda}{2}i_1( E)^2 + \mu i_1( E^2)\\
  {\hat{\hat{\sigma}}}   &= \lambda i_1( E)I + 2\mu E\\
  \mathcal{A} &= \lambda i_1(H)I + \mu (H + H^T)

``Two parameters Mooney-Rivlin law``

Incompressible material.

.. math::

  {W} = c_1(j_1( C) - 3) + c_2(j_2( C)-3)
  \intertext{with the additional constraint:}
  i_3( C) = 1

where :math:`c_1` and :math:`c_2` are given coefficients and

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

  {\hat{\hat{\sigma}}}   &= 2c_1 \frac{\partial j_1}{\partial C}(C) + 2c_2 \frac{\partial j_2}{\partial C}(C)  \\
  \mathcal{B} &= 4 c_1 \frac{\partial^2 j_1}{\partial C^2}(C) + 4c_2 \frac{\partial^2 j_2}{\partial C^2}(C) \\
  \mathcal{A}_{ijkl} &= (\mathcal{B}_{ijkl} + \mathcal{B}_{jikl})/2


The incompressibility constraint :math:`i_3( C) = 1` is handled with a Lagrange multiplier :math:`p` (the pression) 

constraint: :math:`\sigma = -pI \Rightarrow {\hat{\hat{\sigma}}} = -p\nabla\Phi\nabla\Phi^{-T}\det\nabla\Phi`

.. math::

  1 - i_3(\nabla\Phi) &= 0 \\
  -\int_{\Omega_0} (\det\nabla\Phi  -1) q  dx &= 0 ~~~ \forall q


.. math::

  B &= -\int_{\Omega_0} p(\nabla\Phi)^{-T} \det \nabla\Phi : \nabla v  dx \\
  K &= \int_{\Omega_0} \left( p(\nabla\Phi)^{-T}(\nabla h)^{T}(\nabla\Phi)^{-T}\det\nabla\Phi : \nabla v  dx - 
  p(\nabla\Phi)^{-T}(\det \nabla\Phi(\nabla\Phi)^{-T}:\nabla h) : \nabla v \right)  dx\\
  &= \int_{\Omega_0} p(\nabla h^T\nabla\Phi^{-T}):(\nabla\Phi^{-1}\nabla v)\det\nabla\Phi dx - \int_{\Omega_0} p(\nabla\Phi^{-T}:\nabla h)(\nabla\Phi^{-T}:\nabla v)\det\nabla\Phi dx


``Ciarlet-Geymonat law``

.. math::

  {W} &= \gamma_1i_1( E) + \frac{\lambda}{2}i_2( E) + 8ci_3( E) - \frac{\gamma_1}{2} \log \det  C


``Generalized Blatz-Ko law``

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

previous models are valid in volumic domains. Corresponding plane strain 2D models can be obtained by restricting the stress tensor and the fourth order tensor :math:`\mathcal{A}` to their plane components.  



Add an nonlinear elasticity brick to a model
++++++++++++++++++++++++++++++++++++++++++++

This brick represents a large strain elasticity problem. It is defined in the files :file:`getfem/getfem_nonlinear_elasticity.h` and :file:`getfem/getfem_nonlinear_elasticity.cc`. The function adding this brick to a model is ::

  ind = getfem::add_nonlinear_elasticity_brick
    (md, mim, varname, AHL, dataname, region = -1);

where ``AHL`` is an object of type ``getfem::abstract_hyperelastic_law`` which represents the considered hyperelastic law. It has to be chosen between: ::

  getfem::SaintVenant_Kirchhoff_hyperelastic_law AHL;
  getfem::Ciarlet_Geymonat_hyperelastic_law AHL;
  getfem::Mooney_Rivlin_hyperelastic_law AHL;
  getfem::plane_strain_hyperelastic_law AHL(pAHL);
  getfem::generalized_Blatz_Ko_hyperelastic_law AHL;

The Saint-Venant Kirchhoff law is a linearized law defined with the two Lame coefficients, Ciarlet Geymonat law is defined with the two Lame coefficients and an additional coefficient and the Mooney-Rivlin law is defined with two coefficients and is to be used with the large strain incompressibility condition. The plane strain hyperelastic law take a pointer on an hyperelastic law as a parameter and performs a 2D plane strain approximation.

``md`` is the model variable, ``mim`` the integration method, ``varname`` the string being the name of the variable on which the term is added, ``dataname`` the string being the name of the data in the model representing the coefficients of the law (can be constant or decribe on a finite element method) and ``region`` is the region on which the term is considered (by default, all the mesh). 


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

