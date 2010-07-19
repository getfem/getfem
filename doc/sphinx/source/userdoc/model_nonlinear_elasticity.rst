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


Let us also recall that

.. math::

  \frac{\partial (M^{-1})}{\partial M}(M;H) &= -M^{-1}HM^{-1}\\
  &= -H:(M^{-T}M^{-T}) \textrm{\`a confirmer..}

The notation :math:`A:B` denotes the Frobenius product :math:`A:B = \displaystyle\sum_{ij}A_{ij}B_{ij}`. This product has the following properties:

.. math::

  A:B &= \mbox{tr }(A^TB) = \mbox{tr }(AB^T) = \mbox{tr }(BA^T) = \mbox{tr }(B^TA),\\
  A:BC &= B^TA:C,\\
  A:BC &= AC^T:B,\\
  \mbox{tr }(ABC) &= \mbox{tr }(B^TA^TC^T)


Note also that

.. math::

  \frac{\partial i_j}{\partial E}(C;H) = 2 \frac{\partial i_j}{\partial C}(C;H).


Potential elastic energy and its derivative
###########################################

The stress in the reference configuration can be describe by the second Piola-Kirchhoff stress tensor :math:`{\hat{\hat{\sigma}}} = \nabla\Phi^{-1}\sigma\nabla\Phi^{-t}~\det \nabla\Phi` where :math:`\sigma` is the Cauchy stress tensor in the deformed configuration :math:`\Omega_t`. An hyper-elastic constitutive law is given by

.. math::

  {\hat{\hat{\sigma}}} &= -\frac{\partial}{\partial E} {W}(E)

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

where :math:`\mathcal{A}` is the symmetric :math:`3\times3\times3\times3` tensor given by :math:`\mathcal{A}_{ijkl} = ((\frac{\partial W}{\partial E})_{ijkl} + (\frac{\partial W}{\partial E})_{jikl})/2`.

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

  {W} &= c_1(i_1( C) - 3) + c_2(i_2( C)-2)\\
  &= 2c_1i_1( E) + 4c_2(i_2( E)+i_1( E))
  \intertext{with the additional constraint:}
  i_3( C) = 1

where :math:`c_1` and :math:`c_2` are given coefficients.

.. math::

  {\hat{\hat{\sigma}}} &= (2c_1 + 4c_2)I + 4c_2(i_1(L)I - L)\\
  \mathcal{A} &= 2c_2(2i_1(H)I - H - H^T)


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


Add an nonlinear elasticity brick to a model
++++++++++++++++++++++++++++++++++++++++++++

This brick represents a large strain elasticity problem. It is defined in the files :file:`getfem/getfem_nonlinear_elasticity.h` and :file:`getfem/getfem_nonlinear_elasticity.cc`. The function adding this brick to a model is ::

  ind = getfem::add_nonlinear_elasticity_brick
    (md, mim, varname, AHL, dataname, region = -1);

where ``AHL`` is an object of type ``getfem::abstract_hyperelastic_law`` which represents the considered hyperelastic law. It has to be chosen between: ::

  getfem::SaintVenant_Kirchhoff_hyperelastic_law AHL;
  getfem::Ciarlet_Geymonat_hyperelastic_law AHL;
  getfem::Mooney_Rivlin_hyperelastic_law AHL;

The Saint-Venant Kirchhoff law is a linearized law defined with the two Lame coefficients, Ciarlet Geymonat law is defined with the two Lame coefficients and an additional coefficient and the Mooney-Rivlin law is defined with two coefficients and is to be used with the large strain incompressibility condition.

``md`` is the model variable, ``mim`` the integration method, ``varname`` the string being the name of the variable on which the term is added, ``dataname`` the string being the name of the data in the model representing the coefficients of the law (can be constant or decribe on a finite element method) and ``region`` is the region on which the term is considered (by default, all the mesh). 


The program :file:`nonlinear_elastostatic.cc` in :file:`tests` directory and :file:`demo_nonlinear_elasticity.m` in :file:`interface/tests/matlab` directory are some examples of use of this brick with or without an incompressibility condition.


It can be noted that the add of a new hyperelastic constitutive law is rther easy. It is sufficient to add the expression of the strain energy, the stress tensor and the derivative of the stress tensor. See the file  :file:`getfem/getfem_nonlinear_elasticity.h` for more details.


It is also furnished a function which computes the Von Mises or Tresca stress: ::

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

