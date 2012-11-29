.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks, Nitsche's method

.. _ud-model-Nitsche:


Nitsche's method for dirichlet and contact boundary conditions
--------------------------------------------------------------

Nitsche's method is a very attractive methods since it allows to take into
account Dirichlet type boundary conditions or contact with friction boundary conditions in a weak way without the use of Lagrange multipliers.  |gf| provides a generic implementation of Nitche's method. The advantage of Nitsche's method, which is to tranform a Dirichlet boundary condition into weak terms similarly as a Neumann boundary condition, is paid by the fact that the implementation is equation dependent. This method needs the use of an approximation of the corresponding Neumann term. Thus, in order to add a boundary condition with Nitsche's method on a variable of a model, the corresponding brick has to have access to an approximation of the Neumann term of all the partial differential terms applied to this variables. In the following, considering a variable :math:`u`, we will denote by

.. math::
  G(u,p)

the sum of all the Neumann terms on its variables. The additional parameter :math:`p` is introduced here because the Neumann term of a variable :math:`u` may depend on some other variables of the model (this is the case for instance for mixed formulations of incompressible elasticity). This additional parameter is here to describe what happens when the Neumann term depend on a other variable of the model. Of course, in complex situations, the Neumann term may depend on several variables. This is taken into account in the bricks implementing Nitsche's method. The Neumann terms depend also frequently on some parameters of the model (elasticity coefficients ...) but this is assumed to be contained in the the expression of :math:`G(u,p)`.

For instance, if there is a Laplace term (:math:`\Delta u`) is applied on the variable :math:`u`, the Neumann term will be :math:`G(u) = \Frac{\partial u}{\partial n}` where :math:`n` is the outward unit normal on the considered boundary. If :math:`u` represents the displacement of a deformable body, the Neumann term will be :math:`G(u) = \sigma(u)n`, where :math:`\sigma(u)` is the stress tensor depending on the consitutive law. Of course, in that case :math:`G(u)` depends on some body parameters. If additionally a mixed incompressibility brick is added with a variable :math:`p` denoting the pressure, the Neumann term on :math:`u` will depend on :math:`p` in the following way: :math:`G(u,p) = \sigma(u)n - pn`

In order to propose a generic implementation in which the brick proposing Nitsche's method are not dependent on the partial differential terms applied to the concerned variables, each brick adding a partial differential term is asked to give the expression of the corresponding Neumann term. Of course, it makes the building of a brick a little bit more complicated. So, this mechanism is not mandatory to build a new brick, but of course, it is mandatory if a brick implementing Nitsche's method is to be applied. An internal mechanism of |gf| controls this.
 

IMPORTANT: Contrarily to other bricks, the order in which the bricks implementing Nitche's method is important. Nitche's bricks have to be added after all the brick having a Neumann term (i.e. implementing partial differential terms) on the corresponding variable. Here again, an internal mechanism controls this.


Neumann term declaration for a brick
++++++++++++++++++++++++++++++++++++

In orer to write the tangent terms for Nitche's method we need not only the expression of the Neumann terms :math:`G(u,p)` but also the derivative with respect to :math:`u` and  :math:`p` and also ,in the nonlinear cases, the corresponding second derivatives.

In order to declare a Neumann term, a brick has to derive the object `Neumann_elem_term` (defined in `getfem_model.h`) and overload the virtual method `compute_Neumann_term`. The last parameter of the method `compute_Neumann_term` is the local variable number which is 0 for :math:`u`, 1 for :math:`p` if :math:`p` has been declared as a supplementary variable necessary to build the Neumann term, and so on. The supplementary variables have to be stored in the structure in the vector `auxilliary_variables' of `Neumann_elem_term` object. They also have to be declared when the brick is added with the method `add_auxilliary_variables_of_Neumann_terms` of the `model` object.


The first parameter of the method `compute_Neumann_term` is an integer denoting what should be provided as output tensor. For instance if the last parameter is equal to 0, the parameter is equal to:

   - 1 for :math:`output_i = (G(u,p))_i`
   - 2 for :math:`output_{ij} = (D_uG(u,p)[\varphi_i])_j`
   - 3 for :math:`output_{ijk} = D^2_{uu} (G(u,p)[\varphi_i, \varphi_j])_k`

where :math:`\varphi_i` is the finite element shape function for the variable :math:`u`. Now, if the last parameter is equal to 1 and if this corresponds to the variable :math:`p`, the output should be

   - 1 for :math:`output_i = (G(u,p))_i`
   - 2 for :math:`output_{ij} = (D_pG(u,p)[\psi_i])_j`
   - 3 for :math:`output_{ijk} = D^2_{up} (G(u,p)[\varphi_i, \psi_j])_k`

where :math:`\psi_i` is the finite element shape function for the variable :math:`p`.

No assistance is provided for the parameters which can intervene in :math:`G(u,p)`. Generally the structure representing the Neumann term has to store  what is necessary to compute the parameters.

Exemples of Neumann terms can be found in `getfem_models.cc` for the generic elliptic brick, the linearized eleasticity brick and the linear incompressibility brick.

Once the Neumann term is built, it should be added by the method `add_Neumann_term` of the `model` object when the assembly is called.

   

Generic Nitsche's method for a Dirichlet condition 
++++++++++++++++++++++++++++++++++++++++++++++++++

Assume that the variable :math:`u` is considered and that on wants to prescribe the condition

.. math::
  Hu = g

on a part :math:`\Gamma_D`  of the boundary of the considered domain. Here :math:`H` is considered equal to one in the scalar case or can be either the identity matrix in the vectorial case either a singular matrix having only 1 or 0 as eigenvalues. This allow here to prescribe only the normal or tangent component of :math:`u`. For instance if one wants to prescribe only the normal component, :math:`H` will be chosen to be equal to :math:`nn^T` where :math:`n` is the outward unit normal on :math:`\Gamma_D`.

Nitsche's method to prescribe this dirichlet condition consists in adding to the weak formulation of the problem the following term

.. math::
  \int_{\Gamma_D} \Frac{1}{\gamma}(Hu-g-\gamma HG(u,p)).(Hv) - \theta(Hu-g).(HD_uG(u,p)[v])d\Gamma,

where :math:`\gamma` and :math:`\theta` are two parameters of Nitsche's method and :math:`v` is the test function corresponding to :math:`u`. The parameter :math:`\theta` can be chosen positive or negative. :math:`\theta = 1` corresponds to the more standard method which leads to a symmetric tangent term in standard situations, :math:`\theta = 0` corresponds to a non-symmetric method which has the advantage to have a reduced number of terms and especially not to need the second derivatives of :math:`G(u,p)` in the nonlinear case, and :math:`\theta = -1` is a kind of skew-symmetric method which ensure an inconditonal coercivity (which means independent of :math:`\gamma`) at least in standard situations.
The parameter :math:`\gamma` is a kind of penalization parameter (although the method is consistent) which is taken to be :math:`\gamma = \gamma_0 h_T` where :math:`\gamma_0` is taken uniform on the mesh and :math:`h_T` is the diameter of the element :math:`T`. Note that, in standard situations, except for :math:`\theta = -1` the parameter :math:`\gamma_0` has to be taken sufficiently small in order to ensure the convergence of Nitsche's method.

Now, let us derive the tangent term corresponding to Nitsche's method. We will still consider that the Neumann term depends both on the variable  :math:`u` and on an auxilliary variable  :math:`p`. Of course, in practical case, there could be more than one auxilliary variable or zero. The tangent term reads as

.. math::
  &\int_{\Gamma_D} \Frac{1}{\gamma}(H\delta_u-\gamma HD_uG(u,p)[\delta_u]).(Hv) - \theta(H\delta_u).(HD_uG(u,p)[v])d\Gamma \\
  &-\int_{\Gamma_D} \theta(Hu-g).(HD^2_{uu}G(u,p)[v,\delta_u])d\Gamma \\
  &-\int_{\Gamma_D} (HD_pG(u,p)[\delta_p]).(Hv) + \theta(Hu-g).(HD^2_{up}G(u,p)[v,\delta_p])d\Gamma

where :math:`\delta_u` and :math:`\delta_p` are the incremental variable correpsonding to :math:`u` and :math:`p`, respectively.


The bricks adding a Dirichlet condition with Nitsche's method to a model are the following::

  getfem::add_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &gamma0name, size_type region,
      scalar_type theta = scalar_type(1),
      const std::string &dataname = std::string());


This function adds a Dirichlet condition on the variable `varname` and the mesh
region `region`. This region should be a boundary. The Dirichlet
condition is prescribed with Nitsche's method. `dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem; scalar or vector valued, depending on the variable
on which the Dirichlet condition is prescribed. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionnaly coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionnaly coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. Returns the brick index in the model.
CAUTION: This brick has to be added in the model after all the bricks
corresponding to partial differential terms having a Neumann term.
Moreover, This brick can only be applied to bricks declaring their
Neumann terms.
::


  getfem::add_normal_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &gamma0name, size_type region,
      scalar_type theta = scalar_type(1),
      const std::string &dataname = std::string());


This function adds a Dirichlet condition to the normal component of the vector
(or tensor) valued variable `varname` and the mesh region `region`.
This region should be a boundary. The Dirichlet
condition is prescribed with Nitsche's method. `dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionnaly coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionnaly coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. Returns the brick index in the model.
CAUTION: This brick has to be added in the model after all the bricks
corresponding to partial differential terms having a Neumann term.
Moreover, This brick can only be applied to bricks declaring their
Neumann terms. 
(This brick is not fully tested)
::

  getfem::add_generalized_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &gamma0name, size_type region, scalar_type theta,
      const std::string &dataname, const std::string &Hname);



This function adds a Dirichlet condition on the variable `varname` and the mesh
region `region`.
This version is for vector field. It prescribes a condition
:math:`Hu = r` where :math:`H` is a matrix field. The region should be a
boundary. This region should be a boundary.  The Dirichlet
condition is prescribed with Nitsche's method.
CAUTION : the matrix H should have all eigenvalues equal to 1 or 0.
`dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionnaly coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionnaly coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. `Hname` is the data
corresponding to the matrix field `H`. It has to be a constant matrix
or described on a scalar fem. Returns the brick index in the model.
CAUTION: This brick has to be added in the model after all the bricks
corresponding to partial differential terms having a Neumann term.
Moreover, This brick can only be applied to bricks declaring their
Neumann terms.
(This brick is not fully tested)

.. _nitsche_contact_small_def_section:

Generic Nitsche's method for contact with friction condition 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We describe here the use of Nitsch's method to prescribe a contact with Coulomb friction condition in the small deformations framework. This corresponds to a weak integral contact condition which as some similarity with the ones which use Lagrange multipliers describe in the corresponding section, see :ref:`weak_integral_contact_section`

In order to simplify notations, let use denote by :math:`P_{n,\mathscr{F}}` the following map which corresponds to a couple of projections:

.. math::
	P_{n,\mathscr{F}}(x) = -(x.n)_- n + P_{B(0,\mathscr{F}(x.n)_-)}(x - (x.n)n)

This application make the projection of the normal part of :math:`x` on :math:`\Reel_-` and the tangential part on the ball of center :math:`0` and radius :math:`\mathscr{F}(x.n)_-`, where :math:`\mathscr{F}` is the friction coefficient.

Using this, and considering that the sliding velocity is approximated by :math:`\alpha(u_{_T} - w_{_T})` where the expression of :math:`\alpha` and :math:`w_{_T}` depend on the time integration scheme used (see :ref:`weak_integral_contact_section`), Nitsche's term for contact with friction reads as:

.. math::
	&-\int_{\Gamma_C} \theta \gamma G(u,p)\cdot D_u G(u,p)[v] d\Gamma \\
	&+\int_{\Gamma_C} \gamma P_{n,\mathscr{F}}(G(u,p) - \Frac{Au}{\gamma} + \Frac{gap}{\gamma}n + \Frac{\alpha w_{_T}}{\gamma})\cdot(\theta D_u G(u,p)[v] - \Frac{v}{\gamma}) d\Gamma.

where :math:`\Gamma_C` is the contact boundary, :math:`G(u,p)` is the Neumann term which represents here :math:`\sigma n` the stress at the contact boundary and :math:`A` is the :math:`d\times d` matrix

.. math::
	A = \alpha I_d + (1-\alpha)n n^T

The corresponding tangent terms can be written as follows denoting :math:`\zeta(u,p) = G(u,p) - \Frac{Au}{\gamma} + \Frac{gap}{\gamma}n + \Frac{\alpha w_{_T}}{\gamma}`:

.. math::
	&-\int_{\Gamma_C}\theta\gamma(D_uG(u,p)[\delta_u])\cdot D_u G(u,p)[v] d\Gamma \\
	&+\int_{\Gamma_C}\gamma(\nabla P_{n,\mathscr{F}}(\zeta(u,p)))(D_uG(u,p)[\delta_u] - \Frac{A\delta u}{\gamma})\cdot (\theta D_u G(u,p)[v] - \Frac{v}{\gamma}) d\Gamma \\
	&+\int_{\Gamma_C} \theta\gamma\left( P_{n,\mathscr{F}}(\zeta(u,p))-G(u,p)\right)\cdot D^2_{uu} G(u,p)[v,\delta_u] d\Gamma \\
	&+\int_{\Gamma_C}\gamma (\nabla P_{n,\mathscr{F}}(\zeta(u,p))+I_d)(D_pG(u,p)[\delta_p])\cdot (\theta D_u G(u,p)[v] - \Frac{v}{\gamma}) d\Gamma \\
	&+\int_{\Gamma_C} \theta\gamma\left( P_{n,\mathscr{F}}(\zeta(u,p)) - G(u,p)\right)\cdot D^2_{up} G(u,p)[v,\delta_p] d\Gamma,

still considering that the Neumann term depends both on the variable  :math:`u` and on an auxilliary variable :math:`p` and with

.. math::
	\nabla P_{n,\mathscr{F}}(x) = H(-x_n) n n^T + \left\{ \begin{array}{l} (I_d-nn^T) \mbox{ if } \|x_t\| \le \mathscr{F}(x_n)_- \\ \Frac{\mathscr{F}(x_n)_-}{\|x_t\|} (I_d - \Frac{x_tx_t^T}{\|x_t\|^2} - nn^T) - \Frac{\mathscr{F}H(-x_n)}{\|x_t\|} x_t n^T \mbox{ otherwise, } \end{array}\right.

where :math:`x_n = x.n`, :math:`x_t = x - x_n n` and :math:`H(\cdot)` is the Heaviside function :math:`H(x) = 0` for :math:`x < 0` and :math:`H(x) = 1` for :math:`x \ge 0` (for :math:`x \in \R^2`, the term :math:`I_d - \Frac{x_tx_t^T}{\|x_t\|^2} - nn^T` vanishes).
	

Note that for the variant with :math:`\theta=0` a majority of terms vanish.