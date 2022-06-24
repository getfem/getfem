.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. index:: models, model bricks, Nitsche's method

.. _ud-model-Nitsche:


Nitsche's method for dirichlet and contact boundary conditions
--------------------------------------------------------------

|gf| provides a generic implementation of Nitche's method which allows to account for Dirichlet type or contact with friction boundary conditions in a weak sense without the use of Lagrange multipliers.
The method is very attractive because it transforms a Dirichlet boundary condition into a weak term similar to a Neumann boundary condition.
However, this advantage is at the cost that the implementation of Nitche's method is model dependent, since it requires an approximation of the corresponding Neumann term.
In order to add a boundary condition with Nitsche's method on a variable of a model, the corresponding brick needs to have access to an approximation of the Neumann term of all partial differential terms applied to this variable.
In the following, considering a variable :math:`u`, we will denote by

.. math::
  G

the sum of all Neumann terms on this variable.
Note that the Neumann term :math:`G` will often depend on the variable :math:`u` but it may also depend on other variables of the model.
This is the case for instance for mixed formulations of incompressible elasticity.
The Neumann terms depend also frequently on some parameters of the model (elasticity coefficients ...) but this is assumed to be contained in its expression.

For instance, if there is a Laplace term (:math:`\Delta u`), applied on the variable :math:`u`, the Neumann term will be :math:`G = \dfrac{\partial u}{\partial n}` where :math:`n` is the outward unit normal on the considered boundary.
If :math:`u` represents the displacements of a deformable body, the Neumann term will be :math:`G = \sigma(u)n`, where :math:`\sigma(u)` is the stress tensor depending on the constitutive law.
Of course, in that case :math:`G` also depends on some material parameters.
If additionally a mixed incompressibility brick is added with a variable :math:`p` denoting the pressure, the Neumann term on :math:`u` will depend on :math:`p` in the following way:
:math:`G = \sigma(u)n - pn`

In order to allow a generic implementation in which the brick imposing Nitsche's method will work for every partial differential term applied to the concerned variables, each brick adding a partial differential term to a model is required to give its expression via a GWFL (generic weak form language) expression.

These expressions are utilized in a special method of the model object::

  expr = md.Neumann_term(variable, region)

which allows to automatically derive an expression for the sum of all Neumann terms, by scanning the expressions provided by all partial differential term bricks and performing appropriate manipulations.
Of course it is required that all volumic bricks were added to the model prior to the call of this method.
The derivation of the Neumann term works only for second order partial differential equations.
A generic implementation for higher order pde would be more complicated.


Generic Nitsche's method for a Dirichlet condition
++++++++++++++++++++++++++++++++++++++++++++++++++

Assume that the variable :math:`u` is considered and that one wants to prescribe the condition

.. math::
  Hu = g

on a part :math:`\Gamma_D`  of the boundary of the considered domain.
Here :math:`H` is considered equal to one in the scalar case or can be either the identity matrix in the vectorial case either a singular matrix having only 1 or 0 as eigenvalues.
This allow here to prescribe only the normal or tangent component of :math:`u`.
For instance if one wants to prescribe only the normal component, :math:`H` will be chosen to be equal to :math:`nn^T` where :math:`n` is the outward unit normal on :math:`\Gamma_D`.

Nitsche's method for prescribing this Dirichlet condition consists in adding the following term to the weak formulation of the problem

.. math::
  \int_{\Gamma_D} \dfrac{1}{\gamma}(Hu-g-\gamma HG).(Hv) - \theta(Hu-g).(HD_uG[v])d\Gamma,

where :math:`\gamma` and :math:`\theta` are two parameters of Nitsche's method and :math:`v` is the test function corresponding to :math:`u`.
The parameter :math:`\theta` can be chosen positive or negative. :math:`\theta = 1` corresponds to the more standard method which leads to a symmetric tangent term in standard situations, :math:`\theta = 0` corresponds to a non-symmetric method which has the advantage of a reduced number of terms and not requiring the second derivatives of :math:`G` in the nonlinear case, and :math:`\theta = -1` is a kind of skew-symmetric method which ensures an inconditonal coercivity (which means independent of :math:`\gamma`) at least in standard situations.
The parameter :math:`\gamma` is a kind of penalization parameter (although the method is consistent) which is taken to be :math:`\gamma = \gamma_0 h_T` where :math:`\gamma_0` is taken uniform on the mesh and :math:`h_T` is the diameter of the element :math:`T`.
Note that, in standard situations, except for :math:`\theta = -1` the parameter :math:`\gamma_0` has to be taken sufficiently small in order to ensure the convergence of Nitsche's method.

The bricks adding a Dirichlet condition with Nitsche's method to a model are the following::

  getfem::add_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &Neumannterm,
      const std::string &gamma0name, size_type region,
      scalar_type theta = scalar_type(1),
      const std::string &dataname = std::string());


This function adds a Dirichlet condition on the variable `varname` and the mesh
region `region`. This region should be a boundary. `Neumannterm`
is the expression of the Neumann term (obtained by the Green formula)
described as an expression of GWFL. This term can be obtained with
md.Neumann_term(varname, region) once all volumic bricks have
been added to the model. The Dirichlet
condition is prescribed with Nitsche's method. `dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem; scalar or vector valued, depending on the variable
on which the Dirichlet condition is prescribed. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionally coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionally coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. Returns the brick index in the model.
::


  getfem::add_normal_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &Neumannterm,
      const std::string &gamma0name, size_type region,
      scalar_type theta = scalar_type(1),
      const std::string &dataname = std::string());


This function adds a Dirichlet condition to the normal component of the vector
(or tensor) valued variable `varname` and the mesh region `region`.
This region should be a boundary. `Neumannterm`
is the expression of the Neumann term (obtained by the Green formula)
described as an expression of GWFL. This term can be obtained with
md.Neumann_term(varname, region) once all volumic bricks have
been added to the model. The Dirichlet
condition is prescribed with Nitsche's method. `dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionally coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionally coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. Returns the brick index in the model.
(This brick is not fully tested)
::

  getfem::add_generalized_Dirichlet_condition_with_Nitsche_method
     (model &md, const mesh_im &mim, const std::string &varname,
      const std::string &Neumannterm,
      const std::string &gamma0name, size_type region, scalar_type theta,
      const std::string &dataname, const std::string &Hname);



This function adds a Dirichlet condition on the variable `varname` and the mesh
region `region`.
This version is for vector field. It prescribes a condition
:math:`Hu = r` where :math:`H` is a matrix field. The region should be a
boundary. This region should be a boundary. `Neumannterm`
is the expression of the Neumann term (obtained by the Green formula)
described as an expression of GWFL. This term can be obtained with
md.Neumann_term(varname, region) once all volumic bricks have
been added to the model. The Dirichlet
condition is prescribed with Nitsche's method.
CAUTION : the matrix H should have all eigenvalues equal to 1 or 0.
`dataname` is the optional
right hand side of the Dirichlet condition. It could be constant or
described on a fem. `gamma0name` is the
Nitsche's method parameter. `theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionally coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionally coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary
even for nonlinear problems. `Hname` is the data
corresponding to the matrix field `H`. It has to be a constant matrix
or described on a scalar fem. Returns the brick index in the model.
(This brick is not fully tested)

.. _nitsche_contact_small_def_section:

Generic Nitsche's method for contact with friction condition
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We describe here the use of Nitsche's method to prescribe a contact with Coulomb friction condition in the small deformations framework. This corresponds to a weak integral contact condition which as some similarity with the ones which use Lagrange multipliers describe in the corresponding section, see :ref:`weak_integral_contact_section`

In order to simplify notations, let use denote by :math:`P_{n,\mathscr{F}}` the following map which corresponds to a couple of projections:

.. math::
	P_{n,\mathscr{F}}(x) = -(x.n)_- n + P_{B(0,\mathscr{F}(x.n)_-)}(x - (x.n)n)

This application make the projection of the normal part of :math:`x` on :math:`\rm I\hspace{-0.15em}R_-` and the tangential part on the ball of center :math:`0` and radius :math:`\mathscr{F}(x.n)_-`, where :math:`\mathscr{F}` is the friction coefficient.

Using this, and considering that the sliding velocity is approximated by :math:`\alpha(u_{_T} - w_{_T})` where the expression of :math:`\alpha` and :math:`w_{_T}` depend on the time integration scheme used (see :ref:`weak_integral_contact_section`), Nitsche's term for contact with friction reads as:

.. math::
	&-\int_{\Gamma_C} \theta \gamma G\cdot D_u G[v] d\Gamma \\
	&+\int_{\Gamma_C} \gamma P_{n,\mathscr{F}}(G - \dfrac{Au}{\gamma} + \dfrac{gap}{\gamma}n + \dfrac{\alpha w_{_T}}{\gamma})\cdot(\theta D_u G[v] - \dfrac{v}{\gamma}) d\Gamma.

where :math:`\Gamma_C` is the contact boundary, :math:`G` is the Neumann term which represents here :math:`\sigma n` the stress at the contact boundary and :math:`A` is the :math:`d\times d` matrix

.. math::
	A = \alpha I_d + (1-\alpha)n n^T

Note that for the variant with :math:`\theta=0` a majority of terms vanish.




The following function adds a contact condition with or without Coulomb
friction on the variable
`varname_u` and the mesh boundary `region`.  `Neumannterm`
is the expression of the Neumann term (obtained by the Green formula)
described as an expression of GWFL. This term can be obtained with
md.Neumann_term(varname, region) once all volumic bricks have
been added to the model. The contact condition
is prescribed with Nitsche's method. The rigid obstacle should
be described with the data `dataname_obstacle` being a signed distance to
the obstacle (interpolated on a finite element method).
`gamma0name` is the Nitsche's method parameter.
`theta` is a scalar value which can be
positive or negative. `theta = 1` corresponds to the standard symmetric
method which is conditionally coercive for  `gamma0` small.
`theta = -1` corresponds to the skew-symmetric method which is
inconditionally coercive. `theta = 0` is the simplest method
for which the second derivative of the Neumann term is not necessary.
The optional parameter `dataexpr_friction_coeff` is the friction
coefficient which could be any expression of GWFL.
Returns the brick index in the model.::


  getfem::add_Nitsche_contact_with_rigid_obstacle_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &Neumannterm,
   const std::string &expr_obs, const std::string &dataname_gamma0,
   scalar_type theta_,
   std::string dataexpr_friction_coeff,
   const std::string &dataname_alpha,
   const std::string &dataname_wt,
   size_type region);
