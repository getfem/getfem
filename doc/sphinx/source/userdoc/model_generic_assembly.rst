.. $Id: model_generic_assembly.rst 3655 2010-07-19 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-generic-assembly:


Generic assembly bricks
-----------------------


A mean to add a term either on one variable or on several ones is to directly use GWFL, the generic weak form language described in Section :ref:`ud-gasm-high`. The more general way is to use::

   size_type getfem::add_nonlinear_term(md, mim, expr,
                         region = -1, is_sym = false, is_coercive = false);

This adds a brick to the model ``md``, using the integration method ``mim``, the assembly string ``expr`` on the mesh region ``region``. If the result is symmetric, you can specify it on the 5th argument and if it is coercive on the 6th argument. The latter indications of symmetry and coercivness are used to determine the right linear solver. If you are not so sure, it is preferable not to indicate anything.

However, this brick consider that the expression is nonlinear. This brick is especially indicated to obtain nonlinear coupled terms between several variables. This means in particular that the assembly of the term is performed at each call of the assembly of the model and that a Newton algorithm will be used to solve the problem. If the term is indeed linear, you should use instead::

  size_type getfem::add_linear_term(md, mim, expr,
                         region = -1, is_sym = false, is_coercive = false);

with the same arguments. Conversely, this brick alway assume that the term corresponding to ``expr`` is linear and the assembly will be performed only once if the data used do not change. Thus, you have to care that your expression is indeed linear (affine in fact) with respect to each variable. Otherwise, the result is of course not guaranted. Source terms in the expression are taken into account. Still for linear problem, it is possible to perform the assembly of a sole source term thanks to::

  size_type getfem::add_source_term(md, mim, expr, region = -1);

with again the same arguments except the symmetry and coercivness. This brick performs the assembly of the corresponding order 1 term (residual vector) and add it as a right hand side to the problem. The assembly will be performed only once, so the term should not depend on the variables of the model (but could depend of course on the constants).


For instance, if one wants to solve a Poisson problem on a predefined variable ``u`` of the model, one may use the corresponding pre-defined bricks (see below) or simply use::

  getfem::add_nonlinear_term(md, mim, "Grad_u.Grad_Test_u - F*Test_u", -1, true, true);

where ``F`` is a pre-defined constant of the model representing the right hand side. Of course, doing so, Newton's algorithms will be called. So, the more appropriate manner is to use the linear bricks as follows::

  getfem::add_linear_term(md, mim, "Grad_u.Grad_Test_u", -1, true, true);
  getfem::add_source_term(md, mim, "F*Test_u");






Note that for the moment, the use of GWFL is not possible for complex valued problems.
