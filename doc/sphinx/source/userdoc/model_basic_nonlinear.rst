.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-basic-nonlinear:


Basic nonlinear brick
---------------------

This brick represents a weak term of the form

.. math::

   \int_{\Omega} f(u,\lambda)\cdot v\ dx + \ldots

where :math:`f` is a function given by a string and interpreted with muparser. This means in particular that this bricks need that Getfem++ is build with muparser library being installed. Here also, :math:`u` is the unknown and  :math:`\lambda` is an optional real parameter. This brick can be used to add basic nonlinear term such as :math:`u^2` or :math:`e^u`.

The function which adds this brick to a model is::

  ind_brick = add_basic_nonlinear_brick(md, mim, varname, const std::string &f,
		const std::string &dfdu, region = size_type(-1),
		dataname_parameter = "");

where ``varname`` is the name of the variable on which the term will be added, ``f`` is the string containing the expression of the function,  ``dfdu`` is the string containing the expression of the derivative of the function with respect to the variable,  ``region`` is an optional mesh region and  ``dataname_parameter`` is the name of the optional scalar parameter.

Note that in the expression of ``f`` the variable and the parameter should be represented by their respective names. For instance, to add the nonlinear term :math:`\lambda e^u` on a model on a variable ``u`` and a parameter ``lambda``, the command is::

	add_basic_nonlinear_brick(md, mim, "u", "lambda*exp(u)",
	                          "lambda*exp(u)", size_type(-1), "lambda");

