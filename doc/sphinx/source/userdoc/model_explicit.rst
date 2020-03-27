.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. index:: models, model bricks

.. _ud-model-explicit:



Other "explicit" bricks
-----------------------

Two (very simple) bricks allow to add some explicit terms to the tangent system.

The function::

  indbrick = getfem::add_explicit_matrix(md, varname1, varname2, B
                                         issymmetric = false,
                                         iscoercive = false);

adds a brick which just adds the matrix ``B`` to the tangent system relatively to
the variables ``varname1`` and ``varname2``. The given matrix should have as many
rows as the dimension of ``varname1`` and as many columns as the dimension of
``varname2``. If the two variables are different and if ``issymmetric`` is set to
true then the transpose of the matrix is also added to the tangent system (default
is false). Set ``iscoercive`` to true if the term does not affect the coercivity
of the tangent system (default is false). The matrix can be changed by the
command::

  getfem::set_private_data_matrix(md, indbrick, B);

The function::

  getfem::add_explicit_rhs(md, varname, L);

adds a brick which just add the vector ``L`` to the right hand side of the tangent
system relatively to the variable ``varname``. The given vector should have the
same size as the variable ``varname``. The value of the vector can by changed by
the command::

  getfem::set_private_data_rhs(md, indbrick, L);

