.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-solvers:

Predefined solvers
------------------

Although it will be more convenient to build a specific
solver for some problems, a generic solver is available to test your models quickly. It
can also be taken as an example to build your own solver. It is defined in
:file:`src/getfem/getfem_model_solvers.h` and :file:`src/getfem_model_solvers.cc` and the call is::

  getfem::standard_solve(md, iter);

where ``md`` is the model object and ``iter`` is an iteration object from |gmm|.
See also the next section for an example of use.

Note that |sLU| is used as a default linear solver on "small" problems. You can also link |mumps| with |gf| (see section :ref:`ud-linalg`) and use the parallel version. For nonlinear problems, A Newton method (also called Newton-Raphson method) is used.

Note also that it is possible to disable some variables
(with the method md.disable_variable(varname) of the model object) in order to
solve the problem only with respect to a subset of variables (the
disabled variables are the considered as data) for instance to
replace the global Newton strategy with a fixed point one.

Let us recall that a standard initialization for the iter object is the folowwing (see Gmm++ documentation on :ref:`gmm-iter`)::

  gmm::iteration iter(1E-7, 1, 200);

where ``1E-7`` is the relative tolerance for the stopping criterion, `1` is the noisy option and `200` is the maximum number of iterations. The stopping criterion of Newton's method is build as follows. For a relative tolerance :math:`\varepsilon`, the algorithm stops when:

.. math::

  \min\left( \|F(u)\|_1 / \max(L, 10^{-25}) ~, ~~ \|h\|_1 / \max(\|u\|_1, 10^{-25})\right) < \varepsilon

where :math:`F(u)` is the residual vector, :math:`\|\cdot\|_1` is the classical 1-norm in :math:`\R^n`, :math:`h` is the search direction given by Newton's algorithm, :math:`L` is the norm of an estimated external loads (coming from source term and Dirichlet bricks) and :math:`u` is the current state of the searched variable. The maximum taken with :math:`10^{-25}` is to avoid pathological cases when :math:`L` and/or :math:`u` are vanishing.



