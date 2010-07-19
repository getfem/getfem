.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-solvers:

Predefined solvers
------------------

Of course, for many problems, it will be more convenient to make a specific
solver. Even so, one generic solver is available to test your models quickly. It
can also be taken as an example to build your own solvers. It is defined in
:file:`getfem/getfem_model_solvers.h` and the call is::

  getfem::standard_solve(md, iter);

where ``md`` is the model object and ``iter`` is an iteration object from |gmm|.
See also the next section for an example of use.

Note that |sLU| is used by default on "small" problems. You can also link
|mumps| with |gf| (see section :ref:`ud-linalg`) and used the parallel version.
