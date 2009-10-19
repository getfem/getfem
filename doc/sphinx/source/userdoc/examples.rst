.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-examples:

Example: Laplacian program
==========================

The program ``laplacian`` is provided in the directory ``tests`` of |gf|
distribution. This program computes the solution of the Poisson problem in a
parellepiped domain in any dimension with various finite element methods and
elements. This program can be used as a model to build application programs. It is
built when a ``make check`` is done on the root directory of |gf| (or just with
``cd tests; make laplacian``).

Once the program is compiled you can test it executing the command::

  $ cd tests
  $ ./laplacian laplacian.param

The file ``laplacian.param`` is the parameter file. You can edit it and test
various situation. The program prints the :math:`L^2` and :math:`H^1` error from
an exact solution.

The program ``elastostatic`` is built in a same way and compute the solution of
linear elasticity problem. Many more examples can be found in the tests directory.
