.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-superlu:


Interface with SuperLU
============================


It is possible to call SuperLU 3.0 (https://portal.nersc.gov/project/sparse/superlu/superlu_3.0.tar.gz) from |gmm|. The following function defined in the file ``gmm/gmm_superlu_interface.h`` is available::

  SuperLU_solve(A, X, B, condest, permc_spec = 1)

solves the system ``AX = B`` where A is a sparse matrix of base type ``float, double, std::complex<float>, or std::complex<double>``. ``permc_spec`` should be 0, 1 or 2 for respectively use the natural ordering, use minimum degree ordering on structure of ``A'A`` or use minimum degree ordering on structure of ``A'+A`` (1 is the default value), ``condest`` should be a reference on a double, it returns an estimate of the condition number of the matrix ``A``.

To use these functions, you need to install SuperLU and compile your code with the additional options::

  g++ ...  -DGMM_USES_SUPERLU (dir_of_superlu)/superlu.a -lblas -I(dir_of_superlu)

Some other functionalities of SuperLU can be interfaced.

