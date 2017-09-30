.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-denselu:

Dense LU decomposition
========================================


The following procedures are available in the file ``gmm/gmm\_dense\_lu.h`` for dense real and complex matrices (``gmm::dense_matrix<T>``, ``gmm::row_matrix< std::vector<T> >`` and ``gmm::col_matrix< std::vector<T> >``)::

  gmm::lu_factor(M, ipvt) : compute the LU factorization of M in M. ipvt should be
                       an gmm::lapack_ipvt (of size gmm::mat_nrows(M))
                       which will contain the indices of the pivots.

  gmm::lu_solve(LU, ipvt, x, b) : solve the system LUx = b. LU is the LU
                             factorization which has to be computed first.

  gmm::lu_solve(M, x, b) : solve the system Mx=b calling the lu factorization on
                      a copy of M.

  gmm::lu_solve_transposed(LU, ipvt, x, b) : solve the system transposed(LU)x = b.
                                        LU is the LU factorization which
                                        has to be computed first.

  gmm::lu_inverse(LU, ipvt, A) : compute the inverse of LU in A. LU is the LU
                            factorization which has to be computed first

  gmm::lu_inverse(A) : invert A calling the LU factorization and the latter
                  procedure.

  gmm::lu_det(LU, ipvt) : compute the determinant of LU. LU is the LU
                     factorization which has to be computed first

  gmm::lu_det(A) : compute the determinant of A calling the LU factorization
              and the latter function.
