.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-triangular:

Solving triangular systems
========================================


If ``M`` is a triangular matrix (upper or lower) and ``X`` a vector containing the right hand side, the following procedures solve the system :math:`x \leftarrow M^{-1}x`. The vector ``X`` contains the result::

   gmm::upper_tri_solve(M, X, false) // Solving an upper triangular system
   gmm::upper_tri_solve(M, X, true)  // Solving an upper triangular system
                                     // assuming there is 1 on the diagonal
   gmm::lower_tri_solve(M, X, false) // Solving a lower triangular system
   gmm::lower_tri_solve(M, X, true)  // Solving a lower triangular system
                                     // assuming there is 1 on the diagonal

components which are lower the diagonal are ignored by ``gmm::upper_tri_solve`` and components which are upper the diagonal are ignored by ``gmm::lower_tri_solve``.
