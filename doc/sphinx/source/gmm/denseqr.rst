.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-denseqr:

Dense QR factorisation, eigenvalues and eigenvectors
====================================================


The following procedures are available in the file ``gmm/gmm\_dense\_qr.h`` for dense real and complex matrices::


  gmm::qr_factor(M, Q, R) // compute the QR factorization of M in Q and R
                          // (Householder version)

  implicit_qr_algorithm(M, eigval, double tol = 1E-16) // compute the
     // eigenvalues of M using the implicit QR factorisation (Householder and
     // Francis QR step version). eigval should be a vector of appropriate size
     // in which the eigenvalues will be computed. If the matrix have
     // complex eigenvalues, please use a complex vector.

  implicit_qr_algorithm(M, eigval, shvect, double tol = 1E-16) // idem,
     // compute additionally the schur vectors in the matrix shvect.

  symmetric_qr_algorithm(M, eigval, double tol = 1E-16) // idem for symmetric
     // real and hermitian complex matrices (based on Wilkinson QR step)

  symmetric_qr_algorithm(M, eigval, eigvect, double tol = 1E-16) // idem,
     // compute additionally the eigenvectors in the matrix eigvect.



`Remark`: The computation of eigenvectors for non hermitian matrices is not yet implemented. You can use for the moment the functions ``geev_interface_left`` and ``geev_interface_right`` from the LAPACK interface (see ``gmm/gmm_lapack_interface.h``). These LAPACK functions compute right and left eigenvectors.


The following function defined in the file ``gmm/gmm\_condition\_number.h``::

   gmm::condition_number(M)

compute the condition number of a matrix ``M``. This function uses a dense QR algorithm and thus is only usable for dense matrices.
