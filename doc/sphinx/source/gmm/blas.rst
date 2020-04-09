.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-blas:

Basic linear algebra operations
========================================


The same choice has been made as in MTL to provide basic operations as functions not as operators. The advantages are that it is clearer to see where are the linear algebra operations in the program and the programming of optimized basic linear algebra operations is greatly simplified.


scale and scaled
----------------

``gmm::scale`` is used to multiply a vector or a matrix with a scalar factor::

  gmm::scale(V, 10.0);  // V * 10.0 ---> V

If one not needs to multiply the vector but wants to use the multiplied vector in an expression  ``gmm::scaled`` gives a reference to a multiplied vector. This is only a reference, no operation is made until this reference is used somewhere. For instance::

  std::cout << gmm::scaled(V, 10.0) << std::endl;

print to the standard output the vector ``V`` multiplied by ``10.0`` without changing ``V``.

transposition
-------------

``gmm::transposed(M)`` gives a possibility modifiable reference on the transposed matrix of ``M``.

imaginary and real part
-----------------------

For a complex matrix ``M`` or a complex vector ``V``,
``gmm::real_part(M)``, ``gmm::real_part(V)``, ``gmm::imag_part(M)`` or ``gmm::imag_part(V)`` give a possibility modifiable reference on the real or imaginary part of the matrix or vector (for instance ``gmm::clear(gmm::imag_part(M))`` will set to zero the imaginary part of a matrix ``M``). These functions cannot be applied to real matrices or vectors.

conjugate
---------

For a matrix ``M`` or a vector ``V``,
``gmm::conjugated(M)`` and ``gmm::conjugated(V)`` give a constant reference on the conjugated vector or matrix. Of course, for a real vectors this has no effect (and no cost at all). Note : ``gmm::conjugated(M)`` transposes the matrix ``M`` so that this is the hermitian conjugate of ``M``. If you need only the conjugate of each component you have to use both transposition and conjugate with ``gmm::conjugated(gmm::transposed(M))`` or equivalently  ``gmm::transposed(gmm::conjugated(M))``.


add
---

addition of vectors or matrices. It is alway possible to mix different type of vector or matrices in the operations. The following operations are valid::

  std::vector<double> V1(10);
  gmm::wsvector<double> V2(10);
  gmm::clear(V1);
  ...
  gmm::add(V1, V2); // V1 + V2 --> V2
  cout << gmm::vref(V2);

  gmm::add(V1, gmm::scaled(V2, -2.0), V2); // V1 - 2.0 * V2 --> V2
  cout << gmm::vref(V2);

  gmm::row_matrix< std::vector<double> > M1(10, 10);
  gmm::col_matrix< gmm::wsvector<double> > M2(1000, 1000);

  // M1 + (sub matrix of M2) ---> (sub matrix of M2)
  gmm::add(M1, gmm::sub_matrix(M2, gmm::sub_interval(4,10)));


IMPORTANT : all the vectors have to have the same size, no resize will be automatically done. If a vector has not the good size, an error will be thrown.

mult
----

Matrix-vector or matrix-matrix multiplication. Again, all the matrices and vectors have to have the good size. The following operations are valid::

  std::vector<double> V1(10);
  gmm::wsvector<double> V2(10);
  ...
  gmm::row_matrix< std::vector<double> > M1(10, 10);
  ...

  gmm::mult(M1, V2, V1);  // M1 * V2 --> V1

  gmm::mult(M1, V2, V2, V1);  // M1 * V2 + V2 --> V1

  gmm::mult_add(M1, V2, V1);  // M1 * V2 + V1 --> V1

  gmm::mult(M1, gmm::scaled(V2, -1.0), V2, V1);  // M1 * (-V2) + V2 --> V1

  gmm::col_matrix< gmm::wsvector<double> > M2(10, 10);
  gmm::col_matrix< gmm::vsvector<double> > M3(10, 10);
  ...

  gmm::mult(M1, M2, M3); // M1 * M2 ---> M3

  gmm::mult(gmm::sub_matrix(M1, sub_interval(0, 3)),
            gmm::sub_matrix(M2, sub_interval(4, 3)),
            gmm::sub_matrix(M3, sub_interval(2, 3)));



norms
-----

::

  gmm::vect_norm1(V)  // sum of the modulus of the components of vector V.
  gmm::vect_norm2(V)  // Euclidean norm of vector V.
  gmm::vect_dist2(V1, V2)  // Euclidean distance between V1 and V2.
  gmm::vect_norminf(V)    // infinity norm of vector V.
  gmm::mat_euclidean_norm(M) // Euclidean norm of matrix ``M``
                             // (called also Frobenius norm).
  gmm::mat_maxnorm(M) // Max norm (defined as max(|m_ij|; i,j ))
  gmm::mat_norm1(M)   // max(sum(|m_ij|, i), j)
  gmm::mat_norminf(M) // max(sum(|m_ij|, j), i)


trace
-----

``gmm::mat_trace(M)`` gives the trace of matrix ``M``.

scalar product
--------------


  for vectors only, ``gmm::vect_sp(V1, V2)`` gives the scalar product between ``V1`` and ``V2``. For complex vectors, this do not conjugate ``V1``, you can use ``gmm::vect_sp(V1, gmm::conjugated(V2))`` or ``gmm::vect_hp(V1, V2)`` which is equivalent.
