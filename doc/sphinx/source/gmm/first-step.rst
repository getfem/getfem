.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-first-step:


First steps with |gmm|
============================


How can I invert a matrix ?
---------------------------

It is not possible in |gmm| to invert all kind of matrices. For the moment, the only mean to invert a matrix is to use the dense LU decomposition (thus, only for dense matrices). An example::

  gmm::dense_matrix<double> M(3, 3), M2(3,3), M3(3,3);
  gmm::copy(gmm::identity_matrix(), M);  // M = Id.
  gmm::scale(M, 2.0);                    // M = 2 * Id.
  M(1,2) = 1.0;

  gmm::copy(M, M2);

  gmm::lu_inverse(M);

  gmm::mult(M, M2, M3);

  std::cout << M << " times " << M2 << " is equal to " << M3 << endl;

see the section corresponding to dense LU decomposition for more details. The type ``gmm::dense_matrix<double>`` can be replaced by ``gmm::row_matrix< std::vector<double> >`` or ``gmm::col_matrix< std::vector<double> >``.

How can I solve a linear system ?
---------------------------------

You have more than one possibility to solve a linear system. If you have a dense matrix, the best may be to use the LU decomposition. An example::

  gmm::dense_matrix<double> M(3, 3);
  gmm::clear(M);                  // M = 0.
  M(0,0) = M(1,1) = M(2,2) = 2.0; // M = 2 * Id.
  M(1,2) = 1.0;

  std::vector<double> X(3), B(3), Bagain(3);
  B[0] = 1.0; B[1] = 2.0; B[2] = 3.0;  // B = [1 2 3]

  gmm::lu_solve(M, X, B);

  gmm::mult(M, X, Bagain);

  std::cout << M << " times " << gmm::vref(X) << " is equal to " << gmm::vref(Bagain) << endl;


If, now, you have a sparse system coming for example from a pde discretization, you have various iterative solvers, with or without preconditioners. This is an example with a precontionned GMRES::

  int nbdof = 1000; // number of degrees of freedom.
  gmm::row_matrix< gmm::rsvector<double> > M(nbdof, nbdof); // a sparse matrix
  std::vector<double> X(nbdof), B(nbdof); // Unknown and left hand side.

  ... here the assembly of the pde discretization stiffness matrix ...
  ... and left hand side ...


  // computation of a preconditioner (ILUT)
  gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(M, 10, 1e-4);

  gmm::iteration iter(1E-8);  // defines an iteration object, with a max residu of 1E-8

  gmm::gmres(M, X, B, P, 50, iter);  // execute the GMRES algorithm

  std::cout << "The result " << gmm::vref(X) << endl;


How can I transform a vector into a matrix and reshape it ?
-----------------------------------------------------------

In |gmm|, a vector is not considered as a matrix. If you need to use a vector as a (1 by n) row matrix or (n by 1) column matrix in a computation, you have to use::

   gmm::row_vector(V) // gives a reference on V considered as
                      // a (1 by n) row matrix
   gmm::col_vector(V) // gives a reference on V considered as
                      // a (n by 1) col matrix

for instance, you can transform a vector into a dense matrix with::

  std::vector<double> V(50);

  // ... computation of V

  gmm::dense_matrix<double> M(1, gmm::vect_size(V));
  gmm::copy(gmm::row_vector(V), M);


Then you can also reshape matrix ``M`` with::

  gmm::reshape(M, 10, 5);


What is the better way to resize a matrix ?
-------------------------------------------

You can change the dimensions of a matrix, if it is not a reference, using::

  gmm::resize(M, m, n);

This function respects the intersection between the original matrix and the resized matrix, and new components are set to zero. An important thing is that it is based on the resize method of ``std::vector``, thus no memory free is done when the size of the new matrix is smaller than the original one.

If you do not need to keep old values of the components, or if you want to really free the surplus of memory, you can resize a matrix using ``std::swap`` as follows::

  MATRIX_TYPE M(m1, n1);

  ... your code

  { MATRIX_TYPE(m2, n2) M2; std::swap(M, M2); } // resize matrix M.

Of course, this works also for a vector.
