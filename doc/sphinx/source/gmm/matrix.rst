.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-matrix:

Matrix and Vector type provided by |gmm|
========================================


The convention is that any vector or matrix type (except if it is a  reference)
can be instantiated with the constructors::

  Vector V(n);        // build a vector of size n.
  Matrix M(n, m);     // build a matrix with n rows and m columns.

No other constructor is used inside |gmm| and you should not use any other if you want your code
to be compatible with any matrix and vector type.

It is assumed that each vector type interfaced with |gmm| allows to
access to a component with the following syntax::

  a = V[i];    // read the ith component of V.
  V[i] = b;    // write the ith component of V.


The write access being available if the vector is not a constant reference. For a matrix::

  a = M(i, j); // read the component at row i and column j of M.
  M(i, j) = b; //  write the component at row i and column j of M.

Again the write access is available if the matrix is not a const reference. Generally, especially for sparse matrices, this access is not very efficient. Linear algebra procedures access to the components of the vectors and matrices via iterators. (see section  :ref:`gmm-inside`)

It is also not recommended (at all) to use the original copy operator for vectors or matrices. Generally, it will not do the appropriate job. instead, you have to use the method::

  gmm::copy(V, W);  //  W <-- V


which works for all correctly interfaced matrix and vector type, even if ``V`` is not of the same type as ``W`` (``V`` could be sparse and ``W`` dense for instance).

in |gmm|, a vector is not a (n by 1) matrix, it is a one dimensional object. If you need to use a vector as a (n by 1) column matrix or a (1 by n) row matrix, you can do it with::

   gmm::row_vector(V) // gives a reference on V considered as
                      // a (1 by n) row matrix
   gmm::col_vector(V) // gives a reference on V considered as
                      // a (n by 1) col matrix

In the following, the template parameter ``T`` will represent a scalar type like ``double`` or ``std::complex<double>``.


dense vectors
-------------

|gmm| interfaces ``std::vector<T>`` so you can use it as your basic dense vector type.
If you need to interface another type of dense vector you can see in ``gmm/gmm_interface.h``
some examples.

sparse vectors
--------------

|gmm| provides two types of sparse vectors: ``gmm::wsvector<T>`` and ``gmm::rsvector<T>``. ``gmm::wsvector<T>`` is optimized for write operations and ``gmm::rsvector<T>`` is optimized for read operations. It should be appropriate to use ``gmm::wsvector<T>`` for assembling procedures and then to copy the vector in a ``gmm::rsvector<T>`` for the solvers. Those two vector types can be used to create row major or column major matrices (see section  :ref:`gmmracmat`).

skyline vectors
---------------

The type ``gmm::slvector<T>`` defines a skyline vector, in the sense that only an interval of this vector is stored. With this type of vector you can build skyline matrices as ``gmm::row_matrix< gmm::slvector<T> >`` (see next section :ref:`gmmracmat`).

.. _gmmracmat:

generic row and column matrices
-------------------------------

|gmm| provides the two following types of matrices: ``gmm::row_matrix<VECT>`` and ``gmm::col_matrix<VECT>`` where ``VECT`` should be a valid (i.e. interfaced) vector type.
Those two type of matrices store an array of ``VECT`` so the memory is not contiguous. Initializations are::

  gmm::row_matrix< std::vector<double> > M1(10, 10);  // dense row matrix
  gmm::col_matrix< gmm::wsvector<double> > M2(5, 20); // sparse column matrix

Of course ``gmm::row_matrix<VECT>`` is a row matrix and it is impossible to access to a particular column of this matrix.


``gmm::mat_nrows(M)`` gives the number of rows of a matrix and ``gmm::mat_ncols(M)`` the number of columns.

dense matrices
--------------

It is recommended to use the type::

  gmm::dense_matrix<T>

to represent a dense matrix type because it is compatible with the Fortran format (column major) and some operations are interfaced with blas and Lapack (see section  :ref:`gmm-lapack`). It is considered as a column and row matrix (column preferred) which means that you can access both to the columns and rows.

However, matrix types as ``gmm::row_matrix< std::vector<double> >`` or ``gmm::col_matrix< std::vector<double> >`` represent also some dense matrices.

sparse matrices
---------------

Similarly, ``gmm::row_matrix< gmm::wsvector<double> >`` or ``gmm::col_matrix< gmm::rsvector<double> >`` represents some sparse matrices, but |gmm| provides also two types of classical sparse matrix types::

  gmm::csr_matrix<T>
  gmm::csc_matrix<T>

The type ``gmm::csr_matrix<T>`` represents a compressed sparse row matrix and ``gmm::csc_matrix<T>`` a compressed sparse column matrix. The particularity of these two types of matrices is to be read only, in the sense that it is not possible to access at a particular component to write on it (the operation is too expansive). The only write operation permitted is ``gmm::copy``. The right way to use these matrices is first to execute the write operations on another type of matrix like ``gmm::row_matrix< gmm::wsvector<double> >`` then to do a copy::

  gmm::row_matrix< gmm::wsvector<double> > M1;
  ...
  assembly operation on M1
  ...
  M1(i,j) = b;
  ...
  gmm::csc_matrix<double> M2;
  gmm::clean(M1, 1E-12);
  gmm::copy(M1, M2);

Matrices ``gmm::csr_matrix<T>`` and ``gmm::csc_matrix<T>`` have the advantage to have a standard format (interfacable with Fortran code) and to have a compact format (contiguous in memory). To be able to be compatible with Fortran programs a second template parameter exists on these type, you can declare::

  gmm::csc_matrix<double, 1> M1;
  gmm::csr_matrix<double, 1> M2;

The ``1`` means that a shift will be done on all the indices.
