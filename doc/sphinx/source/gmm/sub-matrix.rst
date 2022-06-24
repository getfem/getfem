.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-sub:


sub-vectors and sub-matrices
============================

It is possible to obtain any sub-vector or sub-matrix of a fully interfaced object. There are four types of sub indices::

  gmm::sub_interval(first, length);

represents an interval whose first index is ``first`` and length is ``length`` ( for instance ``gmm::sub_interval(10, 3);`` represents the indices ``{10, 11, 12}``).

::

  gmm::sub_slice(first, length, step);

represents also an interval in which one index over ``step`` is taken. ( for instance ``gmm::sub_slice(10, 3, 2);`` represents the indices ``{10, 12, 14}``)

::

  gmm::sub_index(CONT c);

represents the sub-index which is the collection of index contained in the container ``c``. For instance::

  std::vector<size_t> c(3);
  c[0] = 1; c[1] = 3; c[2] = 16;
  gmm::sub_index(c);


represents the indices ``{1, 3, 16}``.

`VERY IMPORTANT` : the container ``c`` has to be `sorted` from the smaller index to the greater one (i.e. with increasing order) and no repetition is allowed.


For unsorted index such as permutation, a special type of sub index is defined::

  gmm::unsorted_sub_index(CONT c);

Some algorithms are a little bit slower with unsorted sub indices.

Now ``gmm::sub_vector(V, subi)`` gives a reference to a sub-vector::

  gmm::vsvector<double> V(10);
  V[5] = 3.0;
  std::cout << gmm::sub_vector(V, gmm::sub_interval(2, 3)) << std::endl;

prints to the standard output ``V[2], V[3]`` and ``V[4]``.

``gmm::sub_matrix(V, subi1, subi2)`` gives a reference to a sub-matrix. For instance::

  gmm::col_matrix< gmm::wsvector<double> > M(5, 20);
  M(3, 2) = 5.0;
  std::cout << gmm::sub_matrix(M, gmm::sub_interval(2, 3), gmm::sub_interval(2, 3))
            << std::endl;

prints to the output a sub-matrix. If the two sub-indices are equal, it is possible to omit the second. For instance::

  gmm::col_matrix< gmm::wsvector<double> > M(5, 20);
  M(3, 2) = 5.0;
  std::cout << gmm::sub_matrix(V, gmm::sub_interval(2, 3)) << std::endl;

The reference on sub_matrix is writable if the corresponding matrix is writable (so you can copy on a sub_matrix, add sub-matrices ...).

row and column of a matrix
--------------------------

``gmm::mat_row(M, i)`` gives a (possibly writable) reference to the row ``i`` of matrix ``M``, and ``gmm::mat_col(M, i)``  gives a (possibly writable) reference to the column ``i``. It is not possible to access to the rows if ``M`` is a column matrix and to the columns if it is a row matrix. It is possible to use ``gmm::mat_const_row(M, i)`` and ``gmm::mat_const_col(M, i)`` to have constant references.
