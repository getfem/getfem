.. $Id$

.. include:: ../replaces.txt

.. highlight:: none

.. _gmm-inside:


Deeper inside |gmm|
============================


The linalg_traits structure
---------------------------

The major principle of |gmm| is that each vector and matrix type has a corresponding structure (which is never instantiated) named ``linalg_traits`` containing all informations on it. For instance, the component ``linalg_type`` of this structure is set to ``abstract_vector`` or ``abstract_matrix`` if the corresponding type represent a vector or a matrix. If ``V`` is an interfaced type of vector and ``M`` an interface type of matrix, it is possible to access to this component with::

  typename gmm::linalg_traits<V>::linalg_type ...  // should be abstract_vector
  typename gmm::linalg_traits<M>::linalg_type ...  // should be abstract_matrix

The types ``abstract_vector`` and ``abstract_matrix`` are defined in ``gmm/gmm_def.h``. They are void type allowing to specialize generic algorithms.

For a vector type, the following informations are available::

  typename gmm::linalg_traits<V>::value_type     --> type of the components of the
                                                     vector
  typename gmm::linalg_traits<V>::reference      --> type of reference on a component
  typename gmm::linalg_traits<V>::is_reference   --> if the vector is a simple
                                                     reference or an instantiated vector
  typename gmm::linalg_traits<V>::linalg_type    --> should be abstract_vector
  typename gmm::linalg_traits<V>::index_sorted    --> linalg_true or linalg_false
  typename gmm::linalg_traits<V>::const_iterator --> const iterator to iterate on the
                                                     components of the vector in
                                                     order to read them.
  typename gmm::linalg_traits<V>::iterator       --> iterator to iterate on the
                                                     components of the vector in
                                                     order to read or write them.
  typename gmm::linalg_traits<V>::storage_type   --> should be abstract_sparse,
                                                     abstract_skyline or
                                                     abstract_dense

  typename gmm::linalg_traits<V>::origin_type    --> the type of vector itself
                                                     or the type of referenced
                                                     vector for a reference.

  gmm::linalg_traits<V>::size(v)     --> a method which gives the size of the vector.
  gmm::linalg_traits<V>::begin(v)    --> a method which gives an iterator on the
                                         beginning of the vector
  gmm::linalg_traits<V>::end(v)      --> iterator on the end of the vector
  gmm::linalg_traits<V>::origin(v)   --> gives a void pointer allowing to identify
                                         the vector
  gmm::linalg_traits<V>::do_clear(v) --> make a clear on the vector

  gmm::linalg_traits<V>::access(o, it, ite, i) --> return the ith component or a
                                          reference on the ith component. o is a
                                          pointer o type ``origin_type *'' or
                                          ``const origin_type *''.

  gmm::linalg_traits<V>::clear(o, it, ite) --> clear the vector. o is a
                                          pointer o type ``origin_type *'' or
                                          ``const origin_type *''.

and for a matrix type::

  typename gmm::linalg_traits<M>::value_type     --> type of the components of the
                                                     matrix
  typename gmm::linalg_traits<M>::reference      --> type of reference on a component
  typename gmm::linalg_traits<M>::is_reference   --> if the matrix is a simple
                                                     reference or an instantiated matrix
  typename gmm::linalg_traits<M>::linalg_type    --> should be abstract_matrix
  typename gmm::linalg_traits<M>::storage_type   --> should be abstract_sparse,
                                                     abstract_skyline or
                                                     abstract_dense
  typename gmm::linalg_traits<M>::index_sorted    --> linalg_true or linalg_false
  typename gmm::linalg_traits<M>::sub_orientation --> should be row_major, col_major
                                                      row_and_col or col_and_row.
  typename gmm::linalg_traits<M>::sub_col_type      --> type of reference on a column
                                                      (if the matrix is not row_major)
  typename gmm::linalg_traits<M>::const_sub_col_type --> type of const reference on a
                                                       column
  typename gmm::linalg_traits<M>::col_iterator      --> iterator on the columns
  typename gmm::linalg_traits<M>::const_col_iterator --> const iterator on the columns
  typename gmm::linalg_traits<M>::sub_row_type      --> type of reference on a row
                                                      (if the matrix is not col_major)
  typename gmm::linalg_traits<M>::const_sub_row_type --> type of const reference on a
                                                       row
  typename gmm::linalg_traits<M>::const_row_iterator --> const iterator on the rows
  typename gmm::linalg_traits<M>::row_iterator       --> iterator on the rows

  typename gmm::linalg_traits<M>::origin_type    --> the type of vector itself
                                                     or the type of referenced
                                                     vector for a reference.

  gmm::linalg_traits<M>::nrows(m)     --> methods which gives the number of rows of
                                          the matrix
  gmm::linalg_traits<M>::ncols(m)     --> number of columns
  gmm::linalg_traits<M>::row_begin(m) --> iterator on the first row (if not col_major)
  gmm::linalg_traits<M>::row_end(m)   --> iterator on the end of the rows
  gmm::linalg_traits<M>::col_begin(m) --> iterator on the first column
                                          (if not row_major)
  gmm::linalg_traits<M>::col_end(m)   --> iterator on the end of the columns
  gmm::linalg_traits<M>::row(it)      --> gives the reference on a row with an iterator
                                          (if not col_major)
  gmm::linalg_traits<M>::col(it)      --> gives the reference on a column with an
                                          iterator  (if not row_major)
  gmm::linalg_traits<M>::origin(m)    --> gives a void pointer allowing to identify
                                          the matrix
  gmm::linalg_traits<M>::access(it,i) --> return the ith component or a reference
                                          on the ith component of the row or
                                          column pointed by it.
  gmm::linalg_traits<M>::do_clear(m)  --> make a clear on the matrix


This is this structure you have to fill in to interface a new vector or matrix type. You can see some examples in ``gmm/gmm_interface.h`` . Most of the generic algorithms are in ``gmm/gmm_blas.h`` .


How to iterate on the components of a vector
--------------------------------------------

Here is an example which accumulate the components of a vector. It is assumed that ``V`` is a vector type and ``v`` an instantiated vector::

  typename gmm::linalg_traits<V>::value_type r(0); // scalar in which we accumulate
  typename gmm::linalg_traits<V>::const_iterator it = vect_const_begin(v); // beginning of v
  typename gmm::linalg_traits<V>::const_iterator ite = vect_const_end(v); // end of v

  for (; it != ite; ++it)  // loop on the components
    r += *it;              // accumulate the components

This piece of code will work with every kind of interfaced vector.

For sparse or skyline vectors, it is possible to obtain the index of the components pointed by the iterator with ``it.index()``. Here is the example of the scalar product of two sparse or skyline vectors, assuming ``V1`` and ``V2`` are two vector types and ``v1``, ``v2`` two corresponding instantiated vectors::

   typename gmm::linalg_traits<V1>::const_iterator it1 = vect_const_begin(v1),
   typename gmm::linalg_traits<V1>::const_iterator ite1 = vect_const_end(v1);
   typename gmm::linalg_traits<V2>::const_iterator it2 = vect_const_begin(v2),
   typename gmm::linalg_traits<V2>::const_iterator ite2 = vect_const_end(v2);
   typename gmm::linalg_traits<V1>::value_type r(0); // it is assumed that V2 have a
                                                // compatible value_type

   while (it1 != ite1 && it2 != ite2) {  // loops on the components
     if (it1.index() == it2.index()) {
       res += (*it1) * (*it2));          // if the indices are equals accumulate
       ++it1;
       ++it2;
     }
     else if (it1.index() < it2.index())
       ++it1;
     else
       ++it2;
   }

This algorithm use the fact that indices are increasing in a sparse vector. This code will not work for dense vectors because dense vector iterators do not have the method ``it.index()``.

How to iterate on a matrix
--------------------------

You can iterate on the rows of a matrix if it is not a column major matrix and on the columns of a matrix if it is not a row major matrix (the type ``gmm::dense_matrix<T>`` has is sub orientation type as col_and_rox, so you can iterate on both rows and columns).

If you need not to be optimal, you can use a basic loop like that::

  for (size_t i = 0; i < gmm::mat_nrows(m); ++i) {
    typename gmm::linalg_traits<M>::const_sub_row_type row = mat_const_row(M, i);

    ...

    std::cout << "norm of row " << i << " : " << vect_norm2(row) << std::endl;
  }

But you can also use iterators, like that::

  typename gmm::linalg_traits<M>::const_row_iterator it = mat_row_const_begin(m);
  typename gmm::linalg_traits<M>::const_row_iterator ite = mat_row_const_end(m);

  for (; it != ite; ++it) {
    typename gmm::linalg_traits<M>::const_sub_row_type
      row = gmm::linalg_traits<M>::row(it);

    ...

    std::cout << "norm of row " << i << " : " << vect_norm2(row) << std::endl;
  }


How to make your algorithm working on all type of matrices
----------------------------------------------------------

For this, you will generally have to specialize it. For instance, let us take a look at the code for ``gmm::nnz`` which count the number of stored components (in fact, the real ``gmm::nnz`` algorithm is specialized in most of the cases so that it does not count the components one by one)::

  template <class L> inline size_type nnz(const L& l) {
    return nnz(l, typename linalg_traits<L>::linalg_type());
  }

  template <class L> inline size_type nnz(const L& l, abstract_vector) {
    typename linalg_traits<L>::const_iterator it = vect_const_begin(l);
    typename linalg_traits<L>::const_iterator ite = vect_const_end(l);
    size_type res(0);
    for (; it != ite; ++it) ++res;
    return res;
  }

  template <class L> inline size_type nnz(const L& l, abstract_matrix) {
    return nnz(l,  typename principal_orientation_type<typename
                   linalg_traits<L>::sub_orientation>::potype());
  }

  template <class L> inline size_type nnz(const L& l, row_major) {
    size_type res(0);
    for (size_type i = 0; i < mat_nrows(l); ++i)
      res += nnz(mat_const_row(l, i));
    return res;
  }

  template <class L> inline size_type nnz(const L& l, col_major) {
    size_type res(0);
    for (size_type i = 0; i < mat_ncols(l); ++i)
      res += nnz(mat_const_col(l, i));
    return res;
  }


The first function dispatch on the second or the third function respectively if the parameter is a vector or a matrix. The third function dispatch again on the fourth and the fifth function respectively if the matrix is row_major or column major. Of course, as the function are declared ``inline``, at least the two dispatcher functions will not be implemented. Which means that this construction is not costly.

