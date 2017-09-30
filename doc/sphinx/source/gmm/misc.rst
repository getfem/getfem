.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-misc:

Miscellaneous methods
========================================


::

  gmm::vect_size(V); // gives the size of the vector V.

::

  gmm::resize(V, n); // Change the size of the vector V.
                     // Preserve the min(n, vect_size(V)) first components.
                     // Do not work for references.
  gmm::resize(M, m, n); // Change the dimensions of matrix M.
                        // Preserve the
                        // min(m, mat_nrows(M)) x min(n, mat_ncols(M))
                        // first components. Do not work for references.
  gmm::reshape(M, m, n);  // returns the m-by-n matrix whose elements
                          // are taken columnwise from M.
                          // An error results if M does not have m*m
                          // elements. Works only with dense_matrix<T> for
                          // the moment.


::

  gmm::nnz(V); // gives the number of stored components of the vector V.
  gmm::nnz(M); // gives the total number of stored components of the matrix M.



::

  gmm::mat_nrows(M) // gives the number of rows of a matrix M.
  gmm::mat_ncols(M) // gives the number of columns of a matrix M.


::

  gmm::write(o, V); // print the vector V to the output stream o.
  gmm::write(o, M); // print the matrix M to the output stream o.

Most of the time it is more convenient to use::

  std::cout << gmm::vref(V) << std::endl;
  std::cout << M << std::endl;


::

  gmm::clear(V); // set to zero all the components of the vector V;
  gmm::clear(M); // set to zero all the components of the matrix M;


::

  gmm::clean(V, 1E-10); // set to zero all the components of the vector V
                        // whose modulus is less or equal to 1E-10
  gmm::clean(M, 1E-10); // idem for a matrix M.


::

  gmm::fill_random(V); // fill a dense vector V with random number
                       //  between -1 and 1
  gmm::fill_random(V, cfill); // fill a dense or sparse vector with random
                       // numbers. cfill should be between 0.0 qnd 1.0 and
                       // represent the ratio of filled components.
  gmm::fill_random(M); // fill a dense matrix M with random number
  gmm::fill_random(M, cfill); // fill a dense or sparse matrix M with random
                       // numbers.

