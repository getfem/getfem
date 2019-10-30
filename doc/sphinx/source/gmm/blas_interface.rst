.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: none

.. _gmm-lapack:


Interface with BLAS, LAPACK or ATLAS
======================================

For better performance on dense matrices, it is possible to interface some operations of the type ``gmm::dense_matrix<T>`` with ``LAPACK`` (http://www.netlib.org/lapack/) or ``ATLAS`` (http://math-atlas.sourceforge.net/), for ``T = float, double, std::complex<float> or std::complex<double>``. In fact, concerning ``ATLAS`` no specific interface has been made untill now, so the fortran interface of ``ATLAS`` should be used.

to use this interface you have first to define ``GMM_USES_LAPACK`` before including |gmm| \ files::

  #define GMM_USES_LAPACK
  #include <gmm/gmm.h>

  ... your code


or specify -DGMM_USES_LAPACK on the command line of your compiler. Of course, you have also to link ``LAPACK`` or ``ATLAS`` libraries. For example on a standard linux configuration and g++ compiler the adding libraries to link ``LAPACK`` are::

  g++ ...  -llapack -lblas -lgfortanbegin -lgfortran

and to link  ``ATLAS``::

  g++ ... /usr/lib/atlas/liblapack.a /usr/lib/atlas/libblas.a -latlas -lgfortranbegin -lgfortran

The library ``libgfortranbegin`` and ``libgfortran`` are specific to g++ compiler and may vary for other compilers.


Ask your system administrator if this configuration does not work.

The following operations are interfaced::

  vect_norm2(std::vector<T>)

  vect_sp(std::vector<T>, std::vector<T>)
  vect_sp(scaled(std::vector<T>), std::vector<T>)
  vect_sp(std::vector<T>, scaled(std::vector<T>))
  vect_sp(scaled(std::vector<T>), scaled(std::vector<T>))

  vect_hp(std::vector<T>, std::vector<T>)
  vect_hp(scaled(std::vector<T>), std::vector<T>)
  vect_hp(std::vector<T>, scaled(std::vector<T>))
  vect_hp(scaled(std::vector<T>), scaled(std::vector<T>))

  add(std::vector<T>, std::vector<T>)
  add(scaled(std::vector<T>, a), std::vector<T>)

  mult(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)
  mult(transposed(dense_matrix<T>), dense_matrix<T>, dense_matrix<T>)
  mult(dense_matrix<T>, transposed(dense_matrix<T>), dense_matrix<T>)
  mult(transposed(dense_matrix<T>), transposed(dense_matrix<T>),
       dense_matrix<T>)
  mult(conjugated(dense_matrix<T>), dense_matrix<T>, dense_matrix<T>)
  mult(dense_matrix<T>, conjugated(dense_matrix<T>), dense_matrix<T>)
  mult(conjugated(dense_matrix<T>), conjugated(dense_matrix<T>),
       dense_matrix<T>)

  mult(dense_matrix<T>, std::vector<T>, std::vector<T>)
  mult(transposed(dense_matrix<T>), std::vector<T>, std::vector<T>)
  mult(conjugated(dense_matrix<T>), std::vector<T>, std::vector<T>)
  mult(dense_matrix<T>, scaled(std::vector<T>), std::vector<T>)
  mult(transposed(dense_matrix<T>), scaled(std::vector<T>),
       std::vector<T>)
  mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),
       std::vector<T>)

  mult_add(dense_matrix<T>, std::vector<T>, std::vector<T>)
  mult_add(transposed(dense_matrix<T>), std::vector<T>, std::vector<T>)
  mult_add(conjugated(dense_matrix<T>), std::vector<T>, std::vector<T>)
  mult_add(dense_matrix<T>, scaled(std::vector<T>), std::vector<T>)
  mult_add(transposed(dense_matrix<T>), scaled(std::vector<T>),
           std::vector<T>)
  mult_add(conjugated(dense_matrix<T>), scaled(std::vector<T>),
           std::vector<T>)

  mult(dense_matrix<T>, std::vector<T>, std::vector<T>, std::vector<T>)
  mult(transposed(dense_matrix<T>), std::vector<T>, std::vector<T>,
       std::vector<T>)
  mult(conjugated(dense_matrix<T>), std::vector<T>, std::vector<T>,
       std::vector<T>)
  mult(dense_matrix<T>, scaled(std::vector<T>), std::vector<T>,
       std::vector<T>)
  mult(transposed(dense_matrix<T>), scaled(std::vector<T>),
       std::vector<T>, std::vector<T>)
  mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),
       std::vector<T>, std::vector<T>)
  mult(dense_matrix<T>, std::vector<T>, scaled(std::vector<T>),
       std::vector<T>)
  mult(transposed(dense_matrix<T>), std::vector<T>,
       scaled(std::vector<T>), std::vector<T>)
  mult(conjugated(dense_matrix<T>), std::vector<T>,
       scaled(std::vector<T>), std::vector<T>)
  mult(dense_matrix<T>, scaled(std::vector<T>), scaled(std::vector<T>),
    std::vector<T>)
  mult(transposed(dense_matrix<T>), scaled(std::vector<T>),
       scaled(std::vector<T>), std::vector<T>)
  mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),
       scaled(std::vector<T>), std::vector<T>)

  lower_tri_solve(dense_matrix<T>, std::vector<T>, k, b)
  upper_tri_solve(dense_matrix<T>, std::vector<T>, k, b)
  lower_tri_solve(transposed(dense_matrix<T>), std::vector<T>, k, b)
  upper_tri_solve(transposed(dense_matrix<T>), std::vector<T>, k, b)
  lower_tri_solve(conjugated(dense_matrix<T>), std::vector<T>, k, b)
  upper_tri_solve(conjugated(dense_matrix<T>), std::vector<T>, k, b)

  lu_factor(dense_matrix<T>, std::vector<int>)
  lu_solve(dense_matrix<T>, std::vector<T>, std::vector<T>)
  lu_solve(dense_matrix<T>, std::vector<int>, std::vector<T>,
           std::vector<T>)
  lu_solve_transposed(dense_matrix<T>, std::vector<int>, std::vector<T>,
           std::vector<T>)
  lu_inverse(dense_matrix<T>)
  lu_inverse(dense_matrix<T>, std::vector<int>, dense_matrix<T>)

  qr_factor(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)

  implicit_qr_algorithm(dense_matrix<T>, std::vector<T>)
  implicit_qr_algorithm(dense_matrix<T>, std::vector<T>,
                        dense_matrix<T>)
  implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >)
  implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >,
                        dense_matrix<T>)

Of course, it is not difficult to interface another operation if needed.

The following interface does not correspond to an algorithm existing in |gmm|:

The interface to ``gesvd`` (singular value decomposition)::

   svd(dense_matrix<T> &X, dense_matrix<T> &U,
       dense_matrix<T> &Vt, std::vector<T> sigma);
   svd(dense_matrix<std::complex<T> > &X, dense_matrix<std::complex<T> > &U,
       dense_matrix<std::complex<T> > &Vt, std::vector<T> sigma);
