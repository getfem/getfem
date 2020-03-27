.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-export:


Input and output with Harwell-Boeing and Matrix Market formats
==============================================================

Including the file ``gmm/gmm_inoutput.h`` you will be able to load and save matrices with Harwell-Boeing and Matrix Market formats. Concerning the Harwell-Boeing format, only the type ``gmm::csc_matrix<double>`` and ``gmm::csc_matrix<std::complex<double> >`` has been interfaced, so you can execute::

  gmm::Harwell_Boeing_save("filename", A); // save the matrix A .
  gmm::Harwell_Boeing_load("filename", A); // load the matrix A.

If ``A`` is not a  ``gmm::csc_matrix<double>`` or a ``gmm::csc_matrix<std::complex<double> >`` a copy is made.

Concerning the Matrix Market format, it is possible to save a ``gmm::csc_matrix<double>`` or a  ``gmm::csc_matrix<std::complex<double> >`` and to load a ``gmm::row_matrix<VECT>`` or a ``gmm::col_matrix<VECT>``::

  gmm::MatrixMarket_save("filename", A); // save a csc_matrix.
  gmm::MatrixMarket_load("filename", A); // load a row_matrix or a col_matrix

