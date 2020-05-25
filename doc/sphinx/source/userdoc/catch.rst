.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-catch:

Catch errors
============

Errors used in |gf| are defined in the file :file:`gmm/gmm_except.h`. In order to
make easier the error catching all errors derive from the type
``std::logic_error`` defined in the file ``stdexcept`` of the S.T.L.

A standard procedure, ``GMM_STANDARD_CATCH_ERROR``, is defined in
:file:`gmm/gmm_except.h`. This procedure catches all errors and prints the error
message when an error occurs. It can be used in the main procedure of the program
as follows::

  int main(void) {
    try {
      ... main program ...
    } GMM_STANDARD_CATCH_ERROR;
  }
