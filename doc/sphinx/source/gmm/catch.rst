.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-catch:


Catch errors
============================


Errors used in |gmm| are defined in the file ``gmm/gmm\_except.h``. In order to make easier  the error catching all errors derive from the type ``std::logic\_error`` defined in the file `` stdexcept`` of the S.T.L.

A standard procedure, ``GMM\_STANDARD\_CATCH\_ERROR``, is defined in ``gmm/gmm\_except.h``. This procedure catches all errors and print the error message when an error occurs. It can be used in the main procedure of the program as follows::

  int main(void) \{ 
    try \{ 
      ... main program ... 
        \} 
     GMM\_STANDARD\_CATCH\_ERROR;
  \}


It is highly recommended to catch the errors at least in the main function, because if you do not so, you will not be able to see error messages.
