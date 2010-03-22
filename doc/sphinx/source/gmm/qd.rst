.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-qd:


How to use |gmm| with QD type (double-double and quad-double)
===============================================================

The QD library (see http://www.cs.berkeley.edu/\verb\~\yozo or http://www.nersc.gov/\verb\~\dhb/mpdist/mpdist.html) is an efficient library for double-double (32 decimal digits) and quad-double (approx. 64 decimal digits). Once you installed this library on your system you have to link your program with QD library (with -lqd). In your program, include the header files of QD with::

  #include <qd/dd.h>
  #include <qd/qd.h>
  #include <qd/fpu.h>


Then the two type ``dd_real`` and ``qd_real`` will be usable with |gmm|. You will also be able to use ``std::complex<dd_real>`` and ``std::complex<qdreal>``

IMPORTANT : do not forget to initialize QD before using it with the following call::

  unsigned int old_cw;
  fpu_fix_start(&old_cw);

This disables the 80 bits precision of x86 processors which conflicts with QD. Once you finished to use QD you can reactivate it with::

  fpu_fix_end(&old_cw);

(see the QD documentation for more details).

