.. $Id: intro.rst 3327 2009-11-04 16:42:35Z lsaavedr $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _gmm-noverif:


How to disable verifications
============================


On some type of matrices such as ``gmm::dense_matrix`` some verification are made on the range of indices. This could deteriorate  the performance of your code but is satisfactory in the developpment stage. You can disable these verifications adding a ``-dNDEBUG`` to the compiler options.

