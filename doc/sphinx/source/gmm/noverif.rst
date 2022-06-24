.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _gmm-noverif:


How to disable verifications
============================


On some type of matrices such as ``gmm::dense_matrix`` some verification are made on the range of indices. This could deteriorate  the performance of your code but is satisfactory in the development stage. You can disable these verifications adding a ``-dNDEBUG`` to the compiler options.

