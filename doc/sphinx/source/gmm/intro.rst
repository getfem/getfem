.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++


.. _gmm-intro:

Introduction
============

|gmm| provides some basic types of sparse and dense matrices and vectors. It provides some generic operations on them (copy, addition, multiplication, sub-vector and sub-matrices, solvers ... ). The syntax of |gmm| is very close to MTL and ITL (see http://www.mtl4.org/). Especially, the code for most of the iterative solvers has been imported from ITL. The performance of |gmm| is also close to the one of MTL, sometimes better. The difference is that basically |gmm| has been written to be able to interface other libraries and gives an access to sub matrices and sub vectors in all cases. Also some optimizations has been made for matrix-matrix multiplication, usage of reference has been somewhat cleared, const qualifier usage is clarified, and we hope it is somewhat easier to use.


.. include:: ../license.txt
