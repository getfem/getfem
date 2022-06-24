.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc_gmm:

Gmm library
-----------

Description
^^^^^^^^^^^

|gmm| is a template linear algebra library which was originally designed to make an
interface between the need in linear algebra procedures of |gf| and existing free
linear algebra libraries (MTL, Superlu, Blas, Lapack originally). It rapidly
evolves to an independent self-consistent library with its own vector and matrix
types. It is now used as a base linear algebra library by several other projects.

However, it preserves the characteristic to be a potential interface for more
specific packages. Any vector or matrix type having the minimum of compatibility
can be used by generic algorithms of |gmm| writing a ``linalg_traits`` structure.

A |gmm| standalone version is distributed since release 1.5 of |gf|. It is
however developed inside the |gf| project even though since release 3.0 it is
completely independent of any |gf| file.

In addition to the linear algebra procedures, it furnishes also the following
utilities to |gf|.

* Fix some eventual compatibility problems in :file:`gmm_std.h`.

* Error, warning and trace management in :file:`gmm_except.h`.

* Some extended math definitions in :file:`gmm_def.h`.

See :ref:`gmm` documenation for more details.

Files
^^^^^

All files in src/gmm


State
^^^^^

For the moment, |gmm| cover the needs of |gf| concerning the basic linear algebra
procedures.


Perspectives
^^^^^^^^^^^^

There is potentially several points to be improved in |gmm| (partial
introduction of expression template for some base types of matrix and vectors,
think about the way to represent in a more coherent manner sparse sub-vectors
and sub-matrices, introduction of C++ concepts, etc.). However, since |gmm|
globally cover the needs of |gf| and since there exists some other project like
`Glas <http://glas.sourceforge.net/>`_ to build a reference C++ library for
linear algebra, a global change seem to be unnecessary. This part
is considered to be stabilized.

The current vocation of |gmm| is to continue to collect generic algorithms and
interfaces to some other packages (DIFFPACK for instance) in order to cover new needs of the whole
project. The library is now frequently used as a separate package and has also
the vocation to collect the contribution of any person who propose some
improvements, new algorithms or new interfaces.


