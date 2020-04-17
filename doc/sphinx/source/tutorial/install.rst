.. $Id: install.rst 4738 2014-07-27 12:25:54Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-install-tut:

How to install
==============

Since |gf| is developed on linux (Ubuntu), the installation is simpler on linux, especially on Debian-based distributions (Debian/Ubuntu/Mint). However, |gf| can be installed also on other linux distributions, on Mac os X and Windows. In order to compile |gf| from sources, you need a recent C++ complier (supporting C++ 14 standard) and a recent version of python.

The main dependences of |Gf| on other libraries are

* git client, automake, autoconf and libtool if you  want to build binaries
  from git version to get the latest changes.

* Python development files (Python.h etc.) and also the |np| and |sp| packages if
  you want to build the python interface.

* sequential MUMPS package (direct solver for sparse matrices) if you want to use it instead of the SuperLU version distributed along with |gf|.

* Parallel MUMPS, METIS and MPI4PY packages if you want to use the MPI parallelized version of |gf|.

* qhull package for mesh generation and fictitious domain applications

* BLAS and LAPACK packages

|gf| C++ library can be build on its own or together with the Python, Scilab and/or Matlab interface.

You can also install the stable release of Getfem on linux distributions using the corresponding package management system.


More specific information on how to build Getfem C++ library can be found on the `download and install <../download.html>`_ page.
