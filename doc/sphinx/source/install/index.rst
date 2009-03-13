.. include:: ../replaces.txt

.. highlightlang:: none

.. _install-index:

***********************
  Installing GetFEM++
***********************

In a debian/ubuntu system
=========================

If you have a problem installing the packages, please report it as a bug.

Edit ``/etc/apt/sources.list`` and add the following lines::

  deb http://apt.dim.uchile.cl distro main
  deb-src http://apt.dim.uchile.cl distro main

where

distro = `debian` xor `ubuntu`,

then, do a ``aptitude update`` and:

* ``aptitude install libgetfem++-dev`` for install GetFEM++ kernel library.
* ``aptitude install python-getfem`` for install Python interface.
* There isn't a package to install Matlab interface. See below...

In a general unix based systems
===============================

Since we use standard |gnu| tools, the installation of the |py| |gf| is
somewhat standard.

Requirements
------------

It requires the python developpement files (python.h etc.) to be available
(package `python-all-dev` in debian distribution), and also the numpy
package to be installed (package `python-numpy` in debian distribution).
In case of troubles with a non-GNU compiler, gcc/g++ (>= 3.0) should be a
safe solution (package `build-essential` in debian distribution).

If you want to:

* mesh generation, it requires the package qhull installed on your system
  (package `libqhull-dev` in debian distribution).

* build binaries from svn to get the latest changes, improvements, bugfixes,
  new bugs, etc. It requires an svn client, automake, and libtool.

Download sources
----------------

There are two ways to get |gf|, either as a compressed package (stable
release) or via anonymous svn access (unstable releases).

The latest stable release of |gf| is `getfem++-3.1.tar.gz
<http://download.gna.org/getfem/stable/getfem++-3.1.tar.gz>`_.

 * download package::

     wget http://download.gna.org/getfem/stable/getfem++-3.1.tar.gz

 * unpack::

     tar xzf getfem++-3.1.tar.gz

 * and go to the root directory of getfem::
     
     cd getfem++-3.1

The latest unstable releases is:

 * checkout over SVN protocol (TCP 3690)::

     svn co svn://svn.gna.org/svn/getfem/trunk getfem

 * or checkout over HTTP protocol (TCP 80)::

     svn co http://svn.gna.org/svn/getfem/trunk getfem

 * go to the root directory of getfem::

     cd getfem/getfem++
 
 * and run ``autogen.sh`` script::

     bash autogen.sh

Compilling
----------

Configure with::

  ./configure --enable-python=yes

If you want to use a specific **BLAS** library, you may have to supply the
necessary link flags and libs to the configure script, for example with::

  ./configure --enable-python=yes BLAS_LIBS="-L/usr/lib/sse2/atlas -lblas"

More specific instruccions can be found in the README\* files of the
distribution.

.. warning::

   * you should not use a different compiler than the one that was used
     for |gf|.
   * you should have built the |gf| static library (i.e. do not use
     ``./configure --disable-static`` when building |gf|).
   * On linux/x86_64 platforms, a mandatory option when building |gf| (and
     any static library linked to them) is the ``--with-pic`` option of
     their ``./configure`` script.

Then start the compilation with::

  make

and finally install with::

  make install
