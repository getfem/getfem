.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-install:

How to install
==============

Since we use standard |gnu| tools, the installation of the |gf| library is
somewhat standard.

Requirements
------------

* If you want to build binaries from svn to get the latest changes,
  improvements, bugfixes, new bugs, etc. It requires an svn client,
  automake, and libtool.

* If you want to build |py| |gfi|, it requires the python
  developpement files (python.h etc.) to be available (package
  ``python-all-dev`` in debian distribution), and also the |np| and |sp|
  packages to be installed (package ``python-numpy`` and ``python-scipy``
  in debian
  distribution). In case of troubles with a non-|gnu| compiler,
  gcc/g++ (>= 4.1) should be a safe solution (package
  ``build-essential`` in debian distribution).

* If you want mesh generation, it requires the package qhull
  installed on your system (package ``libqhull-dev`` in debian
  distribution).

* If you want use mathematical parser capabilities, it requieres
  the package muParser installed on your system (package
  ``libmuparser-dev`` in debian distribution).

Download sources
----------------

There are two ways to get |gf|, either as a compressed package (stable
release) or via anonymous svn access (unstable releases).

The latest stable release of |gf| is `getfem-4.0.0.tar.gz
<http://download.gna.org/getfem/stable/getfem-4.0.0.tar.gz>`_.

 * download package::

     $ wget http://download.gna.org/getfem/stable/getfem-4.0.0.tar.gz

 * unpack::

     $ tar xzf getfem-4.0.0.tar.gz

 * and go to the root directory of |gf|::

     $ cd getfem-4.0.0/

The latest unstable releases is:

 * checkout over SVN protocol (TCP 3690)::

     $ svn co svn://svn.gna.org/svn/getfem/trunk/getfem++ getfem

 * or checkout over HTTP protocol (TCP 80)::

     $ svn co http://svn.gna.org/svn/getfem/trunk/getfem++ getfem

 * go to the root directory of |gf|::

     $ cd getfem/

 * and run ``autogen.sh`` script::

     $ bash autogen.sh


Compilling
----------

Configure with::

  $ ./configure

then start the compilation with::

  $ make

and finally install with::

  $ make install

Configure Options
^^^^^^^^^^^^^^^^^

* If you want to use a different compiler than the one chosen
  automatically by the ``./configure`` script, just specify its
  name on the command line::

    $ ./configure CXX=mycompiler

* If you want to use a specific **BLAS** library, you may have to
  supply the necessary link flags and libs to the configure script
  with::

    $ ./configure BLAS_LIBS="-L/path/to/lib -lfoo -lbar ....etc"

  for example::

    $ ./configure BLAS_LIBS="-L/usr/lib/sse2/atlas -lblas"

* If you want to set the prefix directory where to install the library
  you can use the ``--prefix`` option (the default prefix directory is
  ``/usr/local``)::

    $ ./configure --prefix=my_dest_dir

* If you want build |py| |gfi|, use ``--enable-python=yes`` option.

  .. warning::

     * you should not use a different compiler than the one that was used
       for |gf|.
     * you should have built the |gf| static library (i.e. do not use
       ``./configure --disable-static`` when building |gf|).
     * On linux/x86_64 platforms, a mandatory option when building |gf|
       (and any static library linked to them) is the ``--with-pic``
       option of their ``./configure`` script.

Note that there are other options to the configure script. A
``./configure --help`` will list them. Most important ones are
``--enable-matlab``, ``--enable-python`` and ``--enable-scilab``
to build the |gfi|.

More specific instructions can be found in the README\* files of the
distribution.
