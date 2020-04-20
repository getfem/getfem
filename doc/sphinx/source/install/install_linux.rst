.. $Id$

.. include:: ../replaces.txt

.. highlight:: none

.. _ud-install-linux:

How to install from sources on Linux
====================================

See `download and install <../download.html>`_ page for general requirements and the download of the last stable version of |gf|.



Download sources
----------------

There are two ways to get |gf|, either as a compressed package for the stable
release (file getfem-|version|.tar.gz downloadable on the page `download and install <../download.html>`_)  or via anonymous git access (current unstable version under development).

You can get the last stable version directly with

 * download package:

     $ wget http://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-|version|.tar.gz

 * unpack:

     $ tar xzf getfem-|version|.tar.gz

 * and go to the root directory of |gf|:

     $ cd getfem-|version|/

The current git version is:

* checkout over GIT protocol::

     $ git clone https://git.savannah.nongnu.org/git/getfem.git

* go to the root directory of |gf|::

     $ cd getfem

* and run ``autogen.sh`` script (you need m4, automake and libtool) ::

     $ bash autogen.sh


Compiling
----------

Configure with::

  $ ./configure

then start the compilation with::

  $ make

and finally install with::

  $ make install

You can find some additional help on how to build the Matlab interface on Ubuntu on the page of `Mirko Windhoff <http://windhoff.net/wiki/how_to/build_getfem_matlab_toolbox_on_ubuntu_linux>`_.



Configure Options
^^^^^^^^^^^^^^^^^

* If you want to use a different compiler than the one chosen
  automatically by the ``./configure`` script, just specify its
  name on the command line::

    $ ./configure CXX=mycompiler

* If you want to build one of the interfaces, use::

    $ ./configure ``--enable-python``
    $ ./configure ``--enable-scilab``
    $ ./configure ``--enable-matlab``

  depending on the interface you want to build. Note that the python interface
  is build by default.

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

* By default, the python interface is built and for python 3 version. You can disable the built of the python interface with::

    $ ./configure --disable-python

Note that there are other options to the configure script. A
``./configure --help`` will list them.


.. warning::

     * On linux/x86_64 platforms, a mandatory option when building |gf|
       (and any static library linked to them) is the ``--with-pic``
       option of their ``./configure`` script.



Scilab interface
^^^^^^^^^^^^^^^^

The installation of the |sci| |gf| toolbox can be somewhat tricky, since it combines a C++ compiler, libraries and |sci| interaction. In case of troubles with a
non-GNU compiler, gcc/g++ (>= 8.0) should be a safe solution.


.. caution::

   * The minimal |sci| release is the 5.2.2.

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| toolbox (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

   * you should have use the ``--enable-scilab`` option to configure the |gf| sources (i.e. ``./configure --enable-scilab`` ...)

You may also use ``--with-scilab-toolbox-dir=toolbox_dir`` to change the default toolbox installation directory (``gfdest_dir/getfem_toolbox``). Use ``./configure --help`` for more options.


With this, since the Scilab interface is contained into the |gf| sources (in the directory interface/src) you can compile both the |gf| library and the Scilab interface by ::

  make

Optionally, you can install it with ::

  make install

If you want to use a different compiler than the one chosen automatically by the ``./configure`` script, just specify its name on the command line: ``./configure CXX=mycompiler``.


Once getfem is compiled:

  - Go to the Scilab GetFEM interface install directory (interface/src/scilab if the installation is not done)

  - launch Scilab

  - load the GetFEM toolbox with:
    ``exec loader.sce;``

  - You can try to launch a demo with:
    ``cd demos;``
    ``exec demo_static_contact.sce;``

Octave interface
^^^^^^^^^^^^^^^^

You have first to install |octv| with the developpement package

.. caution::

   * You have first to install |octv|, minimal release 4.1.1 with the developpement package such
     that the command ``mkoctfile`` is available (liboctave-dev package on Debian, for instance)

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| toolbox (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

   * you should have use the ``--enable-octave`` option to configure the |gf| sources (i.e. ``./configure --enable-octave`` ...)



The last step is to add the path to the toolbox in the octave path:

* you can put ``addpath('toolbox_dir', '-begin')`` to your ``$HOME/.octaverc`` file
* you can simply use the ``addpath`` command in the octave command line. 

Matlab interface
^^^^^^^^^^^^^^^^

The installation of the |gfi| toolbox can be somewhat tricky, since it combines a
C++ compiler, libraries and |Mlab| interaction... In case of troubles with a
non-GNU compiler, gcc/g++ (>= 8.0) should be a safe solution.


.. caution::

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

   * you should have use the --enable-matlab option to configure the |gf| sources (i.e. ./configure --enable-matlab ...)

You may also use ``--with-matlab-toolbox-dir=toolbox_dir`` to change the default toolbox installation directory (``gfdest_dir/getfem_toolbox``). Use ``./configure --help`` for more options.


With this, since the Matlab interface is contained into the |gf| sources (in the directory interface/src) you can compile both the |gf| library and the Matlab interface by ::

  make

An optional step is ``make check`` in order to check the matlab interface (this
sets some environment variables and runs the ``check_all.m`` script which is the ``tests/matlab`` directory of the distribution) and install it (the libraries
will be copied in ``gfdest_dir/lib``, while the MEX-File and M-Files will be
copied in ``toolbox_dir``)::

  make install

If you want to use a different compiler than the one chosen automatically by the ``./configure`` script, just specify its name on the command line: ``./configure CXX=mycompiler``.

When the library is installed, you may have to set the ``LD_LIBRARY_PATH``
environment variable to the directory containing the ``libgetfem.so`` and
``libgetfemint.so``, which is ``gfdest_dir/lib``::

  export LD_LIBRARY_PATH=gfdest_dir/lib # if you use bash

The last step is to add the path to the toolbox in the matlab path:

* you can set the environment variable ``MATLABPATH`` to ``toolbox_dir``
  (``export MATLABPATH=toolbox_dir`` for example).
* you can put ``addpath('toolbox_dir')`` to your ``$HOME/matlab/startup.m``

A very classical problem at this step is the incompatibility of the C and C++ libraries used by Matlab. Matlab is distributed with its own libc and libstdc++ libraries. An error message of the following type occurs when one tries to use a command of the interface::

  /usr/local/matlab14-SP3/bin/glnxa64/../../sys/os/??/libgcc_s.so.1:
  version `GCC_?.?' not found (required by .../gf_matlab.mex??).

In order to fix this problem one has to enforce Matlab to load the C and C++ libraries of the system. There is two possibilities to do this. The most radical is to delete the C and C++ libraries distributed along with Matlab (if you have administrator privileges ...!) for instance with::

  mv /usr/local/matlab14-SP3/sys/os/??/libgcc_s.so.1 libgcc_s.so.1_old
  mv /usr/local/matlab14-SP3/sys/os/??/libstdc++_s.so.6 libstdc++_s.so.6_old
  mv /usr/local/matlab14-SP3/sys/os/??/libgfortran.so.3 libgfortran.so.3_old

The second possibility is to set the variable LDPRELOAD before launching Matlab for instance with (depending on the system)::

  LD_PRELOAD=/usr/lib/libgcc_s.so:/usr/lib/libstdc++.so.6 matlab

More specific instructions can be found in the ``README*`` files of the
distribution.

A few precompiled versions of the Matlab interface are available on the `download and install <../download.html>`_ page of |gf|.


A second problem arising with recent distribution of Matlab (2016a), is the incompatibility of some libraries with ILP64 version of MKL loaded by MATLAB which uses 64 bits integers instead of 32 bits ones contrarily to most system blas/lapack libraries. New releases of |gf| are compatible with both 64 bits and 32 bits integer blas/lapack libraries. However, for instance, Mumps should be recompiled in a 64 bit integer version to be compatible with MKL ILP64. Mumps version on the system is the 32 bits integer version. If problem of this kind are encountered, you can try to force Matlab to load 32 bit blas and lapack libraries with::

  LD_PRELOAD=/usr/lib/libblas.so:/usr/lib/liblapack.so matlab


