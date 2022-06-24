.. $Id$

.. include:: ../replaces.txt

.. highlight:: none

.. _ud-install-mac:

How to install from sources on MacOS X
======================================

See `download and install <../download.html>`_ page for general requirements and the download of the last stable version of |gf|.

First, verify that you have installed the following components on your system:

  - Xcode
  - Xquartz
  - Homebrew

(Xquartz is not strictly necessary but more confortable).

Then, if you download the current git version

   $ brew install m4

   $ brew install automake

   $ brew install libtool

For the sequential mumps,

   $ brew tap brewsci/num

   $ brew install brewsci-mumps --without-mpi

For the parallel one, just forget --without-mpi and install also mpi and metis.

For Qhull

   $ brew install qhull

For Python

   $ pip install numpy

   $ pip install scipy

   $ pip install matplotlib


Download sources
----------------

There are two ways to get |gf|, either as a compressed package for the stable
release (file getfem-|version|.tar.gz downloadable on the page `download and install <../download.html>`_)  or via anonymous git access (current unstable version under development).

You can get the last stable version directly with

 * download package:

     $ curl -# "http://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-|version|.tar.gz" -o "getfem-|version|.tar.gz"

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
  is build by default and for python 3 version.

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

Note that there are other options to the configure script. A
``./configure --help`` will list them.

Octave interface
^^^^^^^^^^^^^^^^

The compilation of the Octave interface is performed with the ``--enable-octave`` option of the ``configure`` script.

First, you need ``octave`` and ``mkoctfile`` commands accessible from your shell prompt (for instance invoking ``brew install octave``).


The last step is to add the path to the toolbox in the octave path:

* you can put ``addpath('toolbox_dir', '-begin')`` to your ``$HOME/.octaverc`` file
* you can simply use the ``addpath`` command in the octave command line. 


Matlab interface
^^^^^^^^^^^^^^^^

The compilation of the Matlab interface (with the ``--enable-matlab`` option of the ``configure`` script) may fail due to a bad configuration of the Matlab compiler `mex`.

First, you need ``matlab`` and ``mex`` commands accessible from your shell prompt. If not, add ``Applications/MATLAB_RXXXX.app/bin`` on your path (for instance with ``export PATH=$PATH:Applications/MATLAB_RXXXX.app/bin`` if your shell is ``bash`` and for ``XXXX`` your Matlab installed version. Alternatively, you can make symbolic links to ``matlab`` and ``mex`` executable in ``/usr/local/bin`` thanks to the command ``sudo ln -s Applications/MATLAB_RXXXX.app/bin/matlab matlab`` and ``sudo ln -s Applications/MATLAB_RXXXX.app/bin/mex mex``.

Then, you will probably have to run

    $ mex -setup

To produce the correct ``mexopts.sh`` file in the ``.matlab/`` directory of your home directory. If it still does not work, then you can try to modify the ``.matlab/mexopts.sh`` or replace it. Some ``mexopts.sh`` specially adapted to macOS X/Xcode are available on the internet (See for instance here for `MATLAB_R2015 <https://gist.github.com/varunagrawal/811e05ee4ca0f6a9952d>`_).




.. caution::

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|).

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

You can then try to execute one of the demo program in ``getfem_root_directory/interface/tests/matlab``.


A few precompiled versions of the Matlab interface are available on the `download and install <../download.html>`_ page of |gf|.









Scilab interface
^^^^^^^^^^^^^^^^

The installation of the |sci| |gf| toolbox can be somewhat tricky, since it combines a C++ compiler, libraries and |sci| interaction. In case of troubles with a
non-GNU compiler, gcc/g++ (>= 4.8) should be a safe solution.


.. caution::

   * The minimal |sci| release is the 5.2.2.

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|).

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

