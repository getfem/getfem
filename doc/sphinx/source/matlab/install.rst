.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _mlab-install:

Installation 
============

The installation of the |gfi| toolbox can be somewhat tricky, since it combines a
C++ compiler, libraries and |Mlab| interaction... In case of troubles with a
non-GNU compiler, gcc/g++ (>= 4.1) should be a safe solution.

.. caution::

   * you should not use a different compiler than the one that was used for the |gf| library.

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

   * you should have use the --enable-matlab option to configure the |gf| sources (i.e. ./configure--enable-matlab ...)

You may also use ``--with-toolbox-dir=toolbox_dir`` to change the default toolbox installation directory (``gfdest_dir/getfem_toolbox``). Use ``./configure --help`` for more options.


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

  export LD_LIBRARY_PATH=gfdest_dir/lib # if you use ksh or bash

The last step is to add the path to the toolbox in the matlab path:

* you can set the environment variable ``MATLABPATH`` to ``toolbox_dir``
  (``export MATLABPATH=toolbox_dir`` for example).
* you can put ``addpath('toolbox_dir')`` to your ``$HOME/matlab/startup.m``

More specific instructions can be found in the ``README*`` files of the
distribution.
