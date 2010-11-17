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

  export LD_LIBRARY_PATH=gfdest_dir/lib # if you use ksh or bash

The last step is to add the path to the toolbox in the matlab path:

* you can set the environment variable ``MATLABPATH`` to ``toolbox_dir``
  (``export MATLABPATH=toolbox_dir`` for example).
* you can put ``addpath('toolbox_dir')`` to your ``$HOME/matlab/startup.m``

A very classical problem at this step is the incompatibility of the C and C++ libraries used by Matlab. Matlab is distributed with its own libc and libstdc++ libraries. An error message of the following type occurs when one tries to use a command of the interface::

  /usr/local/matlab14-SP3/bin/glnxa64/../../sys/os/??/libgcc_s.so.1:
  version `GCC_?.?' not found (required by .../gf_matlab.mex??).

In order to fix this problem one has to enforce Matlab to load the C and C++ libraries of the system. There is two possibilities to do this. The most radical is to delete the C and C++ libraries distributed along with Matlab (if you have administrator privileges ...!) for instance with::

  rm /usr/local/matlab14-SP3/sys/os/??/libgcc_s.so.1
  rm /usr/local/matlab14-SP3/sys/os/??/libstdc++_s.so.6

The second possibility is to set the variable LDPRELOAD before launching Matlab for instance with (depending on the system)::

  LD_PRELOAD=/usr/lib/libgcc_s.so:/usr/lib/libstdc++.so.6 matlab

More specific instructions can be found in the ``README*`` files of the
distribution.

In particular, instruction for the installation on Mac OS can be found here :ref:`mlab-install_mac`.

A few precompiled versions of the Matlab interface are available on the download page of |gf|.
