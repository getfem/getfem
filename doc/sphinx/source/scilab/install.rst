.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _sci-install:

Installation 
============

The installation of the |gfi| toolbox can be somewhat tricky, since it combines a
C++ compiler, libraries and |sci| interaction... In case of troubles with a
non-GNU compiler, gcc/g++ (>= 4.1) should be a safe solution.

.. caution::

   * you should not use a different compiler than the one that was used for the |gf| library.

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

   * you should have use the --enable-scilab option to configure the |gf| sources (i.e. ./configure --enable-matlab ...)

You may also use ``--with-toolbox-dir=toolbox_dir`` to change the default toolbox installation directory (``gfdest_dir/getfem_toolbox``). Use ``./configure --help`` for more options.


With this, since the Scilab interface is contained into the |gf| sources (in the directory interface/src) you can compile both the |gf| library and the Scilab interface by ::

  make

An optional step is ``make check`` in order to check the scilab interface ... and install it ( ... )::

  make install

If you want to use a different compiler than the one chosen automatically by the ``./configure`` script, just specify its name on the command line: ``./configure CXX=mycompiler``.

When the library is installed, 

completer les instructions ...



...





A very classical problem at this step is the incompatibility of the C and C++ libraries used by Scilab. Scilab is distributed with its own libc and libstdc++ libraries. An error message of the following type occurs when one tries to use a command of the interface::

  /usr/local/matlab14-SP3/bin/glnxa64/../../sys/os/??/libgcc_s.so.1:
  version `GCC_?.?' not found (required by .../gf_matlab.mex??).

In order to fix this problem one has to enforce Scilab to load the C and C++ libraries of the system. There is two possibilities to do this. The most radical is to delete the C and C++ libraries distributed along with Matlab (if you have administrator privileges ...!) for instance with::

  rm `a completer`/libgcc_s.so.1
  rm `a completer`/libstdc++_s.so.6

The second possibility is to set the variable LDPRELOAD before launching Matlab for instance with (depending on the system)::

  LD_PRELOAD=/usr/lib/libgcc_s.so:/usr/lib/libstdc++.so.6 scilab

More specific instructions can be found in the ``README*`` files of the
distribution.
