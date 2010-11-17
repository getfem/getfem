.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _sci-install:

Installation 
============

The installation of the |sci| |gf| toolbox can be somewhat tricky, since it combines a
C++ compiler, libraries and |sci| interaction. In case of troubles with a
non-GNU compiler, gcc/g++ (>= 4.1) should be a safe solution.

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

  - Go to the scilab getfem++ interface install directory (interface/src/scilab if the installation is not done)
 
  - launch scilab

  - load the getfem++ toolbox with:
    ``exec loader.sce;``

  - You can try to launch a demo with:
    ``cd demos;``
    ``exec demo_static_contact.sce;``

