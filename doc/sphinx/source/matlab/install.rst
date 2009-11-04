.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _mlab-install:

Installation 
============

The installation of the |gfi| toolbox can be somewhat tricky, since it combines a
C++ compiler, libraries and |Mlab| interaction... In case of troubles with a
non-GNU compiler, gcc/g++ (>= 3.0) should be a safe solution.

.. caution::

   * you should not use a different compiler than the one that was used for |gf|.

   * you should have built the |gf| static library (i.e. do not use ``./configure
     --disable-static`` when building |gf|). On linux/x86_64 platforms, a
     mandatory option when building |gf| and |gfi| (and any static library linked
     to them) is the ``--with-pic`` option of their ``./configure`` script.

Here we assume that |gf| was installed in the directory ``gfdest_dir`` (i.e.  you
ran ``./configure --prefix=gfdest_dir`` before compiling and installing |gf|, the
default value being ``/usr/local``).

Unpack the |gfi| archive and run the configure script::

  gzip -dc getfem-interface-2.0.tar.gz | tar xvf -+
  cd getfem-interface-2.0+

If you did install |gf|, then running ``./configure`` or ``./configure
--prefix=gfdest_dir`` should be sufficient.

Nevertheless, if ``gfdest_dir/bin`` is not in the ``PATH``, then you will have to
to provide the path to the ``getfem-config`` script with
``--with-getfem-config=gfdest_dir/bin/getfem-config``.

You may also use ``--with-toolbox-dir=toolbox_dir`` to change the default toolbox
installation directory (``gfdest_dir/getfem_toolbox``). Use ``./configure
--help`` for more options.

When the ``configure`` is done, you can compile the toolbox (use ``gmake`` if
your default ``make`` is not the GNU one)::

  make

An optional step is ``make check`` in order to check the matlab interface (this
sets some environment variables and runs the ``check_all.m`` script which is the
``tests/matlab`` directory of the distribution) and install it (the libraries
will be copied in ``gfdest_dir/lib``, while the MEX-File and M-Files will be
copied in ``toolbox_dir``)::

  make install

If you want to use a different compiler than the one chosen automatically by the
``./configure`` script, just specify its name on the command line: ``./configure
CXX=mycompiler``.

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
