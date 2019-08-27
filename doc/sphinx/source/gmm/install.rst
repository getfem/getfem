.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _gmm-install:

Installation
============

Since we use standard GNU tools, the installation of the |gmm| library is somewhat standard.

Note that if you use |gf|, you do not have to install |gmm| since |gf| is provided with its own version of |gmm|.

Moreover, as |gmm| is a template library, no compilation is needed to install it. If the |gmm|  archive is on your current directory you can unpack it and enter inside the directory of the distribution  with the commands::

  gunzip -c gmm-x.xx.tar.gz | tar xvf -
  cd  gmm-x.xx

Then you you have to run the configure script just typing::

  ./configure

or if you want to set the prefix directory where to install the library you can use the ``--prefix`` option (the default prefix directory is ``/usr/local``)::

  ./configure --prefix=\textit{dest_dir}

then start the installation with::

  make install

You can also check if your configuration is correct with::

  make check

which compiles random tests.

If you want to use a different compiler than the one chosen
automatically by the ``./configure`` script, just specify its
name on the command line::

  ./configure CXX=mycompiler

More specific instructions can be found in the ``README*`` files of
the distribution.

Now, to use |gmm| in you programs, the simpler manner is to include the file ``gmm/gmm.h`` which includes all the template library. If the compilation time is too important, the minimum to be included is contained is the file ``gmm/gmm\_kernel.h`` (vectors and matrix types, blas, sub vector and sub matrices).

DO NOT FORGET to catch errors messages. See the corresponding section.
