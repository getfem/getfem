.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-install-windows:

How to install |gf| from sources on Windows
===========================================

See the `download and install <../download.html>`_ page for general
requirements and the download of the last stable version of |gf|.

Building |gf| in an environment hostile to open-source like Windows
is naturally more painful than on Linux or MacOS X. Several possibilities
exist to build |gf| with a commercial C++ compiler (see the ``mscv``
directory on the git version of |gf| for some MSCV project files).
However, we describe in the following how to build the sequential version of
|gf| with both Python and Matlab interfaces under
`Mingw <http://www.mingw.org/>`_
and `Msys <http://www.mingw.org/wiki/MSYSwhich>`_ which provide a
minimal GNU environment for Windows. A similar installation should
possible with `Cygwin <https://www.cygwin.com/>`_.

  - The first step is to install Mingw in 64 bits version
    (x86_64, posix, sjlj options,
    see `Mingw-64 web page <https://mingw-w64.org/doku.php/download>`_).
    For this you can usr an installer or simply download Mingw and copy it
    in the directory ``C:\mingw64``

  - Install also Msys 64 bits version the same way (either with an installer
    or directly copy Msys in the directory ``C:\msys``)

  - If it is not done by the installer, add to the Windows path variable the
    directories ``C:\mingw64\bin;C:\msys;C:\msys\bin``. To this aim, go in the
    start menu of Windows, search ``"system"`` and select
    ``"edit the system environement variables"`` then edit the path variable.

  - Open a Windows command console and enter ``msys``. This should open an Msys
    shell console. For who is familiar with sh/bash, you normally feel
    in a more civilized environment. Nothing to be too much happy,
    this environment
    is somehow minimalist (it aims to). Normally, Msys creates you a home
    directory (try the command ``pwd`` to see the Msys path and note that it
    corresponds to ``C:\msys\home\login`` windows path, where ``login`` is
    your home directory name). You have first to create a file in that path
    named ``.profile`` and containing at list the two following lines::

      export CPATH="/usr/local/include"
      export LIBRARY_PATH="$LIBRARY_PATH:/usr/local/lib"

    In theory, this should not be necessary, because ``/usr/local/include``
    and ``/usr/local/lib`` are default directories, but ok, in my config,
    it does not work without. Note that you can use any text editor to create
    this file, for instance ``vim`` of Msys (there is also good versions of
    emacs for windows or you can use the rudimentary blocnote of Windows).

  - You will have to compile and install BLAS and LAPACK versions (it would
    be possible to avoid that, but it is not a difficult step). Download a
    version of tar.gz source package of blas and lapack into you home Msys
    directory. Untar the two archives (for instance with
    ``tar xvzf blas-?.?.?.tar.gz``,  ``tar xvzf lapack-?.?.?.tar.gz``,
    ``?.?.?`` being the version numbers) and then
    enter into the BLAS source directory, compile and install BLAS with::

      $ cd blas-?.?.?
      $ make
      $ cp blas_LINUX.a /usr/local/lib/libblas.a
      $ cd ..

    Do the same for LAPACK with::
   
      $ cd lapack-?.?.?
      $ cp make.inc.example make.inc
      $ cp ../BLAS-?.?.?/blas_LINUX.a librefblas.a
      $ make
      $ cp liblapack.a /usr/local/lib
      $ cd ..

  - You now need an installation of QHULL library in order to have access to
    the meshing and Xfem tools of |gf|
    (see `Qhull <http://www.qhull.org>`_ and
    `Qhull install instructions <http://www.qhull.org/README.txt>`_). 
    Download the sources of Qhull in you Msys home (similarly as what you have
    done for BLAS and LAPACK), untar them and enter into the Qhull source
    directory with Msys. You can compile and install Qhull simply with::

      $ make SO=dll
      $ make install
  
  - Similarly, we will compile and install now a sequential version of
    the sparse linear solver `MUMPS <http://mumps.enseeiht.fr/>`_.
    Again, it is not strictly necessary since a version of
    `Superlu <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`_ is
    included with the |gf| sources, but, especially in 3D, MUMPS can solve
    really larger systems. So, download a MUMPS source package on your
    Msys home directory, untar it and make the following steps::

      $ cd MUMPS_?.?.?
      $ cp Make.inc/Makefile.inc.generic.SEQ Makefile.inc

    Then, edit the ``Makefile.inc`` file and replace in it ``f90`` by
    ``gfortran`` and ``cc`` by ``gcc``. Then compile MUMPS and install it with::

      $ make all
      $ cp lib/*.a /usr/local/lib
      $ cp libseq/*.a /usr/local/lib
      $ cp include/*.h /usr/local/include


  - At this stage, you should be able to make a first compilation of |gf|
    without the interfaces.  Download the last released version of |gf| on
    `download and install <../download.html>`_ (you can either download the
    git current version, however, it necessitates the additional installation
    of m4, autoconf, automake and libtool). Then, untar the package, enter into
    the source directory of |gf| with Msys and compile with::

      $ ./configure --with-blas="-lblas -lgfortran"
      $ cd superlu
      $ make
      $ cp .libs/libsuperlu.a /usr/local/lib
      $ cd ..
      $ ./configure --with-blas="-lblas -lgfortran" --disable-superlu
      $ make
      $ make check

    the ``make check`` is not necessary, but it is to verify that the
    compilation is correctly done. Note the separate compilation of
    the ``superlu`` library is due to some difficulties with Msys
    command ``ar``.

Build with the Python interface
*******************************

Additionnaly to build the Python interface, you will have first to install a 64bits version of Python 2 on your system together with Numpy and Scipy packages. This is not completely simple, but you can follow the following steps 

  - Install a 64 bits Python 2 or 3 version
    (see `Python website <https://www.python.org/downloads/windows/>`_).
    Then, if it is not done by the installer you used, add ``C:\Pythonxx``
    to your Windows path (where ``xx`` is the version number).
    Close you Msys and Windows shell and re-open them to take into
    account the changes.

  - Install Pip (see `Pip <https://pip.pypa.io/en/latest/installing/>`_)

  - Downloads the precompiled packages of numpy and scipy for 64 bits
    and Python `here <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_
    (i.e. for isntance wheel packages
    ``numpy-1.11.3+mkl-cp27-cp27m-win_amd64.whl``
    and ``scipy-0.19.0-cp27-cp27m-win_amd64.whl`` for python 2).

  - Enter into the directory where you downloaded the two
    wheel Python packages and install them with::

      $ python -m pip install numpy-1.11.3+mkl-cp27-cp27m-win_amd64.whl
      $ python -m pip install scipy-0.19.0-cp27-cp27m-win_amd64.whl

  - Go to the Getfem source directory and perform again::

      $ ./configure --with-blas="-lblas -lgfortran" --disable-superlu

    if all is ok, the configure script should detect the Python installation
    and the presence of Numpy and Scipy packages and allow the compilation
    of the Python interface. Then just perform a::

      $ make

    You just remain to add as a Windows system variable ``PYTHONPATH`` with the
    value ``c:\msys\home\login\getfem-5.?\interface\src\python`` where
    ``login`` and ``5.?`` have tobe adapted. You can either copy the
    directory ``interface\src\python`` where the interface has been built
    in a Python 2 directory.

Still some problems ... to be completed soon

Build with the Matlab interface
*******************************

Here follows the additional step to build the Matlab interface. You have first, of course to have installed a (recent) version of Matlab on your system.

  - You also need a installation of Python, because some Python scripts
    are used to build the interface. You can follow the steps described
    in the previous section for Python interface installation.
    However, for the Matlab interface, Numpy and Scipy are not required.

  - Upload `gnumex <http://gnumex.sourceforge.net/>`_ and run it under Matlab
    (see the indications on gnumex website). Make it generate the
    ``mexopts.bat`` et ``mexopts.stp`` files in an accessible directory.

  - In the source top directory of |gf|, add the file ``gnumex.opts``
    with the following lines (to be adapted)::

      #!/bin/sh
      MATLAB_ROOT="c:\Program Files\MATLAB\R20???"
      MATLAB_RELEASE=14
      MATLAB_INC_DIR="$MATLAB_ROOT\extern\include"
      MEXOPTS=c:\path_to_mexopts\mexopts.bat

  - Re-run the configure script and compile with::

      $ ./configure --with-blas="-lblas -lgfortran" --disable-superlu --enable-matlab
      $ make

    The Matlab interface should be compiled without error. If there is some
    linker errors, go to the ``interface\src\matlab`` directory
    and try to adapt the library list in the matlab/mex call by copy/past the
    command building the interface.

  - Once the interface is properly compiled, you can lauch Matlab, add the
    path of the compiled interface (on the matlab command line) with::

      add_path('c:\msys\home\login\getfem-5.?\interface/tests/matlab')
    
    and try the demo matlab programs of the interface in
    ``interface\tests\matlab``. In order not to have to call the ``addpath``
    command each time you open Matlab, you can add a Windows system variable
    ``MATLABPATH`` set to
    ``c:\msys\home\login\getfem-5.?\interface/tests/matlab``.
    You can also move the ``interface\tests\matlab`` directory into your 
    Matlab installation directory.



