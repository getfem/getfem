Building the python interface on Mingw32:

I assume python 2.4 or better is installed, with the numpy extension.

I have installed mingw, and the msys shell.

In msys, go in the getfem++ directory, and run the ./configure script.

Do not try to use the --enable-python , it does not seem to detect python installations on window..
Just go ahead and compile the getfem interface (i.e. just type "make"). When it is built, go in the "interface/src/python" directory,
and compile the getfem_python.c file:



gcc -c -I../src -I../../../src -I(path to python)/include getfem_python.c

then link with libgetfemint, libgetfem and libpython libraries:

g++ -shared -o _getfem.dll ./getfem_python.o ../.libs/libgetfemint.a ../../../src/.libs/libgetfem.a (path to python)/libs/libpython24.a

and then you have your python extension. 

Now, build the getfem.py wrapper:

../../../interface/bin/extract_pydoc python .. > getfem.py

and it's finished. You can place these two files anywhere, just add the path to them in your PYTHONPATH environment variable before running python.



-----------------------
HOW TO INSTALL GetFEM/Matlab toolbox
-----------------------

The installation of the GetFEM Matlab toolbox requires :

- perl installed in /usr/bin/perl

- MATLAB at least version 6 and mex for compiling MATLAB interface. 

- a decent c++ compiler (gnu/g++ >= 3.0 , compaq/cxx >= 6.2 ). There is no
obligation to use the same compiler as the one with which matlab was compiled.

Hence you just have to compile getfem++ with 
./configure CXX=mycompiler --enable-matlab --disable-shared

Remark: the OpenGL renderer of matlab displays some artifacts with getfem 3D graphics.
You should change it for the zbuffer one with set(gcf, 'renderer', 'zbuffer').

Quick Install:

( tar xzvf getfem++-4.x.x.tar.gz && cd getfem++-4.x.x && ./configure --disable-shared --enable-matlab && make check )

----------------------------------------------------------------------
Detailed Installation process (if quick install does not work..) 
----------------------------------------------------------------------

- unpack archive

- ./configure CXX=mycompiler --enable-matlab --disable-shared

- the mex files and .m files will go in the directory specified by the
  --with-matlab-toolbox-dir option. Nothing else will be installed.

- if you already had an old version of the toolbox in this directory, run
  make clean in order to start on a good basis.

- if configure did not complain, just run
     make 

- if the make succeeded, run 
     make install

- run "make check" to perform some basic checks. Try the various demo in the tests directory
(for example demo_tripod.m)

- alternative build procedure: use the --enable-matlab-rpc
option. This will build a small mex file and a 'getfem_server'
executable. The advantage is that getfem and matlab a not intermixed
and communicate via sockets (UNIX or INET). It can be useful to find if a crash is
coming from getfem, or from matlab. It is also necessary with some compilers
(for example it is not possible to build a mex file with the intel c++ compiler).
The getfem_server will be launched automatically by matlab, or you can launch it
"by hand" with ./getfem_server -tcp (but this requires the portmap service running).

The add the toolbox-dir to your MATLABPATH.

You can check that everything is correctly installed by running matlab, and
typing the name of any getfem++/matlab command (for example gf_delete).  

* If matlab does not find it, check your matlabpath.

* If the function complains about missing arguments, then everything is
 fine. You can start reading the user manual!

