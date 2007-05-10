Building the python interface on Mingw32:

I assume python 2.4 or better is installed, with the numarray extension.

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

