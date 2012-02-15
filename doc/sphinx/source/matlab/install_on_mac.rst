
.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: bash

.. _mlab-install_mac:


Installing the matlab interface for getfem 4.0.0 on snow leopard.
=================================================================

The MATLAB version considered here is a recent one (2009b).

This matlab version requires some specific flags to be used when building getfem. These flags are displayed when I run "mex -v"::

 CFLAGS = -fno-common -no-cpp-precomp -arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5  -fexceptions
 CXXFLAGS = -fno-common -no-cpp-precomp -fexceptions -arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5

Those that are important here are the arch one (you need to build a binary for the same architecture than the matlab one (ppc, ppc64, i386, x86_64)). The -isysroot and the -mmacos-min-version are used to linked against the same system library versions than matlab.


If you want to install qhull (in order to use the levelset stuff), 
you need to install it first. This is optional::

 ----------------------------QHULL INSTALL (optional)
 Build qhull: You need to use the same options that are used for building getfem:

 cd qhull-2010.1/src
 make CCOPTS1="-O2 -arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5"

 # Now install it into a standard location so that getfem configure will detect qhull

 sudo mkdir /usr/include/qhull/
 sudo install *.h /usr/include/qhull/
 sudo install libqhull.a /usr/lib/
 sudo mkdir /Developer/SDKs/MacOSX10.5.sdk/usr/include/qhull/
 sudo install *.h /Developer/SDKs/MacOSX10.5.sdk/usr/include/qhull/
 sudo install libqhull.a /Developer/SDKs/MacOSX10.5.sdk/usr/lib/
 ---------------------------END OF QHULL INSTALL


Hence I will pass them to the ./configure script::

 ./configure --enable-matlab CXXFLAGS="-arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5" CFLAGS="-arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5" --with-matlab-toolbox-dir=$HOME/matlab/getfem

Which should end with an encouraging::

 Lapack library found : -llapack
 ---------------------------------------
 Ready to build getfem
   building MATLAB interface: YES
   building PYTHON interface: NO (requires numpy)
 ---------------------------------------

But if you look at the output of the configure script, and if you happen to use the same matlab version than me, you might see::

 checking for mex... mex
 checking for matlab path...  /Applications/MATLAB_R2009b.app
 checking for mex extension...  .mexmaci
 grep: /Applications/MATLAB_R2009b.app/extern/src/mexversion.c: No such file or directory
 Matlab release is : R

Obviously the configure script failed to recognize the matlab version number...

You now need to edit two files in order to be able to build the getfem toolbox without error:

 - Open src/getfem_interpolated_fem.cc , and replace "uint" by "unsigned" on line 260 and 295

 - open interface/src/matlab/gfm_common.h and add
   #define MATLAB_RELEASE 2009
   at the top of the file

Now launch the compilation (I'm putting -j2 because I have a dual-core)::

 make -j2

It will take a long time to complete (20 minutes)
in order to install the toolbox, just create a directory for it, for example in ::

 mkdir -p $HOME/matlab/getfem

and copy all files into it (the "make install" does not work, unfortunately)::

 cp -pr interface/src/matlab/* $HOME/matlab/getfem

remove the assert.m which is useless.

now launch Matlab. In order to be able to use the toolbox, add it to your matlab path::

 >> addpath('~/matlab/getfem')

Test that the mex file loads correctly::

 >> gf_workspace('stats')
 message from [gf_workspace]:
 Workspace 0 [main -- 0 objects]


Go to the getfem test directory for matlab::

 >> cd interface/tests/matlab

And try the various tests::

 >> demo_laplacian
 >> demo_tripod
 etc..


