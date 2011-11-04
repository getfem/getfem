readme.txt of the scilab getfem interface

To compile this interface, you will need to work with the scilab-5.2.2
version or the scilab-master version (the version in the git
repository) because this interface uses the new interface api.

So, to compile this interface:
- compile getfem. I use the following configure script:
/configure --prefix=/home/collette/repositories/install_dir/getfem-dev/ \
           --with-pic \
           --enable-superlu=yes \
           --enable-muparser=yes \
           --enable-qhull=yes \
           --enable-scilab \
           --with-scilab-toolbox-dir=<getfem_scilabinstalldir> \
           --with-scilab-prefix=<scilabinstalldir> \
           --with-optimization=-ggdb \
           BLAS_LIBS="-L/usr/lib64/atlas/ -lblas"
  don't forget to install the following package: qhull and muparser.
  to get better performances, install atlas. If you don't install
  the atlas package, remove the BLAS_LIBS line.
Once getfem is compiled:
- go to the scilab getfem++ interface install directory (getfem_scilabinstalldir here).
- launch scilab
Now load the getfem++ toolbox:
- exec loader.sce;
You can try to launch a demo (be careful, there is a lot of work needed before the interface can be really useable).
- cd demos
- exec demo_static_contact.sce;


* Some hints for the compilation of this toolbox under windows.

For the compilation of the toolbox under windows:
- compile getfem + getfem interface using the visual studio 2010 project.
  copy the lib files (from msvc2010/Release) into scilab/src/win32 ou
  src/win64
  If you plan to add support for qhull and muparser in the windows
  library, you must add:
  - for qhull:
    - GETFEM_HAVE_QHULL_QHULL_H in the preprocessor
    - the path to the include where we can find qhull/qhull.h
  - for muparser:
    - GETFEM_HAVE_MUPARSER_H in the preprocessor
    - the path to the include where we can find muParser.h

- download and compile statically qhull using visual studio 2010.
  copy libqhull.lib into scilab/src/win32 ou src/win64

- download and compile statically muparser using visual studio 2010.
  copy libmuparser.lib into scilab/src/win32 ou src/win64

Now, you can go into the scilab directory.
Launch Scilab and do:

exec builder.sce;

Best regards,

Y. Collette (ycollette dot nospam at free dot fr)
