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

Y. Collette (ycollet at freesurf dot fr)
