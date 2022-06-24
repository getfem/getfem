## Copyright (C) 2009-2020 Yann Colette
## 
##  This scilab getfem interface is a part of GetFEM++
## 
##  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
##  under  the  terms  of the  GNU  Lesser General Public License as published
##  by  the  Free Software Foundation;  either version 3 of the License,  or
##  (at your option) any later version along with the GCC Runtime Library
##  Exception either version 3.1 or (at your option) any later version.
##  This program  is  distributed  in  the  hope  that it will be useful,  but
##  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
##  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
##  License and GCC Runtime Library Exception for more details.
##  You  should  have received a copy of the GNU Lesser General Public License
##  along  with  this program;  if not, write to the Free Software Foundation,
##  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.


readme.txt of the scilab getfem interface

To compile this interface, you will need to work with the scilab-5.2.2
version or the scilab-master version (the version in the git
repository) because this interface uses the new interface api.

So, to compile this interface:
- compile getfem. I use the following configure script:
/configure --prefix=/home/collette/repositories/install_dir/getfem-dev/ \
           --with-pic \
           --enable-superlu=yes \
           --enable-qhull=yes \
           --enable-scilab \
           --with-scilab-toolbox-dir=<getfem_scilabinstalldir> \
           --with-scilab-prefix=<scilabinstalldir> \
           --with-optimization=-ggdb \
           BLAS_LIBS="-L/usr/lib64/atlas/ -lblas"
  don't forget to install the following package: qhull.
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
  copy the lib files (from msvc2010/Release) into scilab/src/win32 or
  src/win64
  If you plan to add support for qhull in the windows
  library, you must add:
  - for qhull:
    - GETFEM_HAVE_QHULL_QHULL_H in the preprocessor
    - the path to the include where we can find qhull/qhull.h

- download and compile statically qhull using visual studio 2010.
  copy libqhull.lib into scilab/src/win32 ou src/win64


Now, you can go into the scilab directory.
Launch Scilab and do:

exec builder.sce;

Best regards,

Y. Collette (ycollette dot nospam at free dot fr)
