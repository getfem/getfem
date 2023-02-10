dnl  Copyright (C) 1999-2020 Yves Renard
dnl
dnl  This file is a part of GetFEM
dnl
dnl  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
dnl  under  the  terms  of the  GNU  Lesser General Public License as published
dnl  by  the  Free Software Foundation;  either version 3 of the License,  or
dnl  (at your option) any later version along with the GCC Runtime Library
dnl  Exception either version 3.1 or (at your option) any later version.
dnl  This program  is  distributed  in  the  hope  that it will be useful,  but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License and GCC Runtime Library Exception for more details.
dnl  You  should  have received a copy of the GNU Lesser General Public License
dnl  along  with  this program;  if not, write to the Free Software Foundation,
dnl  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.



dnl Process this file with autoconf to produce a configure script.
dnl ------------------------------------------------------------------------
dnl initialisation
dnl ------------------------------------------------------------------------

dnl ./configure: sh internal 2K buffer overflow on HP-UX 9.xx
dnl thus, updating cache ./config.cache avoided.
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl

AC_INIT(getfem, 5.4.2)
MAJOR_VERSION="5"
MINOR_VERSION="4"
PATCH_VERSION="2"
AC_DEFINE_UNQUOTED([MAJOR_VERSION],$MAJOR_VERSION,[Major version number])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$MINOR_VERSION,[Minor version number])
AC_DEFINE_UNQUOTED([PATCH_VERSION],$PATCH_VERSION,[Patch version number])

AC_CONFIG_SRCDIR([install-sh])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER(config.h)
AX_PREFIX_CONFIG_H(src/getfem/getfem_arch_config.h,GETFEM) 
AX_PREFIX_CONFIG_H(src/gmm/gmm_arch_config.h,GMM) 
AC_PREREQ(2.61)
AC_ARG_PROGRAM

dnl ------------------------------------------------------------------------
dnl   init automake
dnl ------------------------------------------------------------------------

dnl AM_INIT_AUTOMAKE($PACKAGE,$VERSION)
AM_INIT_AUTOMAKE([1.11 parallel-tests]) 
WARNING_MSG=""

CPPFLAGS="$CPPFLAGS -I/usr/local/include"

dnl --------------------------
dnl set the optimization level
dnl --------------------------

AC_ARG_WITH(optimization,
	    AC_HELP_STRING([--with-optimization=FLAG],[Set the optimization level (-O3 by default)]),
	    [with_optimization=$withval],
	    [with_optimization='-O3']
	    )	    

dnl ---------------------------PARA LEVEL--------------------------
paralevel=0
AC_ARG_ENABLE(paralevel,
   [AS_HELP_STRING([--enable-paralevel[=level]],[enable the parallel version of GetFEM (use MPI and METIS)])],
   [ case $enableval in
        yes | "") paralevel=2;;
        no) ;;
        *) paralevel=$enableval ;;
     esac
])

if test $paralevel -ge 1; then
  CPPFLAGS="$CPPFLAGS -DGETFEM_PARA_LEVEL=$paralevel"
fi;
dnl ---------------------------END OF PARA LEVEL-------------------



dnl -----------------------------------------------
dnl test du c++
dnl -----------------------------------------------

USER_CXXFLAGS="$CXXFLAGS"
USER_CFLAGS="$CFLAGS"
AX_PROG_CXX_MPI([test $paralevel -ge 1],[usempi=yes],[usempi=no])
AX_PROG_CC_MPI([test "x$usempi" = "xyes"],,[usempi=no])
AX_PROG_FC_MPI([test "x$usempi" = "xyes"],[CPPFLAGS="$CPPFLAGS -DGMM_USES_MPI"],[usempi=no])

if  test "x$usempi" = "xyes"; then
AC_CHECK_LIB(mpi_cxx, _init,,)
fi

AC_PROG_CXXCPP
CXXFLAGS="${USER_CXXFLAGS}"
CFLAGS="${USER_CFLAGS}"
SUPLDFLAGS=""
AC_FC_LIBRARY_LDFLAGS

AC_LANG([C++])

if test "x$prefix" = "xNONE"; then
  GFPREFIX=/usr/local;
else
  GFPREFIX="$prefix";
fi;

dnl AC_CXX_FULL_SPECIALIZATION_SYNTAX (c)Luc Maisonobe v 1.1.1.1 (2001/07/26)     0.5.41 
dnl with some modification to test partial specialization
AC_CACHE_CHECK(whether the compiler recognizes the partial specialization syntax,
ac_cv_cxx_partial_specialization_syntax,
[AC_DIAGNOSE([obsolete],[Instead of using `AC_LANG', `AC_LANG_SAVE', and `AC_LANG_RESTORE',
you should use `AC_LANG_PUSH' and `AC_LANG_POP'.])dnl
AC_LANG_SAVE
 AC_LANG([C++])
 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
template<class T> class A        { public : int f () const { return 1; } };
template<class T> class A<T*>    { public : int f () const { return 0; } };]], [[
A<float*> a; return a.f();]])],[ac_cv_cxx_partial_specialization_syntax=yes],[ac_cv_cxx_partial_specialization_syntax=no])
 AC_LANG_POP([])
])
if test "$ac_cv_cxx_partial_specialization_syntax" != yes; then
  echo "Your compiler ($CXX) does not support partial template specialization, trash it"
  exit 1;
fi

AC_CANONICAL_HOST

echo "you are compiling GetFEM on a $host"

case $CXX in
	
 *g++* | c++)
	GCCVER=`$CXX --version | head -1 | cut -d ' ' -f3`
	echo "Using the GNU g++ compiler $GCCVER"
	AC_CHECK_CXX_FLAG([$with_optimization],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-fmessage-length=0],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-fvisibility-inlines-hidden],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-ftemplate-depth-100],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-std=c++14], CXXFLAGS,[AC_CHECK_CXX_FLAG([-std=c++0x],CXXFLAGS,[AC_MSG_ERROR([g++ do not support option -std=c++14. Update g++ to at least release 4.9])] )])
	AC_CHECK_CXX_FLAG([-fPIC],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wall -W -Wextra],CXXFLAGS)
	dnl AC_CHECK_CXX_FLAG([-pedantic],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wshadow],CXXFLAGS)
	dnl AC_CHECK_CXX_FLAG([-Wno-terminate],CXXFLAGS)
	dnl AC_CHECK_CXX_FLAG([-Wno-implicit-fallthrough],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wno-unknown-pragmas],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wno-variadic-macros],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wno-unused-but-set-variable],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wpointer-arith],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wcast-qual],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wwrite-strings],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wconversion],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wredundant-decls],CXXFLAGS)
	dnl -Wno-long-double fixes a warning on Darwin
	dnl AC_CHECK_CXX_FLAG([-Wno-long-double],CXXFLAGS)
	AC_CHECK_CXX_FLAG([-Wno-long-long],CXXFLAGS)
	dnl -rdynamic used for backtraces
	AC_CHECK_CXX_FLAG([-rdynamic],SUPLDFLAGS)

dnl	SUPLDFLAGS="-rdynamic" # -rdynamic for backtraces
        CFLAGS="$CFLAGS $with_optimization -fPIC"
	;;
 icc | icpc)
	echo "Using INTEL icc"
	CXXFLAGS="$CXXFLAGS $with_optimization -Xc -ansi -std=c++14"
        dnl CXXFLAGS="$CXXFLAGS $with_optimization -Xc -ansi"
        CFLAGS="$CFLAGS $with_optimization -Xc -ansi"
	;;
 *)
	echo "Using a unknown compiler"
 	CXXFLAGS="$CXXFLAGS $with_optimization -std=c++14"
        dnl CXXFLAGS="$CXXFLAGS $with_optimization"
        CFLAGS="$CFLAGS $with_optimization"
	;;
esac

AC_SUBST(SUPLDFLAGS)

dnl ------------------------------------------------------------------------
dnl   init libtools for shared libraries
dnl ------------------------------------------------------------------------

dnl option pic-only is not working: a libtool bug ...
LT_INIT([pic-only disable-shared])
AC_SUBST([LIBTOOL_DEPS])


dnl -----------------------------------------------
dnl MATLAB Interface
dnl -----------------------------------------------

# list of pseudo functions
PSEUDO_FUNCTIONS_LOC=`$srcdir/bin/extract_doc $srcdir/interface/src pseudo_loc`
echo $PSEUDO_FUNCTIONS_LOC
PSEUDO_FUNCTIONS=`$srcdir/bin/extract_doc $srcdir/interface/src pseudo_gen`
MATLAB_OBJ_DIRS=`$srcdir/bin/extract_doc $srcdir/interface/src mobj_dirs`
AC_SUBST(PSEUDO_FUNCTIONS)
AC_SUBST(PSEUDO_FUNCTIONS_LOC)
AC_SUBST(MATLAB_OBJ_DIRS)

AC_ARG_ENABLE(matlab,
 [AS_HELP_STRING([--enable-matlab],[turn on/off matlab support])],
 [case "${enableval}" in
   yes) usematlab=YES ;;
   no)  usematlab=NO ;;
   *) AC_MSG_ERROR([bad value ${enableval} for --enable-matlab]) ;;
 esac],[usematlab=NO])

AC_ARG_WITH(matlab-toolbox-dir,
            [AS_HELP_STRING([--with-matlab-toolbox-dir],[directory in which the matlab interface will be installed])],
            TOOLBOXDIR="$withval",TOOLBOXDIR="$GFPREFIX/getfem_toolbox")
AC_SUBST(TOOLBOXDIR)

if test "$usematlab" != NO; then
  AC_CHECK_PROGS(MEX, mex)
  if test x"$MEX" = x""; then
    AC_CHECK_PROGS(MEX, mex.bat)
    if test x"$MEX" = x""; then
      if test x$usematlab = xYES; then
        AC_MSG_ERROR([Impossible to build the matlab interface without mex -- specify its full path with the MEX=/path/to/mex option, or use --enable-matlab=no])
        exit 1
      fi
    else
      MEX=gnumex;
      MATLAB_COM_EXT=".dll";
      echo "You are using Matlab on a windows platform (assuming MingW compiler)";
      if test -f gnumex.opts; then
         echo "sourcing gnumex.opts.."
         source gnumex.opts;         
         echo "MATLAB_ROOT=$MATLAB_ROOT"
         echo "Matlab release is : R$MATLAB_RELEASE"
      elif test x$usematlab = xYES; then
        echo "You need to fill the gnumex.opts file, for example (use MSys-style paths, not DOS-style paths)"
        echo '#!/bin/sh'
        echo 'MATLAB_ROOT="c:\\MATLAB6p5"'
        echo 'MATLAB_RELEASE=13'
        echo 'MATLAB_INC_DIR="$MATLAB_ROOT\\extern\\include"'
        echo 'MEXOPTS=c:\\gnumex\\mexopts.bat'
        echo "when this is done, check that the gnumex script works correctly"
        echo " (i.e. gnumex gnumex.opts -v prints the rights options to use the MinGW gcc)"
        exit 1
      fi
    fi
  else
     dnl thanks to paolo for pointing the 'twin mex' problem
     if $(echo "" | $MEX 2>&1 | grep 'This is .*TeX'); then
	  AC_MSG_ERROR([the mex binary which is in the PATH appears to be part of LaTeX, not matlab !! run ./configure MEX=/path/to/matlab/mex]);
     fi;
    
     # MATLAB_ROOT=`$MEX -v 2>&1 | grep "MATLAB  " | awk '{print $4}'|sed -e '2,$d'`
     MEX_EXE=`which $MEX`
     MEX_EXE=`readlink -e $MEX_EXE`
     MATLAB_ROOT=`echo $MEX_EXE | sed -e 's/bin.*$//'`
     MEXEXT=$MATLAB_ROOT/bin/mexext
     Matlab_INC_DIR=$MATLAB_ROOT/extern/include
     echo "checking for matlab path... " $MATLAB_ROOT
     MATLAB_COM_EXT=.`$MEXEXT`
     echo "checking for mex extension... " $MATLAB_COM_EXT
#    MATLAB_RELEASE=`grep "MATLAB R" $MATLAB_ROOT/extern/src/mexversion.c | awk '{print $4}' | sed -e 's/R//'`
#    MATLAB_RELEASE=`grep "full_ver="$(which $MEX) | sed 's/[[^0-9]]//g'` # double brackets are for escaping reasons.
     MATLAB_RELEASE=`matlab -nosplash -nojvm -r "a=version;display(a);exit" | grep "(R2" | sed -n '1p' | sed 's/^.*R//g' | sed 's/).*$//'`
     echo "Matlab release is : R$MATLAB_RELEASE"
  fi
fi
AM_CONDITIONAL(BUILDMEX, test x$usematlab = xYES)

AC_SUBST(MATLAB_ROOT)
AC_SUBST(MATLAB_INC_DIR)
AC_SUBST(MATLAB_RELEASE)
AC_SUBST(MATLAB_COM_EXT)
AC_SUBST(MEX)

AM_CONDITIONAL(USE_MINGW_MEX, test x"$MATLAB_COM_EXT" = x".dll")



dnl ----------------------------
dnl RPCs -- matlab interface communication with a separated getfem process
dnl useful for debugging..
GETFEM_SERVER="";
use_rpc="no";
AC_ARG_ENABLE(matlab-rpc,
 [AS_HELP_STRING([--enable-matlab-rpc],[enable use of RPCs for matlab interface])],
 [ matlab_rpc="yes"; use_rpc="yes";
   echo "Matlab mex-file will use sun RPCs in order to communicate with the getfem server"],
 [matlab_rpc="no"])

if test x$use_rpc = xyes; then
  GETFEM_SERVER="getfem_server";
  AC_ARG_WITH(rpc-include,
              [AS_HELP_STRING([--with-rpc-include],[directory in which the rpc/rpc.h header can be found])],
              RPC_INC_DIR="-I$withval",RPC_INC_DIR="")
  case $host in
        *alpha*)
                RPC_LIB="-lrpc";
                ;;
	*darwin*)
	        RPC_LIB="";
		;;
        *)
                RPC_LIB="-lnsl";
                ;;
  esac
  AC_ARG_WITH(rpc-lib,
              [AS_HELP_STRING([--with-rpc-lib],[linker flags for the RPC library])],
              RPC_LIB="$withval")
  AC_SUBST(RPC_INC_DIR)
  AC_SUBST(RPC_LIB)
  AC_DEFINE_UNQUOTED(USE_RPC, 1, [Use rpc for getfem communication with matlab])
fi;
AC_SUBST(GETFEM_SERVER)
AM_CONDITIONAL(BUILDMEXRPC, test x$matlab_rpc = xyes)


dnl the pb is that we cannot link the libstdc++.so in the mex-file without horrible problems
dnl with dynamic_casts (with matlab 6.5 -- the pb seems to have disappeared since matlab-7). 
dnl Hence the gf_matlab.mexglx should be linked against the libstdc++.a ..
STDCPP_STATICLIBS=""

if test $usematlab = xYES; then
  dnl ------------------------------------
  dnl COMPILER SETTINGS
  compiler_type=dontcare
  case $CXX in
   *g++* | c++)
	case $host in
	x86_64-*)
	       echo "Compiling on an x86_64 architecture..."
	       ;;
        *-darwin*)
               echo "Compiling on Darwin (MacOS)"
		;;
	*)
		STDCPP_STATICLIBS=$($CXX -print-file-name=libstdc++.a)
		echo "The MEX file will be linked against the static c++ library '$STDCPP_STATICLIBS'"
		;;
	esac
	;;
   *icc | *icpc)
	dnl a small remark: with icpc 8.0, the getfem_server will crash 
	dnl at the first exception throwed (except with -g)
	dnl the fix is to pass the -static flag at the linker
	dnl unfortunately, the lovely libtool assumes that icpc won't
	dnl understand it, and removes it. I hate libtool.
	dnl so I added the -Wl,-static -- it works for now.
	GFSERVERFLAGS="-Wl,-static -static"
	;;
   *)
	;;
  esac
fi
AC_SUBST(GFSERVERFLAGS)
AC_SUBST(STDCPP_STATICLIBS)






dnl -----------------------------------------------
dnl Octave Interface
dnl -----------------------------------------------

# list of pseudo functions
# PSEUDO_FUNCTIONS_LOC=`$srcdir/bin/extract_doc $srcdir/interface/src pseudo_loc`
# echo $PSEUDO_FUNCTIONS_LOC
# PSEUDO_FUNCTIONS=`$srcdir/bin/extract_doc $srcdir/interface/src pseudo_gen`
# MATLAB_OBJ_DIRS=`$srcdir/bin/extract_doc $srcdir/interface/src mobj_dirs`
# AC_SUBST(PSEUDO_FUNCTIONS)
# AC_SUBST(PSEUDO_FUNCTIONS_LOC)
# AC_SUBST(MATLAB_OBJ_DIRS)

AC_ARG_ENABLE(octave,
 [AS_HELP_STRING([--enable-octave],[turn on/off octave support])],
 [case "${enableval}" in
   yes) useoctave=YES ;;
   no)  useoctave=NO ;;
   *) AC_MSG_ERROR([bad value ${enableval} for --enable-octave]) ;;
 esac],[useoctave=NO])

AC_ARG_WITH(octave-toolbox-dir,
            [AS_HELP_STRING([--with-octave-toolbox-dir],[directory in which the octave interface will be installed])],
            OCTAVETOOLBOXDIR="$withval",OCTAVETOOLBOXDIR="$GFPREFIX/getfem_octavetoolbox")
AC_SUBST(OCTAVETOOLBOXDIR)

if test "$useoctave" != NO; then
  AC_CHECK_PROGS(MKOCTFILE, mkoctfile)
  if test x"$MKOCTFILE" = x""; then
    if test x$useoctave = xYES; then
      AC_MSG_ERROR([Impossible to build the octave interface without mkoctfile -- specify its full path with the MKOCTFILE=/path/to/mkoctfile option, or use --enable-octave=no])
        exit 1
    fi
  else
     # MATLAB_ROOT=`$MEX -v 2>&1 | grep "MATLAB  " | awk '{print $4}'|sed -e '2,$d'`
     MKOCTFILE_EXE=`which $MKOCTFILE`
     MKOCTFILE_EXE=`readlink -e $MKOCTFILE_EXE`
     OCTAVE_ROOT=`echo $MKOCTFILE_EXE | sed -e 's/bin.*$//'`
     OCTAVE_COM_EXT=".mex"
     OCTAVE_RELEASE=`octave --version | grep "GNU Octave" | sed -n '1p' | sed 's/^.*version //g'`
     echo "Octave release is $OCTAVE_RELEASE"
  fi
fi
AM_CONDITIONAL(BUILDOCTAVE, test x$useoctave = xYES)

AC_SUBST(OCTAVE_RELEASE)
AC_SUBST(OCTAVE_COM_EXT)
AC_SUBST(MKOCTFILE)

AM_CONDITIONAL(USE_MINGW_OCTAVE, test x"$OCTAVE_COM_EXT" = x".dll")

dnl ----------------------------------------------
dnl python 
dnl ----------------------------------------------

AC_ARG_ENABLE(python,
 [AS_HELP_STRING([--enable-python],[turn on/off python support])],
 [case "${enableval}" in
   yes) usepython=YES ;;
   no)  usepython=NO ;;
   *) AC_MSG_ERROR([bad value ${enableval} for --enable-python]) ;;
 esac],[usepython=YES])

if test x$usepython = xYES; then
  AM_PATH_PYTHON(3.6, usepython=YES, usepython=NO)
  if test x$usepython = xNO; then
    AM_PATH_PYTHON(2.7, usepython=YES, usepython=NO)
  fi
fi

case $host in
    *mingw*)
        PYTHON_CC_ARG="--compiler=mingw32"
	PYTHON_EXTRA_EXT_PARAM="extra_compile_args=[['-D MS_WIN64']]"
        ;;
    *)
        PYTHON_CC_ARG=""
	PYTHON_EXTRA_EXT_PARAM=""
        ;;
esac

AC_SUBST(PYTHON_CC_ARG)
AC_SUBST(PYTHON_EXTRA_EXT_PARAM)

if test "x$usepython" = "xYES"; then
  echo "Building with python ($PYTHON) support (use --enable-python=no to disable it) "
  echo "You will need the python-numpy and python-scipy packages."
dnl  AM_PATH_PYTHON(2.2)
  AX_PYTHON_DEVEL
  ac_python_numpy=`$PYTHON -c 'import numpy; print("YES")' 2>/dev/null`
  if  test "x$ac_python_numpy" != "xYES"; then
    usepython=NO
    MSG="PYTHON DISABLED: numpy not found. You need to install the python-numpy package."
    echo $MSG
    WARNING_MSG="$WARNING_MSG\n$MSG"
  fi
  ac_python_scipy=`$PYTHON -c 'import scipy; print("YES")' 2>/dev/null`
  if  test "x$ac_python_scipy" != "xYES"; then
    usepython=NO
    MSG="PYTHON DISABLED: scipy not found. You need to install the python-scipy package."
    echo $MSG
    WARNING_MSG="$WARNING_MSG\n$MSG"
  fi
  if  test $paralevel -ge 1; then
    ac_python_mpi4py=`$PYTHON -c 'import mpi4py; print("YES")' 2>/dev/null`
    if test "x$ac_python_mpi4py" != "xYES"; then
      usepython=NO
      MSG="PARALLEL PYTHON DISABLED: mpi4py not found. You need to install the python-mpi4py package."
      echo $MSG
      WARNING_MSG="$WARNING_MSG\n$MSG"
    fi
  fi
  PYTHON_SO=`$PYTHON -c "import sysconfig, sys; get = sysconfig.get_config_var; sys.stdout.write(get('EXT_SUFFIX') or get('SO') or '.so');"`
  AC_SUBST(PYTHON_SO)
fi

AM_CONDITIONAL(BUILDPYTHON, test x$usepython = xYES)
AM_CONDITIONAL(BUILDPYTHONPAR, test x$ac_python_mpi4py = xYES)

dnl -----------------------------------------------
dnl SCILAB Interface
dnl -----------------------------------------------

m4_include([m4/scilab.m4])

REQUIRED_SCILAB_MAJOR=6
REQUIRED_SCILAB_MINOR=0
REQUIRED_SCILAB_MICRO=0

AC_CHECK_SCILAB

GETFEM_INTERFACE_PATH="`readlink -f $srcdir`"
GETFEM_BUILD_INTERFACE_PATH="`readlink -f $PWD`"
AC_SUBST(GETFEM_INTERFACE_PATH)
AC_SUBST(GETFEM_BUILD_INTERFACE_PATH)

dnl if the scilab directory doesn't exists, we copy the 
dnl scilab sources into the build directory

if test "x$usescilab" == "xYES"
then
  currentdir=`pwd`
  if test ! -f $currentdir/interface/src/scilab/builder.sce
  then
    echo "Copying Scilab toolbox src in the build directory"
    mkdir -p $currentdir/interface/src/scilab/
    cp -r $srcdir/interface/src/scilab/* $currentdir/interface/src/scilab
  fi
fi

AC_ARG_WITH(scilab-toolbox-dir,
            [AS_HELP_STRING([--with-scilab-toolbox-dir],[directory in which the scilab interface will be installed])],
            SCILAB_TOOLBOX_DIR="$withval",SCILAB_TOOLBOX_DIR="$GFPREFIX/getfem_toolbox")
AC_SUBST(SCILAB_TOOLBOX_DIR)


dnl -------------------------------BLAS----------------------------------

dnl why I hate autoconf: if the code below is put into a separate file,
dnl the generated ./configure will stop if no Fortran compiler is found. always. even
dnl if no AC_FC_FUNC is executed.
acx_blas_ok=no

AC_ARG_WITH(blas,
        [AS_HELP_STRING([--with-blas=<lib>],[use BLAS library <lib>])])
case $with_blas in
        yes | "") ;;
        no) acx_blas_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o| builtin) BLAS_LIBS="$with_blas" ;;
        *) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
if test x"$FC" = "x"; then
  echo "No fortran compiler found, assuming c-name for SGEMM is 'sgemm_'"
  sgemm=sgemm_
  dgemm=dgemm_
else
  AC_FC_FUNC(sgemm)
  AC_FC_FUNC(dgemm)
fi
acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"
echo "BLAS_LIBS=$BLAS_LIBS"
# First, check BLAS_LIBS environment variable
if test "x$BLAS_LIBS" = xbuiltin; then
  echo "Using builtin blas lib";
  BLAS_LIBS=""
else

if test $acx_blas_ok = no; then
  if test "x$BLAS_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
        AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
        AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
        AC_MSG_RESULT($acx_blas_ok)
        LIBS="$save_LIBS"
  fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS"
        AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(atlas, ATL_xerbla,
                [AC_CHECK_LIB(f77blas, $sgemm,
                [AC_CHECK_LIB(cblas, cblas_dgemm,
                        [acx_blas_ok=yes
                         BLAS_LIBS="-lf77blas -latlas $FCLIBS"],
                        [], [-lf77blas -latlas])],
                        [], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm,
                [AC_CHECK_LIB(dgemm, $dgemm,
                [AC_CHECK_LIB(sgemm, $sgemm,
                        [acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
                        [], [-lblas])],
                        [], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
        if test "x$GCC" != xyes; then # only works with Sun CC
                AC_CHECK_LIB(sunmath, acosp,
                        [AC_CHECK_LIB(sunperf, $sgemm,
                                [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
        fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(complib.sgimath, $sgemm,
                     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm,
                [AC_CHECK_LIB(essl, $sgemm,
                        [acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
                        [], [-lblas $FLIBS])])
fi

# Generic BLAS library with fortran dependency?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas $FCLIBS"])
fi

fi # if BLAS_LIBS=builtin

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
	echo "OK, You have working BLAS libs ! Using $BLAS_LIBS" ; HAVE_VENDOR_BLAS=1
else
        echo " *** YOU DONT HAVE BLAS! *** Using a cheap replacement" ; HAVE_VENDOR_BLAS=0
fi

dnl ACX_BLAS([ echo "OK, You have working BLAS libs !"; HAVE_VENDOR_BLAS=1 ], [echo "YOU DONT HAVE BLAS! Using a cheap replacement" ; HAVE_VENDOR_BLAS=0])
LIBS="$LIBS $BLAS_LIBS"
CPPFLAGS="$CPPFLAGS -DGMM_USES_BLAS"

useblas64support=NO
AC_ARG_ENABLE(blas64-support,
        [AS_HELP_STRING([--enable-blas64-support],[enable the 64 bits integer blas and lapack support])],
        [ case $enableval in
          yes | "") useblas64support=YES;;
          no) useblas64support=NO;;
          esac], [useblas64support=NO])

if test x$useblas64support = xNO && test x$usematlab = xYES; then
  # 32 or 64 bits blas interface
  AC_RUN_IFELSE([AC_LANG_SOURCE([[
  #include <stdlib.h>
  int main() {
    int a; long b;
    int *pa; long *pb; void *pc;

    if (sizeof(long) > sizeof(int)) {
      a = 0x3E0024; b = 0;
      pb = &b; pc = (void *)(pb); pa = (int *)(pc); *pa = a;
      if (int(b) != a) exit(1);
      *pa = 0x3E0024;
      if (int(b) != a) exit(1);
      if (b != long(a)) exit(1);
    } 
  } ]])],[useblas64support=YES],[],[])
fi

if test x$useblas64support = xYES; then
  AC_DEFINE_UNQUOTED(USE_BLAS64_INTERFACE, 1, [Use blas with 64 bits integers])
fi

dnl ------------------------------LAPACK TEST--------------------------------

if test x"$acx_blas_ok" = xyes; then
#This name mangling detection is currently not in use
#  if test x"$FC" = "x"; then
#    dgetrf=dgetrf_
#  else
#    AC_FC_FUNC(dgetrf)
#  fi;

  # LAPACK linked to by default?  (is sometimes included in BLAS lib)
  AC_CHECK_FUNC(dgetrf_, [acx_lapack_ok=yes; LAPACK_LIBS=""])
  if test "x$acx_lapack_ok" != xyes; then
    AC_CHECK_LIB(lapack, dgetrf_, [acx_lapack_ok=yes; LAPACK_LIBS="-llapack"])
  fi

  if test x"$acx_lapack_ok" = xyes; then
     CPPFLAGS="$CPPFLAGS -DGMM_USES_LAPACK"
     LIBS="$LAPACK_LIBS $LIBS"
  fi
fi

dnl -----------------------------END OF LAPACK TEST--------------------------


dnl ---------------------------multithread-blas--------------------------

multithread_blas=NO
AC_ARG_ENABLE(multithread-blas,
   [AS_HELP_STRING([--enable-multithread-blas],[enable openblas to be multithreaded])],
   [ case $enableval in
        yes | "") multithread_blas=YES;;
        no) multithread_blas=NO;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-multithread-blas], [multithread_blas=NO]) ;;
     esac
])

if test x$multithread_blas = xNO; then
  AC_CHECK_LIB(dl, dlsym, [acx_dl_ok=yes; DL_LIBS="-ldl"])
  if test x"$acx_dl_ok" = xyes; then
     AC_DEFINE_UNQUOTED([FORCE_SINGLE_THREAD_BLAS],1,[enable openblas to be multithreaded])
     LIBS="$DL_LIBS $LIBS"
  fi
fi;
dnl --------------------------end of multithread-blas--------------------






dnl ---------------------------OPENMP------------------------------
useopenmp=0
AC_ARG_ENABLE(openmp,
   [AS_HELP_STRING([--enable-openmp],[enable the multihreaded version of GetFEM])],
   [ case $enableval in
        yes | "") useopenmp=YES ;;
        no) useopenmp=NO ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
     esac],
   [useopenmp=NO]
)

if test x$useopenmp = xYES; then
  AC_OPENMP
  if test "x$ac_cv_prog_cxx_openmp" != "xunsupported" && test "x$ac_cv_prog_cxx_openmp" != "x"; then
    AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
    CPPFLAGS="$CPPFLAGS -DGETFEM_HAS_OPENMP"
  else
    AC_MSG_ERROR([OpenMP support not found. Use --enable-openmp=no flag to compile GetFEM without OpenMP]);
  fi
fi;
dnl ---------------------------END OF OPENMP-----------------------


dnl ------------------------------SuperLU config-------------------------
AC_ARG_ENABLE(superlu,
 [AS_HELP_STRING([--enable-superlu],[turn on/off SuperLU support])],
 [case "${enableval}" in
   yes) usesuperlu=YES ;;
   no)  usesuperlu=NO ;;
   *) AC_MSG_ERROR([bad value ${enableval} for --enable-superlu]) ;;
 esac],[usesuperlu=YES])

SUPERLU_CPPFLAGS=""
SUPERLU_SRC=""
SUPERLU_LIBS=""
SUPERLU_MAKEFILE=""

if test x$usesuperlu = xYES; then
  echo "Building with SuperLU support (use --enable-superlu=no to disable it)"
  if test x"$FC" = "x"; then
    sgemm="sgemm_"
  else
    AC_FC_FUNC(sgemm)
    echo "FC=$FC"
  fi
  case $sgemm in
    sgemm)
          F77_CALL_C="NOCHANGE";
          ;;
    sgemm_)
          F77_CALL_C="ADD_";
          ;;
    SGEMM)
          F77_CALL_C="UPCASE";
          ;;
    sgemm__)
          F77_CALL_C="ADD__";
          ;;
    *)
          AC_MSG_ERROR(["superlu won't handle this calling convention: sgemm -> $sgemm"])
          ;;
  esac
  SUPERLU_CPPFLAGS="$CPPFLAGS -DUSE_VENDOR_BLAS -DF77_CALL_C=$F77_CALL_C"
  SUPERLU_SRC="superlu"
  case $host in
    *apple*)
        SUPERLU_LIBS="../$SUPERLU_SRC/libsuperlu.la"
        ;;
    *mingw*)
        SUPERLU_LIBS="../$SUPERLU_SRC/.libs/libsuperlu.a"
        ;;
    *)
        SUPERLU_LIBS="`readlink -f .`/$SUPERLU_SRC/libsuperlu.la"
        ;;
  esac
  SUPERLU_MAKEFILE="$SUPERLU_SRC/Makefile"
else
  echo "Building without SuperLU support (use --enable-superlu=yes to enable it)"
  AC_CHECK_LIB([superlu], [dCreate_CompCol_Matrix],[],
               [AC_MSG_ERROR([SuperLU library not found])])

  AC_CHECK_HEADERS(
  [superlu/colamd.h superlu/slu_Cnames.h \
   superlu/slu_cdefs.h superlu/slu_ddefs.h superlu/slu_sdefs.h superlu/slu_zdefs.h \
   superlu/slu_dcomplex.h superlu/slu_scomplex.h],
  [usesuperlu="YES"],
  [
    if test "x$usesuperlu" = "xYES"; then
      AC_MSG_ERROR([header files of superlu not found. Use --enable-superlu=yes flag]);
    fi;
  ])

  SUPERLU_LIBS="-lsuperlu"
  LIBS="$SUPERLU_LIBS $LIBS"
fi

AC_SUBST([SUPERLU_CPPFLAGS])
AC_SUBST([SUPERLU_SRC])
AC_SUBST([SUPERLU_LIBS])
AM_CONDITIONAL(USEBLASLITE, test x$HAVE_VENDOR_BLAS = x0)
echo "Configuration of SuperLU done"


dnl ----------------EXPERIMENTAL PARTS OF THE LIBRARY--------------------
EXPER=""
AC_ARG_ENABLE(experimental,
        [AS_HELP_STRING([--enable-experimental],[compile experimental parts of the library])],
[ if   test "x$enableval" = "xyes" ; then EXPER="-DEXPERIMENTAL_PURPOSE_ONLY"; fi], [EXPER=""])
CPPFLAGS="$CPPFLAGS $EXPER"

dnl -----------------------------QD TESTS--------------------------------
AC_ARG_WITH(qd-lib-dir,
        [AS_HELP_STRING([--with-qd-lib-dir],[directory in which the libqd.a can be found])],
	QDLIB="$withval/libqd.a",QDLIB="$GFPREFIX/lib/libqd.a")
AC_ARG_WITH(qd-include-dir,
        [AS_HELP_STRING([--with-qd-include-dir],[directory in which the qd.h header can be found])],
	QDINC="-I$withval",QDINC="-I$GFPREFIX/include")
AC_ARG_ENABLE(dd,
 [AS_HELP_STRING([--enable-dd],[enable the use of the qd library (some computation will be done with double-double precision, useful for high order FEMs)])],
 [ if   test "x$enableval" = "xyes" ; then useQDlib="yes"; QD_PREC="double"; fi], [useQDlib="no"])
AC_ARG_ENABLE(qd,
 [AS_HELP_STRING([--enable-qd],[enable the use of the qd library (some computation will be done with quad-double precision, useful for high order FEMs)])],
 [ if   test "x$enableval" = "xyes" ; then useQDlib="yes"; QD_PREC="quad"; fi], [if test "x$useQDlib" = "xyes"; then useQDlib="yes"; else useQDlib="no"; fi])
if test "x$useQDlib" = "xyes" ; then  
  LIBS="$LIBS $QDLIB -lm"
  CPPFLAGS="$CPPFLAGS $QDINC"
dnl #define NO_INLINE
  AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <qd/qd_real.h>
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <iostream>
int main() {
  unsigned int old_cw;
  int ok;
  fpu_fix_start(&old_cw);
  qd_real q = 1.0;
  qd_real qq = qd_real("0.01");
  qd_real qqq = "1.010101010101010101010101010101010101010101010101010101010101010E0";
  dd_real d = 1.0;
  dd_real dd = dd_real("0.1");
  dd_real ddd = "1.1111111111111111111111111111111E0";
  for (int i=0; i < 100; ++i) { d += dd; dd *= dd_real("0.1"); }
  for (int i=0; i < 100; ++i) { q += qq; qq *= qd_real("0.01"); }
  std::cerr << "d = " << d << std::endl << "q = " << q << std::endl;
  std::cerr << abs(q - qqq) << std::endl;
  std::cerr << abs(d - ddd) << std::endl;
  if (abs(q - qqq) < 1e-63 && abs(d -ddd) < 1e-31) ok = 1;
  else ok = 0;
  fpu_fix_end(&old_cw); return 1-ok;
}
  ]])],[echo "checking if qd library is working...yes"],[ echo "QD library is not working (check config.log)"; exit 1],[])
  AC_DEFINE_UNQUOTED([HAVE_QDLIB],1,[defined if the qd library was found and is working])
  HAVE_QDLIB=1;
  if test "x$QD_PREC" = "xquad"; then
    AC_DEFINE_UNQUOTED([QDLIB_USE_QUAD],1,[defined if quad-doubles are to be used instead of double-double])
  fi;
fi;
dnl -----------------------------END QD TESTS--------------------------------

dnl ------------------------------QHULL TEST---------------------------------
useQHULL="no"
AC_ARG_ENABLE(qhull,
 [AS_HELP_STRING([--enable-qhull],[enable the use of the qhull library (required for generation of non regular meshes)])],
 [ if   test "x$enableval" = "xyes" ; then useQHULL="yes"; fi], [useQHULL="test"])
QHULL_LIBS=""
save_LIBS="$LIBS";

if test "x$useQHULL" = "xno"; then
  echo "Building with libqhull explicitly disabled";
else
  AC_CHECK_LIB(qhull_r, qh_new_qhull, [QHULL_LIBS="-lqhull_r"],
  [
    AC_CHECK_LIB(qhullstatic_r, qh_new_qhull,[QHULL_LIBS="-lqhullstatic_r"],[QHULL_LIBS=""])
  ])
  AC_CHECK_HEADERS(libqhull_r/qhull_ra.h,[useQHULL="yes"],
  [
    if test "x$useQHULL" = "xyes"; then
      AC_MSG_ERROR([header files libqhull_r/qhull_ra.h not found. Use --enable-qhull=no flag]);
      useQHULL="no"
    fi;
  ])
  echo "Building with libqhull (use --enable-qhull=no to disable it)"
fi;
AM_CONDITIONAL(QHULL, test x$useQHULL = xyes)

LIBS="$QHULL_LIBS $save_LIBS"
AC_SUBST([QHULL_LIBS])
echo "Configuration of qhull done"
dnl -----------------------------END OF QHULL TEST---------------------------

dnl ------------------------------MUMPS TEST------------------------------
MUMPSINC=""
AC_ARG_WITH(mumps-include-dir,
 [AS_HELP_STRING([--with-mumps-include-dir],[directory in which the dmumps.h header can be found])],
 [case $withval in
   -I* ) MUMPSINC="$withval";;
   * ) MUMPSINC="-I$withval";;
  esac],
 [MUMPSINC="-I$GFPREFIX/include"]
)
CPPFLAGS="$CPPFLAGS $MUMPSINC"

MUMPS_LIBS=""
case $host in
    *mingw*)
        MUMPS_SEQ_LIBS="-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord"
        ;;
    *apple*)
        MUMPS_SEQ_LIBS="-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord -lgomp"
        ;;
    *)
        MUMPS_SEQ_LIBS="-lsmumps_seq -ldmumps_seq -lcmumps_seq -lzmumps_seq"
        ;;
esac
acx_mumps_ok="no"
usemumps="no"
AC_ARG_ENABLE(mumps,
 [AS_HELP_STRING([--enable-mumps],[enable the use of the (sequential) MUMPS library. A direct solver for large sparse linear systems.])],
 [case $enableval in
   yes | "") usemumps="yes"; acx_mumps_ok="yes"; MUMPS_LIBS="$MUMPS_SEQ_LIBS";;
   no) usemumps="no";;
  esac],
 [usemumps="test"; acx_mumps_ok="test"; MUMPS_LIBS="$MUMPS_SEQ_LIBS"]
)

AC_ARG_ENABLE(par-mumps,
 [AS_HELP_STRING([--enable-par-mumps],[enable the use of the parrallel MUMPS library. A direct solver for large sparse linear systems.])],
 [case $enableval in
   yes | "") usemumps="yes"; MUMPS_LIBS="-lsmumps -ldmumps -lcmumps -lzmumps";;
   no) usemumps="no";;
  esac],
 [if test $paralevel -ge 1; then
    usemumps="test"; acx_mumps_ok="test"; MUMPS_LIBS="-lsmumps -ldmumps -lcmumps -lzmumps"
  fi;]
)

AC_ARG_WITH(mumps,
 [AS_HELP_STRING([--with-mumps=<lib>],[use MUMPS library <lib>])],
 [case $with_mumps in
   yes | "") usemumps="yes";;
   no) acx_mumps_ok="no" ;;
   -* | */* | *.a | *.so | *.so.* | *.o| builtin) MUMPS_LIBS="$with_mumps"; acx_mumps_ok="yes" ;;
   *) MUMPS_LIBS=`echo $with_mumps | sed -e 's/^/-l/g;s/ \+/ -l/g'` ; usemumps="yes";;
  esac]
)

save_LIBS="$LIBS";
if test "x$usemumps" = "xno" -o "x$acx_mumps_ok" = "xno"; then
  echo "Building with MUMPS explicitly disabled";
else
 AC_SEARCH_LIBS(smumps_c, [`echo $MUMPS_LIBS | sed -e 's/^-l//g;s/ -l/ /g'`],
   [usemumps="yes"],
   [if test "x$acx_mumps_ok" = "xyes"; then
     AC_MSG_ERROR([The function smumps_c couldn't be found in the provided MUMPS libraries.]);
    fi;
    usemumps="no"]
 )
 AC_SEARCH_LIBS(dmumps_c, [`echo $MUMPS_LIBS | sed -e 's/^-l//g;s/ -l/ /g'`],
   [usemumps="yes"],
   [if test "x$acx_mumps_ok" = "xyes"; then
     AC_MSG_ERROR([The function dmumps_c couldn't be found in the provided MUMPS libraries.]);
    fi;
    usemumps="no"]
 )
 AC_SEARCH_LIBS(cmumps_c, [`echo $MUMPS_LIBS | sed -e 's/^-l//g;s/ -l/ /g'`],
   [usemumps="yes"],
   [if test "x$acx_mumps_ok" = "xyes"; then
     AC_MSG_ERROR([The function cmumps_c couldn't be found in the provided MUMPS libraries.]);
    fi;
    usemumps="no"]
 )
 AC_SEARCH_LIBS(zmumps_c, [`echo $MUMPS_LIBS | sed -e 's/^-l//g;s/ -l/ /g'`],
   [usemumps="yes"],
   [if test "x$acx_mumps_ok" = "xyes"; then
     AC_MSG_ERROR([The function zmumps_c couldn't be found in the provided MUMPS libraries.]);
    fi;
    usemumps="no"]
 )
 AC_CHECK_HEADERS([smumps_c.h dmumps_c.h cmumps_c.h zmumps_c.h],
   [usemumps="yes"],
   [if test "x$acx_mumps_ok" = "xyes"; then
     AC_MSG_ERROR([header file dmumps_c.h not found.]);
    fi;
    usemumps="no"]
 )

 if test "x$usemumps" = "xyes"; then
   echo "Building with MUMPS (use --enable-mumps=no to disable it)"
   LIBS="$MUMPS_LIBS $save_LIBS"
 else
   MUMPS_LIBS=""
   LIBS="$save_LIBS"
 fi;
fi;

AM_CONDITIONAL(MUMPS, test x$usemumps = xyes)
AC_SUBST([MUMPS_LIBS])
echo "Configuration of MUMPS done"
dnl ---------------------------END OF MUMPS TEST--------------------------

dnl ---------------------------METIS--------------------------
METIS_LIBS=""
AC_ARG_ENABLE(metis,
 [AS_HELP_STRING([--enable-metis],[enable the use of the METIS library.])],
 [case $enableval in
   yes | "") usemetis="yes" ;;
   no) usemetis="no"; METIS_LIBS="" ;;
  esac],
 [usemetis="test"]
)

if test $paralevel -ge 1 -a "x$usemetis" = "xno"; then
  echo "Parallel getfem requires the METIS library, --enable-metis=no will be ignored";
  usemetis="yes"
fi;

if test "x$usemetis" = "xno"; then
  echo "Building without METIS";
else
  if test $paralevel -ge 1; then
    AC_CHECK_LIB(metis, METIS_PartGraphRecursive,
                 [usemetis="yes"],
                 [usemetis="no";
                  AC_MSG_ERROR([METIS library required for parallel getfem was not found])
                 ])

    METIS_LIBS="-lmetis"
    LIBS="$LIBS $METIS_LIBS"
    AC_DEFINE_UNQUOTED([HAVE_METIS],1,[defined if the Metis library was found and is working])
    AC_CHECK_LIB(metis, METIS_SetDefaultOptions, [usemetisnew="yes"],
                 [AC_DEFINE_UNQUOTED([HAVE_METIS_OLD_API],1,
                                     [defined if the Metis library found is older than version 4])
                 ])
    echo "Building with METIS (use --enable-metis=no to disable it)"
    if test "x$usemetisnew" = "xyes"; then
      AC_CHECK_HEADERS(metis.h,
                       [usemetis="yes"],
                       [usemetis="no";
                        AC_MSG_ERROR([metis.h header required for parallel getfem was not found])
                       ])
    fi;
  fi;
fi;

AM_CONDITIONAL(METIS, test x$usemetis = xyes)
AC_SUBST([METIS_LIBS])


dnl ---------------------------END OF METIS--------------------------


AC_CHECK_HEADERS(sys/times.h,[],[SUPERLU_CPPFLAGS="$SUPERLU_CPPFLAGS -DNO_TIMER"])
AC_CHECK_HEADERS(cxxabi.h)
dnl ---------------------------- CHECK FOR __PRETTY_FUNCTION__ MACRO --------
AC_CACHE_CHECK([for __PRETTY_FUNCTION__], ac_cv_have_pretty_function, [
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [
                [ const char *s = __PRETTY_FUNCTION__; ]])],
                [ ac_cv_have_pretty_function="yes" ],
                [ ac_cv_have_pretty_function=="no"  ])])
if test "x$ac_cv_have_pretty_function" = "xyes"; then
        AC_DEFINE_UNQUOTED(HAVE_PRETTY_FUNCTION,1,[gcc style __PRETTY_FUNCTION__ macro])
fi;     


dnl ---------------------------- CHECK FOR GLIBC BACKTRACE availability -----
AC_CACHE_CHECK([for execinfo.h and backtrace], ac_cv_have_backtrace, [
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                [[ #include <execinfo.h>  ]],
                [[ void* trace[256]; int n = backtrace(trace, 256); ]])],
                [ ac_cv_have_backtrace="yes" ],
                [ ac_cv_have_backtrace="no"  ])])
if test "x$ac_cv_have_backtrace" = "xyes"; then
        AC_DEFINE_UNQUOTED(HAVE_BACKTRACE,1,[glibc backtrace function])
fi;     

dnl ---------------------------- CHECK FOR feenableexcept -----
AC_CACHE_CHECK([for fenv.h and feenableexcept], ac_cv_have_feenableexcept, [
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                [[ #include <fenv.h>           ]], 
                [[ feenableexcept(FE_DIVBYZERO | FE_INVALID); ]])],
                [ ac_cv_have_feenableexcept="yes" ],
                [ ac_cv_have_feenableexcept="no"  ])])
if test "x$ac_cv_have_feenableexcept" = "xyes"; then
        AC_DEFINE_UNQUOTED(HAVE_FEENABLEEXCEPT,1,[glibc floating point exceptions control])
fi;

BUILDER=`whoami`
AC_SUBST(BUILDER)
BUILDDATE=`date +%D,%H:%M:%S`
AC_SUBST(BUILDDATE)
CONFIGURE_ARGS=$ac_configure_args
AC_SUBST(CONFIGURE_ARGS)
LIBTOOL_VERSION_INFO="-version-info ${MAJOR_VERSION}:${MINOR_VERSION}:${PATCH_VERSION}"
AC_SUBST(LIBTOOL_VERSION_INFO)

dnl AC_CHECK_PROGS(RANLIB, ranlib)


dnl ------------ for distclean of meshes ---------------------
j="tests/meshes/disc_P2_h4.mesh"
if test -L $j || test ! -f $j; then
  DISTCLEANMESH="";
else
  DISTCLEANMESH="#";
fi;
AC_SUBST(DISTCLEANMESH)



dnl -----------------------------------------------
dnl Outputs
dnl -----------------------------------------------

IM_METHODS=`$srcdir/bin/extract_doc $srcdir/interface/src cubature`
IM_METHODS_LOC=`$srcdir/bin/extract_doc $srcdir/interface/src cubature_loc`
AC_SUBST(IM_METHODS)
AC_SUBST(IM_METHODS_LOC)


AC_CONFIG_FILES(                                                        \
Makefile                                                                \
m4/Makefile                                                             \
cubature/Makefile                                                       \
$SUPERLU_MAKEFILE                                                       \
doc/Makefile                                                            \
doc/sphinx/Makefile                                                     \
src/Makefile                                                            \
tests/Makefile                                                          \
contrib/Makefile                                                        \
contrib/icare/Makefile                                                  \
contrib/delaminated_crack/Makefile                                      \
contrib/bimaterial_crack_test/Makefile                                  \
contrib/xfem_stab_unilat_contact/Makefile                               \
contrib/mixed_elastostatic/Makefile                                     \
contrib/xfem_contact/Makefile                                           \
contrib/crack_plate/Makefile                                            \
contrib/aposteriori/Makefile                                            \
contrib/level_set_contact/Makefile                                      \
contrib/static_contact_gears/Makefile                                   \
contrib/test_plasticity/Makefile                                        \
contrib/opt_assembly/Makefile                                           \
contrib/continuum_mechanics/Makefile                                    \
bin/Makefile                                                            \
interface/Makefile                                                      \
interface/src/Makefile                                                  \
interface/src/matlab/Makefile                                           \
interface/src/matlab/private/Makefile                                   \
interface/src/octave/Makefile                                           \
interface/src/octave/private/Makefile                                   \
interface/src/python/Makefile                                           \
interface/src/python/setup.py                                           \
interface/src/scilab/Makefile                                           \
interface/src/scilab/sci_gateway/c/builder_gateway_c.sce                \
interface/tests/Makefile                                                \
interface/tests/meshes/Makefile                                         \
interface/tests/matlab-octave/Makefile                                  \
interface/tests/matlab-octave/private/Makefile                          \
interface/tests/python/Makefile                                         \
getfem-config                                                           \
getfem-config-notinstalled                                              \
gmm-config)
AC_OUTPUT
chmod a+x getfem-config-notinstalled
chmod a+x getfem-config
chmod a+x gmm-config

dnl -----------------------------------------------
dnl Symbolic links for the meshes in tests/meshes
dnl -----------------------------------------------

if test -z ""`echo $srcdir | grep "^/"`; then
  addpathm="../"
else
  addpathm=""
fi

if test ! -d tests/meshes; then
  ln -s $addpathm$srcdir/tests/meshes tests/meshes
fi;


dnl configuration sum-up

echo
echo "------------------------------------------------------------------------------"
echo
echo "Libraries Used:"
echo "---------------"
echo



if test "x$useQDlib" = "xyes" ; then  
  echo "- QD library found. High precision (${QD_PREC}-double precision) polynomials"
  echo "  and integration methods are enabled.";
else
  echo "- QD library not found (don't worry, this library is only recommended for very specific uses)."
fi;

if test "x$useQHULL" = "xyes"; then
  echo "- Qhull found. Using the Qhull library for delaunay triangulations."
else
  echo "- Qhull not found. Mesh generation will be disabled."
fi;

if test "x$usemumps" = "xyes"; then
  echo "- Mumps found. A direct solver for large sparse linear systems."
else
  echo "- Mumps not found. Not using the MUMPS library for large sparse linear systems."
fi;

if test x"$acx_lapack_ok" = xyes; then
  echo "- Lapack library found: $LAPACK_LIBS"
else
  echo "- Lapack library not found: generic (less effective) algorithms will be used"
fi

if test "x$HAVE_VENDOR_BLAS" = "x0"; then
  echo "- *** No usable blas library was found ***"
  echo "  A generic BLAS implementation will be used, however you should "
  echo "  consider installing a faster BLAS, such as ATLAS"
else
  echo "- BLAS library found. Link options: $BLAS_LIBS"
fi;
echo "  You can give the location of your prefered blas library with either"
echo "  the --with-blas=<lib> option, or the BLAS_LIBS environment variable"
echo '  for example: ./configure BLAS_LIBS="-L/usr/lib/sse2/atlas/ -lblas"'
if test x$useblas64support = xYES; then
  echo "  Use blas with 64 bits integers or 32/64 bits compatibility mode"
else
  echo "  Use blas with 32 bits integers"
fi
echo -e "\n\n"

echo "-----------------------------------------------------------------------"
echo "Ready to build getfem"
echo "  building MATLAB interface: $usematlab"
echo "  building OCTAVE interface: $useoctave"
echo "  building PYTHON interface: $usepython (requires numpy, scipy and also mpi4py for the parallel version)"
echo "  building SCILAB interface: $usescilab"
echo "  If you want to build the shared library of GetFEM, use --enable-shared"
echo "  (by default, only the static one will be built)"
echo "-----------------------------------------------------------------------"

case $host in
  x86_64-*)
	if test $usematlab = "YES" -o $usepython = "YES"; then
          if test $pic_mode != "yes"; then
            echo "!!!!!"
            echo "!!!!! Your build will fail because you did not use the --with-pic option"
            echo "!!!!! This is required for the getfem interfaces on x86_64"
            echo ""
          fi
        fi      
  ;;
esac

if test "x$MSG" != "x"; then
  echo -e "\n\nWARNINGS during the configure:\n$MSG\n\n"
fi
