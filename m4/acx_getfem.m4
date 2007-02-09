AC_DEFUN([ACX_GETFEM],
[
AC_ARG_WITH(getfem-config, [  --with-getfem-config : path and name of the getfem-config script (needed if it is not in the PATH)], getfem_config="$withval", getfem_config="getfem-config")

if test -x $getfem_config; then
  $getfem_config; ret="$?";
  if test "$ret" -eq 0; then
    echo "the configuration script $getfem_config seems to be ok.."
  else  
    AC_MSG_ERROR(["there is a problem with getfem-config (returned $ret). Aborting"])
  fi;
else
  AC_CHECK_PROG([gfconfig],[$getfem_config],1,0)
  if test $gfconfig -ne 1; then
    AC_MSG_ERROR(
[The getfem-config script was not found. Either the getfem library
is not installed, either the script is not in the PATH.
You can specify the full path an name of the script with --with-getfem-config])
  fi;
fi;
AC_MSG_CHECKING([for getfem++ configuration flags])
GETFEM_CPPFLAGS=`$getfem_config --cflags`
GETFEM_LIBS=`$getfem_config --libs`
GETFEM_LIBS_LA=`$getfem_config --libs-la`
GETFEM_STATICLIBS=`$getfem_config --static-libs`
GETFEM_CXX=`$getfem_config --cxx`
GETFEM_VERSION=`$getfem_config --version`
GETFEM_BUILD=`$getfem_config --build`
GETFEM_LIB_LA=`$getfem_config --libs-la | sed -e 's/\.la .*$/\.la/'`

AC_SUBST(GETFEM_LIBS)
AC_SUBST(GETFEM_STATICLIBS)
AC_SUBST(GETFEM_LIBS_LA)
AC_SUBST(GETFEM_LIB_LA)
AC_SUBST(GETFEM_CPPFLAGS)
AC_SUBST(GETFEM_CXX)
AC_MSG_RESULT([done])
AC_MSG_NOTICE(["Found GetFem++ version $GETFEM_VERSION, build: $GETFEM_CXX,$GETFEM_BUILD"])
])dnl ACX_GETFEM