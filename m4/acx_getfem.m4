dnl Copyright (C) 2004-2020 Julien Pommier
dnl 
dnl This file is  free software;  you  can  redistribute  it  and/or modify it
dnl under  the  terms  of the  GNU  Lesser General Public License as published
dnl by  the  Free Software Foundation;  either version 3 of the License,  or
dnl (at your option) any later version along with the GCC Runtime Library
dnl Exception either version 3.1 or (at your option) any later version.
dnl This program  is  distributed  in  the  hope  that it will be useful,  but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl License and GCC Runtime Library Exception for more details.
dnl You  should  have received a copy of the GNU Lesser General Public License
dnl along  with  this program;  if not, write to the Free Software Foundation,
dnl Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.



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
