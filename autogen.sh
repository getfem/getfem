#!/bin/bash
#
# Copyright (C) 2001-2009 Yves Renard
#
# This file is a part of GETFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

function die {
      echo "ERROR: $1";
          exit 1
}

aclocal -I ./m4 || die "aclocal failed";
libtoolize -f || glibtoolize -f || die "libtoolize failed";
autoheader || die "autoheader failed";
autoreconf
autoconf || die "autoconf failed";
#pas de ./ devant les noms des makefiles !!!
automake -a --gnu `find . -name Makefile.am | sed -e 's@\./\(.*\)\.am@\1@g'` || die "automake failed";
echo "autogen.sh is ok, you can run the ./configure script"
