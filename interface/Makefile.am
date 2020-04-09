#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM
#
#  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

#http://sources.redhat.com/automake/automake.html#Local-Macros
#ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests

install-toolbox: all
	cd src && $(MAKE) install-toolbox

#EXTRA_DIST=BUGS README.alpha_cxx README.sgiCC README.darwin gnumex autogen.sh

#DISTCHECK_CONFIGURE_FLAGS = ...PLOPLOPLOPAREMPLIR
host_canonical=@host_cpu@-@host_vendor@-@host_os@

bindist: all	
	mkdir bdist
	cp @srcdir@/src/*.m bdist
	cp -a @srcdir@/src/private bdist
	cp -a @srcdir@/src/@gf* bdist
	cp src/*.m bdist
	cp src/gf_matlab@MATLAB_COM_EXT@ bdist
	mv bdist getfem-matlab-$(host_canonical) && tar czvf getfem-matlab-$(host_canonical).tar.gz getfem-matlab-$(host_canonical)
