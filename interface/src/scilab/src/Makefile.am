## SIVP - Scilab Image and Video Processing toolbox
## Copyright (C) 2009 Yann COLLETTE
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

libsp_getdir=@SCILAB_TOOLBOX_INSTALL_DIR@/lib

libsp_get_DATA = c/libsp_get.so

all: c/libsp_get.so

c/libs_get.so: builder.sce
	cd c
	$(SCILAB_EXE) -nw -nb -f builder_c.sce

clean:
	rm -f *.a *.o *.lo *.so

distclean:
	rm -fr *.a *.o *.lo *.so builder.sce .libs
