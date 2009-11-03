## Scilab Getfem Interface
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

top_scilabdir=`pwd`

scilabbasedir=$(top_scilabdir)/builder.sce \
	$(top_scilabdir)/DESCRIPTION-FUNCTIONS \
	$(top_scilabdir)/loader.sce \
	$(top_scilabdir)/readme.txt \
	$(top_scilabdir)/changelog.txt \
	$(top_scilabdir)/cleaner.sce \
	$(top_scilabdir)/Makefile.am \
	$(top_scilabdir)/sci_getfem.iss \
	$(top_scilabdir)/DESCRIPTION \
	$(top_scilabdir)/license.txt \
	$(top_scilabdir)/makefile_builder.sce

scilabsrccdir=$(top_scilabdir)/src/c/copyright \
	$(top_scilabdir)/src/c/makefile \
	$(top_scilabdir)/src/c/loader.sce \
	$(top_scilabdir)/src/c/makefile.in \
	$(top_scilabdir)/src/c/builder_c.sce \
	$(top_scilabdir)/src/c/ls.dat \
	$(top_scilabdir)/src/c/README \
	$(top_scilabdir)/src/c/rk4.dat \
	$(top_scilabdir)/src/c/cleaner.sce \
	$(top_scilabdir)/src/c/configure \
	$(top_scilabdir)/src/c/configure.in \
	$(top_scilabdir)/src/c/FILELIST \
	$(top_scilabdir)/src/c/machine.h.in  

toolboxdir=@SCILAB_TOOLBOX_DIR@

all:
	@SCILAB_EXE@ -nw -f @srcdir@/makefile_builder.sce

install-toolbox:
	$(mkinstalldirs) $(toolboxdir)/demos
	$(mkinstalldirs) $(toolboxdir)/etc
	$(mkinstalldirs) $(toolboxdir)/help
	$(mkinstalldirs) $(toolboxdir)/help/en_US
	$(mkinstalldirs) $(toolboxdir)/help/en_US/examples
	$(mkinstalldirs) $(toolboxdir)/help/fr_FR
	$(mkinstalldirs) $(toolboxdir)/help/fig
	$(mkinstalldirs) $(toolboxdir)/help/mml
	$(mkinstalldirs) $(toolboxdir)/jar
	$(mkinstalldirs) $(toolboxdir)/macros
	$(mkinstalldirs) $(toolboxdir)/sci_gateway
	$(mkinstalldirs) $(toolboxdir)/sci_gateway/c
	$(mkinstalldirs) $(toolboxdir)/src
	$(mkinstalldirs) $(toolboxdir)/src/c
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/Cray
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/GCC
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/Linux
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/RS6000
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/SGI
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/SPARC
	$(mkinstalldirs) $(toolboxdir)/src/c/DOC
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos $(top_scilabdir)/demos/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/etc $(top_scilabdir)/etc/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help $(top_scilabdir)/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US $(top_scilabdir)/help/en_US/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US $(top_scilabdir)/help/en_US/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples $(top_scilabdir)/help/en_US/examples/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples $(top_scilabdir)/help/en_US/examples/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fig $(top_scilabdir)/help/fig/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/mml $(top_scilabdir)/help/mml/*.mml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/jar $(top_scilabdir)/jar/*.jar
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(top_scilabdir)/macros/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(top_scilabdir)/macros/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(top_scilabdir)/macros/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(top_scilabdir)/macros/lib
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway $(top_scilabdir)/sci_gateway/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c $(top_scilabdir)/sci_gateway/c/*.[ch]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c $(top_scilabdir)/sci_gateway/c/*.sce
	@INSTALL@ -D -m 744 -t $(toolboxdir)/sci_gateway/c $(top_scilabdir)/sci_gateway/c/*.so
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src $(top_scilabdir)/src/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c $(top_scilabdir)/src/c/*.[ch]
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c $(scilabsrccdir)
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Cray $(top_scilabdir)/src/c/MACHINES/Cray/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/GCC $(top_scilabdir)/src/c/MACHINES/GCC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Linux $(top_scilabdir)/src/c/MACHINES/Linux/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/RS6000 $(top_scilabdir)/src/c/MACHINES/RS6000/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SGI $(top_scilabdir)/src/c/MACHINES/SGI/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SPARC $(top_scilabdir)/src/c/MACHINES/SPARC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/DOC $(top_scilabdir)/src/c/DOC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ $(scilabbasedir)
