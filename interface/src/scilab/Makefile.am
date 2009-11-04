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

scilabbuilddir=$(top_srcdir)/interface/src/scilab_build

scilabbasedir=$(srcdir)/builder.sce \
	$(srcdir)/DESCRIPTION-FUNCTIONS \
	$(srcdir)/loader.sce \
	$(srcdir)/readme.txt \
	$(srcdir)/changelog.txt \
	$(srcdir)/cleaner.sce \
	$(srcdir)/Makefile.am \
	$(srcdir)/sci_getfem.iss \
	$(srcdir)/DESCRIPTION \
	$(srcdir)/license.txt \
	$(srcdir)/makefile_builder.sce

scilabsrccdir=$(srcdir)/src/c/copyright \
	$(srcdir)/src/c/makefile \
	$(srcdir)/src/c/loader.sce \
	$(srcdir)/src/c/makefile.in \
	$(srcdir)/src/c/builder_c.sce \
	$(srcdir)/src/c/ls.dat \
	$(srcdir)/src/c/README \
	$(srcdir)/src/c/rk4.dat \
	$(srcdir)/src/c/cleaner.sce \
	$(srcdir)/src/c/configure \
	$(srcdir)/src/c/configure.in \
	$(srcdir)/src/c/FILELIST \
	$(srcdir)/src/c/machine.h.in  

toolboxdir=@SCILAB_TOOLBOX_DIR@

all:
	$(mkinstalldirs) $(top_srcdir)/interface/src/scilab_build
	@cp -r $(srcdir) $(top_srcdir)/interface/src/scilab_build
	@SCILAB_EXE@ -nw -f $(top_srcdir)/interface/src/scilab_build/makefile_builder.sce

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
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos $(scilabbuilddir)/demos/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/etc $(scilabbuilddir)/etc/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help $(scilabbuilddir)/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US $(scilabbuilddir)/help/en_US/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US $(scilabbuilddir)/help/en_US/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples $(scilabbuilddir)/help/en_US/examples/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples $(scilabbuilddir)/help/en_US/examples/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fig $(scilabbuilddir)/help/fig/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/mml $(scilabbuilddir)/help/mml/*.mml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/jar $(scilabbuilddir)/jar/*.jar
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(scilabbuilddir)/macros/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(scilabbuilddir)/macros/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(scilabbuilddir)/macros/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros $(scilabbuilddir)/macros/lib
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway $(scilabbuilddir)/sci_gateway/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c $(scilabbuilddir)/sci_gateway/c/*.[ch]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c $(scilabbuilddir)/sci_gateway/c/*.sce
	@INSTALL@ -D -m 744 -t $(toolboxdir)/sci_gateway/c $(scilabbuilddir)/sci_gateway/c/*.so
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src $(scilabbuilddir)/src/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c $(scilabbuilddir)/src/c/*.[ch]
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c $(scilabsrccdir)
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Cray $(scilabbuilddir)/src/c/MACHINES/Cray/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/GCC $(scilabbuilddir)/src/c/MACHINES/GCC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Linux $(scilabbuilddir)/src/c/MACHINES/Linux/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/RS6000 $(scilabbuilddir)/src/c/MACHINES/RS6000/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SGI $(scilabbuilddir)/src/c/MACHINES/SGI/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SPARC $(scilabbuilddir)/src/c/MACHINES/SPARC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/DOC $(scilabbuilddir)/src/c/DOC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ $(scilabbasedir)
