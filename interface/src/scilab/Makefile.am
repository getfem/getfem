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

scilabbuilddir=@GETFEM_BUILD_INTERFACE_PATH@/interface/src/scilab

scilabbasedir=builder.sce \
	loader.sce \
	readme.txt \
	changelog.txt \
	cleaner.sce \
	Makefile.am \
	sci_getfem.iss \
	license.txt \
	makefile_builder.sce

scilabsrccdir=src/c/copyright \
	src/c/makefile \
	src/c/loader.sce \
	src/c/makefile.in \
	src/c/builder_c.sce \
	src/c/ls.dat \
	src/c/README \
	src/c/rk4.dat \
	src/c/cleaner.sce \
	src/c/configure \
	src/c/configure.in \
	src/c/FILELIST \
	src/c/machine.h.in  

EXTRA_DIST=$(toolboxdir)/demos/*.sc[ie] \
	$(toolboxdir)/etc/* \
	$(toolboxdir)/help/*.sce \
	$(toolboxdir)/help/en_US/*.sce \
	$(toolboxdir)/help/en_US/*.xml \
	$(toolboxdir)/help/en_US/examples/CHAPTER \
	$(toolboxdir)/help/en_US/examples/*.xml \
	$(toolboxdir)/help/fig/* \
	$(toolboxdir)/help/mml/*.mml \
	$(toolboxdir)/jar/*.jar \
	$(toolboxdir)/macros/*.bin \
	$(toolboxdir)/macros/*.sc[ie] \
	$(toolboxdir)/macros/names \
	$(toolboxdir)/macros/lib \
	$(toolboxdir)/sci_gateway/*.sce \
	$(toolboxdir)/sci_gateway/c/*.[ch] \
	$(toolboxdir)/sci_gateway/c/*.sce \
	$(toolboxdir)/sci_gateway/c/*.so \
	$(toolboxdir)/src/*.sce \
	$(toolboxdir)/src/c/*.[ch] \
	$(toolboxdir)/src/c/*.so \
	$(toolboxdir)/src/c/MACHINES/Cray/* \
	$(toolboxdir)/src/c/MACHINES/GCC/* \
	$(toolboxdir)/src/c/MACHINES/Linux/* \
	$(toolboxdir)/src/c/MACHINES/RS6000/* \
	$(toolboxdir)/src/c/MACHINES/SGI/* \
	$(toolboxdir)/src/c/MACHINES/SPARC/* \
	$(toolboxdir)/src/c/DOC/*

toolboxdir=@SCILAB_TOOLBOX_DIR@

.NOTPARALLEL: *

all:
	@SCILAB_EXE@ -nw -nb -f $(scilabbuilddir)/makefile_builder.sce

install:
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
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos demos/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/etc etc/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help *.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fig help/fig/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/mml help/mml/*.mml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/jar jar/*.jar
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/lib
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway sci_gateway/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.[ch]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.sce
	@INSTALL@ -D -m 744 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.so
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src src/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c src/c/*.[ch]
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c src/c/*.so
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c $(scilabsrccdir)
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Cray src/c/MACHINES/Cray/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/GCC src/c/MACHINES/GCC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Linux src/c/MACHINES/Linux/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/RS6000 src/c/MACHINES/RS6000/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SGI src/c/MACHINES/SGI/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SPARC src/c/MACHINES/SPARC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/DOC src/c/DOC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ $(scilabbasedir)

clean:
	@SCILAB_EXE@ -nw -nb -f $(scilabbuilddir)/makefile_cleaner.sce
