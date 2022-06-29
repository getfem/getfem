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

GETFEM_LIB_LA = ../../../src/libgetfem.la

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I@srcdir@ -I@srcdir@/..

gfpythondir=$(pythondir)/getfem
gfpyexecdir=$(pyexecdir)/getfem

gfpython_PYTHON = getfem.py __init__.py
nodist_gfpyexec_PYTHON = _getfem$(PYTHON_SO)

EXTRA_DIST = getfem_python.c

# $(warning PSEUDO_FUNCTIONS= $(PSEUDO_FUNCTIONS))

CLEANFILES = getfem.py _getfem*.so *.o getfem_python_c.o getfem_python_c.c getfem.pyc __pycache__/*

if BUILDPYTHONPAR

pythoncom=python-com-par

else

pythoncom=python-com

endif

if BUILDPYTHON

getfem.py: ../libgetfemint.la $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc $(srcdir)/.. $(pythoncom) > getfem.py || ( rm getfem.py ; /bin/false )


PYTHON_CC=$(CC)

getfem_python_c.c : getfem_python.c
	cp $(srcdir)/getfem_python.c getfem_python_c.c

# ARCHFLAGS is set to empty to disable universal binaries with python 2.5 on macos 10.5
_getfem$(PYTHON_SO): getfem_python_c.c ../libgetfemint.la $(GETFEM_LIB_LA)
	touch _getfem$(PYTHON_SO) && rm _getfem$(PYTHON_SO)
	ARCHFLAGS="" CC="$(CC)" CFLAGS="$(CFLAGS)" $(PYTHON) setup.py -v build @PYTHON_CC_ARG@ --build-temp . --build-base . --build-lib . --force
#LDSHARED="$(CXX) -shared" ARCHFLAGS="" python setup.py -v build --build-temp . --build-base . --build-lib . --force


all: _getfem$(PYTHON_SO) getfem.py
#pyexec_LTLIBRARIES = libgfpython.la
#libgfpython_la_LIBADD = ../.libs/libgetfemint.a @GETFEM_STATICLIBS@
#libgfpython_la_SOURCES = \
#	getfem_python.c

endif BUILDPYTHON
