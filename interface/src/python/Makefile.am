GETFEM_LIB_LA = ../../../src/libgetfem.la

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I@srcdir@ -I@srcdir@/..

gfpythondir=$(pythondir)/getfem
gfpyexecdir=$(pyexecdir)/getfem

gfpython_PYTHON = getfem.py __init__.py
nodist_gfpyexec_PYTHON = _getfem.so

EXTRA_DIST = getfem_python.c

# $(warning PSEUDO_FUNCTIONS= $(PSEUDO_FUNCTIONS))

getfem.py: @PSEUDO_FUNCTIONS@ $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc $(srcdir)/.. python-com > getfem.py || ( rm getfem.py ; /bin/false )

CLEANFILES = getfem.py _getfem.so *.o getfem_python_c.o getfem_python_c.c

if BUILDPYTHON

PYTHON_CC=$(CC)

getfem_python_c.c : getfem_python.c
	cp $(srcdir)/getfem_python.c getfem_python_c.c

# ARCHFLAGS is set to empty to disable universal binaries with python 2.5 on macos 10.5
_getfem.so: getfem_python_c.c ../libgetfemint.la $(GETFEM_LIB_LA)
	touch _getfem.so && rm _getfem.so
	ARCHFLAGS="" python setup.py -v build --build-temp . --build-base . --build-lib . --force
#LDSHARED="$(CXX) -shared" ARCHFLAGS="" python setup.py -v build --build-temp . --build-base . --build-lib . --force


# getfem_python_reference.html: getfem.py _getfem.so
#	cp getfem.py getfem_python_reference.py
#	(export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):../../../src/.libs && pydoc -w getfem_python_reference) &&  if test -d $(top_srcdir)/interface/doc; then cp getfem_python_reference.html $(top_srcdir)/interface/doc/getfem_python_reference.html; fi;
#	rm -f getfem_python_reference.py

all: _getfem.so getfem.py
#pyexec_LTLIBRARIES = libgfpython.la
#libgfpython_la_LIBADD = ../.libs/libgetfemint.a @GETFEM_STATICLIBS@
#libgfpython_la_SOURCES = \
#	getfem_python.c

endif BUILDPYTHON
