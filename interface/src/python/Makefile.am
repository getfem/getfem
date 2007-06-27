if BUILDPYTHON

GETFEM_LIB_LA = ../../../src/libgetfem.la @LIBS@

INCLUDES = $(PYTHON_CPPFLAGS) -I@srcdir@ -I@srcdir@/..

gfpythondir=$(pythondir)/getfem
gfpyexecdir=$(pyexecdir)/getfem

gfpython_PYTHON = getfem.py __init__.py
nodist_gfpyexec_PYTHON = _getfem.so

EXTRA_DIST = getfem.base.py getfem.py getfem_python.c

PSEUDO_MFUNCTIONS = \
	$(srcdir)/../gf_util.cc \
	$(srcdir)/../gf_cvstruct_get.cc \
	$(srcdir)/../gf_geotrans.cc \
	$(srcdir)/../gf_geotrans_get.cc \
	$(srcdir)/../gf_compute.cc \
	$(srcdir)/../gf_mesh_fem.cc \
	$(srcdir)/../gf_mesh_fem_set.cc \
	$(srcdir)/../gf_mesh_fem_get.cc \
	$(srcdir)/../gf_mesh_im.cc \
	$(srcdir)/../gf_mesh_im_set.cc \
	$(srcdir)/../gf_mesh_im_get.cc \
	$(srcdir)/../gf_eltm.cc \
	$(srcdir)/../gf_mesh.cc \
	$(srcdir)/../gf_mesh_set.cc \
	$(srcdir)/../gf_mesh_get.cc \
	$(srcdir)/../gf_slice.cc \
	$(srcdir)/../gf_slice_get.cc \
	$(srcdir)/../gf_slice_set.cc \
	$(srcdir)/../gf_precond.cc \
	$(srcdir)/../gf_precond_get.cc \
	$(srcdir)/../gf_spmat.cc \
	$(srcdir)/../gf_spmat_set.cc \
	$(srcdir)/../gf_spmat_get.cc \
	$(srcdir)/../gf_linsolve.cc \
	$(srcdir)/../gf_asm.cc \
	$(srcdir)/../gf_fem.cc \
	$(srcdir)/../gf_fem_get.cc \
	$(srcdir)/../gf_integ.cc \
	$(srcdir)/../gf_integ_get.cc \
	$(srcdir)/../gf_workspace.cc \
	$(srcdir)/../gf_delete.cc

getfem.py: getfem.base.py $(PSEUDO_MFUNCTIONS)
	$(top_srcdir)/interface/bin/extract_pydoc python $(srcdir)/.. > getfem.py || ( rm getfem.py ; /bin/false )

PYTHON_CC=$(CC)

_getfem.so: getfem_python.c ../libgetfemint.la $(GETFEM_LIB_LA)
	touch _getfem.so && rm _getfem.so
	LDSHARED="$(CXX) -shared" python setup.py -v build --build-temp . --build-base . --build-lib . --force


doc: getfem.py
	pydoc -w getfem && mv getfem.html $(top_srcdir)/interface/doc/getfem_python_reference.html

all: _getfem.so getfem.py
#pyexec_LTLIBRARIES = libgfpython.la
#libgfpython_la_LIBADD = ../.libs/libgetfemint.a @GETFEM_STATICLIBS@
#libgfpython_la_SOURCES = \
#	getfem_python.c

endif BUILDPYTHON