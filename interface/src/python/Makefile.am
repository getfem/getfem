GETFEM_LIB_LA = ../../../src/libgetfem.la

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
	$(srcdir)/../gf_eltm.cc \
	$(srcdir)/../gf_mesh.cc \
	$(srcdir)/../gf_mesh_get.cc \
	$(srcdir)/../gf_mesh_set.cc \
	$(srcdir)/../gf_slice.cc \
	$(srcdir)/../gf_slice_get.cc \
	$(srcdir)/../gf_slice_set.cc \
	$(srcdir)/../gf_mesh_fem.cc \
	$(srcdir)/../gf_mesh_fem_get.cc \
	$(srcdir)/../gf_mesh_fem_set.cc \
	$(srcdir)/../gf_mesh_im.cc \
	$(srcdir)/../gf_mesh_im_get.cc \
	$(srcdir)/../gf_mesh_im_set.cc \
	$(srcdir)/../gf_precond.cc \
	$(srcdir)/../gf_precond_get.cc \
	$(srcdir)/../gf_spmat.cc \
	$(srcdir)/../gf_spmat_get.cc \
	$(srcdir)/../gf_spmat_set.cc \
	$(srcdir)/../gf_linsolve.cc \
	$(srcdir)/../gf_asm.cc \
	$(srcdir)/../gf_fem.cc \
	$(srcdir)/../gf_fem_get.cc \
	$(srcdir)/../gf_integ.cc \
	$(srcdir)/../gf_integ_get.cc \
	$(srcdir)/../gf_workspace.cc \
	$(srcdir)/../gf_delete.cc \
	$(srcdir)/../gf_mdbrick.cc \
	$(srcdir)/../gf_mdbrick_get.cc \
	$(srcdir)/../gf_mdbrick_set.cc \
	$(srcdir)/../gf_mdstate.cc \
	$(srcdir)/../gf_mdstate_get.cc \
	$(srcdir)/../gf_mdstate_set.cc \
	$(srcdir)/../gf_model.cc \
	$(srcdir)/../gf_model_get.cc \
	$(srcdir)/../gf_model_set.cc \
	$(srcdir)/../gf_levelset.cc \
	$(srcdir)/../gf_levelset_get.cc \
	$(srcdir)/../gf_levelset_set.cc \
	$(srcdir)/../gf_mesh_levelset.cc \
	$(srcdir)/../gf_mesh_levelset_get.cc \
	$(srcdir)/../gf_mesh_levelset_set.cc

getfem.py: getfem.base.py $(PSEUDO_MFUNCTIONS) ../../bin/extract_pydoc
	$(top_srcdir)/interface/bin/extract_pydoc python $(srcdir)/.. > getfem.py || ( rm getfem.py ; /bin/false )

clean: clean-am
	rm -f getfem.py _getfem.so getfem.html
	rm -f $(top_srcdir)/interface/doc/getfem_python_reference.html
	rm -f *.o

if BUILDPYTHON

PYTHON_CC=$(CC)


# ARCHFLAGS is set to empty to disable universal binaries with python 2.5 on macos 10.5
_getfem.so: getfem_python.c ../libgetfemint.la $(GETFEM_LIB_LA)
	touch _getfem.so && rm _getfem.so
	ARCHFLAGS="" python setup.py -v build --build-temp . --build-base . --build-lib . --force
#LDSHARED="$(CXX) -shared" ARCHFLAGS="" python setup.py -v build --build-temp . --build-base . --build-lib . --force


getfem_python_reference.html: getfem.py _getfem.so
	cp getfem.py getfem_python_reference.py
	(export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):../../../src/.libs && pydoc -w getfem_python_reference) &&  if test -d $(top_srcdir)/interface/doc; then cp getfem_python_reference.html $(top_srcdir)/interface/doc/getfem_python_reference.html; fi;
	rm -f getfem_python_reference.py

all: _getfem.so getfem.py getfem_python_reference.html
#pyexec_LTLIBRARIES = libgfpython.la
#libgfpython_la_LIBADD = ../.libs/libgetfemint.a @GETFEM_STATICLIBS@
#libgfpython_la_SOURCES = \
#	getfem_python.c

endif BUILDPYTHON
