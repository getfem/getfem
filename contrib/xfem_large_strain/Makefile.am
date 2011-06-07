# SUBDIRS = 

check_PROGRAMS = xfem_large_strain nonlinear_incomp_xfem vertex_large_strain test_xfem_large_strain 
CLEANFILES = 

xfem_large_strain_SOURCES = xfem_large_strain.cc
nonlinear_incomp_xfem_SOURCES = nonlinear_incomp_xfem.cc
test_xfem_large_strain_SOURCES = test_xfem_large_strain.cc
vertex_large_strain_SOURCES = vertex_large_strain.cc \
                              getfem_nonlinear_elastoptim.h \
                              compressible_getfem_nonlinear_elastoptim.h
SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/xfem_large_strain/xfem_large_strain.pl
EXTRA_DIST = \
	xfem_large_strain.pl            \
	vertex_large_strain.param       \
	xfem_large_strain.param  	\
	nonlinear_incomp_xfem.param        \
	vertex_large_strain.param

TESTS_ENVIRONMENT = perl
