# SUBDIRS = 

check_PROGRAMS = xfem_large_strain nonlinear_incomp_xfem vertex_large_strain linear_incomp_xfem
CLEANFILES = 

xfem_large_strain_SOURCES = xfem_large_strain.cc
linear_incomp_xfem_SOURCES = linear_incomp_xfem.cc
nonlinear_incomp_xfem_SOURCES = nonlinear_incomp_xfem.cc
vertex_large_strain_SOURCES = vertex_large_strain.cc \
                              getfem_nonlinear_elastoptim.h \
                              compressible_getfem_nonlinear_elastoptim.h

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = xfem_large_strain.pl

EXTRA_DIST = \
	xfem_large_strain.pl            \
	vertex_large_strain.param       \
	xfem_large_strain.param  	\
	nonlinear_incomp_xfem.param        \
	vertex_large_strain.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
