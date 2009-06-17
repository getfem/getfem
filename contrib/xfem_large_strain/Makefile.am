# SUBDIRS = 

check_PROGRAMS = xfem_large_strain linear_incomp_xfem

CLEANFILES = 

xfem_large_strain_SOURCES = xfem_large_strain.cc
linear_incomp_xfem_SOURCES = linear_incomp_xfem.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/xfem_large_strain/xfem_large_strain.pl
EXTRA_DIST = \
	xfem_large_strain.pl                  \
	xfem_large_strain.param  	      \
	linear_incomp_xfem.param

TESTS_ENVIRONMENT = perl
