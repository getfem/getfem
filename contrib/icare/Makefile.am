# SUBDIRS = 

check_PROGRAMS =                 \
	icare

CLEANFILES = icare.data

icare_SOURCES = icare.cc icare.h

MUMPS_LIBS = @MUMPS_LIBS@
# MUMPS_CFLAGS = @MUMPS_CFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src # $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = icare.pl

EXTRA_DIST = \
	icare.pl                  \
	icare.param		\
	navier_stokes_cylinder1.mesh \
	navier_stokes_cylinder2.mesh

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
