# SUBDIRS = 

check_PROGRAMS =                 \
	navier_stokes

CLEANFILES = navier_stokes.data

navier_stokes_SOURCES = navier_stokes.cc


INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = $(top_srcdir)/tests/navier_stokes.pl

EXTRA_DIST = \
	navier_stokes.pl                  \
	navier_stokes.param


TESTS_ENVIRONMENT = perl
