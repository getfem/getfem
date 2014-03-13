# SUBDIRS = 

check_PROGRAMS =                 \
	mixed_elastostatic

CLEANFILES = 

mixed_elastostatic_SOURCES = mixed_elastostatic.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = mixed_elastostatic.pl

EXTRA_DIST = \
	mixed_elastostatic.pl                  \
	mixed_elastostatic.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
