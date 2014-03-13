# SUBDIRS = 

check_PROGRAMS =                 \
	bimaterial_crack_test

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = bimaterial_crack_test.pl

EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
