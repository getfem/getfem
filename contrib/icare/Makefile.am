# SUBDIRS = 

check_PROGRAMS =                 \
	icare

CLEANFILES = icare.data

icare_SOURCES = icare.cc icare.h


INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = $(top_srcdir)/contrib/icare/icare.pl

EXTRA_DIST = \
	icare.pl                  \
	icare.param


TESTS_ENVIRONMENT = perl
