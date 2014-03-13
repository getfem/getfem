# SUBDIRS = 

check_PROGRAMS =                 \
        test_contact

CLEANFILES = 

test_contact_SOURCES = contact_problem.cpp test_contact.cpp contact_problem.h

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = contact_problem.pl

EXTRA_DIST = \
	contact_problem.pl              \
	contact_2D.param		\
	contact_3D.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
