# SUBDIRS = 

check_PROGRAMS =                 \
        test_contact

CLEANFILES = 

test_contact_SOURCES = contact_problem.ccp test_contact.cpp

INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = contact_problem.o ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = $(top_srcdir)/contrib/level_set_contact/contact_problem.pl

EXTRA_DIST = \
	contact_problem.pl              \
	contact_2D.param		\
	contact_3D.param


TESTS_ENVIRONMENT = perl
