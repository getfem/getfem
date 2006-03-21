# SUBDIRS = 

check_PROGRAMS =                 \
	static_friction

CLEANFILES = normal_stress tangential_stress

static_friction_SOURCES = static_friction.cc

INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = $(top_srcdir)/contrib/static_friction/static_friction.pl

EXTRA_DIST = \
	static_friction.pl                  \
	static_friction.param


TESTS_ENVIRONMENT = perl
