# SUBDIRS = 

check_PROGRAMS =                 \
	static_friction

CLEANFILES = normal_stress tangential_stress

static_friction_SOURCES = static_friction.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = $(abs_top_srcdir)/contrib/static_friction/static_friction.pl

EXTRA_DIST = \
	static_friction.pl                  \
	static_friction.param

LOG_COMPILER = perl
