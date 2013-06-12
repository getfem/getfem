# SUBDIRS = 

check_PROGRAMS =                 \
	mixed_elastostatic

CLEANFILES = 

mixed_elastostatic_SOURCES = mixed_elastostatic.cc

SUPLDFLAGS = @SUPLDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm  $(SUPLDFLAGS)

TESTS = $(abs_top_srcdir)/contrib/mixed_elastostatic/mixed_elastostatic.pl

EXTRA_DIST = \
	mixed_elastostatic.pl                  \
	mixed_elastostatic.param


TESTS_ENVIRONMENT = perl
