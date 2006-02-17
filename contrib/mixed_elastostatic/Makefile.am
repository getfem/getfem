# SUBDIRS = 

check_PROGRAMS =                 \
	mixed_elastostatic

CLEANFILES = 

mixed_elastostatic_SOURCES = mixed_elastostatic.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/mixed_elastostatic/mixed_elastostatic.pl

EXTRA_DIST = \
	mixed_elastostatic.pl                  \
	mixed_elastostatic.param


TESTS_ENVIRONMENT = perl
