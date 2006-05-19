# SUBDIRS = 

check_PROGRAMS =                 \
	crack_mindlin crack_bilaplacian

CLEANFILES = 

crack_mindlin_SOURCES = crack_mindlin.cc

crack_bilaplacian_SOURCES = crack_bilaplacian.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/crack_mindlin/crack_mindlin.pl

EXTRA_DIST = \
	crack_mindlin.pl                  \
	crack_mindlin.param \
	crack_bilaplacian.param


TESTS_ENVIRONMENT = perl
