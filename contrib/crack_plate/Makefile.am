# SUBDIRS = 

if QHULL
check_PROGRAMS =                 \
	crack_mindlin crack_bilaplacian
else
check_PROGRAMS =
endif

CLEANFILES = 

if QHULL
crack_mindlin_SOURCES = crack_mindlin.cc
crack_bilaplacian_SOURCES = crack_bilaplacian.cc
endif

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

if QHULL
TESTS = $(top_srcdir)/contrib/crack_mindlin/crack_mindlin.pl
else
TESTS =
endif

EXTRA_DIST = \
	crack_mindlin.pl                  \
	crack_mindlin.param \
	crack_bilaplacian.param


TESTS_ENVIRONMENT = perl
