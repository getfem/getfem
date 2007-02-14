# SUBDIRS = 

check_PROGRAMS = aposteriori

CLEANFILES = 

aposteriori_SOURCES = aposteriori.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/aposteriori/aposteriori.pl 

EXTRA_DIST = \
	aposteriori.pl                  \
	aposteriori.param


TESTS_ENVIRONMENT = perl
