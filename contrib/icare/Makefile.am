# SUBDIRS = 

check_PROGRAMS =                 \
	icare

CLEANFILES = icare.data

icare_SOURCES = icare.cc icare.h

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/icare/icare.pl

EXTRA_DIST = \
	icare.pl                  \
	icare.param


TESTS_ENVIRONMENT = perl
