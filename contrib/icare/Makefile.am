# SUBDIRS = 

check_PROGRAMS =                 \
	icare

CLEANFILES = icare.data

icare_SOURCES = icare.cc icare.h

SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm  $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/icare/icare.pl

EXTRA_DIST = \
	icare.pl                  \
	icare.param		\
	navier_stokes_cylinder1.mesh \
	navier_stokes_cylinder2.mesh

TESTS_ENVIRONMENT = perl
