# SUBDIRS = 

check_PROGRAMS = signorini_schwarz

CLEANFILES = 

signorini_schwarz_SOURCES = signorini_schwarz.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/signorini_schwarz/signorini_schwarz.pl

EXTRA_DIST = \
	signorini_schwarz.param      \
	signorini_schwarz.pl          


TESTS_ENVIRONMENT = perl
