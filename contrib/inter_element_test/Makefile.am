# SUBDIRS = 

check_PROGRAMS = inter_element

CLEANFILES = 

inter_element_SOURCES = inter_element.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/inter_element_test/inter_element.pl

EXTRA_DIST = inter_element.pl

TESTS_ENVIRONMENT = perl
