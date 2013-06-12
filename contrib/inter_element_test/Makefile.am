# SUBDIRS = 

check_PROGRAMS = inter_element

CLEANFILES = 

inter_element_SOURCES = inter_element.cc

SUPLDFLAGS = @SUPLDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm  $(SUPLDFLAGS)

TESTS = $(abs_top_srcdir)/contrib/inter_element_test/inter_element.pl

EXTRA_DIST = inter_element.pl

TESTS_ENVIRONMENT = perl
