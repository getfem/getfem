# SUBDIRS = 

check_PROGRAMS = inter_element

CLEANFILES = 

inter_element_SOURCES = inter_element.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = inter_element.pl

EXTRA_DIST = inter_element.pl

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
