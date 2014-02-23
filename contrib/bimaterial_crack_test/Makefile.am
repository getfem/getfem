# SUBDIRS = 

check_PROGRAMS =                 \
	bimaterial_crack_test    \
	crack

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc
crack_SOURCES = crack.cc crack_exact_solution.cc crack_exact_solution.h

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = $(abs_top_srcdir)/contrib/bimaterial_crack_test/bimaterial_crack_test.pl \
        $(abs_top_srcdir)/contrib/bimaterial_crack_test/crack.pl

EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param          \
	crack.pl                  \
	crack.param 

LOG_COMPILER = perl
