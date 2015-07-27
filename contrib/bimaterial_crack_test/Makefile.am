# SUBDIRS = 

check_PROGRAMS =                 \
	bimaterial_crack_test    \
	crack

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc getfem_Xfem.cc
crack_SOURCES = crack.cc crack_exact_solution.cc crack_exact_solution.h getfem_Xfem.cc getfem_spider_fem.h getfem_Xfem.h


AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = bimaterial_crack_test.pl \
        crack.pl

EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param          \
	crack.pl                  \
	crack.param 

LOG_COMPILER = perl
