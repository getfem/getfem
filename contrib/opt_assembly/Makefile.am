check_PROGRAMS =  opt_assembly

CLEANFILES =

opt_assembly_SOURCES = opt_assembly.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS =  opt_assembly.pl

EXTRA_DIST = \
	opt_assembly.pl

LOG_COMPILER = perl
