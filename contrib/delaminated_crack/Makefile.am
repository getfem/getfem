# SUBDIRS = 

check_PROGRAMS =                 \
	delaminated_crack

CLEANFILES = 

delaminated_crack_SOURCES = delaminated_crack.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src 
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = delaminated_crack.pl

EXTRA_DIST = \
	delaminated_crack.pl                  \
	delaminated_crack.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl