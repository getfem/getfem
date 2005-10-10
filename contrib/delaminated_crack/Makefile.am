# SUBDIRS = 

check_PROGRAMS =                 \
	delaminated_crack

CLEANFILES = 

delaminated_crack_SOURCES = delaminated_crack.cc

INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = $(top_srcdir)/contrib/icare/delaminated_crack.pl

EXTRA_DIST = \
	delaminated_crack.pl                  \
	delaminated_crack.param


TESTS_ENVIRONMENT = perl
