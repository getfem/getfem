# SUBDIRS = 

check_PROGRAMS =                 \
	bimaterial_crack_test

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc
SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/bimaterial_crack_test/bimaterial_crack_test.pl \
        $(top_srcdir)/contrib/bimaterial_crack_test/crack.pl
EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param


TESTS_ENVIRONMENT = perl
