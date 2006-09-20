# SUBDIRS = 

check_PROGRAMS =                 \
	bimaterial_crack_test		\
	crack

CLEANFILES = 

bimaterial_crack_test_SOURCES = bimaterial_crack_test.cc
crack_SOURCES = crack.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/bimaterial_crack_test/bimaterial_crack_test.pl
TESTS = $(top_srcdir)/contrib/bimaterial_crack_test/crack.pl
EXTRA_DIST = \
	bimaterial_crack_test.pl                  \
	bimaterial_crack_test.param          \
	crack.pl                  \
	crack.param 


TESTS_ENVIRONMENT = perl
