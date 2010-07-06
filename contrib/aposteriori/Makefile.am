# SUBDIRS = 

check_PROGRAMS = aposteriori aposteriori_laplacian

CLEANFILES = 

aposteriori_SOURCES = aposteriori.cc
aposteriori_laplacian_SOURCES = aposteriori_laplacian.cc
SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/aposteriori/aposteriori.pl \
        $(top_srcdir)/contrib/aposteriori/aposteriori_laplacian.pl 

EXTRA_DIST = \
	aposteriori_laplacian.pl      	 \
	aposteriori_laplacian.param      \
	aposteriori.pl                   \
	aposteriori.param


TESTS_ENVIRONMENT = perl
