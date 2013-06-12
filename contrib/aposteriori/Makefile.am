# SUBDIRS = 

check_PROGRAMS = aposteriori aposteriori_laplacian

CLEANFILES = 

aposteriori_SOURCES = aposteriori.cc
aposteriori_laplacian_SOURCES = aposteriori_laplacian.cc
SUPLDFLAGS = @SUPLDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(abs_top_srcdir)/contrib/aposteriori/aposteriori.pl \
        $(abs_top_srcdir)/contrib/aposteriori/aposteriori_laplacian.pl 

EXTRA_DIST = \
	aposteriori_laplacian.pl      	 \
	aposteriori_laplacian.param      \
	aposteriori.pl                   \
	aposteriori.param


TESTS_ENVIRONMENT = perl
