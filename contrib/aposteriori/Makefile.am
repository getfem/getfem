# SUBDIRS = 

check_PROGRAMS = aposteriori aposteriori_laplacian

CLEANFILES = 

aposteriori_SOURCES = aposteriori.cc
aposteriori_laplacian_SOURCES = aposteriori_laplacian.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/aposteriori/aposteriori.pl \
        $(top_srcdir)/contrib/aposteriori/aposteriori_laplacian.pl 

EXTRA_DIST = \
	aposteriori_laplacian.pl      	 \
	aposteriori_laplacian.param      \
	aposteriori.pl                   \
	aposteriori.param


TESTS_ENVIRONMENT = perl
