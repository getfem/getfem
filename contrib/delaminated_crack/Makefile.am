# SUBDIRS = 

check_PROGRAMS =                 \
	delaminated_crack

CLEANFILES = 

delaminated_crack_SOURCES = delaminated_crack.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/delaminated_crack/delaminated_crack.pl

EXTRA_DIST = \
	delaminated_crack.pl                  \
	delaminated_crack.param


TESTS_ENVIRONMENT = perl
