# SUBDIRS = 

check_PROGRAMS = xfem_contact xfem_dirichlet

CLEANFILES = 

xfem_contact_SOURCES = xfem_contact.cc
xfem_dirichlet_SOURCES = xfem_dirichlet.cc

SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/xfem_contact/xfem_contact.pl

EXTRA_DIST = \
	xfem_contact.pl                  \
	xfem_contact.param               \
	xfem_dirichlet.param


TESTS_ENVIRONMENT = perl
