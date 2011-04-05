# SUBDIRS = 

check_PROGRAMS = xfem_contact xfem_dirichlet xfem_stokes

CLEANFILES = 

xfem_contact_SOURCES = xfem_contact.cc
xfem_dirichlet_SOURCES = xfem_dirichlet.cc
xfem_stokes_SOURCES = xfem_stokes.cc

SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/xfem_contact/xfem_contact.pl

EXTRA_DIST = \
	xfem_contact.pl                  \
	xfem_contact.param               \
	xfem_dirichlet.param             \
	xfem_stokes.param


TESTS_ENVIRONMENT = perl
