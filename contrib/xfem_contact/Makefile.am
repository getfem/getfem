if METIS
optprogs = xfem_dirichlet
else
optprogs = 
endif

check_PROGRAMS = xfem_contact xfem_stokes $(optprogs)

CLEANFILES = 

xfem_contact_SOURCES = xfem_contact.cc
if METIS
xfem_dirichlet_SOURCES = xfem_dirichlet.cc
endif
xfem_stokes_SOURCES = xfem_stokes.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = xfem_contact.pl

EXTRA_DIST = \
	xfem_contact.pl                  \
	xfem_contact.param               \
	xfem_dirichlet.param             \
	xfem_stokes.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
