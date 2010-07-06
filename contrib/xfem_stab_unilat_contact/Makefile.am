# SUBDIRS = 

check_PROGRAMS = xfem_stab_unilat_contact

CLEANFILES = 

xfem_stab_unilat_contact_SOURCES = xfem_stab_unilat_contact.cc


SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/xfem_stab_unilat_contact/xfem_stab_unilat_contact.pl

EXTRA_DIST = \
	xfem_stab_unilat_contact.param               

TESTS_ENVIRONMENT = perl
