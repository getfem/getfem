# SUBDIRS = 

check_PROGRAMS = xfem_stab_unilat_contact

CLEANFILES = 

xfem_stab_unilat_contact_SOURCES = xfem_stab_unilat_contact.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = xfem_stab_unilat_contact.pl

EXTRA_DIST = 				\
	xfem_stab_unilat_contact.param	\
	xfem_stab_unilat_contact.pl               

LOG_COMPILER = perl
