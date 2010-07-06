# SUBDIRS = 

check_PROGRAMS =  contact contact_continuation

CLEANFILES = 

contact_SOURCES = contact.cc
contact_continuation_SOURCES = contact_continuation.cc my_getfem_Coulomb_friction.h
SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/contact_grd_trans/contact.pl
EXTRA_DIST = \
	contact.pl                  \
	contact.param 
	contact_continuation.param



TESTS_ENVIRONMENT = perl
