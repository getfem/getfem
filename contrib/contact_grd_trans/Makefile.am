# SUBDIRS = 

check_PROGRAMS =  contact contact_continuation_load contact_continuation_time

CLEANFILES = 

contact_SOURCES = contact.cc
contact_continuation_load_SOURCES = contact_continuation_load.cc
contact_continuation_time_SOURCES = contact_continuation_time.cc my_getfem_Coulomb_friction.h
SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/contact_grd_trans/contact.pl
EXTRA_DIST = \
	contact.pl                  \
	contact.param 
	contact_continuation_load.param 
	contact_continuation_time.param 



TESTS_ENVIRONMENT = perl
