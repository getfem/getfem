# SUBDIRS = 

check_PROGRAMS =  contact  contact_continuation_load  contact_continuation_time

CLEANFILES = 

contact_SOURCES = contact.cc
contact_continuation_load_SOURCES = contact_continuation_load.cc
contact_continuation_time_SOURCES = contact_continuation_time.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/contact_grd_trans/contact.pl
EXTRA_DIST = \
	contact.pl                  \
	contact.param 
	contact_continuation_load.param 
	contact_continuation_time.param 



TESTS_ENVIRONMENT = perl
