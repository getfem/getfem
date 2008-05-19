# SUBDIRS = 

check_PROGRAMS =  contact  

CLEANFILES = 

contact_SOURCES = contact.cc
SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/contact_grd_trans/contact.pl
EXTRA_DIST = \
	contact.pl                  \
	contact.param 


TESTS_ENVIRONMENT = perl
