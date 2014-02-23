# SUBDIRS = 

check_PROGRAMS =  contact contact_continuation_load contact_continuation_time

CLEANFILES = 

contact_SOURCES = contact.cc
contact_continuation_load_SOURCES = contact_continuation_load.cc
contact_continuation_time_SOURCES = contact_continuation_time.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = $(abs_top_srcdir)/contrib/contact_grd_trans/contact.pl

EXTRA_DIST = \
	contact.pl                  \
	contact.param 
	contact_continuation_load.param 
	contact_continuation_time.param 

LOG_COMPILER = perl
