# SUBDIRS = 

check_PROGRAMS = static_contact_gears static_contact_gears_u1_u2

CLEANFILES = 

static_contact_gears_SOURCES = static_contact_gears.cc
static_contact_gears_u1_u2_SOURCES = static_contact_gears_u1_u2.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = 

EXTRA_DIST = 


TESTS_ENVIRONMENT = perl
