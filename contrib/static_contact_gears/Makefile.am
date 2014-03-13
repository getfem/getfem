# SUBDIRS = 

check_PROGRAMS = static_contact_gears static_contact_gears_u1_u2

CLEANFILES = 

static_contact_gears_SOURCES = static_contact_gears.cc
static_contact_gears_u1_u2_SOURCES = static_contact_gears_u1_u2.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = 

EXTRA_DIST = \
	gear1.msh \
	gear1_2D.msh \
	gear1_2teeth.msh \
	gear2.msh \
	gear2_2D.msh \
	gear2_2teeth.msh \
	static_contact_planetary_1.msh \
	static_contact_planetary_2.msh \
	static_contact_planetary_3.msh \
	static_contact_planetary_4.msh \
	static_contact_planetary_5.msh \
	static_contact_planetary.py \
	static_contact_gears.param \
	static_contact_gears_2D.param \
	static_contact_gears_2teeth.param

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
