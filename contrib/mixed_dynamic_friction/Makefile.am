# SUBDIRS = 

check_PROGRAMS =  mixed_dynamic_friction  mixed_scalar_hyperbolic

CLEANFILES = 

mixed_dynamic_friction_SOURCES = mixed_dynamic_friction.cc
mixed_scalar_hyperbolic_SOURCES = mixed_scalar_hyperbolic.cc

SUPLDFLAGS = @SUPLDFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

TESTS = $(abs_top_srcdir)/contrib/mixed_dynamic_friction/mixed_dynamic_friction.pl
EXTRA_DIST = \
	mixed_dynamic_friction.pl                  \
	mixed_dynamic_friction.param               \
	mixed_scalar_hyperbolic.param              \
	mixed_scalar_hyperbolic.net

TESTS_ENVIRONMENT = perl
