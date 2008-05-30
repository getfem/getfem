# SUBDIRS = 

check_PROGRAMS =  mixed_dynamic_friction  mixed_scalar_hyperbolic

CLEANFILES = 

mixed_dynamic_friction_SOURCES = mixed_dynamic_friction.cc
mixed_scalar_hyperbolic_SOURCES = mixed_scalar_hyperbolic.cc

SUPLDFLAGS = @SUPLDFLAGS@
MUMPS_LIBS = @MUMPS_LIBS@
MUMPS_CFLAGS = @MUMPS_CFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src $(MUMPS_CFLAGS)
LDADD    = ../../src/libgetfem.la -lm $(MUMPS_LIBS) $(SUPLDFLAGS)

TESTS = $(top_srcdir)/contrib/mixed_dynamic_friction/mixed_dynamic_friction.pl
EXTRA_DIST = \
	mixed_dynamic_friction.pl                  \
	mixed_dynamic_friction.param               \
	mixed_scalar_hyperbolic.param

TESTS_ENVIRONMENT = perl
