# SUBDIRS = 

if QHULL
check_PROGRAMS = crack_mindlin crack_bilaplacian mortar_bilaplacian 
else
check_PROGRAMS =
endif

CLEANFILES = 

if QHULL
crack_mindlin_SOURCES = crack_mindlin.cc
mortar_bilaplacian_SOURCES = mortar_bilaplacian.cc mortar_bilaplacian.h
crack_bilaplacian_SOURCES = crack_bilaplacian.cc crack_bilaplacian_singularities.cc crack_bilaplacian_problem.cc crack_bilaplacian.h crack_bilaplacian_moment.cc crack_bilaplacian_tools.cc crack_bilaplacian_sif.cc
endif

SUPLDFLAGS = @SUPLDFLAGS@
INCLUDES = -I$(top_srcdir)/src -I../../src
LDADD    = ../../src/libgetfem.la -lm $(SUPLDFLAGS)

if QHULL
TESTS = $(top_srcdir)/contrib/crack_plate/crack_mindlin.pl
else
TESTS =
endif

EXTRA_DIST =		 		\
	crack_mindlin.pl		\
	crack_mindlin.param 		\
	crack_bilaplacian.param 	\
	mortar_bilaplacian.param 	\
	crack_bilaplacian_mixed.param

#	serie.param \
#	serie_bis.param \
#

TESTS_ENVIRONMENT = perl
