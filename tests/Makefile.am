check_PROGRAMS = \
	dynamic_array \
	dynamic_tas \
	laplacien \
	elasto_statique
dynamic_array_SOURCES = dynamic_array.C 
dynamic_tas_SOURCES = dynamic_tas.C 
laplacien_SOURCES = laplacien.C
elasto_statique_SOURCES = elasto_statique.C

INCLUDES = -I$(top_srcdir)/src
LDADD    = ../src/.libs/libgetfem.a -lm
CXXFLAGS = @APPSCXXFLAGS@

TESTS             = check.pl
TESTS_ENVIRONMENT = perl


clean-local:
	rm -rf cxx_repository
