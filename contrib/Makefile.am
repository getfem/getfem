bin_PROGRAMS = \
	compare_solutions \
	interpole_solution
compare_solutions_SOURCES = compare_solutions.C
interpole_solution_SOURCES = interpole_solution.C


CXXFLAGS = @APPSCXXFLAGS@

INCLUDES = -I$(top_srcdir)/src 
LDADD    = $(top_srcdir)/src/.libs/libgetfem.so -lm


clean-local:
	rm -rf cxx_repository
