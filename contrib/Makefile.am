bin_PROGRAMS = \
	compare_solutions \
	interpole_solution
compare_solutions_SOURCES = compare_solutions.C
interpole_solution_SOURCES = interpole_solution.C


CXXFLAGS = @APPSCXXFLAGS@

INCLUDES = -I$(top_srcdir)/src 
LDADD    = $(top_srcdir)/src/.libs/libgetfem.a -lm -L/usr/X11R6/lib -lX11


clean-local:
	rm -rf cxx_repository
