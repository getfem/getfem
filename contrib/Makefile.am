bin_PROGRAMS = \
	compare_solutions \
	interpole_solution
compare_solutions_SOURCES = compare_solutions.C
interpole_solution_SOURCES = interpole_solution.C


GETFEMPATH = @GETFEMPATH@
CXXFLAGS = @APPSCXXFLAGS@

# a arranger, pour MTL et misc
INCLUDES = -I$(top_srcdir)/misc -I$(top_srcdir)/../$(GETFEMPATH)/src -I/opt/mtl/include -I/usr/local/mtl/include
LDADD    = $(top_srcdir)/../$(GETFEMPATH)/src/.libs/libgetfem.a $(top_srcdir)/misc/.libs/libgetfemmisc.a -lm -L/usr/X11R6/lib -lX11


clean-local:
	rm -rf cxx_repository
