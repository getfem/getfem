bin_PROGRAMS = \
	getfem_import_export \
	compare_solutions \
	interpole_solution \
	draw_mesh_doedges

bin_SCRIPTS = draw_mesh

getfem_import_export_SOURCES = getfem_import_export.C
compare_solutions_SOURCES = compare_solutions.C
interpole_solution_SOURCES = interpole_solution.C
draw_mesh_doedges_SOURCES = draw_mesh_doedges.C

EXTRA_DIST = \
	$(bin_SCRIPTS) \
	draw_mesh_example1.mesh \
	draw_mesh_example2.mesh \
	draw_mesh_example3.mesh \
	draw_mesh_example4.mesh \
	draw_mesh_example5.mesh \
	draw_mesh_example6.mesh

INCLUDES = -I$(top_srcdir)/src -I../src
LDADD    = ../src/.libs/libgetfem.so -lm @SUPLDFLAGS@

