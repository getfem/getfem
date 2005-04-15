SUBDIRS = icare 

bin_PROGRAMS = \
	draw_mesh_doedges

bin_SCRIPTS = draw_mesh

draw_mesh_doedges_SOURCES = draw_mesh_doedges.cc

EXTRA_DIST = \
	$(bin_SCRIPTS) \
	draw_mesh_example1.mesh \
	draw_mesh_example2.mesh \
	draw_mesh_example3.mesh \
	draw_mesh_example4.mesh \
	draw_mesh_example5.mesh \
	draw_mesh_example6.mesh

INCLUDES = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@

CLEANFILES = ii_files/*