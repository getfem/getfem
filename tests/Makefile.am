check_PROGRAMS = \
	dynamic_array \
	dynamic_tas \
	elasto_statique  \
	geo_trans \
	geo_trans_inv \
	integration \
	laplacian \
	linalg \
	poly \
	speed_test_for_gmm \
	test_assembly \
	test_int_set \
	test_link_fem \
	test_mat_elem \
	test_mesh \
	test_tree_sorted

CLEANFILES = \
	laplacian.res laplacian.mesh laplacian.dataelt \
	elasto_statique.mesh test_mesh.mesh toto.mat test_mat_elem.mesh \
	ii_files/* auto_gmm*


dynamic_array_SOURCES = dynamic_array.C 
dynamic_tas_SOURCES = dynamic_tas.C 
test_assembly_SOURCES = test_assembly.C
laplacian_SOURCES = laplacian.C
elasto_statique_SOURCES = elasto_statique.C
integration_SOURCES = integration.C
poly_SOURCES = poly.C
linalg_SOURCES = linalg.C
test_mesh_SOURCES = test_mesh.C
geo_trans_SOURCES = geo_trans.C
geo_trans_inv_SOURCES = geo_trans_inv.C
speed_test_for_gmm_SOURCES = speed_test_for_gmm.C
test_int_set_SOURCES = test_int_set.C
test_link_fem_SOURCES = test_link_fem.C
test_tree_sorted_SOURCES = test_tree_sorted.C
test_mat_elem_SOURCES = test_mat_elem.C

INCLUDES = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@

TESTS = \
	$(top_srcdir)/tests/dynamic_array.pl \
	$(top_srcdir)/tests/dynamic_tas.pl \
	$(top_srcdir)/tests/test_mat_elem.pl \
	$(top_srcdir)/tests/test_assembly.pl \
	$(top_srcdir)/tests/laplacian.pl \
	$(top_srcdir)/tests/elasto_statique.pl  \
	$(top_srcdir)/tests/integration.pl \
	$(top_srcdir)/tests/linalg.pl \
	$(top_srcdir)/tests/poly.pl \
	$(top_srcdir)/tests/test_mesh.pl \
	$(top_srcdir)/tests/geo_trans.pl \
	$(top_srcdir)/tests/geo_trans_inv.pl \
	$(top_srcdir)/tests/speed_test_for_gmm.pl \
	$(top_srcdir)/tests/test_int_set.pl \
	$(top_srcdir)/tests/test_link_fem.pl \
	$(top_srcdir)/tests/make_gmm_test.pl \
	$(top_srcdir)/tests/test_tree_sorted.pl

EXTRA_DIST = \
	dynamic_array.pl \
	dynamic_tas.pl \
	test_assembly.pl \
	laplacian.pl \
	elasto_statique.pl  \
	integration.pl \
	linalg.pl \
	poly.pl \
	test_mesh.pl \
	geo_trans.pl \
	geo_trans_inv.pl \
	speed_test_for_gmm.pl \
	test_int_set.pl \
	test_link_fem.pl \
	test_tree_sorted.pl \
	test_mat_elem.pl \
	laplacian.param \
	elasto_statique.param \
	test_link_fem.param \
	geo_trans_inv.param \
	gmm_test_sub.C \
	gmm_test_mult.C \
	gmm_test_mat_mult.C \
	gmm_test_lusolve.C

TESTS_ENVIRONMENT = perl
