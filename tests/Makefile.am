check_PROGRAMS = \
	dynamic_array \
	dynamic_tas \
	laplacien \
	elasto_statique  \
	integration \
	poly \
	test_mesh \
	geo_trans \
	test_int_set \
	test_tree_sorted

EXTRA_DIST = \
	dynamic_array.pl \
	dynamic_tas.pl \
	laplacien.pl \
	elasto_statique.pl  \
	integration.pl \
	poly.pl \
	test_mesh.pl \
	geo_trans.pl \
	test_int_set.pl \
	test_tree_sorted.pl \
	laplacien.param \
	elasto_statique.param

dynamic_array_SOURCES = dynamic_array.C 
dynamic_tas_SOURCES = dynamic_tas.C 
laplacien_SOURCES = laplacien.C
elasto_statique_SOURCES = elasto_statique.C
integration_SOURCES = integration.C
poly_SOURCES = poly.C
test_mesh_SOURCES = test_mesh.C
geo_trans_SOURCES = geo_trans.C
test_int_set_SOURCES = test_int_set.C
test_tree_sorted_SOURCES = test_tree_sorted.C

INCLUDES = -I$(top_srcdir)/src
LDADD    = ../src/.libs/libgetfem.a -lm
CXXFLAGS = @APPSCXXFLAGS@

TESTS = \
	dynamic_array.pl \
	dynamic_tas.pl \
	laplacien.pl \
	elasto_statique.pl  \
	integration.pl \
	poly.pl \
	test_mesh.pl \
	geo_trans.pl \
	test_int_set.pl \
	test_tree_sorted.pl

TESTS_ENVIRONMENT = perl


clean-local:
	rm -rf cxx_repository
