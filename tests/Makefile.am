check_PROGRAMS = \
	dynamic_array      \
	dynamic_tas        \
	test_small_vector  \
	test_kdtree	   \
	test_rtree	   \
	test_interpolation \
	elastostatic       \
	nonlinear_elastostatic   \
	helmholtz          \
	geo_trans          \
	geo_trans_inv      \
	integration        \
	laplacian          \
	poly               \
	schwarz_additive    \
	test_assembly      \
	test_int_set       \
	test_link_fem      \
	test_mat_elem      \
	test_slice         \
	test_norm	   \
	test_mesh          \
	test_tree_sorted   \
	plasticity

CLEANFILES = \
	laplacian.res laplacian.mesh laplacian.dataelt \
	elastostatic.mesh test_mesh.mesh toto.mat test_mat_elem.mesh \
	ii_files/* auto_gmm*


dynamic_array_SOURCES = dynamic_array.C 
dynamic_tas_SOURCES = dynamic_tas.C 
test_small_vector_SOURCES = test_small_vector.C
test_kdtree_SOURCES = test_kdtree.C
test_rtree_SOURCES = test_rtree.C
test_assembly_SOURCES = test_assembly.C
laplacian_SOURCES = laplacian.C
test_interpolation_SOURCES = test_interpolation.C
elastostatic_SOURCES = elastostatic.C
nonlinear_elastostatic_SOURCES = nonlinear_elastostatic.C
helmholtz_SOURCES = helmholtz.C
integration_SOURCES = integration.C
poly_SOURCES = poly.C
test_mesh_SOURCES = test_mesh.C
geo_trans_SOURCES = geo_trans.C
geo_trans_inv_SOURCES = geo_trans_inv.C
test_int_set_SOURCES = test_int_set.C
test_link_fem_SOURCES = test_link_fem.C
test_tree_sorted_SOURCES = test_tree_sorted.C
test_mat_elem_SOURCES = test_mat_elem.C
test_slice_SOURCES = test_slice.C
test_norm_SOURCES = test_norm.C
schwarz_additive_SOURCES = schwarz_additive.C
plasticity_SOURCES = plasticity.C
#test_fem_SOURCES = test_fem.C

INCLUDES = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@

plasticity_LDADD = ../src/libgetfem.la -lm @SUPLDFLAGS@ $(HOME)/source++/SuperLU/superlu.a -lblas -lg2c
plasticity_INCLUDES = $(INCLUDES) -I$(HOME)/source++/SuperLU/

TESTS = \
	$(top_srcdir)/tests/dynamic_array.pl      \
	$(top_srcdir)/tests/dynamic_tas.pl        \
	$(top_srcdir)/tests/test_small_vector.pl  \
	$(top_srcdir)/tests/test_kdtree.pl        \
	$(top_srcdir)/tests/test_rtree.pl         \
	$(top_srcdir)/tests/test_interpolation.pl \
	$(top_srcdir)/tests/test_mat_elem.pl      \
	$(top_srcdir)/tests/test_assembly.pl      \
	$(top_srcdir)/tests/laplacian.pl          \
	$(top_srcdir)/tests/elastostatic.pl       \
	$(top_srcdir)/tests/helmholtz.pl          \
	$(top_srcdir)/tests/integration.pl        \
	$(top_srcdir)/tests/poly.pl               \
	$(top_srcdir)/tests/test_mesh.pl          \
	$(top_srcdir)/tests/geo_trans.pl          \
	$(top_srcdir)/tests/geo_trans_inv.pl      \
	$(top_srcdir)/tests/schwarz_additive.pl   \
	$(top_srcdir)/tests/test_int_set.pl       \
	$(top_srcdir)/tests/test_link_fem.pl      \
	$(top_srcdir)/tests/test_slice.pl         \
	$(top_srcdir)/tests/test_tree_sorted.pl   \
	$(top_srcdir)/tests/make_gmm_test.pl     

EXTRA_DIST = \
	dynamic_array.pl                  \
	dynamic_tas.pl                    \
	test_small_vector.pl		  \
	test_kdtree.pl                    \
	test_rtree.pl                     \
	test_interpolation.pl             \
	test_assembly.pl                  \
	laplacian.pl                      \
	elastostatic.pl                   \
	helmholtz.pl                      \
	integration.pl                    \
	poly.pl                           \
	test_mesh.pl                      \
	geo_trans.pl                      \
	geo_trans_inv.pl                  \
	test_int_set.pl                   \
	test_link_fem.pl                  \
	test_slice.pl			  \
	test_tree_sorted.pl               \
	test_mat_elem.pl                  \
	laplacian.param                   \
	schwarz_additive.pl               \
	schwarz_additive.param            \
	elastostatic.param                \
	nonlinear_elastostatic.param      \
	test_link_fem.param               \
	geo_trans_inv.param               \
	make_gmm_test.pl                  \
	gmm_torture01_lusolve.C           \
	gmm_torture05_mult.C              \
	gmm_torture06_mat_mult.C          \
	gmm_torture10_qr.C                \
	gmm_torture15_sub.C               \
	gmm_torture20_iterative_solvers.C

TESTS_ENVIRONMENT = perl
