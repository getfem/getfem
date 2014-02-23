if QHULL
optprogs = test_mesh_im_level_set
optpl = test_mesh_im_level_set.pl
else
optprogs =
optpl = 
endif

check_PROGRAMS =                 \
	dynamic_array            \
	dynamic_tas              \
	test_int_set             \
	test_tree_sorted         \
	poly                     \
	test_small_vector        \
	test_kdtree	         \
	test_rtree	         \
	test_mesh                \
	test_slice               \
	integration              \
	geo_trans_inv            \
	test_mat_elem            \
	test_interpolation       \
	test_assembly            \
	test_norm	         \
	test_interpolated_fem    \
	laplacian                \
	elastostatic             \
	stokes                   \
	helmholtz                \
	plate                    \
	nonlinear_elastostatic   \
	schwarz_additive         \
	$(optprogs)              \
	plasticity               \
	dynamic_friction         \
	bilaplacian

CLEANFILES = \
	laplacian.res laplacian.mesh laplacian.dataelt \
	elasto_static.mesh test_mesh.mesh toto.mat test_mat_elem.mesh \
	helmholtz.vtk plate.mesh plate.vtk nonlinear_elastostatic.mesh \
	nonlinear_elastostatic.mf nonlinear_elastostatic.mfd \
	nonlinear_elastostatic.dx nonlinear_elastostatic.vtk \
	plasticity.mesh plasticity.U plasticity.sigmabar plasticity.vtk \
	plasticity.meshfem plasticity.coef elasto_static_err.vtk \
	ii_files/* auto_gmm* dyn*.txt

dynamic_array_SOURCES = dynamic_array.cc 
dynamic_tas_SOURCES = dynamic_tas.cc 
test_small_vector_SOURCES = test_small_vector.cc
test_kdtree_SOURCES = test_kdtree.cc
test_rtree_SOURCES = test_rtree.cc
test_assembly_SOURCES = test_assembly.cc
laplacian_SOURCES = laplacian.cc
plate_SOURCES = plate.cc
test_interpolation_SOURCES = test_interpolation.cc
elastostatic_SOURCES = elastostatic.cc
stokes_SOURCES = stokes.cc
nonlinear_elastostatic_SOURCES = nonlinear_elastostatic.cc
helmholtz_SOURCES = helmholtz.cc
integration_SOURCES = integration.cc
poly_SOURCES = poly.cc
test_mesh_SOURCES = test_mesh.cc
geo_trans_inv_SOURCES = geo_trans_inv.cc
test_int_set_SOURCES = test_int_set.cc
test_interpolated_fem_SOURCES = test_interpolated_fem.cc
test_tree_sorted_SOURCES = test_tree_sorted.cc
test_mat_elem_SOURCES = test_mat_elem.cc
test_slice_SOURCES = test_slice.cc
test_norm_SOURCES = test_norm.cc
schwarz_additive_SOURCES = schwarz_additive.cc
plasticity_SOURCES = plasticity.cc
if QHULL
test_mesh_im_level_set_SOURCES = test_mesh_im_level_set.cc 
endif
dynamic_friction_SOURCES = dynamic_friction.cc
bilaplacian_SOURCES = bilaplacian.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

#plasticity_LDADD = ../src/libgetfem.la -lm @SUPLDFLAGS@ $(HOME)/source++/SuperLU/superlu.a -lblas -lg2c
#plasticity_INCLUDES = $(INCLUDES) -I$(HOME)/source++/SuperLU/

TESTS = \
	dynamic_array.pl           \
	dynamic_tas.pl             \
	test_int_set.pl            \
	test_tree_sorted.pl        \
	poly.pl                    \
	test_small_vector.pl       \
	test_kdtree.pl             \
	test_rtree.pl              \
	geo_trans_inv.pl           \
	test_norm.pl               \
	test_mesh.pl               \
	test_interpolation.pl      \
	test_mat_elem.pl           \
	test_slice.pl              \
	integration.pl             \
	test_assembly.pl           \
	test_interpolated_fem.pl   \
	laplacian.pl               \
	elastostatic.pl            \
	stokes.pl                  \
	plate.pl                   \
	$(optpl)                   \
	nonlinear_elastostatic.pl  \
	dynamic_friction.pl        \
	plasticity.pl              \
	helmholtz.pl               \
	schwarz_additive.pl        \
	bilaplacian.pl

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
	stokes.pl                         \
	helmholtz.pl                      \
	integration.pl                    \
	poly.pl                           \
	test_mesh.pl                      \
	geo_trans_inv.pl                  \
	test_int_set.pl                   \
	test_interpolated_fem.pl          \
	test_slice.pl			  \
	test_norm.pl			  \
	test_mesh_im_level_set.pl         \
	test_tree_sorted.pl               \
	test_mat_elem.pl                  \
	laplacian.param                   \
	bilaplacian.param                 \
	bilaplacian.pl                    \
	plate.param                       \
	plate.pl                          \
	schwarz_additive.pl               \
	schwarz_additive.param            \
	elastostatic.param                \
	stokes.param                      \
	nonlinear_elastostatic.pl         \
	dynamic_friction.pl               \
	dynamic_friction.param            \
	plasticity.pl                     \
	nonlinear_elastostatic.param      \
	test_interpolated_fem.param       \
	geo_trans_inv.param

LOG_COMPILER = perl
