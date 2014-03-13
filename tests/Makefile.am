if QHULL
optprogs = test_mesh_generation test_mesh_im_level_set crack
optpl = test_mesh_im_level_set.pl \
        crack.pl
else
optprogs =
optpl = 
endif

check_PROGRAMS =                   \
	dynamic_array              \
	dynamic_tas                \
	test_int_set               \
	test_tree_sorted           \
	poly                       \
	test_small_vector          \
	test_kdtree	           \
	test_rtree	           \
	test_mesh                  \
	test_slice                 \
	integration                \
	geo_trans_inv              \
	test_mat_elem              \
	test_interpolation         \
	test_assembly              \
	test_norm	           \
	test_interpolated_fem      \
	test_range_basis           \
	laplacian                  \
	laplacian_with_bricks      \
	elastostatic               \
	stokes                     \
	helmholtz                  \
	plate                      \
	nonlinear_elastostatic     \
	nonlinear_membrane         \
	schwarz_additive           \
	$(optprogs)                \
	plasticity                 \
	plasticity_old_brick       \
	dynamic_friction           \
	bilaplacian		   \
	heat_equation              \
	wave_equation 		   \
	cyl_slicer		   \
	test_large_sliding_contact \
	test_continuation

CLEANFILES = \
	laplacian.res laplacian.mesh laplacian.dataelt 			    \
	elasto_static.mesh test_mesh.mesh toto.mat test_mat_elem.mesh       \
	helmholtz.vtk plate.mesh plate.vtk nonlinear_elastostatic.mesh      \
	nonlinear_elastostatic.mf nonlinear_elastostatic.mfd                \
	nonlinear_elastostatic.dx nonlinear_elastostatic.vtk                \
	plasticity.mesh plasticity.U plasticity.sigmabar plasticity.vtk     \
	plasticity.meshfem plasticity.coef elasto_static_err.vtk            \
	ii_files/* auto_gmm* dyn*.txt dynamic_friction.data                 \
	dynamic_friction.U dynamic_friction.mfu  *.sl time FN0              \
	nonlinear_elastostatic.U crack.mesh cut.mesh nonlinear_membrane.mfd \
	nonlinear_membrane.mesh test_range_basis.mesh nonlinear_membrane.mf \
	nonlinear_membrane.vtk feminterpolation.vtk circleslice.vtk

dynamic_array_SOURCES = dynamic_array.cc 
dynamic_tas_SOURCES = dynamic_tas.cc 
test_small_vector_SOURCES = test_small_vector.cc
test_kdtree_SOURCES = test_kdtree.cc
test_rtree_SOURCES = test_rtree.cc
test_assembly_SOURCES = test_assembly.cc
laplacian_SOURCES = laplacian.cc
laplacian_with_bricks_SOURCES = laplacian_with_bricks.cc
plate_SOURCES = plate.cc
test_interpolation_SOURCES = test_interpolation.cc
elastostatic_SOURCES = elastostatic.cc
stokes_SOURCES = stokes.cc
nonlinear_elastostatic_SOURCES = nonlinear_elastostatic.cc
nonlinear_membrane_SOURCES = nonlinear_membrane.cc
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
test_range_basis_SOURCES = test_range_basis.cc
test_norm_SOURCES = test_norm.cc
schwarz_additive_SOURCES = schwarz_additive.cc
plasticity_SOURCES = plasticity.cc
plasticity_old_brick_SOURCES = plasticity_old_brick.cc
if QHULL
test_mesh_generation_SOURCES = test_mesh_generation.cc 
test_mesh_im_level_set_SOURCES = test_mesh_im_level_set.cc 
crack_SOURCES = crack.cc
endif
dynamic_friction_SOURCES = dynamic_friction.cc
bilaplacian_SOURCES = bilaplacian.cc
heat_equation_SOURCES = heat_equation.cc
wave_equation_SOURCES = wave_equation.cc
cyl_slicer_SOURCES = cyl_slicer.cc
test_continuation_SOURCES = test_continuation.cc
test_large_sliding_contact_SOURCES = test_large_sliding_contact.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@ @BOOST_LIBS@

TESTS = \
	dynamic_array.pl              \
	dynamic_tas.pl                \
	test_int_set.pl               \
	test_tree_sorted.pl           \
	poly.pl                       \
	test_small_vector.pl          \
	test_kdtree.pl                \
	test_rtree.pl                 \
	geo_trans_inv.pl              \
	test_norm.pl                  \
	test_mesh.pl                  \
	test_interpolation.pl         \
	test_mat_elem.pl              \
	test_slice.pl                 \
	integration.pl                \
	test_assembly.pl              \
	test_interpolated_fem.pl      \
	test_range_basis.pl           \
	laplacian.pl                  \
	laplacian_with_bricks.pl      \
	elastostatic.pl               \
	stokes.pl                     \
	plate.pl                      \
	$(optpl)                      \
	nonlinear_elastostatic.pl     \
	nonlinear_membrane.pl         \
	dynamic_friction.pl           \
	plasticity.pl                 \
	plasticity_old_brick.pl       \
	helmholtz.pl                  \
	schwarz_additive.pl           \
	bilaplacian.pl    	          \
	heat_equation.pl              \
	wave_equation.pl   	          \
	test_large_sliding_contact.pl \
	cyl_slicer.pl	              \
	make_gmm_test.pl

EXTRA_DIST =                               			\
	dynamic_array.pl                   			\
	dynamic_tas.pl                     			\
	test_small_vector.pl		   			\
	test_kdtree.pl                     			\
	test_rtree.pl                      			\
	test_interpolation.pl              			\
	test_assembly.pl                   			\
	laplacian.pl                       			\
	laplacian_with_bricks.pl 	   			\
	elastostatic.pl                    			\
	stokes.pl                          			\
	helmholtz.pl                       			\
	integration.pl                     			\
	poly.pl                            			\
	test_mesh.pl                       			\
	geo_trans_inv.pl                   			\
	test_int_set.pl                    			\
	test_interpolated_fem.pl           			\
	test_slice.pl			   			\
	test_norm.pl			   			\
	test_mesh_im_level_set.pl          			\
	test_tree_sorted.pl                			\
	test_mat_elem.pl                   			\
	crack.pl                           			\
	crack.param                        			\
	laplacian.param                    			\
	laplacian_with_bricks.param        			\
	test_range_basis.param             			\
	test_range_basis.pl                			\
	bilaplacian.param                  			\
	bilaplacian.pl                     			\
	plate.param                        			\
	plate.pl                           			\
	schwarz_additive.pl                			\
	schwarz_additive.param             			\
	elastostatic.param                 			\
	stokes.param                       			\
	nonlinear_elastostatic.pl          			\
	nonlinear_membrane.pl              			\
	nonlinear_membrane.param           			\
	dynamic_friction.pl                			\
	dynamic_friction.param             			\
	dynamic_friction.net               			\
	plasticity.pl                      			\
	plasticity_old_brick.pl                			\
	plasticity.param                      			\
	plasticity_old_brick.param             			\
	nonlinear_elastostatic.param       			\
	test_interpolated_fem.param        			\
	test_large_sliding_contact.pl                           \
	test_large_sliding_contact.param                        \
	geo_trans_inv.param                			\
	heat_equation.pl                   			\
	heat_equation.param                			\
	wave_equation.pl                   			\
	wave_equation.param                			\
	cyl_slicer.pl                           		\
	test_continuation.param                                 \
	make_gmm_test.pl                   			\
	gmm_torture01_lusolve.cc           			\
	gmm_torture05_mult.cc              			\
	gmm_torture06_mat_mult.cc          			\
	gmm_torture10_qr.cc                			\
	gmm_torture15_sub.cc               			\
	gmm_torture20_iterative_solvers.cc 			\
        meshes/disc_2D_degree3.mesh        			\
	meshes/donut_regulier_72_elements.mesh    		\
	meshes/disc_P2_h0_3.mesh                 		\
        meshes/donut_regulier_8_elements_288ddl.mesh		\
	meshes/disc_P2_h0_5.mesh				\
	meshes/disc_P2_h1.mesh					\
	meshes/disc_P2_h2.mesh					\
	meshes/disc_P2_h4.mesh					\
	meshes/sphere_with_quadratic_tetra_16000_elts.mesh	\
	meshes/disc_P2_h6.mesh					\
	meshes/sphere_with_quadratic_tetra_2000_elts.mesh	\
	meshes/disc_P2_h8.mesh					\
	meshes/sphere_with_quadratic_tetra_400_elts.mesh	\
	meshes/disc_P2_h10.mesh					\
	meshes/sphere_with_quadratic_tetra_80_elts.mesh		\
	meshes/donut_regulier_32_elements.mesh			\
	meshes/sphere_with_quadratic_tetra_8_elts.mesh		\
	meshes/disc_with_a_hole.mesh				\
	meshes/punch2D_1.mesh					\
	meshes/punch2D_2.mesh					\
	meshes/multi_body.mesh					\
	meshes/donut_regulier_512_elements.mesh

LOG_COMPILER = perl

TESTS_ENVIRONMENT = perl
