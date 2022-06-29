#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM
#
#  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

if QHULL
optprogs = test_mesh_generation test_mesh_im_level_set crack thermo_elasticity_electrical_coupling
optpl = test_mesh_im_level_set.pl \
        crack.pl                  \
	thermo_elasticity_electrical_coupling.pl \
	test_mesh_generation.pl
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
	test_assembly_assignment   \
	test_interpolated_fem      \
	test_internal_variables    \
	test_condensation          \
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
	bilaplacian		   \
	heat_equation              \
	wave_equation 		   \
	cyl_slicer		   \
	test_continuation          \
	test_gmm_matrix_functions

CLEANFILES = \
	laplacian.res laplacian.mesh laplacian.dataelt 			    \
	elasto_static.mesh test_mesh.mesh toto.mat test_mat_elem.mesh       \
	helmholtz.vtk helmholtz.vtu plate.mesh plate.vtk 		    \
	nonlinear_elastostatic.mesh					    \
	nonlinear_elastostatic.mf nonlinear_elastostatic.mfd                \
	nonlinear_elastostatic.dx plasticity.mesh plasticity.U              \
        plasticity.sigmabar plasticity.meshfem plasticity.coef              \
	ii_files/* auto_gmm* dyn*.txt *.sl time FN0 *.vtk                   \
	nonlinear_elastostatic.U crack.mesh cut.mesh nonlinear_membrane.mfd \
	nonlinear_membrane.mesh test_range_basis.mesh nonlinear_membrane.mf \
	Q2_incomplete.pos Q2_incomplete.msh

dynamic_array_SOURCES = dynamic_array.cc 
dynamic_tas_SOURCES = dynamic_tas.cc 
test_small_vector_SOURCES = test_small_vector.cc
test_kdtree_SOURCES = test_kdtree.cc
test_rtree_SOURCES = test_rtree.cc
test_assembly_SOURCES = test_assembly.cc
test_assembly_assignment_SOURCES = test_assembly_assignment.cc
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
test_internal_variables_SOURCES = test_internal_variables.cc
test_condensation_SOURCES = test_condensation.cc
test_tree_sorted_SOURCES = test_tree_sorted.cc
test_mat_elem_SOURCES = test_mat_elem.cc
test_slice_SOURCES = test_slice.cc
test_range_basis_SOURCES = test_range_basis.cc
schwarz_additive_SOURCES = schwarz_additive.cc
plasticity_SOURCES = plasticity.cc
if QHULL
test_mesh_generation_SOURCES = test_mesh_generation.cc 
test_mesh_im_level_set_SOURCES = test_mesh_im_level_set.cc 
crack_SOURCES = crack.cc
thermo_elasticity_electrical_coupling_SOURCES = thermo_elasticity_electrical_coupling.cc
endif
bilaplacian_SOURCES = bilaplacian.cc
heat_equation_SOURCES = heat_equation.cc
wave_equation_SOURCES = wave_equation.cc
cyl_slicer_SOURCES = cyl_slicer.cc
test_continuation_SOURCES = test_continuation.cc
test_gmm_matrix_functions_SOURCES = test_gmm_matrix_functions.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -I../src
LDADD    = ../src/libgetfem.la -lm @SUPLDFLAGS@

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
	test_mesh.pl                  \
	test_interpolation.pl         \
	test_mat_elem.pl              \
	test_slice.pl                 \
	integration.pl                \
	test_assembly.pl              \
	test_assembly_assignment.pl   \
	test_interpolated_fem.pl      \
	test_internal_variables.pl    \
	test_condensation.pl          \
	test_range_basis.pl           \
	laplacian.pl                  \
	laplacian_with_bricks.pl      \
	elastostatic.pl               \
	stokes.pl                     \
	plate.pl                      \
	$(optpl)                      \
	nonlinear_elastostatic.pl     \
	nonlinear_membrane.pl         \
	test_continuation.pl   	      \
        plasticity.pl                 \
	helmholtz.pl                  \
	schwarz_additive.pl           \
	bilaplacian.pl    	      \
	heat_equation.pl              \
	wave_equation.pl   	      \
	test_gmm_matrix_functions.pl  \
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
	test_assembly_assignment.pl        			\
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
	test_internal_variables.pl         			\
	test_condensation.pl                                    \
	test_slice.pl			   			\
	test_mesh_im_level_set.pl          			\
	thermo_elasticity_electrical_coupling.pl		\
	thermo_elasticity_electrical_coupling.param		\
	test_tree_sorted.pl                			\
	test_mat_elem.pl                   			\
	test_mesh_generation.pl                                 \
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
	plasticity.pl                      			\
	plasticity.param                      			\
	nonlinear_elastostatic.param       			\
	test_interpolated_fem.param        			\
	test_gmm_matrix_functions.pl              		\
	geo_trans_inv.param                			\
	heat_equation.pl                   			\
	heat_equation.param                			\
	wave_equation.pl                   			\
	wave_equation.param                			\
	cyl_slicer.pl                           		\
	test_continuation.param                                 \
	test_continuation.pl                                    \
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
