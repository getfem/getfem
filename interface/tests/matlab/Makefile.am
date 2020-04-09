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

SUBDIRS = private

if BUILDMEX
TESTS = $(abs_srcdir)/check_all.sh
else
TESTS = 
endif


EXTRA_DIST=						\
	axrot_matrix.m 					\
	check_all.m 					\
	check_all.sh 					\
	check_asm.m 					\
	check_fem.m 					\
	check_geotrans.m 				\
	check_integ.m 					\
	check_mesh_fem.m 				\
	check_oo.m 					\
	check_plot.m 					\
	check_slices.m 					\
	check_spmat.m 					\
	check_mitc.m 					\
	check_workspace.m 				\
	check_plasticity.m 				\
	demo_elasticity.m 				\
	demo_bilaplacian.m 				\
	demo_laplacian_DG.m 				\
	demo_Mindlin_Reissner_plate.m 			\
	demo_laplacian.m 				\
	demo_laplacian_pyramid.m			\
	demo_periodic_laplacian.m 			\
	demo_nonlinear_elasticity.m 			\
	demo_nonlinear_elasticity_anim.m 		\
	demo_plasticity.m 				\
	demo_dynamic_plasticity.m 			\
	demo_continuation.m 				\
	demo_crack.m 					\
	demo_mortar.m 					\
	demo_mesh_generation.m 				\
	demo_topological_optimization.m 		\
	demo_step_by_step.m 				\
	demo_Nitsche_contact_by_hand.m 			\
	demo_structural_optimization.m 			\
	demo_thermo_elasticity_electrical_coupling.m 	\
	demo_refine.m 					\
	demo_stokes_3D_tank.m 				\
	demo_stokes_3D_tank_draw.m 			\
	demo_tripod.m 					\
	demo_Navier_Stokes.m 				\
	demo_convection_rotating_cavity.m 		\
	demo_tripod_alt.m 				\
	demo_tripod_slice_anim.m 			\
	demo_fictitious_domains.m			\
	demo_fictitious_domains_laplacian.m 		\
	demo_contact_fictitious_domain_nitsche.m 	\
	demo_static_contact.m 				\
	demo_large_sliding_contact.m 			\
	demo_wheel_contact.m 				\
	demo_wave2D.m 					\
	demo_wave2D_alt.m 				\
	demo_wave_equation.m 				\
	test_argyris.m 					\
	tripod_anim.m 					\
	tutorial1.m 					\
							\
	private/gfassert.m 				\
	private/asserterr.m

TESTS_ENVIRONMENT = \
	MATLABPATH=@abs_top_builddir@/interface/src:@abs_top_builddir@/interface/tests/matlab:@abs_top_srcdir@/interface/src:@abs_top_srcdir@/interface/tests/matlab \
	MATLAB_ROOT=@MATLAB_ROOT@
	PATH=@abs_top_builddir@/src:$(PATH)
#	LD_LIBRARY_PATH=$(prefix)/lib:@GETFEMLIB@:@abs_top_builddir@/src:$(LD_LIBRARY_PATH) 

toolboxdir=@TOOLBOXDIR@/tests
toolbox_SCRIPTS=$(EXTRA_DIST)
