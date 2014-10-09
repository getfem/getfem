SUBDIRS = private

if BUILDMEX
TESTS = $(abs_srcdir)/check_all.sh
else
TESTS = 
endif


EXTRA_DIST= \
	axrot_matrix.m \
	check_all.m \
	check_all.sh \
	check_asm.m \
	check_bricks.m \
	check_fem.m \
	check_geotrans.m \
	check_integ.m \
	check_mesh_fem.m \
	check_oo.m \
	check_plot.m \
	check_slices.m \
	check_spmat.m \
	check_workspace.m \
	demo_bilaplacian.m \
	demo_laplacian.m \
	demo_periodic_laplacian.m \
	demo_laplacian1D.m \
	demo_nonlinear_elasticity.m \
	demo_nonlinear_elasticity_anim.m \
	demo_plasticity.m \
	demo_continuation.m \
	demo_crack.m \
	demo_mortar.m \
	demo_mesh_generation.m \
	demo_topological_optimization.m \
	demo_structural_optimization.m \
	demo_refine.m \
	demo_stokes_3D_tank.m \
	demo_stokes_3D_tank_draw.m \
	demo_tripod.m \
	demo_Navier_Stokes.m \
	demo_convection_rotating_cavity.m \
	demo_tripod_alt.m \
	demo_tripod_slice_anim.m \
	demo_fictitious_domains.m \
	demo_fictitious_domains_laplacian.m \
	demo_contact_fictitious_domain_nitsche.m \
	demo_static_contact.m \
	demo_large_sliding_contact.m \
	demo_wave2D.m \
	demo_wave2D_alt.m \
	demo_wave_equation.m \
	test_argyris.m \
	tripod_anim.m \
	tutorial1.m \
	\
	private/gfassert.m \
	private/asserterr.m

TESTS_ENVIRONMENT = \
	MATLABPATH=@abs_top_builddir@/interface/src:@abs_top_builddir@/interface/tests/matlab:@abs_top_srcdir@/interface/src:@abs_top_srcdir@/interface/tests/matlab \
	MATLAB_ROOT=@MATLAB_ROOT@
	PATH=@abs_top_builddir@/src:$(PATH)
#	LD_LIBRARY_PATH=$(prefix)/lib:@GETFEMLIB@:@abs_top_builddir@/src:$(LD_LIBRARY_PATH) 

toolboxdir=@TOOLBOXDIR@/tests
toolbox_SCRIPTS=$(EXTRA_DIST)
