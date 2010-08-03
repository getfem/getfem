## Scilab Getfem Interface
## Copyright (C) 2009 Yann COLLETTE
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

toolboxdir=@SCILAB_TOOLBOX_DIR@

scilabbuilddir=@GETFEM_BUILD_INTERFACE_PATH@/interface/src/scilab

scilabbasedir=builder.sce \
	loader.sce \
	readme.txt \
	changelog.txt \
	cleaner.sce \
	Makefile.am \
	sci_getfem.iss \
	license.txt \
	makefile_builder.sce

scilabsrccdir=src/c/copyright \
	src/c/makefile \
	src/c/loader.sce \
	src/c/makefile.in \
	src/c/builder_c.sce \
	src/c/ls.dat \
	src/c/README \
	src/c/rk4.dat \
	src/c/cleaner.sce \
	src/c/configure \
	src/c/configure.in \
	src/c/FILELIST \
	src/c/machine.h.in  

EXTRA_DIST=$(toolboxdir)/builddoc.sce \
        $(toolboxdir)/builder.sce \
        $(toolboxdir)/changelog.txt \
        $(toolboxdir)/cleaner.sce \
        $(toolboxdir)/license.txt \
        $(toolboxdir)/loader.sce \
        $(toolboxdir)/readme.txt \
        $(toolboxdir)/demos/axrot_matrix.sci \
        $(toolboxdir)/demos/check_all.sce \
        $(toolboxdir)/demos/check_asm.sce \
        $(toolboxdir)/demos/check_bricks.sce \
        $(toolboxdir)/demos/check_fem.sce \
        $(toolboxdir)/demos/check_geotrans.sce \
        $(toolboxdir)/demos/check_integ.sce \
        $(toolboxdir)/demos/check_interpolated_fem.sce \
        $(toolboxdir)/demos/check_levelset.sce \
        $(toolboxdir)/demos/check_mesh_fem.sce \
        $(toolboxdir)/demos/check_oo.sce \
        $(toolboxdir)/demos/check_plot.sce \
        $(toolboxdir)/demos/check_slices.sce \
        $(toolboxdir)/demos/check_solve.sce \
        $(toolboxdir)/demos/check_spmat.sce \
        $(toolboxdir)/demos/check_workspace.sce \
        $(toolboxdir)/demos/demo_bilaplacian.sce \
        $(toolboxdir)/demos/demo_convection_rotating_cavity.sce \
        $(toolboxdir)/demos/demo_crack.sce \
        $(toolboxdir)/demos/demo_fictitious_domains.sce \
        $(toolboxdir)/demos/demo_finallaplacian.sce \
        $(toolboxdir)/demos/demo_laplacian1D.sce \
        $(toolboxdir)/demos/demo_laplacian.sce \
        $(toolboxdir)/demos/demo_mortar.sce \
        $(toolboxdir)/demos/demo_nonlinear_elasticity_anim.sce \
        $(toolboxdir)/demos/demo_nonlinear_elasticity.sce \
        $(toolboxdir)/demos/demo_nonlinear_elasticity_U.mat \
        $(toolboxdir)/demos/demo_plasticity.sce \
        $(toolboxdir)/demos/demo_plate.sce \
        $(toolboxdir)/demos/demo_refine.sce \
        $(toolboxdir)/demos/demo_static_contact.sce \
        $(toolboxdir)/demos/demo_step_by_step.sce \
        $(toolboxdir)/demos/demo_stokes_2D_poiseuille_arc.sce \
        $(toolboxdir)/demos/demo_stokes_2D_poiseuille.sce \
        $(toolboxdir)/demos/demo_stokes_2D_tube.sce \
        $(toolboxdir)/demos/demo_stokes_3D_tank_alt.sce \
        $(toolboxdir)/demos/demo_stokes_3D_tank_draw.sce \
        $(toolboxdir)/demos/demo_stokes_3D_tank.sce \
        $(toolboxdir)/demos/demo_stokes_3D_tank_UP.mat \
        $(toolboxdir)/demos/demo_structural_optimization.sce \
        $(toolboxdir)/demos/demo_topological_optimization.sce \
        $(toolboxdir)/demos/demo_tripod_alt.sce \
        $(toolboxdir)/demos/demo_tripod_anim.sce \
        $(toolboxdir)/demos/demo_tripod.sce \
        $(toolboxdir)/demos/demo_tripod_slice_anim.sce \
        $(toolboxdir)/demos/demo_wave2D_alt.sce \
        $(toolboxdir)/demos/demo_wave2D_animate.sce \
        $(toolboxdir)/demos/demo_wave2D.sce \
        $(toolboxdir)/demos/demo_wave_equation.sce \
        $(toolboxdir)/demos/sci_getfem.dem.gateway.sce \
        $(toolboxdir)/demos/test_argyris.sce \
        $(toolboxdir)/demos/test_plasticity_new_brick.sce \
        $(toolboxdir)/demos/tutorial1.sce \
        $(toolboxdir)/demos/data/disc_P2_h0.5.mesh \
        $(toolboxdir)/demos/data/disc_P2_h1.mesh \
        $(toolboxdir)/demos/data/disc_P2_h2.mesh \
        $(toolboxdir)/demos/data/donut_regulier.mesh \
        $(toolboxdir)/demos/data/donut_with_quadratic_tetra_1100_elements.msh \
        $(toolboxdir)/demos/data/holed_disc_with_quadratic_2D_triangles.msh \
        $(toolboxdir)/demos/data/quad.geo \
        $(toolboxdir)/demos/data/quad.msh \
        $(toolboxdir)/demos/data/tank_quadratic_2500.GiD.msh \
        $(toolboxdir)/demos/data/tripod.GiD.msh \
        $(toolboxdir)/demos/data/tripod.mesh \
        $(toolboxdir)/demos/data/tube_2D_spline.GiD.msh \
        $(toolboxdir)/etc/sci_getfem.quit \
        $(toolboxdir)/etc/sci_getfem.start \
        $(toolboxdir)/help/builder_help.sce \
        $(toolboxdir)/help/en_US/ \
        $(toolboxdir)/help/fig/ \
        $(toolboxdir)/help/fr_FR/ \
        $(toolboxdir)/help/latex/ \
        $(toolboxdir)/help/mml/ \
        $(toolboxdir)/help/tmp/ \
        $(toolboxdir)/help/en_US/build_help.sce \
        $(toolboxdir)/help/en_US/gf_asm_pdetoolbc.xml \
        $(toolboxdir)/help/en_US/gf_asm.xml \
        $(toolboxdir)/help/en_US/gf_colormap.xml \
        $(toolboxdir)/help/en_US/gf_compute_Q1grid_interp.xml \
        $(toolboxdir)/help/en_US/gf_compute.xml \
        $(toolboxdir)/help/en_US/gf_cvstruct_get.xml \
        $(toolboxdir)/help/en_US/gf_delete.xml \
        $(toolboxdir)/help/en_US/gf_eltm.xml \
        $(toolboxdir)/help/en_US/gf_fem_get.xml \
        $(toolboxdir)/help/en_US/gf_fem.xml \
        $(toolboxdir)/help/en_US/gf_geotrans_get.xml \
        $(toolboxdir)/help/en_US/gf_geotrans.xml \
        $(toolboxdir)/help/en_US/gf_global_function_get.xml \
        $(toolboxdir)/help/en_US/gf_global_function.xml \
        $(toolboxdir)/help/en_US/gf_integ_get.xml \
        $(toolboxdir)/help/en_US/gf_integ.xml \
        $(toolboxdir)/help/en_US/gf_interpolate_on_grid.xml \
        $(toolboxdir)/help/en_US/gf_levelset_get.xml \
        $(toolboxdir)/help/en_US/gf_levelset_set.xml \
        $(toolboxdir)/help/en_US/gf_levelset.xml \
        $(toolboxdir)/help/en_US/gf_linsolve.xml \
        $(toolboxdir)/help/en_US/gf_mdbrick_get.xml \
        $(toolboxdir)/help/en_US/gf_mdbrick_set.xml \
        $(toolboxdir)/help/en_US/gf_mdbrick.xml \
        $(toolboxdir)/help/en_US/gf_mdstate_get.xml \
        $(toolboxdir)/help/en_US/gf_mdstate_set.xml \
        $(toolboxdir)/help/en_US/gf_mdstate.xml \
        $(toolboxdir)/help/en_US/gf_mesh_fem_get_eval.xml \
        $(toolboxdir)/help/en_US/gf_mesh_fem_get.xml \
        $(toolboxdir)/help/en_US/gf_mesh_fem_set.xml \
        $(toolboxdir)/help/en_US/gf_mesh_fem.xml \
        $(toolboxdir)/help/en_US/gf_mesh_get.xml \
        $(toolboxdir)/help/en_US/gf_mesh_im_get.xml \
        $(toolboxdir)/help/en_US/gf_mesh_im_set.xml \
        $(toolboxdir)/help/en_US/gf_mesh_im.xml \
        $(toolboxdir)/help/en_US/gf_mesh_levelset_get.xml \
        $(toolboxdir)/help/en_US/gf_mesh_levelset_set.xml \
        $(toolboxdir)/help/en_US/gf_mesh_levelset.xml \
        $(toolboxdir)/help/en_US/gf_mesh_set.xml \
        $(toolboxdir)/help/en_US/gf_mesh.xml \
        $(toolboxdir)/help/en_US/gf_model_get.xml \
        $(toolboxdir)/help/en_US/gf_model_set.xml \
        $(toolboxdir)/help/en_US/gf_model.xml \
        $(toolboxdir)/help/en_US/gf_plot_1D.xml \
        $(toolboxdir)/help/en_US/gf_plot_mesh.xml \
        $(toolboxdir)/help/en_US/gf_plot_slice.xml \
        $(toolboxdir)/help/en_US/gf_plot.xml \
        $(toolboxdir)/help/en_US/gf_poly.xml \
        $(toolboxdir)/help/en_US/gf_precond_get.xml \
        $(toolboxdir)/help/en_US/gf_precond.xml \
        $(toolboxdir)/help/en_US/gf_slice_get.xml \
        $(toolboxdir)/help/en_US/gf_slice_set.xml \
        $(toolboxdir)/help/en_US/gf_slice.xml \
        $(toolboxdir)/help/en_US/gf_solve.xml \
        $(toolboxdir)/help/en_US/gf_spmat_get.xml \
        $(toolboxdir)/help/en_US/gf_spmat_set.xml \
        $(toolboxdir)/help/en_US/gf_spmat.xml \
        $(toolboxdir)/help/en_US/gf_typeof.xml \
        $(toolboxdir)/help/en_US/gf_util.xml \
        $(toolboxdir)/help/en_US/gf_workspace.xml \
        $(toolboxdir)/help/en_US/objects.xml \
        $(toolboxdir)/help/en_US/preliminary.xml \
        $(toolboxdir)/help/en_US/examples/another_laplacian.xml \
        $(toolboxdir)/help/en_US/examples/avoiding_bricks.xml \
        $(toolboxdir)/help/en_US/examples/CHAPTER \
        $(toolboxdir)/help/en_US/examples/linear_nonlinear_elast.xml \
        $(toolboxdir)/help/en_US/examples/step_by_step_example.xml \
        $(toolboxdir)/help/en_US/sparses/CHAPTER \
        $(toolboxdir)/help/en_US/sparses/sp_cgne.xml \
        $(toolboxdir)/help/en_US/sparses/sp_cgs.xml \
        $(toolboxdir)/help/en_US/sparses/sp_cholinc.xml \
        $(toolboxdir)/help/en_US/sparses/sp_chol.xml \
        $(toolboxdir)/help/en_US/sparses/sp_chsolve.xml \
        $(toolboxdir)/help/en_US/sparses/sp_gmres.xml \
        $(toolboxdir)/help/en_US/sparses/sp_luinc.xml \
        $(toolboxdir)/help/en_US/sparses/sp_lusolve.xml \
        $(toolboxdir)/help/en_US/sparses/sp_lu.xml \
        $(toolboxdir)/help/en_US/sparses/sp_mgcr.xml \
        $(toolboxdir)/help/fig/gf_fem_get_fig_1.png \
        $(toolboxdir)/help/fig/gf_plot_mesh_fig_1.png \
        $(toolboxdir)/help/fig/hierarchy.png \
        $(toolboxdir)/help/fig/tripodvonmiseswithmesh_small.png \
        $(toolboxdir)/help/fr_FR/build_help.sce \
        $(toolboxdir)/help/mml/avoiding_eq1.mml \
        $(toolboxdir)/help/mml/gf_asm_eq1.mml \
        $(toolboxdir)/help/mml/gf_asm_eq2.mml \
        $(toolboxdir)/help/mml/gf_asm_eq3.mml \
        $(toolboxdir)/help/mml/gf_asm_eq4.mml \
        $(toolboxdir)/help/mml/gf_asm_eq5.mml \
        $(toolboxdir)/help/mml/gf_asm_eq6.mml \
        $(toolboxdir)/help/mml/gf_model_set_eq1.mml \
        $(toolboxdir)/help/mml/gf_model_set_eq2.mml \
        $(toolboxdir)/help/mml/gf_model_set_eq3.mml \
        $(toolboxdir)/help/mml/gf_model_set_eq4.mml \
        $(toolboxdir)/help/mml/preliminary_eq1.mml \
        $(toolboxdir)/help/mml/step_by_step_eq1.mml \
        $(toolboxdir)/help/mml/step_by_step_eq2.mml \
        $(toolboxdir)/help/mml/step_by_step_eq3.mml \
        $(toolboxdir)/help/mml/step_by_step_eq4.mml \
        $(toolboxdir)/help/mml/step_by_step_eq5.mml \
        $(toolboxdir)/jar/scilab_en_US_help.jar \
        $(toolboxdir)/jar/scilab_fr_FR_help.jar \
        $(toolboxdir)/macros/add_empty_bound.bin \
        $(toolboxdir)/macros/add_empty_bound.sci \
        $(toolboxdir)/macros/assert.bin \
        $(toolboxdir)/macros/asserterr.bin \
        $(toolboxdir)/macros/asserterr.sci \
        $(toolboxdir)/macros/assert_field.bin \
        $(toolboxdir)/macros/assert_field.sci \
        $(toolboxdir)/macros/assert.sci \
        $(toolboxdir)/macros/buildmacros.sce \
        $(toolboxdir)/macros/build_options_list.bin \
        $(toolboxdir)/macros/build_options_list.sci \
        $(toolboxdir)/macros/cart2pol.bin \
        $(toolboxdir)/macros/cart2pol.sci \
        $(toolboxdir)/macros/champ3.bin \
        $(toolboxdir)/macros/champ3.sci \
        $(toolboxdir)/macros/cross.bin \
        $(toolboxdir)/macros/cross.sci \
        $(toolboxdir)/macros/dot.bin \
        $(toolboxdir)/macros/dot.sci \
        $(toolboxdir)/macros/gf_asm_pdetoolbc.bin \
        $(toolboxdir)/macros/gf_asm_pdetoolbc.sci \
        $(toolboxdir)/macros/gfassert.bin \
        $(toolboxdir)/macros/gfassert.sci \
        $(toolboxdir)/macros/gf_colormap.bin \
        $(toolboxdir)/macros/gf_colormap.sci \
        $(toolboxdir)/macros/gf_compute_Q1grid_interp.bin \
        $(toolboxdir)/macros/gf_compute_Q1grid_interp.sci \
        $(toolboxdir)/macros/gf_interpolate_on_grid.bin \
        $(toolboxdir)/macros/gf_interpolate_on_grid.sci \
        $(toolboxdir)/macros/gf_mesh_fem_get_eval.bin \
        $(toolboxdir)/macros/gf_mesh_fem_get_eval.sci \
        $(toolboxdir)/macros/gf_plot_1D.bin \
        $(toolboxdir)/macros/gf_plot_1D.sci \
        $(toolboxdir)/macros/gf_plot.bin \
        $(toolboxdir)/macros/gf_plot_mesh.bin \
        $(toolboxdir)/macros/gf_plot_mesh.sci \
        $(toolboxdir)/macros/gf_plot.sci \
        $(toolboxdir)/macros/gf_plot_slice.bin \
        $(toolboxdir)/macros/gf_plot_slice.sci \
        $(toolboxdir)/macros/gf_solve.bin \
        $(toolboxdir)/macros/gf_solve.sci \
        $(toolboxdir)/macros/has_field.bin \
        $(toolboxdir)/macros/has_field.sci \
        $(toolboxdir)/macros/init_pde.bin \
        $(toolboxdir)/macros/init_pde.sci \
        $(toolboxdir)/macros/isauto.bin \
        $(toolboxdir)/macros/isauto.sci \
        $(toolboxdir)/macros/isnumeric.bin \
        $(toolboxdir)/macros/isnumeric.sci \
        $(toolboxdir)/macros/ison.bin \
        $(toolboxdir)/macros/ison.sci \
        $(toolboxdir)/macros/isscalar.bin \
        $(toolboxdir)/macros/isscalar.sci \
        $(toolboxdir)/macros/lib \
        $(toolboxdir)/macros/names \
        $(toolboxdir)/macros/null_space.bin \
        $(toolboxdir)/macros/null_space.sci \
        $(toolboxdir)/macros/repmat.bin \
        $(toolboxdir)/macros/repmat.sci \
        $(toolboxdir)/macros/_setdiff.bin \
        $(toolboxdir)/macros/_setdiff.sci \
        $(toolboxdir)/macros/spdiags.bin \
        $(toolboxdir)/macros/spdiags.sci \
        $(toolboxdir)/macros/surfnorm.bin \
        $(toolboxdir)/macros/surfnorm.sci \
        $(toolboxdir)/macros/test_champ3.sce \
        $(toolboxdir)/macros/test_spdiags.sce \
        $(toolboxdir)/macros/overload/gf_typeof.bin \
        $(toolboxdir)/macros/overload/gf_typeof.sci \
        $(toolboxdir)/macros/overload/init_gf_types.sce \
        $(toolboxdir)/macros/overload/lib \
        $(toolboxdir)/macros/overload/names \
        $(toolboxdir)/macros/overload/%objid_e.bin \
        $(toolboxdir)/macros/overload/%objid_e.sci \
        $(toolboxdir)/macros/overload/%objid_get.bin \
        $(toolboxdir)/macros/overload/%objid_get.sci \
        $(toolboxdir)/macros/overload/%objid_set.bin \
        $(toolboxdir)/macros/overload/%objid_set.sci \
        $(toolboxdir)/sci_gateway/builder_gateway.sce \
        $(toolboxdir)/sci_gateway/cleaner_gateway.sce \
        $(toolboxdir)/sci_gateway/loader_gateway.sce \
        $(toolboxdir)/sci_gateway/c/builder_gateway_c.sce \
        $(toolboxdir)/sci_gateway/c/builder_gateway_c.sce.in \
        $(toolboxdir)/sci_gateway/c/cleaner.sce \
        $(toolboxdir)/sci_gateway/c/gfm_common.c \
        $(toolboxdir)/sci_gateway/c/gfm_common.h \
        $(toolboxdir)/sci_gateway/c/gfm_scilab.c \
        $(toolboxdir)/sci_gateway/c/libscigetfem_c.c \
        $(toolboxdir)/sci_gateway/c/libscigetfem_c.so* \
        $(toolboxdir)/sci_gateway/c/loader.sce \
        $(toolboxdir)/sci_gateway/c/sci_cgne.c \
        $(toolboxdir)/sci_gateway/c/sci_cgs.c \
        $(toolboxdir)/sci_gateway/c/sci_gmres.c \
        $(toolboxdir)/sci_gateway/c/sci_mgcr.c \
        $(toolboxdir)/sci_gateway/c/sci_spchol.c \
        $(toolboxdir)/sci_gateway/c/sci_spcholinc.c \
        $(toolboxdir)/sci_gateway/c/sci_spchsolve.c \
        $(toolboxdir)/sci_gateway/c/sci_splu.c \
        $(toolboxdir)/sci_gateway/c/sci_spluinc.c \
        $(toolboxdir)/sci_gateway/c/sci_splusolve.c \
        $(toolboxdir)/src/builder_src.sce \
        $(toolboxdir)/src/c/bdfactor.c \
        $(toolboxdir)/src/c/bkpfacto.c \
        $(toolboxdir)/src/c/builder_c.sce \
        $(toolboxdir)/src/c/chfactor.c \
        $(toolboxdir)/src/c/cleaner.sce \
        $(toolboxdir)/src/c/configure* \
        $(toolboxdir)/src/c/configure.in \
        $(toolboxdir)/src/c/copy.c \
        $(toolboxdir)/src/c/copyright \
        $(toolboxdir)/src/c/err.c \
        $(toolboxdir)/src/c/err.h \
        $(toolboxdir)/src/c/extras.c \
        $(toolboxdir)/src/c/fft.c \
        $(toolboxdir)/src/c/FILELIST \
        $(toolboxdir)/src/c/givens.c \
        $(toolboxdir)/src/c/hessen.c \
        $(toolboxdir)/src/c/hsehldr.c \
        $(toolboxdir)/src/c/init.c \
        $(toolboxdir)/src/c/iter0.c \
        $(toolboxdir)/src/c/iter.h \
        $(toolboxdir)/src/c/iternsym.c \
        $(toolboxdir)/src/c/itersym.c \
        $(toolboxdir)/src/c/ivecop.c \
        $(toolboxdir)/src/c/libsp_get.so* \
        $(toolboxdir)/src/c/loader.sce \
        $(toolboxdir)/src/c/ls.dat \
        $(toolboxdir)/src/c/lufactor.c \
        $(toolboxdir)/src/c/machine.c \
        $(toolboxdir)/src/c/machine.h \
        $(toolboxdir)/src/c/machine.h.in \
        $(toolboxdir)/src/c/makefile \
        $(toolboxdir)/src/c/makefile.in \
        $(toolboxdir)/src/c/matlab.c \
        $(toolboxdir)/src/c/matlab.h \
        $(toolboxdir)/src/c/matop.c \
        $(toolboxdir)/src/c/matrix2.h \
        $(toolboxdir)/src/c/matrix.h \
        $(toolboxdir)/src/c/matrixio.c \
        $(toolboxdir)/src/c/meminfo.c \
        $(toolboxdir)/src/c/meminfo.h \
        $(toolboxdir)/src/c/memory.c \
        $(toolboxdir)/src/c/memstat.c \
        $(toolboxdir)/src/c/mfunc.c \
        $(toolboxdir)/src/c/norm.c \
        $(toolboxdir)/src/c/oldnames.h \
        $(toolboxdir)/src/c/otherio.c \
        $(toolboxdir)/src/c/pxop.c \
        $(toolboxdir)/src/c/qrfactor.c \
        $(toolboxdir)/src/c/README \
        $(toolboxdir)/src/c/rk4.dat \
        $(toolboxdir)/src/c/schur.c \
        $(toolboxdir)/src/c/solve.c \
        $(toolboxdir)/src/c/sparse2.h \
        $(toolboxdir)/src/c/sparse.c \
        $(toolboxdir)/src/c/sparse.h \
        $(toolboxdir)/src/c/sparseio.c \
        $(toolboxdir)/src/c/spbkp.c \
        $(toolboxdir)/src/c/spchfctr.c \
        $(toolboxdir)/src/c/splufctr.c \
        $(toolboxdir)/src/c/sprow.c \
        $(toolboxdir)/src/c/spswap.c \
        $(toolboxdir)/src/c/submat.c \
        $(toolboxdir)/src/c/svd.c \
        $(toolboxdir)/src/c/symmeig.c \
        $(toolboxdir)/src/c/update.c \
        $(toolboxdir)/src/c/vecop.c \
        $(toolboxdir)/src/c/version.c \
        $(toolboxdir)/src/c/zcopy.c \
        $(toolboxdir)/src/c/zfunc.c \
        $(toolboxdir)/src/c/zgivens.c \
        $(toolboxdir)/src/c/zhessen.c \
        $(toolboxdir)/src/c/zhsehldr.c \
        $(toolboxdir)/src/c/zlufctr.c \
        $(toolboxdir)/src/c/zmachine.c \
        $(toolboxdir)/src/c/zmatio.c \
        $(toolboxdir)/src/c/zmatlab.c \
        $(toolboxdir)/src/c/zmatop.c \
        $(toolboxdir)/src/c/zmatrix2.h \
        $(toolboxdir)/src/c/zmatrix.h \
        $(toolboxdir)/src/c/zmemory.c \
        $(toolboxdir)/src/c/znorm.c \
        $(toolboxdir)/src/c/zqrfctr.c \
        $(toolboxdir)/src/c/zschur.c \
        $(toolboxdir)/src/c/zsolve.c \
        $(toolboxdir)/src/c/zvecop.c \
        $(toolboxdir)/src/c/DOC/fnindex.txt \
        $(toolboxdir)/src/c/DOC/tutorial.txt \
        $(toolboxdir)/src/c/MACHINES/Cray/machine.h \
        $(toolboxdir)/src/c/MACHINES/Cray/makefile \
        $(toolboxdir)/src/c/MACHINES/Cray/patch.1 \
        $(toolboxdir)/src/c/MACHINES/Cray/patch.2 \
        $(toolboxdir)/src/c/MACHINES/Cray/patch.3 \
        $(toolboxdir)/src/c/MACHINES/GCC/machine.h \
        $(toolboxdir)/src/c/MACHINES/GCC/makefile \
        $(toolboxdir)/src/c/MACHINES/Linux/machine.h \
        $(toolboxdir)/src/c/MACHINES/Linux/makefile \
        $(toolboxdir)/src/c/MACHINES/RS6000/machine.c \
        $(toolboxdir)/src/c/MACHINES/RS6000/machine.h \
        $(toolboxdir)/src/c/MACHINES/RS6000/makefile \
        $(toolboxdir)/src/c/MACHINES/SGI/machine.h \
        $(toolboxdir)/src/c/MACHINES/SGI/makefile \
        $(toolboxdir)/src/c/MACHINES/SPARC/machine.h \
        $(toolboxdir)/src/c/MACHINES/SPARC/makefile \
        $(toolboxdir)/tests \
        $(toolboxdir)/tests/nonreg_tests \
        $(toolboxdir)/tests/unit_tests 

.NOTPARALLEL: *

all:
	@SCILAB_EXE@ -nw -nb -f $(scilabbuilddir)/makefile_builder.sce

install:
	$(mkinstalldirs) $(toolboxdir)/demos
	$(mkinstalldirs) $(toolboxdir)/etc
	$(mkinstalldirs) $(toolboxdir)/help
	$(mkinstalldirs) $(toolboxdir)/help/en_US
	$(mkinstalldirs) $(toolboxdir)/help/en_US/examples
	$(mkinstalldirs) $(toolboxdir)/help/fr_FR
	$(mkinstalldirs) $(toolboxdir)/help/fig
	$(mkinstalldirs) $(toolboxdir)/help/mml
	$(mkinstalldirs) $(toolboxdir)/jar
	$(mkinstalldirs) $(toolboxdir)/macros
	$(mkinstalldirs) $(toolboxdir)/sci_gateway
	$(mkinstalldirs) $(toolboxdir)/sci_gateway/c
	$(mkinstalldirs) $(toolboxdir)/src
	$(mkinstalldirs) $(toolboxdir)/src/c
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/Cray
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/GCC
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/Linux
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/RS6000
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/SGI
	$(mkinstalldirs) $(toolboxdir)/src/c/MACHINES/SPARC
	$(mkinstalldirs) $(toolboxdir)/src/c/DOC
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos demos/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/etc etc/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help *.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fig help/fig/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/mml help/mml/*.mml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/jar jar/*.jar
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/lib
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway sci_gateway/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.[ch]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.sce
	@INSTALL@ -D -m 744 -t $(toolboxdir)/sci_gateway/c sci_gateway/c/*.so
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src src/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c src/c/*.[ch]
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c src/c/*.so
	@INSTALL@ -D -m 744 -t $(toolboxdir)/src/c $(scilabsrccdir)
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Cray src/c/MACHINES/Cray/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/GCC src/c/MACHINES/GCC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/Linux src/c/MACHINES/Linux/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/RS6000 src/c/MACHINES/RS6000/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SGI src/c/MACHINES/SGI/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/MACHINES/SPARC src/c/MACHINES/SPARC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/src/c/DOC src/c/DOC/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ $(scilabbasedir)

clean:
	@SCILAB_EXE@ -nw -nb -f $(scilabbuilddir)/makefile_cleaner.sce
