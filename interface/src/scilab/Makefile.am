## Scilab Getfem Interface
## Copyright (C) 2009 Yann COLLETTE
##
## This scilab getfem interface is a part of GetFEM++
##
## GefFEM++ is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
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
	builddoc.sce \
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

EXTRA_DIST=builddoc.sce \
        builder.sce \
        changelog.txt \
        cleaner.sce \
        license.txt \
        loader.sce \
        readme.txt \
	sci_getfem.iss \
	makefile_builder.sce \
        demos/axrot_matrix.sci \
        tests/unit_tests/check_all.sce \
        tests/unit_tests/check_asm.sce \
        tests/unit_tests/check_fem.sce \
        tests/unit_tests/check_geotrans.sce \
        tests/unit_tests/check_integ.sce \
        tests/unit_tests/check_interpolated_fem.sce \
        tests/unit_tests/check_levelset.sce \
        tests/unit_tests/check_mesh_fem.sce \
        tests/unit_tests/check_oo.sce \
        tests/unit_tests/check_plot.sce \
        tests/unit_tests/check_slices.sce \
        tests/unit_tests/check_spmat.sce \
        tests/unit_tests/check_workspace.sce \
        demos/demo_bilaplacian.sce \
	demos/demo_continuation_block.sce \
        demos/demo_continuation.sce \
	demos/demo_continuation_vee.sce \
        demos/demo_convection_rotating_cavity.sce \
        demos/demo_crack.sce \
        demos/demo_fictitious_domains.sce \
        demos/demo_laplacian.sce \
        demos/demo_mesh_generation.sce \
        demos/demo_mortar.sce \
        demos/demo_thermo_elasticity_electrical_coupling.sce \
        demos/demo_nonlinear_elasticity_anim.sce \
        demos/demo_nonlinear_elasticity.sce \
        demos/demo_plasticity.sce \
        demos/demo_plate.sce \
        demos/demo_refine.sce \
        demos/demo_static_contact.sce \
        demos/demo_step_by_step.sce \
        demos/demo_stokes_3D_tank_draw.sce \
        demos/demo_stokes_3D_tank.sce \
        demos/demo_structural_optimization.sce \
        demos/demo_topological_optimization.sce \
        demos/demo_tripod_alt.sce \
        demos/demo_tripod_anim.sce \
        demos/demo_tripod.sce \
        demos/demo_tripod_slice_anim.sce \
        demos/demo_wave2D_alt.sce \
        demos/demo_wave2D_animate.sce \
        demos/demo_wave2D.sce \
        demos/demo_wave_equation.sce \
        demos/sci_getfem.dem.gateway.sce \
        tests/unit_tests/test_argyris.sce \
        tests/unit_tests/test_plasticity_new_brick.sce \
        demos/tutorial1.sce \
        demos/data/disc_P2_h0.5.mesh \
        demos/data/disc_P2_h1.mesh \
        demos/data/disc_P2_h2.mesh \
        demos/data/donut_regulier.mesh \
        demos/data/donut_with_quadratic_tetra_1100_elements.msh \
        demos/data/holed_disc_with_quadratic_2D_triangles.msh \
        demos/data/quad.geo \
        demos/data/quad.msh \
        demos/data/tank_quadratic_2500.GiD.msh \
        demos/data/tripod.GiD.msh \
        demos/data/tripod.mesh \
        demos/data/tube_2D_spline.GiD.msh \
	demos/data/vee_h_0.03.mesh \
        etc/sci_getfem.quit \
        etc/sci_getfem.start \
        help/builder_help.sce \
        help/en_US/build_help.sce \
        help/en_US/gf_asm_pdetoolbc.xml \
        help/en_US/gf_asm.xml \
        help/en_US/gf_colormap.xml \
        help/en_US/gf_compute_Q1grid_interp.xml \
        help/en_US/gf_compute.xml \
        help/en_US/gf_cvstruct_get.xml \
        help/en_US/gf_delete.xml \
        help/en_US/gf_eltm.xml \
        help/en_US/gf_fem_get.xml \
        help/en_US/gf_fem.xml \
        help/en_US/gf_geotrans_get.xml \
        help/en_US/gf_geotrans.xml \
        help/en_US/gf_global_function_get.xml \
        help/en_US/gf_global_function.xml \
        help/en_US/gf_integ_get.xml \
        help/en_US/gf_integ.xml \
        help/en_US/gf_interpolate_on_grid.xml \
        help/en_US/gf_levelset_get.xml \
        help/en_US/gf_levelset_set.xml \
        help/en_US/gf_levelset.xml \
        help/en_US/gf_linsolve.xml \
        help/en_US/gf_mesh_fem_get_eval.xml \
        help/en_US/gf_mesh_fem_get.xml \
        help/en_US/gf_mesh_fem_set.xml \
        help/en_US/gf_mesh_fem.xml \
        help/en_US/gf_mesh_get.xml \
        help/en_US/gf_mesh_im_get.xml \
        help/en_US/gf_mesh_im_set.xml \
        help/en_US/gf_mesh_im.xml \
        help/en_US/gf_mesh_levelset_get.xml \
        help/en_US/gf_mesh_levelset_set.xml \
        help/en_US/gf_mesh_levelset.xml \
        help/en_US/gf_mesh_set.xml \
        help/en_US/gf_mesh.xml \
        help/en_US/gf_model_get.xml \
        help/en_US/gf_model_set.xml \
        help/en_US/gf_model.xml \
        help/en_US/gf_plot_1D.xml \
        help/en_US/gf_plot_mesh.xml \
        help/en_US/gf_plot_slice.xml \
        help/en_US/gf_plot.xml \
        help/en_US/gf_poly.xml \
        help/en_US/gf_precond_get.xml \
        help/en_US/gf_precond.xml \
        help/en_US/gf_slice_get.xml \
        help/en_US/gf_slice_set.xml \
        help/en_US/gf_slice.xml \
        help/en_US/gf_solve.xml \
        help/en_US/gf_spmat_get.xml \
        help/en_US/gf_spmat_set.xml \
        help/en_US/gf_spmat.xml \
        help/en_US/gf_typeof.xml \
        help/en_US/gf_util.xml \
        help/en_US/gf_workspace.xml \
        help/en_US/objects.xml \
        help/en_US/preliminary.xml \
        help/en_US/examples/another_laplacian.xml \
        help/en_US/examples/avoiding_bricks.xml \
        help/en_US/examples/CHAPTER \
        help/en_US/examples/linear_nonlinear_elast.xml \
        help/en_US/examples/step_by_step_example.xml \
        help/en_US/sparses/CHAPTER \
        help/en_US/sparses/sp_cgne.xml \
        help/en_US/sparses/sp_cgs.xml \
        help/en_US/sparses/sp_cholinc.xml \
        help/en_US/sparses/sp_chol.xml \
        help/en_US/sparses/sp_chsolve.xml \
        help/en_US/sparses/sp_gmres.xml \
        help/en_US/sparses/sp_luinc.xml \
        help/en_US/sparses/sp_lusolve.xml \
        help/en_US/sparses/sp_lu.xml \
        help/en_US/sparses/sp_mgcr.xml \
        help/fig/gf_fem_get_fig_1.png \
        help/fig/gf_plot_mesh_fig_1.png \
        help/fig/hierarchy.png \
        help/fig/tripodvonmiseswithmesh_small.png \
        help/fr_FR/build_help.sce \
        help/mml/avoiding_eq1.mml \
        help/mml/gf_asm_eq1.mml \
        help/mml/gf_asm_eq2.mml \
        help/mml/gf_asm_eq3.mml \
        help/mml/gf_asm_eq4.mml \
        help/mml/gf_asm_eq5.mml \
        help/mml/gf_asm_eq6.mml \
        help/mml/gf_model_set_eq1.mml \
        help/mml/gf_model_set_eq2.mml \
        help/mml/gf_model_set_eq3.mml \
        help/mml/gf_model_set_eq4.mml \
        help/mml/preliminary_eq1.mml \
        help/mml/step_by_step_eq1.mml \
        help/mml/step_by_step_eq2.mml \
        help/mml/step_by_step_eq3.mml \
        help/mml/step_by_step_eq4.mml \
        help/mml/step_by_step_eq5.mml \
        jar/scilab_en_US_help.jar \
        macros/add_empty_bound.sci \
        macros/asserterr.sci \
        macros/assert_field.sci \
        macros/assert.sci \
        macros/buildmacros.sce \
        macros/build_options_list.sci \
        macros/cart2pol.sci \
        macros/champ3.sci \
        macros/cross.sci \
        macros/dot.sci \
        macros/gf_asm_pdetoolbc.sci \
        macros/gfassert.sci \
        macros/gf_colormap.sci \
        macros/gf_compute_Q1grid_interp.sci \
        macros/gf_interpolate_on_grid.sci \
        macros/gf_mesh_fem_get_eval.sci \
        macros/gf_plot_1D.sci \
        macros/gf_plot_mesh.sci \
        macros/gf_plot.sci \
        macros/gf_plot_slice.sci \
        macros/gf_solve.sci \
        macros/has_field.sci \
        macros/init_pde.sci \
        macros/isauto.sci \
        macros/isnumeric.sci \
        macros/ison.sci \
        macros/isscalar.sci \
        macros/lib \
        macros/names \
        macros/null_space.sci \
        macros/repmat.sci \
        macros/_setdiff.sci \
        macros/spdiags.sci \
        macros/surfnorm.sci \
        macros/test_champ3.sce \
        macros/test_spdiags.sce \
        macros/overload/gf_typeof.sci \
        macros/overload/init_gf_types.sce \
        macros/overload/lib \
        macros/overload/names \
        macros/overload/%objid_e.sci \
        macros/overload/%objid_get.sci \
        macros/overload/%objid_set.sci \
        sci_gateway/builder_gateway.sce \
        sci_gateway/cleaner_gateway.sce \
        sci_gateway/loader_gateway.sce \
        sci_gateway/c/builder_gateway_c.sce \
        sci_gateway/c/cleaner.sce \
        sci_gateway/c/gfm_common.c \
        sci_gateway/c/gfm_common.h \
        sci_gateway/c/gfm_scilab.cpp \
        sci_gateway/c/libscigetfem_c.c \
        sci_gateway/c/libscigetfem_c.so \
        sci_gateway/c/loader.sce \
        sci_gateway/c/sci_cgne.c \
        sci_gateway/c/sci_cgs.c \
        sci_gateway/c/sci_gmres.c \
        sci_gateway/c/sci_mgcr.c \
        sci_gateway/c/sci_spchol.c \
        sci_gateway/c/sci_spcholinc.c \
        sci_gateway/c/sci_spchsolve.c \
        sci_gateway/c/sci_splu.c \
        sci_gateway/c/sci_spluinc.c \
        sci_gateway/c/sci_splusolve.c \
        src/builder_src.sce \
        src/c/bdfactor.c \
        src/c/bkpfacto.c \
        src/c/builder_c.sce \
        src/c/chfactor.c \
        src/c/cleaner.sce \
        src/c/configure \
        src/c/configure.in \
        src/c/copy.c \
        src/c/copyright \
        src/c/err.c \
        src/c/err.h \
        src/c/extras.c \
        src/c/fft.c \
        src/c/FILELIST \
        src/c/givens.c \
        src/c/hessen.c \
        src/c/hsehldr.c \
        src/c/init.c \
        src/c/iter0.c \
        src/c/iter.h \
        src/c/iternsym.c \
        src/c/itersym.c \
        src/c/ivecop.c \
        src/c/libsp_get.so \
        src/c/loader.sce \
        src/c/ls.dat \
        src/c/lufactor.c \
        src/c/machine.c \
        src/c/machine.h \
        src/c/machine.h.in \
        src/c/makefile \
        src/c/makefile.in \
        src/c/matlab.c \
        src/c/matlab.h \
        src/c/matop.c \
        src/c/matrix2.h \
        src/c/matrix.h \
        src/c/matrixio.c \
        src/c/meminfo.c \
        src/c/meminfo.h \
        src/c/memory.c \
        src/c/memstat.c \
        src/c/mfunc.c \
        src/c/norm.c \
        src/c/oldnames.h \
        src/c/otherio.c \
        src/c/pxop.c \
        src/c/qrfactor.c \
        src/c/README \
        src/c/rk4.dat \
        src/c/schur.c \
        src/c/solve.c \
        src/c/sparse2.h \
        src/c/sparse.c \
        src/c/sparse.h \
        src/c/sparseio.c \
        src/c/spbkp.c \
        src/c/spchfctr.c \
        src/c/splufctr.c \
        src/c/sprow.c \
        src/c/spswap.c \
        src/c/submat.c \
        src/c/svd.c \
        src/c/symmeig.c \
        src/c/update.c \
        src/c/vecop.c \
        src/c/version.c \
        src/c/zcopy.c \
        src/c/zfunc.c \
        src/c/zgivens.c \
        src/c/zhessen.c \
        src/c/zhsehldr.c \
        src/c/zlufctr.c \
        src/c/zmachine.c \
        src/c/zmatio.c \
        src/c/zmatlab.c \
        src/c/zmatop.c \
        src/c/zmatrix2.h \
        src/c/zmatrix.h \
        src/c/zmemory.c \
        src/c/znorm.c \
        src/c/zqrfctr.c \
        src/c/zschur.c \
        src/c/zsolve.c \
        src/c/zvecop.c \
        src/c/DOC/fnindex.txt \
        src/c/DOC/tutorial.txt \
        src/c/MACHINES/Cray/machine.h \
        src/c/MACHINES/Cray/makefile \
        src/c/MACHINES/Cray/patch.1 \
        src/c/MACHINES/Cray/patch.2 \
        src/c/MACHINES/Cray/patch.3 \
        src/c/MACHINES/GCC/machine.h \
        src/c/MACHINES/GCC/makefile \
        src/c/MACHINES/Linux/machine.h \
        src/c/MACHINES/Linux/makefile \
        src/c/MACHINES/RS6000/machine.c \
        src/c/MACHINES/RS6000/machine.h \
        src/c/MACHINES/RS6000/makefile \
        src/c/MACHINES/SGI/machine.h \
        src/c/MACHINES/SGI/makefile \
        src/c/MACHINES/SPARC/machine.h \
        src/c/MACHINES/SPARC/makefile

.NOTPARALLEL: *

all:
	@SCILAB_EXE@ -nw -nb -f $(scilabbuilddir)/makefile_builder.sce

install:
	$(mkinstalldirs) $(toolboxdir)/demos
	$(mkinstalldirs) $(toolboxdir)/demos/data
	$(mkinstalldirs) $(toolboxdir)/etc
	$(mkinstalldirs) $(toolboxdir)/help
	$(mkinstalldirs) $(toolboxdir)/help/en_US
	$(mkinstalldirs) $(toolboxdir)/help/en_US/examples
	$(mkinstalldirs) $(toolboxdir)/help/en_US/sparses
	$(mkinstalldirs) $(toolboxdir)/help/fr_FR
	$(mkinstalldirs) $(toolboxdir)/help/fig
	$(mkinstalldirs) $(toolboxdir)/help/mml
	$(mkinstalldirs) $(toolboxdir)/jar
	$(mkinstalldirs) $(toolboxdir)/macros
	$(mkinstalldirs) $(toolboxdir)/macros/overload
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
	$(mkinstalldirs) $(toolboxdir)/tests
	$(mkinstalldirs) $(toolboxdir)/tests/unit_tests
	$(mkinstalldirs) $(toolboxdir)/tests/nonreg_tests
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos demos/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/demos/data demos/data/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/etc etc/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help help/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US help/en_US/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/examples help/en_US/examples/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/sparses help/en_US/sparses/CHAPTER
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/en_US/sparses help/en_US/sparses/*.xml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fr_FR help/fr_FR/*.sce
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/fig help/fig/*
	@INSTALL@ -D -m 644 -t $(toolboxdir)/help/mml help/mml/*.mml
	@INSTALL@ -D -m 644 -t $(toolboxdir)/jar jar/*.jar
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros macros/lib
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros/overload macros/overload/*.bin
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros/overload macros/overload/*.sc[ie]
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros/overload macros/overload/names
	@INSTALL@ -D -m 644 -t $(toolboxdir)/macros/overload macros/overload/lib
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
