SUBDIRS = @gfCvStruct @gfEltm @gfFem @gfGeoTrans @gfInteg @gfMesh @gfMeshFem @gfMeshIm @gfSpmat @gfPrecond @gfSlice @gfMdBrick @gfMdState @gfModel @gfLevelSet @gfMeshLevelSet private

RPC_INC_DIR=@RPC_INC_DIR@

PSEUDO_MFUNCTIONS = \
	gf_util.cc \
	gf_cvstruct_get.cc \
	gf_geotrans.cc \
	gf_geotrans_get.cc \
	gf_compute.cc \
	gf_mesh_fem.cc \
	gf_mesh_fem_set.cc \
	gf_mesh_fem_get.cc \
	gf_mesh_im.cc \
	gf_mesh_im_set.cc \
	gf_mesh_im_get.cc \
	gf_eltm.cc \
	gf_mesh.cc \
	gf_mesh_set.cc \
	gf_mesh_get.cc \
	gf_mdbrick.cc \
	gf_mdbrick_get.cc \
	gf_mdbrick_set.cc \
	gf_mdstate.cc \
	gf_mdstate_get.cc \
	gf_mdstate_set.cc \
	gf_model.cc \
	gf_model_get.cc \
	gf_model_set.cc \
	gf_slice.cc \
	gf_slice_get.cc \
	gf_slice_set.cc \
	gf_levelset.cc \
	gf_levelset_get.cc \
	gf_levelset_set.cc \
	gf_mesh_levelset.cc \
	gf_mesh_levelset_get.cc \
	gf_mesh_levelset_set.cc \
	gf_precond.cc \
	gf_precond_get.cc \
	gf_linsolve.cc \
	gf_spmat.cc \
	gf_spmat_set.cc \
	gf_spmat_get.cc \
	gf_asm.cc \
	gf_fem.cc \
	gf_fem_get.cc \
	gf_integ.cc \
	gf_integ_get.cc \
	gf_workspace.cc \
	gf_delete.cc

AUTO_M_FILES = $(PSEUDO_MFUNCTIONS:.cc=.m)
$(AUTO_M_FILES): @srcdir@/../../bin/extract_mdoc

EXTRA_DIST = gfm_rpc_mexint.c gfm_mex.c gfm_common.c gfm_common.h $(M_FILES) $(AUTO_M_FILES)

all: $(AUTO_M_FILES) gf_matlab@MATLAB_COM_EXT@

MEX=@MEX@
RPC_LIB = @RPC_LIB@
GNUMEX=@top_srcdir@/gnumex
GNUMEXOPTS=@top_srcdir@/gnumex.opts

#gfi_array.c: ../gfi_array.c
#	ln -s ../gfi_array.c .

GETFEM_LIB_LA = ../../../src/libgetfem.la
GETFEM_STATIC_LIB = ../../../src/.libs/libgetfem.a @LIBS@

if BUILDMEX
if USE_MINGW_MEX

#command extremely sensitive to any modification! fragile! keep the order of the files
# (gfm_mex.c must be first, libstdc++.a must be last)
#virer le -DMATLAB_RELEASE qui marche pas..
#-DMATLAB_RELEASE=@MATLAB_RELEASE@ \ the great windows mex does not understand -D ...
gf_matlab@MATLAB_COM_EXT@: gfm_mex.c gfm_common.c ../libgetfemint.la ../gfi_array.c $(GETFEM_LIB_LA)
	$(GNUMEX) $(GNUMEXOPTS) -output gf_matlab -g @srcdir@/gfm_mex.c \
	@srcdir@/gfm_common.c -I@srcdir@ \
	@srcdir@/../gfi_array.c ../.libs/libgetfemint.a $(GETFEM_STATIC_LIB) @STDCPP_STATICLIBS@
#        /c/MinGW/lib/libstdc++.a
#	cmd /c "$mexbat -v -f c:/gnumex/mexopts.bat gfm_mex.c -output gfm_rpc_mexint gfi*.o gf_*.o matlabint*.o c:\\msys\\1.0\\home\\j\\getfem++-1.5\\src\\.libs\\libgetfem.a getfem_matlab.o c:\\mingw\\lib\\libstdc++.a -Ic:\\msys\\1.0\\home\\j\\mingw_liboncrpc-4.0"
else !USE_MINGW_MEX
if BUILDMEXRPC
gf_matlab@MATLAB_COM_EXT@: ../gfi_rpc_clnt.c gfm_rpc_mexint.c gfm_common.c ../gfi_rpc_xdr.c ../gfi_array.c
	$(MEX) -largeArrayDims -output gf_matlab -g CDEBUGFLAGS="$(CFLAGS)" $(RPC_LIB) \
	-I@srcdir@ -I@srcdir@/.. -DMATLAB_RELEASE=@MATLAB_RELEASE@ -DUSE_RPC \
	@srcdir@/gfm_rpc_mexint.c @srcdir@/gfm_common.c @srcdir@/../gfi_rpc_clnt.c \
        @srcdir@/../gfi_rpc_xdr.c @srcdir@/../gfi_array.c || (rm $@; false)
else !BUILDMEXRPC
# 2006/02/06 I remove the @STDCPP_STATICLIBS@ at the end (added to
# avoid crashes in exception throw code when parts of getfem where
# compiled with ifc (i.e. mumps)) of the command line, as it 
# break the linking with g++-3.3 and matlab R14/R13 on debian (at least) ..
gf_matlab@MATLAB_COM_EXT@: gfm_mex.c gfm_common.c ../libgetfemint.la ../gfi_array.c $(GETFEM_LIB_LA)
	$(MEX) -largeArrayDims -output gf_matlab -g CDEBUGFLAGS="$(CFLAGS)" LD="$(CXX)" \
	-I@srcdir@ -I@srcdir@/.. -DMATLAB_RELEASE=@MATLAB_RELEASE@ \
	@srcdir@/gfm_mex.c @srcdir@/gfm_common.c @srcdir@/../gfi_array.c \
	../.libs/libgetfemint.a $(GETFEM_STATIC_LIB) @STDCPP_STATICLIBS@ || (rm $@; false)
endif !BUILDMEXRPC
endif !USE_MINGW_MEX
endif 

# -largeArrayDims

M_FILES = \
	Contents.m \
	gf_asm_pdetoolbc.m \
	gf_compute_Q1grid_interp.m \
	gf_mesh_fem_get_eval.m \
	gf_plot.m \
	gf_plot_1D.m \
	gf_plot_mesh.m \
	gf_plot_slice.m \
	gf_solve.m \
	gf_colormap.m \
	gfObject.m

#$(AUTO_M_FILES): $(PSEUDO_MFUNCTIONS:%.cc=../%.cc)

.NOTPARALLEL: $(AUTO_M_FILES) $(M_FILES)

%.m: ../%.cc
	@@srcdir@/../../bin/extract_mdoc $<

clean-m-files:
	@echo "cleaning generated m-files"
	rm -f $(AUTO_M_FILES)

clean-local: #clean-m-files
	rm -f gf_matlab@MATLAB_COM_EXT@

toolboxdir=@TOOLBOXDIR@
toolbox_SCRIPTS=gf_matlab@MATLAB_COM_EXT@ $(AUTO_M_FILES) $(M_FILES)
