SUBDIRS = private

RPC_INC_DIR=@RPC_INC_DIR@

PSEUDO_MFUNCTIONS=$(shell ls $(srcdir)/../gf_*.cc)
PSEUDO_MFUNCTIONS_LOC=$(shell cd $(srcdir)/.. && ls gf_*.cc)
AUTO_M_FILES = $(PSEUDO_MFUNCTIONS_LOC:.cc=.m)

# $(warning PSEUDO_MFUNCTIONS= $(PSEUDO_MFUNCTIONS))



gf_mesh.m : $(PSEUDO_MFUNCTIONS) @srcdir@/../../bin/extract_doc
	@srcdir@/../../bin/extract_doc @srcdir@/.. matlab-com || (rm -f gf_mesh.m; /bin/false )


EXTRA_DIST = gfm_rpc_mexint.c gfm_mex.c gfm_common.c gfm_common.h $(M_FILES)

all: gf_mesh.m gf_matlab@MATLAB_COM_EXT@

MEX=@MEX@
RPC_LIB = @RPC_LIB@
GNUMEX=@top_srcdir@/gnumex
GNUMEXOPTS=@top_srcdir@/gnumex.opts

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


.NOTPARALLEL: $(M_FILES)

#%.m: ../%.cc
#	@@srcdir@/../../bin/extract_mdoc $<

clean-m-files:
	@echo "cleaning auto generated m-files and directories"
	rm -f $(AUTO_M_FILES)
	rm -fr \@gf* 

clean-local: clean-m-files
	rm -f gf_matlab@MATLAB_COM_EXT@

toolboxdir=@TOOLBOXDIR@

install:
	$(mkinstalldirs) $(toolboxdir)/
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ *.m
	@INSTALL@ -D -m 644 -t $(toolboxdir)/ $(srcdir)/*.m
	@list='$(shell ls -d @gf*)'; for p in $$list; do \
	  $(mkinstalldirs) $(toolboxdir)/$$p; \
	  @INSTALL@ -D -m 644 -t $(toolboxdir)/$$p $$p/*.m; \
	done
	@INSTALL@ -D -m 744 -t $(toolboxdir)/ gf_matlab@MATLAB_COM_EXT@

uninstall:
	rm -fr $(toolboxdir)