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

RPC_INC_DIR=@RPC_INC_DIR@


AUTO_M_FILES = $(PSEUDO_FUNCTIONS_LOC:.cc=.m)

# $(warning PSEUDO_FUNCTIONS= $(PSEUDO_FUNCTIONS))



gf_mesh.m : ../libgetfemint.la $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc @srcdir@/.. matlab-com || (rm -f gf_mesh.m; /bin/false )


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

gf_matlab@MATLAB_COM_EXT@: gfm_mex.c gfm_common.c ../libgetfemint.la ../gfi_array.c $(GETFEM_LIB_LA)
	matlab -nodesktop -nosplash -nojvm -r "mex -v -output gf_matlab -g ./gfm_mex.c ./gfm_common.c -largeArrayDims -I. -I.. ./../gfi_array.c ../.libs/libgetfemint.a ../../../src/.libs/libgetfem.a ../../../superlu/.libs/libsuperlu.a /msys/local/lib/libsmumps.a /msys/local/lib/libcmumps.a /msys/local/lib/libdmumps.a /msys/local/lib/libzmumps.a /msys/local/lib/libmumps_common.a /msys/local/lib/libmpiseq.a /msys/local/lib/libpord.a /msys/local/lib/liblapack.a /msys/local/lib/libblas.a /msys/local/lib/libqhull.a /Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libgfortran.a c:/Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libstdc++.a c:/Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libquadmath.a; exit(0);"


#	$(GNUMEX) $(GNUMEXOPTS) -output gf_matlab -g @srcdir@/gfm_mex.c \
#	@srcdir@/gfm_common.c -I@srcdir@ \
#	@srcdir@/../gfi_array.c ../.libs/libgetfemint.a $(GETFEM_STATIC_LIB) @STDCPP_STATICLIBS@


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
	gf_colormap.m \
	gfObject.m


.NOTPARALLEL: $(M_FILES)

clean-m-files:
	@echo "cleaning auto generated m-files and directories"
	rm -f $(AUTO_M_FILES)
	rm -fr \@gf* 

clean-local: clean-m-files
	rm -f gf_matlab@MATLAB_COM_EXT@

toolboxdir=@TOOLBOXDIR@

install:
	$(mkinstalldirs) $(toolboxdir)/private
	@INSTALL@ -m 644 -t $(toolboxdir)/ *.m
	@INSTALL@ -m 644 -t $(toolboxdir)/ $(srcdir)/*.m
	@INSTALL@ -m 644 -t $(toolboxdir)/private/ $(srcdir)/private/*.m
	@list='$(MATLAB_OBJ_DIRS)'; for p in $$list; do \
	  $(mkinstalldirs) $(toolboxdir)/$$p; \
	  @INSTALL@ -m 644 -t $(toolboxdir)/$$p $$p/*.m; \
	done
	@INSTALL@ -m 744 -t $(toolboxdir)/ gf_matlab@MATLAB_COM_EXT@

uninstall:
	rm -fr $(toolboxdir)
