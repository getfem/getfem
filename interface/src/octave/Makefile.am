#  Copyright (C) 1999-2017 Yves Renard
#
#  This file is a part of GetFEM++
#
#  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

GETFEM_LIB_LA = ../../../src/libgetfem.la
GETFEM_STATIC_LIB = ../../../src/.libs/libgetfem.a @LIBS@

SUBDIRS = private

AUTO_M_FILES = $(PSEUDO_FUNCTIONS_LOC:.cc=.m)

# $(warning PSEUDO_FUNCTIONS= $(PSEUDO_FUNCTIONS))



gf_mesh.m : ../libgetfemint.la $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc @srcdir@/.. octave-com || (rm -f gf_mesh.m; /bin/false )


EXTRA_DIST = gfm_rpc_mexint.c gfm_mex.c gfm_common.c gfm_common.h $(M_FILES)

all: gf_mesh.m gf_octave@OCTAVE_COM_EXT@

if BUILDOCTAVE
if USE_MINGW_OCTAVE

# Command extremely sensitive to any modification! fragile!
# keep the order of the files
# (gfm_mex.c must be first, libstdc++.a must be last)

gf_octave@OCTAVE_COM_EXT@: gfm_mex.c gfm_common.c ../libgetfemint.la ../gfi_array.c $(GETFEM_LIB_LA)
	octave -nodesktop -nosplash -nojvm -r "mex -v -output gf_octave -g ./gfm_mex.c ./gfm_common.c -largeArrayDims -I. -I.. ./../gfi_array.c ../.libs/libgetfemint.a ../../../src/.libs/libgetfem.a ../../../superlu/.libs/libsuperlu.a /msys/local/lib/libsmumps.a /msys/local/lib/libcmumps.a /msys/local/lib/libdmumps.a /msys/local/lib/libzmumps.a /msys/local/lib/libmumps_common.a /msys/local/lib/libmpiseq.a /msys/local/lib/libpord.a /msys/local/lib/liblapack.a /msys/local/lib/libblas.a /msys/local/lib/libqhull.a /Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libgfortran.a c:/Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libstdc++.a c:/Mingw64/lib/gcc/x86_64-w64-mingw32/4.9.1/libquadmath.a; exit(0);"

#	$(GNUMEX) $(GNUMEXOPTS) -output gf_octave -g @srcdir@/gfm_mex.c \
#	@srcdir@/gfm_common.c -I@srcdir@ \
#	@srcdir@/../gfi_array.c ../.libs/libgetfemint.a $(GETFEM_STATIC_LIB) @STDCPP_STATICLIBS@

else !USE_MINGW_OCTAVE
gf_octave@OCTAVE_COM_EXT@: gfm_mex.c gfm_common.c ../libgetfemint.la ../gfi_array.c $(GETFEM_LIB_LA)
	DL_LD=$(CXX) mkoctfile -v --mex -output gf_octave -g $(LDFLAGS) \
	-I@srcdir@ -I@srcdir@/.. -DOCTAVE_RELEASE=@OCTAVE_RELEASE@ \
	@srcdir@/gfm_mex.c @srcdir@/gfm_common.c @srcdir@/../gfi_array.c \
	../.libs/libgetfemint.a $(GETFEM_STATIC_LIB) @STDCPP_STATICLIBS@ || (rm $@; false)
endif !USE_MINGW_OCTAVE
endif

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
	rm -f gf_octave@OCTAVE_COM_EXT@

toolboxdir=@OCTAVETOOLBOXDIR@

install:
	$(mkinstalldirs) $(toolboxdir)/private
	@INSTALL@ -m 644 -t $(toolboxdir)/ *.m
	@INSTALL@ -m 644 -t $(toolboxdir)/ $(srcdir)/*.m
	@INSTALL@ -m 644 -t $(toolboxdir)/private/ $(srcdir)/private/*.m
	@list='$(MATLAB_OBJ_DIRS)'; for p in $$list; do \
	  $(mkinstalldirs) $(toolboxdir)/$$p; \
	  @INSTALL@ -m 644 -t $(toolboxdir)/$$p $$p/*.m; \
	done
	@INSTALL@ -m 744 -t $(toolboxdir)/ gf_octave@OCTAVE_COM_EXT@

uninstall:
	rm -fr $(toolboxdir)
