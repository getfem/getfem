#  Copyright (C) 1999-2020 Yves Renard
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

if BUILDMEX
subdirMATLAB=matlab
endif

if BUILDOCTAVE
subdirOCTAVE=octave
endif

if BUILDPYTHON
subdirPYTHON=python
endif

if BUILDSCILAB
subdirSCILAB=scilab
endif

SUBDIRS = . $(subdirMATLAB) $(subdirOCTAVE) $(subdirPYTHON) $(subdirSCILAB)

RPC_INC_DIR=@RPC_INC_DIR@

PSEUDO_MFUNCTIONS = 			\
	gf_spmat.cc 			\
	gf_spmat_set.cc 		\
	gf_spmat_get.cc 		\
	gf_linsolve.cc 			\
	gf_util.cc 			\
	gf_cont_struct.cc 		\
	gf_cont_struct_get.cc 		\
	gf_cvstruct_get.cc 		\
	gf_geotrans.cc 			\
	gf_geotrans_get.cc 		\
	gf_compute.cc 			\
	gf_mesh_fem.cc 			\
	gf_mesh_fem_set.cc 		\
	gf_mesh_fem_get.cc 		\
	gf_mesh_im.cc 			\
	gf_mesh_im_set.cc 		\
	gf_mesh_im_get.cc 		\
	gf_mesh_im_data.cc 		\
	gf_mesh_im_data_set.cc 		\
	gf_mesh_im_data_get.cc 		\
	gf_model.cc 			\
	gf_model_get.cc 		\
	gf_model_set.cc 		\
	gf_eltm.cc 			\
	gf_mesher_object.cc 		\
	gf_mesher_object_get.cc 	\
	gf_mesh.cc 			\
	gf_mesh_set.cc 			\
	gf_mesh_get.cc 			\
	gf_slice.cc 			\
	gf_slice_get.cc 		\
	gf_slice_set.cc 		\
	gf_levelset.cc 			\
	gf_levelset_get.cc 		\
	gf_levelset_set.cc 		\
	gf_mesh_levelset.cc 		\
	gf_mesh_levelset_get.cc 	\
	gf_mesh_levelset_set.cc 	\
	gf_precond.cc 			\
	gf_precond_get.cc 		\
	gf_asm.cc 			\
	gf_fem.cc 			\
	gf_fem_get.cc 			\
	gf_integ.cc 			\
	gf_integ_get.cc 		\
	gf_global_function.cc 		\
	gf_global_function_get.cc 	\
	gf_workspace.cc 		\
	gf_delete.cc

EXTRA_DIST = gfi_rpc_clnt.c gfi_rpc_xdr.c gfi_array.c

noinst_LTLIBRARIES = libgetfemint.la
#libgetfemint_a_FLAGS=-D__USE_XOPEN
libgetfemint_la_SOURCES =              	\
	getfem_interface.h  		\
	getfem_interface.cc            	\
	gfi_array.h 			\
	gfi_array.c 			\
	gfi_rpc.h 			\
	getfemint_std.h 		\
	getfemint.h 			\
	getfemint.cc 			\
	getfemint_misc.h 		\
	getfemint_misc.cc 		\
	$(PSEUDO_MFUNCTIONS) 		\
	getfemint_workspace.h 		\
	getfemint_workspace.cc 		\
	getfemint_precond.h 		\
	getfemint_levelset.h 		\
	getfemint_levelset.cc 		\
	getfemint_gsparse.h 		\
	getfemint_gsparse.cc

#libgetfemint_la_INCLUDES = @GETFEM_CPPFLAGS@ #fails with automake 1.6 on macos x tiger
AM_CPPFLAGS = -I$(top_srcdir)/src -I../../src

# -rdynamic for backtraces
#AM_LDFLAGS = -rdynamic

getfem_server_SOURCES = 	\
	gfi_rpc_server.c 	\
	gfi_rpc_svc.c 		\
	gfi_rpc_xdr.c
# for backtrace, need to link with g++ and -rdynamic
getfem_server_LINK=$(CXXLINK)
getfem_server_LDADD = libgetfemint.la ../../src/libgetfem.la @RPC_LIB@ -lm #-lstdc++ <- a corrgier pour icc
getfem_server_LDFLAGS = @GFSERVERFLAGS@ @SUPLDFLAGS@
getfem_server_INCLUDES = -I$(RPC_INC_DIR) -I$(top_srcdir)/src -I../../src
getfem_server_LIBS = libgetfemint.la

EXTRA_PROGRAMS = getfem_server
bin_PROGRAMS = @GETFEM_SERVER@

RPC_LIB = @RPC_LIB@
