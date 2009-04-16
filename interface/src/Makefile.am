if BUILDMEX
subdirMATLAB=matlab
endif

if BUILDPYTHON
subdirPYTHON=python
endif

SUBDIRS = . $(subdirMATLAB) $(subdirPYTHON)

RPC_INC_DIR=@RPC_INC_DIR@

PSEUDO_MFUNCTIONS = \
	gf_spmat.cc \
	gf_spmat_set.cc \
	gf_spmat_get.cc \
	gf_linsolve.cc \
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
	gf_mdbrick.cc \
	gf_mdbrick_get.cc \
	gf_mdbrick_set.cc \
	gf_mdstate.cc \
	gf_mdstate_get.cc \
	gf_mdstate_set.cc \
	gf_model.cc \
	gf_model_get.cc \
	gf_model_set.cc \
	gf_eltm.cc \
	gf_mesh.cc \
	gf_mesh_set.cc \
	gf_mesh_get.cc \
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
	gf_asm.cc \
	gf_fem.cc \
	gf_fem_get.cc \
	gf_integ.cc \
	gf_integ_get.cc \
	gf_workspace.cc \
	gf_delete.cc

EXTRA_DIST = gfi_rpc_clnt.c gfi_rpc_xdr.c gfi_array.c

noinst_LTLIBRARIES = libgetfemint.la
#libgetfemint_a_FLAGS=-D__USE_XOPEN
libgetfemint_la_SOURCES = \
	getfem_interface.cc \
	getfem_interface.h \
	gfi_array.c \
	getfemint.cc \
	getfemint_misc.cc \
	$(PSEUDO_MFUNCTIONS) \
	getfemint_workspace.cc \
	getfemint_mesh.cc \
	getfemint_mesh_fem.cc \
	getfemint_mesh_im.cc \
	getfemint_matelem.cc \
	getfemint_matelemtype.cc \
	getfemint_pfem.cc \
	getfemint_integ.cc \
	getfemint_pgt.cc \
	getfemint_convex_structure.cc \
	gfi_array.h \
	getfemint_convex_structure.h \
	getfemint_mesh.h \
	getfemint_misc.h \
	getfemint_integ.h \
	getfemint_std.h \
	gfi_rpc.h \
	getfemint_matelem.h \
	getfemint_mesh_fem.h \
	getfemint_mesh_im.h \
	getfemint_mdbrick.h \
	getfemint_mdbrick.cc \
	getfemint_mdstate.h \
	getfemint_models.h \
	getfemint_object.h \
	getfemint_pgt.h \
	getfemint_workspace.h \
	getfemint.h \
	getfemint_matelemtype.h \
	getfemint_mesh_slice.h \
	getfemint_precond.h \
	getfemint_levelset.h \
	getfemint_levelset.cc \
	getfemint_mesh_levelset.h \
	getfemint_mesh_levelset.cc \
	getfemint_gsparse.h \
	getfemint_gsparse.cc \
	getfemint_gsparse_misc.h \
	getfemint_gsparse_misc.cc \
	getfemint_pfem.h \
	getfemint_poly.h

#libgetfemint_la_INCLUDES = @GETFEM_CPPFLAGS@ #fails with automake 1.6 on macos x tiger
INCLUDES = -I$(top_srcdir)/src -I../../src

# -rdynamic for backtraces
#AM_LDFLAGS = -rdynamic

getfem_server_SOURCES = gfi_rpc_server.c \
	gfi_rpc_svc.c \
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
