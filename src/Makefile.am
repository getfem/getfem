include_HEADERS = \
	bgeot_approx_integration.h \
	bgeot_config.h \
	bgeot_convex.h \
	bgeot_convex_hull.h \
	bgeot_convex_intersect.h \
	bgeot_convex_ref.h \
	bgeot_convex_structure.h \
	bgeot_geometric_trans.h \
	bgeot_geotrans_inv.h \
	bgeot_matrix.h \
	bgeot_mesh.h \
	bgeot_mesh_structure.h \
	bgeot_poly.h \
	bgeot_poly_integration.h \
	bgeot_simplexify.h \
	bgeot_tensor.h \
	bgeot_vector.h \
	dal_algobase.h \
	dal_alloc.h \
	dal_basic.h \
	dal_bit_vector.h \
	dal_fonc_tables.h \
	dal_ref.h \
	dal_std.h \
	dal_tas.h \
	dal_tree_sorted.h \
	ftool.h \
	gensolv_function.h \
	gensolv_gauss_seidel.h \
	gensolv_line_search.h \
	gensolv_mtl.h \
	gensolv_quasi_newton.h \
	gensolv_secant.h \
	getfem_assembling.h \
	getfem_config.h \
	getfem_export.h \
	getfem_fem.h \
	getfem_mat_elem.h \
	getfem_mat_elem_type.h \
	getfem_mesh.h \
	getfem_mesh_fem.h \
	getfem_regular_meshes.h \
	linkmsg.h

SRC = \
	bgeot_approx_integration.C \
	bgeot_convex_ref.C \
	bgeot_convex_structure.C \
	bgeot_geometric_trans.C \
	bgeot_mesh_structure.C \
	bgeot_poly.C \
	bgeot_poly_integration.C \
	dal_bit_vector.C \
	ftool.C \
	getfem_export.C \
	getfem_fem.C \
	getfem_mat_elem.C \
	getfem_mat_elem_type.C \
	getfem_mesh.C \
	getfem_mesh_fem.C \
	getfem_regular_meshes.C

lib_LTLIBRARIES = libgetfem.la
libgetfm_la_SOURCES = $(SRC)
libgetfm_la_LDFLAGS = ${LIBTOOL_VERSION_INFO}

clean-local:
	rm -rf cxx_repository
