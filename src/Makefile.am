include_HEADERS = \
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
	bgeot_simplexify.h \
	bgeot_tensor.h \
	bgeot_vector.h \
	bgeot_smatrix.h \
	bgeot_abstract_linalg.h \
	bgeot_linalg_interface.h \
	bgeot_generic_solvers.h \
	dal_algobase.h \
	dal_alloc.h \
	dal_basic.h \
	dal_bit_vector.h \
	dal_fonc_tables.h \
	dal_ref.h \
	dal_std.h \
	dal_except.h \
	dal_tas.h \
	dal_tree_sorted.h \
	ftool.h \
	ftool_naming.h \
	getfem_integration.h \
	getfem_poly_composite.h \
	getfem_assembling.h \
	getfem_config.h \
	getfem_export.h \
	getfem_precomp.h \
	getfem_norm.h \
	getfem_derivatives.h \
	getfem_fem.h \
	getfem_mat_elem.h \
	getfem_mat_elem_type.h \
	getfem_mesh.h \
	getfem_mesh_fem.h \
	getfem_regular_meshes.h \
	linkmsg.h

SRC = \
        dal_except.C \
	dal_bit_vector.C \
	bgeot_convex_ref.C \
	bgeot_convex_structure.C \
	bgeot_geometric_trans.C \
	bgeot_geotrans_inv.C \
	bgeot_mesh_structure.C \
	bgeot_poly.C \
	ftool.C \
	getfem_mesh_fem.C \
	getfem_integration.C \
	getfem_poly_composite.C \
	getfem_integration_composite.C \
	getfem_fem.C \
	getfem_fem_virtual_link.C \
	getfem_fem_composite.C \
	getfem_mat_elem.C \
	getfem_mat_elem_type.C \
	getfem_precomp.C \
	getfem_mesh.C \
	getfem_regular_meshes.C \
	getfem_export.C

lib_LTLIBRARIES = libgetfem.la
libgetfem_la_SOURCES = $(SRC)
libgetfem_la_LDFLAGS = ${LIBTOOL_VERSION_INFO}
