include_HEADERS =                      \
	getfem_arch_config.h	       \
	bgeot_config.h                 \
	bgeot_convex.h                 \
	bgeot_convex_hull.h            \
	bgeot_convex_intersect.h       \
	bgeot_convex_ref.h             \
	bgeot_convex_structure.h       \
	bgeot_geometric_trans.h        \
	bgeot_geotrans_inv.h           \
	bgeot_mesh.h                   \
	bgeot_mesh_structure.h         \
	bgeot_poly.h                   \
	bgeot_simplexify.h             \
	bgeot_tensor.h                 \
	bgeot_sparse_tensors.h         \
	bgeot_vector.h                 \
	bgeot_comma_init.h	       \
	dal_algobase.h                 \
	dal_alloc.h                    \
	dal_basic.h                    \
	dal_bit_vector.h               \
	dal_fonc_tables.h              \
	dal_ref.h                      \
	dal_std.h                      \
	dal_except.h                   \
	dal_tas.h                      \
	dal_tree_sorted.h              \
	ftool.h                        \
	ftool_naming.h                 \
	getfem_integration.h           \
	getfem_poly_composite.h        \
	getfem_assembling.h            \
	getfem_assembling_tensors.h    \
	getfem_external_data_fem.h     \
	getfem_config.h                \
	getfem_export.h                \
	getfem_import.h	               \
	getfem_precomp.h               \
	getfem_norm.h                  \
	getfem_derivatives.h           \
	getfem_fem.h                   \
	getfem_Xfem.h                  \
	getfem_im_list.h               \
	getfem_mat_elem.h              \
	getfem_mat_elem_type.h         \
	getfem_mesh.h                  \
	getfem_mesh_fem.h              \
	getfem_mesh_slice.h	       \
	getfem_regular_meshes.h        \
	gmm.h                          \
	gmm_matrix.h                   \
	gmm_solvers.h                  \
	gmm_iter.h                     \
	gmm_inoutput.h                 \
	gmm_vector.h                   \
	gmm_transposed.h               \
	gmm_scaled.h                   \
	gmm_def.h                      \
	gmm_sub_index.h                \
	gmm_sub_vector.h               \
	gmm_sub_matrix.h               \
	gmm_interface.h                \
	gmm_interface_bgeot.h          \
	gmm_solver_cg.h                \
	gmm_solver_constrained_cg.h    \
	gmm_modified_gram_schmidt.h    \
	gmm_dense_lu.h                 \
	gmm_tri_solve.h                \
	gmm_solver_gmres.h             \
	gmm_solver_qmr.h               \
	gmm_solver_cheby.h             \
	gmm_solver_bicgstab.h          \
	gmm_solver_Schwarz_additive.h  \
	gmm_precond_cholesky.h         \
	gmm_precond_choleskyt.h        \
	gmm_precond_mr_approx_inverse.h\
	gmm_precond_diagonal.h         \
	gmm_precond_ilu.h              \
	gmm_precond_ilut.h             \
	gmm_given_rotation.h           \
	gmm_blas.h                     \
	gmm_condest.h                  \
	gmm_opt.h                      \
	linkmsg.h

SRC =                                  \
	dal_bit_vector.C               \
	bgeot_convex_ref.C             \
	bgeot_convex_structure.C       \
	bgeot_geometric_trans.C        \
	bgeot_geotrans_inv.C           \
	bgeot_mesh_structure.C         \
	bgeot_sparse_tensors.C         \
	bgeot_poly.C                   \
	ftool.C                        \
	getfem_mesh_fem.C              \
	getfem_integration.C           \
	getfem_poly_composite.C        \
	getfem_integration_composite.C \
	getfem_fem.C                   \
	getfem_Xfem.C                  \
	getfem_fem_virtual_link.C      \
	getfem_fem_composite.C         \
	getfem_mat_elem.C              \
	getfem_mat_elem_type.C         \
	getfem_precomp.C               \
	getfem_mesh.C                  \
	getfem_mesh_slice.C	       \
	getfem_regular_meshes.C        \
	getfem_export.C                \
	getfem_assembling_tensors.C

lib_LTLIBRARIES = libgetfem.la
libgetfem_la_SOURCES = $(SRC)
libgetfem_la_LDFLAGS = ${LIBTOOL_VERSION_INFO}

CLEANFILES = getfem_arch_config.h ii_files/* *.o.d