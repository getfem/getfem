cmake_minimum_required(VERSION 3.23)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)
include(CheckLibraryExists)

string(ASCII 27 ESC)
message("${ESC}[33m"
        "The cmake build system for GetFEM is an alternative to the "
        "autotools build. It is a work in progress and it does not "
        "support all options supported by the autotools system."
        "${ESC}[0m")
project(GetFEM)

# Set the project version
set(VERSION_MAJOR 5)
set(VERSION_MINOR 4)
set(VERSION_PATCH 2)
# applies to both GetFEM and GMM
set(GETFEM_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
set(GMM_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
# mimic autotools package variables
set(GETFEM_PACKAGE_NAME "getfem")
set(GETFEM_PACKAGE_STRING "getfem ${GETFEM_VERSION}")
set(GETFEM_PACKAGE_TARNAME "getfem")

# Configure options for enabling/disabling Python, Octave, and MATLAB support
option(ENABLE_PYTHON "Enable Python support" OFF)
option(ENABLE_OCTAVE "Enable Octave support" OFF)
option(ENABLE_MATLAB "Enable MATLAB support" OFF)
# Configure option for enabling/disabling OpenMP support
option(ENABLE_OPENMP "Enable OpenMP support" OFF)
# Configure option to enable Qhull support
option(ENABLE_QHULL "Enable Qhull support" ON)
# Configure options for enabling/disabling linear solvers (at least one is required)
option(ENABLE_SUPERLU "Enable SuperLU support" ON) # might be turned off by cmake if SuperLU is not found
option(ENABLE_MUMPS "Enable MUMPS support" ON)     # might be turned off by cmake if MUMPS is not found
# Configure option for enabling/disabling multithreaded BLAS (requires the dl library)
option(ENABLE_MULTITHREADED_BLAS "Enable multithreaded blas support" OFF)
option(GENERATE_GETFEM_IM_LIST_H "Run perl script that (re-)generates header file with integration methods" ON)

option(BUILD_SHARED_LIBS "Build libraries as SHARED, equivalent to BUILD_LIBRARY_TYPE=SHARED" ON)
set(BUILD_LIBRARY_TYPE "SHARED" CACHE STRING
    "Type of library to build, choose among SHARED, STATIC, STATIC_BUNDLE_DEPS")
set_property(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS
             "SHARED" "STATIC" "STATIC_BUNDLE_DEPS")
if(NOT BUILD_LIBRARY_TYPE STREQUAL "SHARED" AND BUILD_SHARED_LIBS)
  message(STATUS "BUILD_LIBRARY_TYPE=${BUILD_LIBRARY_TYPE} option overrides BUILD_SHARED_LIBS=ON")
  set(BUILD_SHARED_LIBS OFF)
elseif(NOT BUILD_SHARED_LIBS AND BUILD_LIBRARY_TYPE STREQUAL "SHARED")
  message(STATUS "BUILD_SHARED_LIBS=OFF option overrides BUILD_LIBRARY_TYPE=SHARED")
  set(BUILD_LIBRARY_TYPE "STATIC")
endif()

add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
set(WITH_OPTIMIZATION "-O2" CACHE STRING "Set the optimization level (default: -O2)")
set(GETFEM_PARA_LEVEL 0 CACHE STRING "Set the GETFEM_PARA_LEVEL option (default: 0)")

# CMake variable for installation path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")


# General tests and configurations
set(CMAKE_CXX_STANDARD 14 CACHE STRING "Set C++ standard version (default: 14))
set(CMAKE_CXX_STANDARD_REQUIRED ON)
check_include_file_cxx("cxxabi.h" GETFEM_HAVE_CXXABI_H)
check_cxx_source_compiles(
  "#include <fenv.h>\nint main() {feenableexcept(FE_DIVBYZERO); return 0;}"
  GETFEM_HAVE_FEENABLEEXCEPT)

# Tests for availability of linear solvers
if(ENABLE_SUPERLU)
  if(NOT BUILD_SHARED_LIBS)
    find_library(SUPERLU_STATIC_LIBS NAMES libsuperlu.a PATHS ${SUPERLU_LIB_DIR})
  endif()
  if(SUPERLU_STATIC_LIBS)
    set(SUPERLU_LIBS SUPERLU_STATIC_LIBS)
  else()
    find_library(SUPERLU_LIBS NAMES superlu PATHS ${SUPERLU_LIB_DIR})
  endif()
  if(SUPERLU_LIBS)
    find_path(SUPERLU_INCLUDE_PATH NAMES "superlu/slu_ddefs.h" HINTS ${SUPERLU_INC_DIR})
    set(CMAKE_REQUIRED_LIBRARIES ${SUPERLU_LIBS})
      check_library_exists(superlu dCreate_CompCol_Matrix "" HAVE_SUPERLU)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()

  if(SUPERLU_INCLUDE_PATH)
    message(STATUS "SUPERLU_INCLUDE_PATH = ${SUPERLU_INCLUDE_PATH}")
  else()
    message(WARNING "SuperLU headers NOT found : you can define SUPERLU_INC_DIR to help cmake find them")
  endif()

  if(HAVE_SUPERLU AND SUPERLU_INCLUDE_PATH)
    message(STATUS "SUPERLU_LIBS = ${SUPERLU_LIBS}")
    set(GMM_USES_SUPERLU 1)
    message(STATUS "Building with SuperLU")
  else()
    set(ENABLE_SUPERLU OFF)
    message(WARNING "Building without SuperLU")
  endif()
else()
  message("Building with SuperLU explicitly disabled")
endif()

if(ENABLE_MUMPS)
  set(MUMPS_LIBS "")
  if(GETFEM_PARA_LEVEL LESS_EQUAL 1) # try to find "?mumps_seq" named libs first
    if(NOT BUILD_SHARED_LIBS)
      message("NOT BUILD SHARED_LIBS")
      find_library(SMUMPS_STATIC_LIB NAMES libsmumps_seq.a libsmumps.a smumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(DMUMPS_STATIC_LIB NAMES libdmumps_seq.a libdmumps.a dmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(CMUMPS_STATIC_LIB NAMES libcmumps_seq.a libcmumps.a cmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(ZMUMPS_STATIC_LIB NAMES libzmumps_seq.a libzmumps.a zmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(MUMPS_COMMON_STATIC_LIB NAMES libmumps_common.a  mumps_common.lib PATHS ${MUMPS_LIB_DIR})
      find_library(PORD_STATIC_LIB NAMES libpord.a pord.lib PATHS ${MUMPS_LIB_DIR})
      find_library(MPISEQ_STATIC_LIB NAMES libmpiseq.a libmpiseq_seq.a mpiseq.lib PATHS ${MUMPS_LIB_DIR})
    endif()
    if(SMUMPS_STATIC_LIB AND DMUMPS_STATIC_LIB AND CMUMPS_STATIC_LIB AND ZMUMPS_STATIC_LIB
       AND MUMPS_COMMON_STATIC_LIB AND PORD_STATIC_LIB AND MPISEQ_STATIC_LIB)
      set(MUMPS_LIBS ${SMUMPS_STATIC_LIB} ${DMUMPS_STATIC_LIB} ${CMUMPS_STATIC_LIB} ${ZMUMPS_STATIC_LIB}
                     ${MUMPS_COMMON_STATIC_LIB} ${PORD_STATIC_LIB} ${MPISEQ_STATIC_LIB})
    else()
      message("NOT BUILD STATIC_LIBS")
      find_library(SMUMPS_LIB NAMES smumps_seq smumps PATHS ${MUMPS_LIB_DIR})
      find_library(DMUMPS_LIB NAMES dmumps_seq dmumps PATHS ${MUMPS_LIB_DIR})
      find_library(CMUMPS_LIB NAMES cmumps_seq cmumps PATHS ${MUMPS_LIB_DIR})
      find_library(ZMUMPS_LIB NAMES zmumps_seq zmumps PATHS ${MUMPS_LIB_DIR})
      find_library(MUMPS_COMMON_LIB NAMES mumps_common PATHS ${MUMPS_LIB_DIR})
      find_library(PORD_LIB NAMES pord PATHS ${MUMPS_LIB_DIR})
      find_library(MPISEQ_LIB NAMES mpiseq PATHS ${MUMPS_LIB_DIR})

      if (SMUMPS_LIB AND DMUMPS_LIB AND CMUMPS_LIB AND ZMUMPS_LIB
          AND MUMPS_COMMON_LIB AND PORD_LIB)
        set(MUMPS_LIBS ${SMUMPS_LIB} ${DMUMPS_LIB} ${CMUMPS_LIB} ${ZMUMPS_LIB}
                       ${MUMPS_COMMON_LIB} ${PORD_LIB})
      endif()

      if(MPISEQ_LIB)
          list(APPEND MUMPS_LIBS ${MPISEQ_LIB})
      else()
          message(WARNING "The library mpiseq was not found and will be ignored.")
      endif()
      message(STATUS "MUMPS_LIBS = ${MUMPS_LIBS}")

    endif()
  else() # GETFEM_PARA_LEVEL = 2
    if(NOT BUILD_SHARED_LIBS)
      find_library(SMUMPS_STATIC_LIB NAMES libsmumps.a smumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(DMUMPS_STATIC_LIB NAMES libdmumps.a dmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(CMUMPS_STATIC_LIB NAMES libcmumps.a cmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(ZMUMPS_STATIC_LIB NAMES libzmumps.a zmumps.lib PATHS ${MUMPS_LIB_DIR})
      find_library(MUMPS_COMMON_STATIC_LIB NAMES libmumps_common.a mumps_common.lib PATHS ${MUMPS_LIB_DIR})
      find_library(PORD_STATIC_LIB NAMES libpord.a pord.lib PATHS ${MUMPS_LIB_DIR})
    endif()
    if(SMUMPS_STATIC_LIB AND DMUMPS_STATIC_LIB AND CMUMPS_STATIC_LIB AND ZMUMPS_STATIC_LIB
       AND MUMPS_COMMON_STATIC_LIB AND PORD_STATIC_LIB)
      set(MUMPS_LIBS ${SMUMPS_STATIC_LIB} ${DMUMPS_STATIC_LIB} ${CMUMPS_STATIC_LIB} ${ZMUMPS_STATIC_LIB}
                     ${MUMPS_COMMON_STATIC_LIB} ${PORD_STATIC_LIB})
    else()
      find_library(SMUMPS_LIB NAMES smumps PATHS ${MUMPS_LIB_DIR})
      find_library(DMUMPS_LIB NAMES dmumps PATHS ${MUMPS_LIB_DIR})
      find_library(CMUMPS_LIB NAMES cmumps PATHS ${MUMPS_LIB_DIR})
      find_library(ZMUMPS_LIB NAMES zmumps PATHS ${MUMPS_LIB_DIR})
      find_library(MUMPS_COMMON_LIB NAMES mumps_common PATHS ${MUMPS_LIB_DIR})
      find_library(PORD_LIB NAMES pord PATHS ${MUMPS_LIB_DIR})
      if (SMUMPS_LIB AND DMUMPS_LIB AND CMUMPS_LIB AND ZMUMPS_LIB
          AND MUMPS_COMMON_LIB AND PORD_LIB)
        set(MUMPS_LIBS ${SMUMPS_LIB} ${DMUMPS_LIB} ${CMUMPS_LIB} ${ZMUMPS_LIB}
                       ${MUMPS_COMMON_LIB} ${PORD_LIB})
      endif()
    endif()
  endif()
#  find_library(MUMPS_COMMON_LIBRARY NAMES mumps_common HINTS ${MUMPS_LIB_DIR})
#  find_library(MUMPS_PORD_LIBRARY NAMES pord HINTS ${MUMPS_LIB_DIR})
#  find_library(MPI_SEQ NAMES mpiseq HINTS ${MUMPS_LIB_DIR})
#  set(MUMPS_LIBRARIES ${DMUMPS_LIB}
#                      ${MUMPS_COMMON_LIBRARY}
#                      ${MUMPS_PORD_LIBRARY}
#                      ${MPI_SEQ})
  #  find_package(MUMPS REQUIRED)
  find_path(MUMPS_INCLUDE_PATH NAMES dmumps_struc.h HINTS ${MUMPS_INC_DIR})
  set(CMAKE_REQUIRED_INCLUDES ${MUMPS_INCLUDE_PATH})
    check_include_file(smumps_c.h HAVE_SMUMPS_C)
    check_include_file(dmumps_c.h HAVE_DMUMPS_C)
    check_include_file(cmumps_c.h HAVE_CMUMPS_C)
    check_include_file(zmumps_c.h HAVE_ZMUMPS_C)
  set(CMAKE_REQUIRED_INCLUDES)

  if(MUMPS_LIBS)
    message(STATUS "MUMPS_LIBS = ${MUMPS_LIBS}")
  else()
    message(WARNING "MUMPS LIB NOT found : you can define MUMPS_LIB_DIR to help cmake find it")
  endif()

  if(MUMPS_INCLUDE_PATH)
    message(STATUS "MUMPS_INCLUDE_PATH = ${MUMPS_INCLUDE_PATH}")
  else()
    message(WARNING "MUMPS headers NOT found : you can define MUMPS_INC_DIR to help cmake find them")
  endif()

  if(MUMPS_LIBS AND MUMPS_INCLUDE_PATH)
    set(GMM_USES_MUMPS 1)
    message(STATUS "Building with MUMPS support")
  else()
    set(ENABLE_MUMPS OFF)
    message(WARNING "Building without MUMPS support")
  endif()
else()
  message("Building with MUMPS explicitly disabled")
endif()

if(NOT ENABLE_MUMPS AND NOT ENABLE_SUPERLU)
  message(FATAL_ERROR
          "Neither MUMPS nor SuperLU activated."
          "At least one direct linear solver is required."
          "Use MUMPS_LIB_DIR, MUMPS_INC_DIR, SUPERLU_LIB_DIR and SUPERLU_INC_DIR")
endif()

# define source files and main build target
set(SOURCES
    src/bgeot_convex_ref.cc
    src/bgeot_convex_ref_simplexified.cc
    src/bgeot_convex_structure.cc
    src/bgeot_ftool.cc
    src/bgeot_geometric_trans.cc
    src/bgeot_geotrans_inv.cc
    src/bgeot_kdtree.cc
    src/bgeot_mesh_structure.cc
    src/bgeot_node_tab.cc
    src/bgeot_poly.cc
    src/bgeot_poly_composite.cc
    src/bgeot_rtree.cc
    src/bgeot_small_vector.cc
    src/bgeot_sparse_tensors.cc
    src/bgeot_torus.cc
    src/dal_backtrace.cc
    src/dal_bit_vector.cc
    src/dal_singleton.cc
    src/dal_static_stored_objects.cc
    src/getfem_assembling_tensors.cc
    src/getfem_contact_and_friction_common.cc
    src/getfem_contact_and_friction_integral.cc
    src/getfem_contact_and_friction_large_sliding.cc
    src/getfem_contact_and_friction_nodal.cc
    src/getfem_context.cc
    src/getfem_continuation.cc
    src/getfem_enumeration_dof_para.cc
    src/getfem_error_estimate.cc
    src/getfem_export.cc
    src/getfem_fem.cc
    src/getfem_fem_composite.cc
    src/getfem_fem_global_function.cc
    src/getfem_fem_level_set.cc
    src/getfem_fourth_order.cc
    src/getfem_generic_assembly_compile_and_exec.cc
    src/getfem_generic_assembly_functions_and_operators.cc
    src/getfem_generic_assembly_interpolation.cc
    src/getfem_generic_assembly_semantic.cc
    src/getfem_generic_assembly_tree.cc
    src/getfem_generic_assembly_workspace.cc
    src/getfem_global_function.cc
    src/getfem_HHO.cc
    src/getfem_im_data.cc
    src/getfem_import.cc
    src/getfem_integration.cc
    src/getfem_integration_composite.cc
    src/getfem_interpolated_fem.cc
    src/getfem_interpolation.cc
    src/getfem_interpolation_on_deformed_domains.cc
    src/getfem_level_set.cc
    src/getfem_level_set_contact.cc
    src/getfem_linearized_plates.cc
    src/getfem_locale.cc
    src/getfem_mat_elem.cc
    src/getfem_mat_elem_type.cc
    src/getfem_mesh.cc
    src/getfem_mesher.cc
    src/getfem_mesh_fem.cc
    src/getfem_mesh_fem_global_function.cc
    src/getfem_mesh_fem_level_set.cc
    src/getfem_mesh_fem_product.cc
    src/getfem_mesh_fem_sum.cc
    src/getfem_mesh_im.cc
    src/getfem_mesh_im_level_set.cc
    src/getfem_mesh_level_set.cc
    src/getfem_mesh_region.cc
    src/getfem_mesh_slice.cc
    src/getfem_mesh_slicers.cc
    src/getfem_models.cc
    src/getfem_model_solvers.cc
    src/getfem_nonlinear_elasticity.cc
    src/getfem_omp.cc
    src/getfem_partial_mesh_fem.cc
    src/getfem_plasticity.cc
    src/getfem_projected_fem.cc
    src/getfem_regular_meshes.cc
    src/getfem_torus.cc)

set(HEADERS
    src/gmm/gmm_algobase.h
    src/gmm/gmm_blas.h
    src/gmm/gmm_blas_interface.h
    src/gmm/gmm_condition_number.h
    src/gmm/gmm_conjugated.h
    src/gmm/gmm_def.h
    src/gmm/gmm_dense_Householder.h
    src/gmm/gmm_dense_lu.h
    src/gmm/gmm_dense_matrix_functions.h
    src/gmm/gmm_dense_qr.h
    src/gmm/gmm_dense_sylvester.h
    src/gmm/gmm_domain_decomp.h
    src/gmm/gmm_except.h
    src/gmm/gmm_feedback_management.h
    src/gmm/gmm.h
    src/gmm/gmm_inoutput.h
    src/gmm/gmm_interface_bgeot.h
    src/gmm/gmm_interface.h
    src/gmm/gmm_iter.h
    src/gmm/gmm_iter_solvers.h
    src/gmm/gmm_kernel.h
    src/gmm/gmm_lapack_interface.h
    src/gmm/gmm_least_squares_cg.h
    src/gmm/gmm_matrix.h
    src/gmm/gmm_modified_gram_schmidt.h
    src/gmm/gmm_MUMPS_interface.h
    src/gmm/gmm_opt.h
    src/gmm/gmm_precond_diagonal.h
    src/gmm/gmm_precond.h
    src/gmm/gmm_precond_ildlt.h
    src/gmm/gmm_precond_ildltt.h
    src/gmm/gmm_precond_ilu.h
    src/gmm/gmm_precond_ilut.h
    src/gmm/gmm_precond_ilutp.h
    src/gmm/gmm_precond_mr_approx_inverse.h
    src/gmm/gmm_range_basis.h
    src/gmm/gmm_real_part.h
    src/gmm/gmm_ref.h
    src/gmm/gmm_scaled.h
    src/gmm/gmm_solver_bfgs.h
    src/gmm/gmm_solver_bicgstab.h
    src/gmm/gmm_solver_cg.h
    src/gmm/gmm_solver_constrained_cg.h
    src/gmm/gmm_solver_gmres.h
    src/gmm/gmm_solver_idgmres.h
    src/gmm/gmm_solver_qmr.h
    src/gmm/gmm_solver_Schwarz_additive.h
    src/gmm/gmm_std.h
    src/gmm/gmm_sub_index.h
    src/gmm/gmm_sub_matrix.h
    src/gmm/gmm_sub_vector.h
    src/gmm/gmm_superlu_interface.h
    src/gmm/gmm_transposed.h
    src/gmm/gmm_tri_solve.h
    src/gmm/gmm_vector.h
    src/gmm/gmm_vector_to_matrix.h
    src/getfem/bgeot_comma_init.h
    src/getfem/bgeot_config.h
    src/getfem/bgeot_convex.h
    src/getfem/bgeot_convex_ref.h
    src/getfem/bgeot_convex_structure.h
    src/getfem/bgeot_ftool.h
    src/getfem/bgeot_geometric_trans.h
    src/getfem/bgeot_geotrans_inv.h
    src/getfem/bgeot_kdtree.h
    src/getfem/bgeot_mesh.h
    src/getfem/bgeot_mesh_structure.h
    src/getfem/bgeot_node_tab.h
    src/getfem/bgeot_permutations.h
    src/getfem/bgeot_poly_composite.h
    src/getfem/bgeot_poly.h
    src/getfem/bgeot_rtree.h
    src/getfem/bgeot_small_vector.h
    src/getfem/bgeot_sparse_tensors.h
    src/getfem/bgeot_tensor.h
    src/getfem/bgeot_torus.h
    src/getfem/dal_backtrace.h
    src/getfem/dal_basic.h
    src/getfem/dal_bit_vector.h
    src/getfem/dal_config.h
    src/getfem/dal_naming_system.h
    src/getfem/dal_singleton.h
    src/getfem/dal_static_stored_objects.h
    src/getfem/dal_tas.h
    src/getfem/dal_tree_sorted.h
    src/getfem/getfem_accumulated_distro.h
    src/getfem/getfem_assembling.h
    src/getfem/getfem_assembling_tensors.h
    src/getfem/getfem_config.h
    src/getfem/getfem_contact_and_friction_common.h
    src/getfem/getfem_contact_and_friction_integral.h
    src/getfem/getfem_contact_and_friction_large_sliding.h
    src/getfem/getfem_contact_and_friction_nodal.h
    src/getfem/getfem_context.h
    src/getfem/getfem_continuation.h
    src/getfem/getfem_convect.h
    src/getfem/getfem_copyable_ptr.h
    src/getfem/getfem_crack_sif.h
    src/getfem/getfem_deformable_mesh.h
    src/getfem/getfem_derivatives.h
    src/getfem/getfem_error_estimate.h
    src/getfem/getfem_export.h
    src/getfem/getfem_fem_global_function.h
    src/getfem/getfem_fem.h
    src/getfem/getfem_fem_level_set.h
    src/getfem/getfem_fourth_order.h
    src/getfem/getfem_generic_assembly_compile_and_exec.h
    src/getfem/getfem_generic_assembly_functions_and_operators.h
    src/getfem/getfem_generic_assembly.h
    src/getfem/getfem_generic_assembly_semantic.h
    src/getfem/getfem_generic_assembly_tree.h
    src/getfem/getfem_global_function.h
    src/getfem/getfem_HHO.h
    src/getfem/getfem_im_data.h
    src/getfem/getfem_im_list.h
    src/getfem/getfem_import.h
    src/getfem/getfem_integration.h
    src/getfem/getfem_interpolated_fem.h
    src/getfem/getfem_interpolation.h
    src/getfem/getfem_level_set_contact.h
    src/getfem/getfem_level_set.h
    src/getfem/getfem_linearized_plates.h
    src/getfem/getfem_locale.h
    src/getfem/getfem_mat_elem.h
    src/getfem/getfem_mat_elem_type.h
    src/getfem/getfem_mesher.h
    src/getfem/getfem_mesh_fem_global_function.h
    src/getfem/getfem_mesh_fem.h
    src/getfem/getfem_mesh_fem_level_set.h
    src/getfem/getfem_mesh_fem_product.h
    src/getfem/getfem_mesh_fem_sum.h
    src/getfem/getfem_mesh.h
    src/getfem/getfem_mesh_im.h
    src/getfem/getfem_mesh_im_level_set.h
    src/getfem/getfem_mesh_level_set.h
    src/getfem/getfem_mesh_region.h
    src/getfem/getfem_mesh_slice.h
    src/getfem/getfem_mesh_slicers.h
    src/getfem/getfem_models.h
    src/getfem/getfem_model_solvers.h
    src/getfem/getfem_Navier_Stokes.h
    src/getfem/getfem_nonlinear_elasticity.h
    src/getfem/getfem_omp.h
    src/getfem/getfem_partial_mesh_fem.h
    src/getfem/getfem_plasticity.h
    src/getfem/getfem_projected_fem.h
    src/getfem/getfem_regular_meshes.h
    src/getfem/getfem_torus.h)

# Create the library target
if(BUILD_SHARED_LIBS)
  add_library(libgetfem SHARED ${SOURCES})
else()
  add_library(libgetfem STATIC ${SOURCES})
endif()
set_target_properties(libgetfem PROPERTIES OUTPUT_NAME "getfem")
set_target_properties(libgetfem PROPERTIES VERSION ${GETFEM_VERSION})

if(MSVC)
  target_compile_options(libgetfem PRIVATE /FIiso646.h)
endif()

target_sources(libgetfem
               PUBLIC FILE_SET header_set1 TYPE HEADERS
                                           BASE_DIRS ${CMAKE_SOURCE_DIR}/src
                                           FILES ${HEADERS}
               PUBLIC FILE_SET header_set2 TYPE HEADERS
                                           BASE_DIRS ${CMAKE_BINARY_DIR}
                                           FILES "${CMAKE_BINARY_DIR}/gmm/gmm_arch_config.h"
                                                 "${CMAKE_BINARY_DIR}/getfem/getfem_arch_config.h")

target_include_directories(libgetfem PRIVATE ${CMAKE_BINARY_DIR})
target_include_directories(libgetfem PRIVATE ${CMAKE_SOURCE_DIR}/src)
if(ENABLE_SUPERLU)
  target_include_directories(libgetfem PRIVATE ${SUPERLU_INCLUDE_PATH})
  target_link_libraries(libgetfem PRIVATE ${SUPERLU_LIBS})
endif()
if(ENABLE_MUMPS)
  target_include_directories(libgetfem PRIVATE ${MUMPS_INCLUDE_PATH})
  target_link_libraries(libgetfem PUBLIC "${MUMPS_LIBS}")
endif()

if(NOT GETFEM_PARA_LEVEL MATCHES "^(0|1|2)$")
  message(FATAL_ERROR "GETFEM_PARA_LEVEL must be 0, 1, or 2")
else()
  if(GETFEM_PARA_LEVEL GREATER_EQUAL 1)
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
      target_link_libraries(libgetfem PRIVATE MPI::MPI_CXX)
    else()
      message(FATAL_ERROR "MPI-enabled C++ compiler not found.")
    endif()
    #if(GETFEM_PARA_LEVEL EQUAL 2)
    #TODO check parallel MUMPS availability
    #endif()
  endif()
endif()


if(BUILD_SHARED_LIBS)
  set_property(TARGET libgetfem PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(libgetfem PRIVATE ${WITH_OPTIMIZATION})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  target_compile_options(libgetfem PRIVATE ${WITH_OPTIMIZATION} -Xc -ansi)
endif()


# Build python interface
if(ENABLE_PYTHON)
  if(NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "Building the python interface requires shared libraries")
  else()
    set(INTERFACE_SOURCES
        interface/src/getfem_interface.h
        interface/src/getfem_interface.cc
        interface/src/gfi_array.h
        interface/src/gfi_array.c
        interface/src/gfi_rpc.h
        interface/src/getfemint_std.h
        interface/src/getfemint.h
        interface/src/getfemint.cc
        interface/src/getfemint_misc.h
        interface/src/getfemint_misc.cc
        interface/src/gf_spmat.cc
        interface/src/gf_spmat_set.cc
        interface/src/gf_spmat_get.cc
        interface/src/gf_linsolve.cc
        interface/src/gf_util.cc
        interface/src/gf_cont_struct.cc
        interface/src/gf_cont_struct_get.cc
        interface/src/gf_cvstruct_get.cc
        interface/src/gf_geotrans.cc
        interface/src/gf_geotrans_get.cc
        interface/src/gf_compute.cc
        interface/src/gf_mesh_fem.cc
        interface/src/gf_mesh_fem_set.cc
        interface/src/gf_mesh_fem_get.cc
        interface/src/gf_mesh_im.cc
        interface/src/gf_mesh_im_set.cc
        interface/src/gf_mesh_im_get.cc
        interface/src/gf_mesh_im_data.cc
        interface/src/gf_mesh_im_data_set.cc
        interface/src/gf_mesh_im_data_get.cc
        interface/src/gf_model.cc
        interface/src/gf_model_get.cc
        interface/src/gf_model_set.cc
        interface/src/gf_mumps_context.cc
        interface/src/gf_mumps_context_set.cc
        interface/src/gf_mumps_context_get.cc
        interface/src/gf_eltm.cc
        interface/src/gf_mesher_object.cc
        interface/src/gf_mesher_object_get.cc
        interface/src/gf_mesh.cc
        interface/src/gf_mesh_set.cc
        interface/src/gf_mesh_get.cc
        interface/src/gf_slice.cc
        interface/src/gf_slice_get.cc
        interface/src/gf_slice_set.cc
        interface/src/gf_levelset.cc
        interface/src/gf_levelset_get.cc
        interface/src/gf_levelset_set.cc
        interface/src/gf_mesh_levelset.cc
        interface/src/gf_mesh_levelset_get.cc
        interface/src/gf_mesh_levelset_set.cc
        interface/src/gf_precond.cc
        interface/src/gf_precond_get.cc
        interface/src/gf_asm.cc
        interface/src/gf_fem.cc
        interface/src/gf_fem_get.cc
        interface/src/gf_integ.cc
        interface/src/gf_integ_get.cc
        interface/src/gf_global_function.cc
        interface/src/gf_global_function_get.cc
        interface/src/gf_workspace.cc
        interface/src/gf_delete.cc
        interface/src/getfemint_workspace.h
        interface/src/getfemint_workspace.cc
        interface/src/getfemint_gprecond.h
        interface/src/getfemint_gmumps.h
        interface/src/getfemint_levelset.h
        interface/src/getfemint_levelset.cc
        interface/src/getfemint_gsparse.h
        interface/src/getfemint_gsparse.cc)

    add_library(libgetfemint OBJECT ${INTERFACE_SOURCES})
    set_property(TARGET libgetfemint PROPERTY POSITION_INDEPENDENT_CODE ON)
    target_include_directories(libgetfemint PRIVATE ${CMAKE_BINARY_DIR})
    target_include_directories(libgetfemint PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_include_directories(libgetfemint PRIVATE ${CMAKE_SOURCE_DIR}/interface/src)
    if(ENABLE_MUMPS)
      target_include_directories(libgetfemint PRIVATE ${MUMPS_INCLUDE_PATH})
    endif()

    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

    # create a folder to build the python extension from
    set(PYTHON_MODULE_PATH ${CMAKE_BINARY_DIR}/python_module)
    file(MAKE_DIRECTORY ${PYTHON_MODULE_PATH})
    message("PYTHON_MODULE_PATH = ${PYTHON_MODULE_PATH}")
    add_custom_target(
      getfem.py ALL
      DEPENDS libgetfemint
      WORKING_DIRECTORY ${PYTHON_MODULE_PATH}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/bin/extract_doc ${CMAKE_SOURCE_DIR}/interface/src python-com > getfem.py
      COMMENT "Generating getfem.py")

    Python3_add_library(_getfem MODULE WITH_SOABI
                        ${PROJECT_SOURCE_DIR}/interface/src/python/getfem_python.c
                        $<TARGET_OBJECTS:libgetfemint>)
    target_include_directories(_getfem PRIVATE ${CMAKE_BINARY_DIR}/getfem)
    target_include_directories(_getfem PRIVATE ${CMAKE_SOURCE_DIR}/interface/src)
    target_include_directories(_getfem PRIVATE ${Python3_NumPy_INCLUDE_DIRS})
    target_link_libraries(_getfem PRIVATE libgetfem)

    # find the installation directory for python libraries (ported from autotools)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "\
import sysconfig; scheme=sysconfig.get_default_scheme().replace('posix_local','posix_prefix'); \
print(sysconfig.get_path('purelib', scheme, vars={'base': '${CMAKE_INSTALL_PREFIX}'}))"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("Python3_SITEARCH = ${Python3_SITEARCH} (not used)")
    message("Python3_SITELIB = ${Python3_SITELIB} (not used)")
    message("PYTHON_SITE_PACKAGES = ${PYTHON_SITE_PACKAGES} (used)")
    message("Python3_NumPy_INCLUDE_DIRS = ${Python3_NumPy_INCLUDE_DIRS} (used)")
  endif()
endif()

if(ENABLE_OCTAVE)
  message(FATAL_ERROR "Building the Octave interface not supported in the cmake configuration")
endif()

if(ENABLE_MATLAB)
  message(FATAL_ERROR "Building the Matlab interface not supported in the cmake configuration")
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    set(GETFEM_HAS_OPENMP 1)
    target_link_libraries(libgetfem PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "Building with OpenMP support")
  else()
    message(FATAL_ERROR "Could not find OpenMP requested by the user")
  endif()
endif()

# in the cmake build, BLAS and Lapack are hard requirements
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(libgetfem PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set(GMM_USES_BLAS 1)
set(GMM_USES_LAPACK 1)

if(ENABLE_QHULL)
#  find_library(QHULL_R_LIB qhull_r)
#  message(STATUS "QHULL_R_LIB = ${QHULL_R_LIB}")
#  target_link_libraries(libgetfem PRIVATE ${QHULL_R_LIB})
  find_package(Qhull REQUIRED COMPONENTS qhull_r)
  if(Qhull_FOUND)
    set(GETFEM_HAVE_LIBQHULL_R_QHULL_RA_H 1)
    target_link_libraries(libgetfem PRIVATE Qhull::qhull_r)
    message(STATUS "Building with Qhull support")
  else()
    set(ENABLE_QHULL OFF)
    message(WARNING "Building without Qhull support")
  endif()
else()
  message("Building with Qhull explicitly disabled")
endif()

if(NOT ENABLE_MULTITHREADED_BLAS)
  if(WIN32)
    # Sur Windows, pas besoin de dl et dlsym
    set(GETFEM_FORCE_SINGLE_THREAD_BLAS 0)
  else()
    find_library(DL_LIB NAMES dl)
    set(CMAKE_REQUIRED_LIBRARIES ${DL_LIB})
    check_library_exists(dl dlsym "" HAVE_DLSYM)
    set(CMAKE_REQUIRED_LIBRARIES ${DL_LIB})
    if(HAVE_DLSYM)
      target_link_libraries(libgetfem PRIVATE ${DL_LIB})
    else()
      message(FATAL_ERROR
              "Could not find dlsym function required for forcing singlethreaded BLAS")
    endif()
    set(GETFEM_FORCE_SINGLE_THREAD_BLAS 1)
  endif()
endif()


# Print build options
message(STATUS "Build options:")
message(STATUS "  WITH_OPTIMIZATION: ${WITH_OPTIMIZATION}")
message(STATUS "  GETFEM_PARA_LEVEL: ${GETFEM_PARA_LEVEL}")
message(STATUS "  ENABLE_PYTHON: ${ENABLE_PYTHON}")
message(STATUS "  ENABLE_OCTAVE: ${ENABLE_OCTAVE}")
message(STATUS "  ENABLE_MATLAB: ${ENABLE_MATLAB}")
message(STATUS "  ENABLE_OPENMP: ${ENABLE_OPENMP}")
message(STATUS "  ENABLE_SUPERLU: ${ENABLE_SUPERLU}")
message(STATUS "  ENABLE_MUMPS: ${ENABLE_MUMPS}")
message(STATUS "  ENABLE_QHULL: ${ENABLE_QHULL}")
message(STATUS "  ENABLE_MULTITHREADED_BLAS: ${ENABLE_MULTITHREADED_BLAS}")
message(STATUS "GetFEM version ${GETFEM_VERSION}")

# Generate configuration header files for gmm and getfem
configure_file(cmake/gmm_arch_config.h.in ${CMAKE_BINARY_DIR}/gmm/gmm_arch_config.h)
configure_file(cmake/getfem_arch_config.h.in ${CMAKE_BINARY_DIR}/getfem/getfem_arch_config.h)

if(GENERATE_GETFEM_IM_LIST_H)
find_program(Perl_EXECUTABLE perl REQUIRED)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/cubature/make_getfem_im_list" MGF_SCRIPT)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/cubature" MGF_ARG)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/src/getfem" MGF_WORKDIR)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/src/getfem/getfem_im_list.h" MGF_OUTPUT)

message(STATUS "Perl_EXECUTABLE = ${Perl_EXECUTABLE}")
message(STATUS "MGF_SCRIPT = ${MGF_SCRIPT}")
message(STATUS "MGF_ARG = ${MGF_ARG}")
message(STATUS "MGF_WORKDIR = ${MGF_WORKDIR}")
message(STATUS "MGF_OUTPUT = ${MGF_OUTPUT}")

add_custom_command(
  OUTPUT MGF_OUTPUT
  WORKING_DIRECTORY MGF_WORKDIR 
  COMMAND ${Perl_EXECUTABLE} ${MGF_SCRIPT} ${MGF_ARG}
  COMMENT "Generating getfem_im_list.h")
endif()

install(TARGETS libgetfem
        FILE_SET header_set1
        FILE_SET header_set2)
if(ENABLE_PYTHON)
  install(TARGETS _getfem
          DESTINATION ${PYTHON_SITE_PACKAGES}/getfem)
  install(FILES "${PYTHON_MODULE_PATH}/getfem.py" "${CMAKE_SOURCE_DIR}/interface/src/python/__init__.py"
          DESTINATION ${PYTHON_SITE_PACKAGES}/getfem)
endif()

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
