/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_laplack_interface.h : specialization of operations for   */
/*                             dense matrices calling lapack.              */
/*     									   */
/* Date : October 7, 2003.                                                 */
/* Authors : Caroline Lecalvez, Caroline.Lecalvez@gmm.insa-tlse.fr         */
/*           Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GMM++                                            */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#if defined(GMM_USES_LAPACK) || defined(GMM_USES_ATLAS)

#ifndef GMM_LAPACK_INTERFACE_H
#define GMM_LAPACK_INTERFACE_H

namespace gmm {

  /* ********************************************************************* */
  /* Operations interfaced for T = float, double, std::complex<float>      */
  /*    or std::complex<double> :                                          */
  /*                                                                       */
  /* mult(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)               */
  /* mult(transposed(dense_matrix<T>), dense_matrix<T>, dense_matrix<T>)   */
  /* mult(dense_matrix<T>, transposed(dense_matrix<T>), dense_matrix<T>)   */
  /* mult(transposed(dense_matrix<T>), transposed(dense_matrix<T>),        */
  /*      dense_matrix<T>)                                                 */
  /* mult(conjugated(dense_matrix<T>), dense_matrix<T>, dense_matrix<T>)   */
  /* mult(dense_matrix<T>, conjugated(dense_matrix<T>), dense_matrix<T>)   */
  /* mult(conjugated(dense_matrix<T>), conjugated(dense_matrix<T>),        */
  /*      dense_matrix<T>)                                                 */
  /*                                                                       */
  /* mult(dense_matrix<T>, std::vector<T>, std::vector<T>)                 */
  /* mult(transposed(dense_matrix<T>), std::vector<T>, std::vector<T>)     */
  /* mult(conjugated(dense_matrix<T>), std::vector<T>, std::vector<T>)     */
  /* mult(dense_matrix<T>, scaled(std::vector<T>), std::vector<T>)         */
  /* mult(transposed(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      std::vector<T>)                                                  */
  /* mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      std::vector<T>)                                                  */
  /*                                                                       */
  /* mult(dense_matrix<T>, std::vector<T>, std::vector<T>, std::vector<T>) */
  /* mult(transposed(dense_matrix<T>), std::vector<T>, std::vector<T>,     */
  /*      std::vector<T>)                                                  */
  /* mult(conjugated(dense_matrix<T>), std::vector<T>, std::vector<T>,     */
  /*      std::vector<T>)                                                  */
  /* mult(dense_matrix<T>, scaled(std::vector<T>), std::vector<T>,         */
  /*      std::vector<T>)                                                  */
  /* mult(transposed(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      std::vector<T>, std::vector<T>)                                  */
  /* mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      std::vector<T>, std::vector<T>)                                  */
  /* mult(dense_matrix<T>, std::vector<T>, scaled(std::vector<T>),         */
  /*      std::vector<T>)                                                  */
  /* mult(transposed(dense_matrix<T>), std::vector<T>,                     */
  /*      scaled(std::vector<T>), std::vector<T>)                          */
  /* mult(conjugated(dense_matrix<T>), std::vector<T>,                     */
  /*      scaled(std::vector<T>), std::vector<T>)                          */
  /* mult(dense_matrix<T>, scaled(std::vector<T>), scaled(std::vector<T>), */
  /*   std::vector<T>)                                                     */
  /* mult(transposed(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      scaled(std::vector<T>), std::vector<T>)                          */
  /* mult(conjugated(dense_matrix<T>), scaled(std::vector<T>),             */
  /*      scaled(std::vector<T>), std::vector<T>)                          */
  /*                                                                       */
  /* lower_tri_solve(dense_matrix<T>, std::vector<T>, k, b)                */
  /* upper_tri_solve(dense_matrix<T>, std::vector<T>, k, b)                */
  /* lower_tri_solve(transposed(dense_matrix<T>), std::vector<T>, k, b)    */
  /* upper_tri_solve(transposed(dense_matrix<T>), std::vector<T>, k, b)    */
  /* lower_tri_solve(conjugated(dense_matrix<T>), std::vector<T>, k, b)    */
  /* upper_tri_solve(conjugated(dense_matrix<T>), std::vector<T>, k, b)    */
  /*                                                                       */
  /* lu_factor(dense_matrix<T>, std::vector<int>)                          */
  /* lu_solve(dense_matrix<T>, std::vector<T>, std::vector<T>)             */
  /* lu_solve(dense_matrix<T>, std::vector<int>, std::vector<T>,           */
  /*          std::vector<T>)                                              */
  /* lu_solve_transposed(dense_matrix<T>, std::vector<int>, std::vector<T>,*/
  /*          std::vector<T>)                                              */
  /* lu_inverse(dense_matrix<T>)                                           */
  /* lu_inverse(dense_matrix<T>, std::vector<int>, dense_matrix<T>)        */
  /*                                                                       */
  /* qr_factor(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)          */
  /*                                                                       */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>)                */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>,                */
  /*                       dense_matrix<T>)                                */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >) */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >, */
  /*                       dense_matrix<T>)                                */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /* Basic defines.                                                        */
  /* ********************************************************************* */

# define BLAS_S float
# define BLAS_D double
# define BLAS_C std::complex<float>
# define BLAS_Z std::complex<double>

  /* ********************************************************************* */
  /* BLAS functions used.                                                  */
  /* ********************************************************************* */
  extern "C" {
    void sgemm_(...); void dgemm_(...); void cgemm_(...); void zgemm_(...);
    void sgemv_(...); void dgemv_(...); void cgemv_(...); void zgemv_(...);
    void strsv_(...); void dtrsv_(...); void ctrsv_(...); void ztrsv_(...);
  }

  /* ********************************************************************* */
  /* LAPACK functions used.                                                */
  /* ********************************************************************* */

  extern "C" {
    void sgetrf_(...); void dgetrf_(...); void cgetrf_(...); void zgetrf_(...);
    void sgetrs_(...); void dgetrs_(...); void cgetrs_(...); void zgetrs_(...);
    void sgetri_(...); void dgetri_(...); void cgetri_(...); void zgetri_(...);
    void sgeqrf_(...); void dgeqrf_(...); void cgeqrf_(...); void zgeqrf_(...);
    void sorgqr_(...); void dorgqr_(...); void cungqr_(...); void zungqr_(...);
    void  sgees_(...); void  dgees_(...); void  cgees_(...);  void zgees_(...);
    void  sgeev_(...); void  dgeev_(...); void  cgeev_(...);  void zgeev_(...);
  }


  /* ********************************************************************* */
  /* mult(A, x, y, z).                                                     */
  /* ********************************************************************* */
  
# define gemv_interface(param1, trans1, param2, trans2, param3, trans3,    \
                        blas_name, base_type, orien)                       \
  inline void mult_spec(param1(base_type), param2(base_type),              \
              param3(base_type), std::vector<base_type > &z, orien) {      \
    trans1(base_type); trans2(base_type); trans3(base_type);               \
    int m(mat_nrows(A)), lda(m), n(mat_ncols(A)), inc(1); gmm::copy(y, z); \
    blas_name(&t, &m, &n, &alpha, &A(0,0), &lda, &x[0], &inc, &beta,       \
              &z[0], &inc);                                                \
  }

  // First parameter
# define gem_p1_n(base_type)  const dense_matrix<base_type > &A
# define gem_trans1_n(base_type) const char t = 'N'
# define gem_p1_t(base_type)                                               \
         const transposed_col_ref<dense_matrix<base_type > *> &_A
# define gem_trans1_t(base_type) dense_matrix<base_type > &A =             \
         *(dense_matrix<base_type > *)(linalg_origin(_A)); const char t = 'T'
# define gem_p1_tc(base_type)                                              \
         const transposed_col_ref<const dense_matrix<base_type > *> &_A
# define gem_p1_c(base_type)                                               \
         const conjugated_col_matrix_const_ref<dense_matrix<base_type > > &_A
# define gem_trans1_c(base_type) dense_matrix<base_type > &A =             \
         *(dense_matrix<base_type > *)(linalg_origin(_A)); const char t = 'C'


  // second parameter 
# define gemv_p2_n(base_type)  const std::vector<base_type > &x
# define gemv_trans2_n(base_type) base_type alpha(1)
# define gemv_p2_s(base_type)                                              \
         const scaled_vector_const_ref<std::vector<base_type > > &_x
# define gemv_trans2_s(base_type) std::vector<base_type > &x =             \
         *(std::vector<base_type > *)(linalg_origin(_x));                  \
         base_type alpha(_x.r)

  // third parameter
# define gemv_p3_n(base_type)  const std::vector<base_type > &y
# define gemv_trans3_n(base_type) base_type beta(1)
# define gemv_p3_s(base_type)                                              \
         const scaled_vector_const_ref<std::vector<base_type > > &_y
# define gemv_trans3_s(base_type) std::vector<base_type > &y =             \
         *(std::vector<base_type > *)(linalg_origin(_y));                  \
         base_type beta(_y.r)  


  // Z <- AX + Y.
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, col_major);

  // Z <- transposed(A)X + Y.
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);
  
  // Z <- transposed(const A)X + Y.
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);
  
  // Z <- conjugated(A)X + Y.
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);

  // Z <- A scaled(X) + Y.
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, col_major);

  // Z <- transposed(A) scaled(X) + Y.
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);
  
  // Z <- transposed(const A) scaled(X) + Y.
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);
  
  // Z <- conjugated(A) scaled(X) + Y.
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_n,
		 gemv_trans3_n, zgemv_, BLAS_Z, row_major);

  // Z <- AX + scaled(Y).
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, col_major);

  // Z <- transposed(A)X + scaled(Y).
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);
  
  // Z <- transposed(const A)X + scaled(Y).
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);
  
  // Z <- conjugated(A)X + scaled(Y).
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);

  // Z <- A scaled(X) + scaled(Y).
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, col_major);
  gemv_interface(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, col_major);

  // Z <- transposed(A) scaled(X) + scaled(Y).
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);
  
  // Z <- transposed(const A) scaled(X) + scaled(Y).
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);
  
  // Z <- conjugated(A) scaled(X) + scaled(Y).
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, sgemv_, BLAS_S, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, dgemv_, BLAS_D, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, cgemv_, BLAS_C, row_major);
  gemv_interface(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, gemv_p3_s,
		 gemv_trans3_s, zgemv_, BLAS_Z, row_major);

  /* ********************************************************************* */
  /* mult(A, x, y).                                                        */
  /* ********************************************************************* */
  
# define gemv_interface2(param1, trans1, param2, trans2, blas_name,        \
                         base_type, orien)                                 \
  inline void mult_spec(param1(base_type), param2(base_type),              \
              std::vector<base_type > &z, orien) {                         \
    trans1(base_type); trans2(base_type); base_type beta(0);               \
    int m(mat_nrows(A)), lda(m), n(mat_ncols(A)), inc(1);                  \
    blas_name(&t, &m, &n, &alpha, &A(0,0), &lda, &x[0], &inc, &beta,       \
              &z[0], &inc);                                                \
  }

  // Y <- AX.
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, sgemv_,
		  BLAS_S, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, dgemv_,
		  BLAS_D, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, cgemv_,
		  BLAS_C, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_n, gemv_trans2_n, zgemv_,
		  BLAS_Z, col_major);

  // Y <- transposed(A)X.
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_n, gemv_trans2_n, zgemv_,
		  BLAS_Z, row_major);
  
  // Y <- transposed(const A)X.
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_n, gemv_trans2_n, zgemv_,
		  BLAS_Z, row_major);
  
  // Y <- conjugated(A)X.
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_n, gemv_trans2_n, zgemv_,
		  BLAS_Z, row_major);

  // Y <- A scaled(X).
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, sgemv_,
		  BLAS_S, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, dgemv_,
		  BLAS_D, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, cgemv_,
		  BLAS_C, col_major);
  gemv_interface2(gem_p1_n, gem_trans1_n, gemv_p2_s, gemv_trans2_s, zgemv_,
		  BLAS_Z, col_major);

  // Y <- transposed(A) scaled(X).
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_t, gem_trans1_t, gemv_p2_s, gemv_trans2_s, zgemv_,
		  BLAS_Z, row_major);
  
  // Y <- transposed(const A) scaled(X).
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_tc, gem_trans1_t, gemv_p2_s, gemv_trans2_s, zgemv_,
		  BLAS_Z, row_major);
  
  // Y <- conjugated(A) scaled(X).
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, sgemv_,
		  BLAS_S, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, dgemv_,
		  BLAS_D, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, cgemv_,
		  BLAS_C, row_major);
  gemv_interface2(gem_p1_c, gem_trans1_c, gemv_p2_s, gemv_trans2_s, zgemv_,
		  BLAS_Z, row_major);

  /* ********************************************************************* */
  /* dense matrix x dense matrix multiplication.                           */
  /* ********************************************************************* */

# define gemm_interface_nn(blas_name, base_type)                           \
  inline void mult_spec(const dense_matrix<base_type > &A,                 \
            const dense_matrix<base_type > &B,                             \
            dense_matrix<base_type > &C, c_mult) {                         \
    const char t = 'N';                                                    \
    int m = mat_nrows(A), lda = m, k = mat_ncols(A), n = mat_ncols(B);     \
    int ldb = k, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &t, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_nn(sgemm_, BLAS_S);
  gemm_interface_nn(dgemm_, BLAS_D);
  gemm_interface_nn(cgemm_, BLAS_C);
  gemm_interface_nn(zgemm_, BLAS_Z);
  
  /* ********************************************************************* */
  /* transposed(dense matrix) x dense matrix multiplication.               */
  /* ********************************************************************* */

# define gemm_interface_tn(blas_name, base_type, is_const)                 \
  inline void mult_spec(                                                   \
         const transposed_col_ref<is_const dense_matrix<base_type > *> &_A,\
         const dense_matrix<base_type > &B,                                \
         dense_matrix<base_type > &C, rcmult) {                            \
    dense_matrix<base_type > &A =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_A));\
    const char t = 'T', u = 'N';                                           \
    int m = mat_ncols(A), k = mat_nrows(A), n = mat_ncols(B), lda = k;     \
    int ldb = k, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_tn(sgemm_, BLAS_S,);
  gemm_interface_tn(dgemm_, BLAS_D,);
  gemm_interface_tn(cgemm_, BLAS_C,);
  gemm_interface_tn(zgemm_, BLAS_Z,);
  gemm_interface_tn(sgemm_, BLAS_S, const);
  gemm_interface_tn(dgemm_, BLAS_D, const);
  gemm_interface_tn(cgemm_, BLAS_C, const);
  gemm_interface_tn(zgemm_, BLAS_Z, const);

  /* ********************************************************************* */
  /* dense matrix x transposed(dense matrix) multiplication.               */
  /* ********************************************************************* */

# define gemm_interface_nt(blas_name, base_type, is_const)                 \
  inline void mult_spec(const dense_matrix<base_type > &A,                 \
         const transposed_col_ref<is_const dense_matrix<base_type > *> &_B,\
         dense_matrix<base_type > &C, c_mult) {                            \
    dense_matrix<base_type > &B =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_B));\
    const char t = 'N', u = 'T';                                           \
    int m = mat_nrows(A), lda = m, k = mat_ncols(A), n = mat_nrows(B);     \
    int ldb = n, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_nt(sgemm_, BLAS_S,);
  gemm_interface_nt(dgemm_, BLAS_D,);
  gemm_interface_nt(cgemm_, BLAS_C,);
  gemm_interface_nt(zgemm_, BLAS_Z,);
  gemm_interface_nt(sgemm_, BLAS_S, const);
  gemm_interface_nt(dgemm_, BLAS_D, const);
  gemm_interface_nt(cgemm_, BLAS_C, const);
  gemm_interface_nt(zgemm_, BLAS_Z, const);

  /* ********************************************************************* */
  /* transposed(dense matrix) x transposed(dense matrix) multiplication.   */
  /* ********************************************************************* */

# define gemm_interface_tt(blas_name, base_type, isA_const, isB_const)     \
  inline void mult_spec(                                                   \
        const transposed_col_ref<isA_const dense_matrix<base_type > *> &_A,\
        const transposed_col_ref<isB_const dense_matrix<base_type > *> &_B,\
        dense_matrix<base_type > &C, r_mult) {                             \
    dense_matrix<base_type > &A =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_A));\
    dense_matrix<base_type > &B =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_B));\
    const char t = 'T', u = 'T';                                           \
    int m = mat_ncols(A), k = mat_nrows(A), n = mat_nrows(B), lda = k;     \
    int ldb = n, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_tt(sgemm_, BLAS_S,,);
  gemm_interface_tt(dgemm_, BLAS_D,,);
  gemm_interface_tt(cgemm_, BLAS_C,,);
  gemm_interface_tt(zgemm_, BLAS_Z,,);
  gemm_interface_tt(sgemm_, BLAS_S, const,);
  gemm_interface_tt(dgemm_, BLAS_D, const,);
  gemm_interface_tt(cgemm_, BLAS_C, const,);
  gemm_interface_tt(zgemm_, BLAS_Z, const,);
  gemm_interface_tt(sgemm_, BLAS_S,, const);
  gemm_interface_tt(dgemm_, BLAS_D,, const);
  gemm_interface_tt(cgemm_, BLAS_C,, const);
  gemm_interface_tt(zgemm_, BLAS_Z,, const);
  gemm_interface_tt(sgemm_, BLAS_S, const, const);
  gemm_interface_tt(dgemm_, BLAS_D, const, const);
  gemm_interface_tt(cgemm_, BLAS_C, const, const);
  gemm_interface_tt(zgemm_, BLAS_Z, const, const);


  /* ********************************************************************* */
  /* conjugated(dense matrix) x dense matrix multiplication.               */
  /* ********************************************************************* */

# define gemm_interface_cn(blas_name, base_type)                           \
  inline void mult_spec(                                                   \
      const conjugated_col_matrix_const_ref<dense_matrix<base_type > > &_A,\
      const dense_matrix<base_type > &B,                                   \
      dense_matrix<base_type > &C, rcmult) {                               \
    dense_matrix<base_type > &A =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_A));\
    const char t = 'C', u = 'N';                                           \
    int m = mat_ncols(A), k = mat_nrows(A), n = mat_ncols(B), lda = k;     \
    int ldb = k, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_cn(sgemm_, BLAS_S);
  gemm_interface_cn(dgemm_, BLAS_D);
  gemm_interface_cn(cgemm_, BLAS_C);
  gemm_interface_cn(zgemm_, BLAS_Z);

  /* ********************************************************************* */
  /* dense matrix x conjugated(dense matrix) multiplication.               */
  /* ********************************************************************* */

# define gemm_interface_nc(blas_name, base_type)                           \
  inline void mult_spec(const dense_matrix<base_type > &A,                 \
      const conjugated_col_matrix_const_ref<dense_matrix<base_type > > &_B,\
      dense_matrix<base_type > &C, c_mult) {                               \
    dense_matrix<base_type > &B =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_B));\
    const char t = 'N', u = 'C';                                           \
    int m = mat_nrows(A), lda = m, k = mat_ncols(A), n = mat_nrows(B);     \
    int ldb = n, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_nc(sgemm_, BLAS_S);
  gemm_interface_nc(dgemm_, BLAS_D);
  gemm_interface_nc(cgemm_, BLAS_C);
  gemm_interface_nc(zgemm_, BLAS_Z);

  /* ********************************************************************* */
  /* conjugated(dense matrix) x conjugated(dense matrix) multiplication.   */
  /* ********************************************************************* */

# define gemm_interface_cc(blas_name, base_type)                           \
  inline void mult_spec(                                                   \
      const conjugated_col_matrix_const_ref<dense_matrix<base_type > > &_A,\
      const conjugated_col_matrix_const_ref<dense_matrix<base_type > > &_B,\
      dense_matrix<base_type > &C, r_mult) {                               \
    dense_matrix<base_type > &A =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_A));\
    dense_matrix<base_type > &B =                                          \
                          *(dense_matrix<base_type > *)(linalg_origin(_B));\
    const char t = 'C', u = 'C';                                           \
    int m = mat_ncols(A), k = mat_nrows(A), lda = k, n = mat_nrows(B);     \
    int ldb = n, ldc = m;                                                  \
    base_type alpha(1), beta(0);                                           \
    blas_name(&t, &u, &m, &n, &k, &alpha,                                  \
	        &A(0,0), &lda, &B(0,0), &ldb, &beta, &C(0,0), &ldc);       \
  }

  gemm_interface_cc(sgemm_, BLAS_S);
  gemm_interface_cc(dgemm_, BLAS_D);
  gemm_interface_cc(cgemm_, BLAS_C);
  gemm_interface_cc(zgemm_, BLAS_Z);
   
  /* ********************************************************************* */
  /* Tri solve.                                                            */
  /* ********************************************************************* */

# define trsv_interface(f_name, loru, param1, trans1, blas_name, base_type)\
  inline void f_name(param1(base_type), std::vector<base_type > &x,        \
                              size_type k, bool is_unit) {                 \
    loru; trans1(base_type); char d = is_unit ? 'U' : 'N';                 \
    int lda(mat_nrows(A)), inc(1), n(k);                                   \
    blas_name(&l, &t, &d, &n, &A(0,0), &lda, &x[0], &inc);                 \
  }

# define trsv_upper const char l = 'U'
# define trsv_lower const char l = 'L'

  // X <- LOWER(A)^{-1}X.
  trsv_interface(lower_tri_solve, trsv_lower, gem_p1_n, gem_trans1_n,
		 strsv_, BLAS_S);
  trsv_interface(lower_tri_solve, trsv_lower, gem_p1_n, gem_trans1_n,
		 dtrsv_, BLAS_D); 
  trsv_interface(lower_tri_solve, trsv_lower, gem_p1_n, gem_trans1_n,
		 ctrsv_, BLAS_C); 
  trsv_interface(lower_tri_solve, trsv_lower, gem_p1_n, gem_trans1_n,
		 ztrsv_, BLAS_Z);
  
  // X <- UPPER(A)^{-1}X.
  trsv_interface(upper_tri_solve, trsv_upper, gem_p1_n, gem_trans1_n,
		 strsv_, BLAS_S);
  trsv_interface(upper_tri_solve, trsv_upper, gem_p1_n, gem_trans1_n,
		 dtrsv_, BLAS_D); 
  trsv_interface(upper_tri_solve, trsv_upper, gem_p1_n, gem_trans1_n,
		 ctrsv_, BLAS_C); 
  trsv_interface(upper_tri_solve, trsv_upper, gem_p1_n, gem_trans1_n,
		 ztrsv_, BLAS_Z);
  
  // X <- LOWER(transposed(A))^{-1}X.
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_t, gem_trans1_t,
		 strsv_, BLAS_S);
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_t, gem_trans1_t,
		 dtrsv_, BLAS_D); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_t, gem_trans1_t,
		 ctrsv_, BLAS_C); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_t, gem_trans1_t,
		 ztrsv_, BLAS_Z);
  
  // X <- UPPER(transposed(A))^{-1}X.
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_t, gem_trans1_t,
		 strsv_, BLAS_S);
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_t, gem_trans1_t,
		 dtrsv_, BLAS_D); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_t, gem_trans1_t,
		 ctrsv_, BLAS_C); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_t, gem_trans1_t,
		 ztrsv_, BLAS_Z);

  // X <- LOWER(transposed(const A))^{-1}X.
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_tc, gem_trans1_t,
		 strsv_, BLAS_S);
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_tc, gem_trans1_t,
		 dtrsv_, BLAS_D); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_tc, gem_trans1_t,
		 ctrsv_, BLAS_C); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_tc, gem_trans1_t,
		 ztrsv_, BLAS_Z);
  
  // X <- UPPER(transposed(const A))^{-1}X.
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_tc, gem_trans1_t,
		 strsv_, BLAS_S);
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_tc, gem_trans1_t,
		 dtrsv_, BLAS_D); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_tc, gem_trans1_t,
		 ctrsv_, BLAS_C); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_tc, gem_trans1_t,
		 ztrsv_, BLAS_Z);

  // X <- LOWER(conjugated(A))^{-1}X.
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_c, gem_trans1_c,
		 strsv_, BLAS_S);
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_c, gem_trans1_c,
		 dtrsv_, BLAS_D); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_c, gem_trans1_c,
		 ctrsv_, BLAS_C); 
  trsv_interface(lower_tri_solve, trsv_upper, gem_p1_c, gem_trans1_c,
		 ztrsv_, BLAS_Z);
  
  // X <- UPPER(conjugated(A))^{-1}X.
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_c, gem_trans1_c,
		 strsv_, BLAS_S);
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_c, gem_trans1_c,
		 dtrsv_, BLAS_D); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_c, gem_trans1_c,
		 ctrsv_, BLAS_C); 
  trsv_interface(upper_tri_solve, trsv_lower, gem_p1_c, gem_trans1_c,
		 ztrsv_, BLAS_Z);
  

  /* ********************************************************************* */
  /* LU decomposition.                                                     */
  /* ********************************************************************* */

# define getrf_interface(lapack_name, base_type) inline                    \
  size_type lu_factor(dense_matrix<base_type > &A, std::vector<int> &ipvt){\
    int m(mat_nrows(A)), n(mat_ncols(A)), lda(m), info;                    \
    lapack_name(&m, &n, &A(0,0), &lda, &ipvt[0], &info);                   \
    return size_type(info);                                                \
  }

  getrf_interface(sgetrf_, BLAS_S);
  getrf_interface(dgetrf_, BLAS_D);
  getrf_interface(cgetrf_, BLAS_C);
  getrf_interface(zgetrf_, BLAS_Z);

  /* ********************************************************************* */
  /* LU solve.                                                             */
  /* ********************************************************************* */

# define getrs_interface(f_name, trans1, lapack_name, base_type) inline    \
  void f_name(const dense_matrix<base_type > &A,                           \
	      const std::vector<int> &ipvt, std::vector<base_type > &x,    \
	      const std::vector<base_type > &b) {                          \
    int n(mat_nrows(A)), info, nrhs(1);                                    \
    gmm::copy(b, x); trans1;                                               \
    lapack_name(&t, &n, &nrhs, &(A(0,0)), &n, &ipvt[0], &x[0], &n, &info); \
  }
  
# define getrs_trans_n const char t = 'N'
# define getrs_trans_t const char t = 'T'

  getrs_interface(lu_solve, getrs_trans_n, sgetrs_, BLAS_S);
  getrs_interface(lu_solve, getrs_trans_n, dgetrs_, BLAS_D);
  getrs_interface(lu_solve, getrs_trans_n, cgetrs_, BLAS_C);
  getrs_interface(lu_solve, getrs_trans_n, zgetrs_, BLAS_Z);
  getrs_interface(lu_solve_transposed, getrs_trans_t, sgetrs_, BLAS_S);
  getrs_interface(lu_solve_transposed, getrs_trans_t, dgetrs_, BLAS_D);
  getrs_interface(lu_solve_transposed, getrs_trans_t, cgetrs_, BLAS_C);
  getrs_interface(lu_solve_transposed, getrs_trans_t, zgetrs_, BLAS_Z);

  /* ********************************************************************* */
  /* LU inverse.                                                           */
  /* ********************************************************************* */

# define getri_interface(lapack_name, base_type) inline                    \
  void lu_inverse(const dense_matrix<base_type > &LU,                      \
       std::vector<int> &ipvt, const dense_matrix<base_type > &A_) {       \
    dense_matrix<base_type >&                                              \
    A = const_cast<dense_matrix<base_type > &>(A_);                        \
    int n(mat_nrows(A)), info, lwork(-1); base_type work1;                 \
    gmm::copy(LU, A);                                                      \
    lapack_name(&n, &A(0,0), &n, &ipvt[0], &work1, &lwork, &info);         \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name(&n, &A(0,0), &n, &ipvt[0], &work[0], &lwork, &info);       \
  }

  getri_interface(sgetri_, BLAS_S);
  getri_interface(dgetri_, BLAS_D);
  getri_interface(cgetri_, BLAS_C);
  getri_interface(zgetri_, BLAS_Z);


  /* ********************************************************************* */
  /* QR factorization.                                                     */
  /* ********************************************************************* */
  
# define geqrf_interface(lapack_name1, base_type) inline                   \
  void qr_factor(dense_matrix<base_type > &A){                             \
    int m(mat_nrows(A)), n(mat_ncols(A)), info, lwork(-1); base_type work1;\
    std::vector<base_type > tau(n);                                        \
    lapack_name1(&m, &n, &A(0,0), &m, &tau[0], &work1  , &lwork, &info);   \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name1(&m, &n, &A(0,0), &m, &tau[0], &work[0], &lwork, &info);   \
    if (info) DAL_THROW(failure_error, "QR factorization failed");         \
  }

  geqrf_interface(sgeqrf_, BLAS_S);
  geqrf_interface(dgeqrf_, BLAS_D);
  geqrf_interface(cgeqrf_, BLAS_C);
  geqrf_interface(zgeqrf_, BLAS_Z);


# define geqrf_interface1(lapack_name1, lapack_name2, base_type) inline    \
  void qr_factor(dense_matrix<base_type > &A, dense_matrix<base_type > &Q){\
    int m(mat_nrows(A)), n(mat_ncols(A)), info, lwork(-1); base_type work1;\
    std::vector<base_type > tau(n);                                        \
    lapack_name1(&m, &n, &A(0,0), &m, &tau[0], &work1  , &lwork, &info);   \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name1(&m, &n, &A(0,0), &m, &tau[0], &work[0], &lwork, &info);   \
    if (info) DAL_THROW(failure_error, "QR factorization failed");         \
    gmm::copy(A, Q);                                                       \
    lapack_name2(&m, &n, &n, &Q(0,0), &m, &tau[0], &work[0],&lwork,&info); \
  }

  geqrf_interface1(sgeqrf_, sorgqr_, BLAS_S);
  geqrf_interface1(dgeqrf_, dorgqr_, BLAS_D);
  geqrf_interface1(cgeqrf_, cungqr_, BLAS_C);
  geqrf_interface1(zgeqrf_, zungqr_, BLAS_Z);



# define geqrf_interface2(lapack_name1, lapack_name2, base_type) inline    \
  void qr_factor(const dense_matrix<base_type > &A,                        \
       dense_matrix<base_type > &Q, dense_matrix<base_type > &R) {         \
    int m(mat_nrows(A)), n(mat_ncols(A)), info, lwork(-1); base_type work1;\
    gmm::copy(A, Q);                                                       \
    std::vector<base_type > tau(n);                                        \
    lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work1  , &lwork, &info);   \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work[0], &lwork, &info);   \
    if (info) DAL_THROW(failure_error, "QR factorization failed");         \
    base_type *p = &R(0,0), *q = &Q(0,0);                                  \
    for (int j = 0; j < n; ++j, q += m-n)                                  \
      for (int i = 0; i < n; ++i, ++p, ++q)                                \
        *p = (j < i) ? base_type(0) : *q;                                  \
    lapack_name2(&m, &n, &n, &Q(0,0), &m, &tau[0], &work[0],&lwork,&info); \
  }

  geqrf_interface2(sgeqrf_, sorgqr_, BLAS_S);
  geqrf_interface2(dgeqrf_, dorgqr_, BLAS_D);
  geqrf_interface2(cgeqrf_, cungqr_, BLAS_C);
  geqrf_interface2(zgeqrf_, zungqr_, BLAS_Z);

  /* ********************************************************************* */
  /* QR algorithm for eigenvalues search.                                  */
  /* ********************************************************************* */

# define gees_interface(lapack_name, base_type)                            \
  template <class VECT> inline void implicit_qr_algorithm(                 \
         const dense_matrix<base_type > &A,  const VECT &eigval_,          \
         dense_matrix<base_type > &Q,                                      \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    int n(mat_nrows(A)), info, lwork(-1), sdim; base_type work1;           \
    dense_matrix<base_type > H(n,n); gmm::copy(A, H);                      \
    char jobvs = (compvect ? 'V' : 'N'), sort = 'N';                       \
    std::vector<double> rwork(n), eigv1(n), eigv2(n);                      \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
                &eigv2[0], &Q(0,0), &n, &work1, &lwork, &rwork[0], &info); \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
		&eigv2[0], &Q(0,0), &n, &work[0], &lwork, &rwork[0],&info);\
    if (info) DAL_THROW(failure_error, "QR algorithm failed");             \
    extract_eig(H, const_cast<VECT &>(eigval_), tol);                      \
  }

# define gees_interface2(lapack_name, base_type)                           \
  template <class VECT> inline void implicit_qr_algorithm(                 \
         const dense_matrix<base_type > &A,  const VECT &eigval_,          \
         dense_matrix<base_type > &Q,                                      \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    int n(mat_nrows(A)), info, lwork(-1), sdim; base_type work1;           \
    dense_matrix<base_type > H(n,n); gmm::copy(A, H);                      \
    char jobvs = (compvect ? 'V' : 'N'), sort = 'N';                       \
    std::vector<double> rwork(n), eigvv(n*2);                              \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work1, &lwork, &rwork[0], &rwork[0], &info); \
    lwork = int(dal::real(work1));                                         \
    std::vector<base_type > work(lwork);                                   \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work[0], &lwork, &rwork[0], &rwork[0],&info);\
    if (info) DAL_THROW(failure_error, "QR algorithm failed");             \
    extract_eig(H, const_cast<VECT &>(eigval_), tol);                      \
  }

  gees_interface(sgees_, BLAS_S);
  gees_interface(dgees_, BLAS_D);
  gees_interface2(cgees_, BLAS_C);
  gees_interface2(zgees_, BLAS_Z);

}

#endif // GMM_LAPACK_INTERFACE_H

#endif // GMM_USES_LAPACK || GMM_USES_ATLAS
