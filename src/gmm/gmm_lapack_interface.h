/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2003-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file gmm_lapack_interface.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date October 7, 2003.
   @brief gmm interface for LAPACK
*/

#ifndef GMM_LAPACK_INTERFACE_H
#define GMM_LAPACK_INTERFACE_H

#include "gmm_blas_interface.h"
#include "gmm_dense_lu.h"
#include "gmm_dense_qr.h"

#if defined(GMM_USES_LAPACK) && !defined(GMM_MATLAB_INTERFACE)

namespace gmm {

  /* ********************************************************************** */
  /* Operations interfaced for T = float, double, std::complex<float>       */
  /*    or std::complex<double> :                                           */
  /*                                                                        */
  /* lu_factor(dense_matrix<T>, std::vector<long>)                          */
  /* lu_solve(dense_matrix<T>, std::vector<T>, std::vector<T>)              */
  /* lu_solve(dense_matrix<T>, std::vector<long>, std::vector<T>,           */
  /*          std::vector<T>)                                               */
  /* lu_solve_transposed(dense_matrix<T>, std::vector<long>, std::vector<T>,*/
  /*          std::vector<T>)                                               */
  /* lu_inverse(dense_matrix<T>)                                            */
  /* lu_inverse(dense_matrix<T>, std::vector<long>, dense_matrix<T>)        */
  /*                                                                        */
  /* qr_factor(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)           */
  /*                                                                        */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>)                 */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<T>,                 */
  /*                       dense_matrix<T>)                                 */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >)  */
  /* implicit_qr_algorithm(dense_matrix<T>, std::vector<std::complex<T> >,  */
  /*                       dense_matrix<T>)                                 */
  /*                                                                        */
  /* geev_interface_right                                                   */
  /* geev_interface_left                                                    */
  /*                                                                        */
  /* schur(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>)               */
  /*                                                                        */
  /* svd(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>, std::vector<T>) */
  /* svd(dense_matrix<T>, dense_matrix<T>, dense_matrix<T>,                 */
  /*     std::vector<std::complex<T> >)                                     */
  /*                                                                        */
  /* ********************************************************************** */

  /* ********************************************************************** */
  /* LAPACK functions used.                                                 */
  /* ********************************************************************** */

  extern "C" {
    void sgetrf_(...); void dgetrf_(...); void cgetrf_(...); void zgetrf_(...);
    void sgetrs_(...); void dgetrs_(...); void cgetrs_(...); void zgetrs_(...);
    void sgetri_(...); void dgetri_(...); void cgetri_(...); void zgetri_(...);
    void sgeqrf_(...); void dgeqrf_(...); void cgeqrf_(...); void zgeqrf_(...);
    void sorgqr_(...); void dorgqr_(...); void cungqr_(...); void zungqr_(...);
    void sormqr_(...); void dormqr_(...); void cunmqr_(...); void zunmqr_(...);
    void sgees_ (...); void dgees_ (...); void cgees_ (...); void zgees_ (...);
    void sgeev_ (...); void dgeev_ (...); void cgeev_ (...); void zgeev_ (...);
    void sgeesx_(...); void dgeesx_(...); void cgeesx_(...); void zgeesx_(...);
    void sgesvd_(...); void dgesvd_(...); void cgesvd_(...); void zgesvd_(...);
  }

  /* ********************************************************************** */
  /* LU decomposition.                                                      */
  /* ********************************************************************** */

# define getrf_interface(lapack_name, base_type)                              \
  inline size_type                                                            \
  lu_factor(dense_matrix<base_type> &A, lapack_ipvt &ipvt) {                  \
    GMMLAPACK_TRACE("getrf_interface");                                       \
    const BLAS_INT m=BLAS_INT(mat_nrows(A)), n=BLAS_INT(mat_ncols(A)), lda(m);\
    BLAS_INT info(-1);                                                        \
    if (m && n) lapack_name(&m, &n, &A(0,0), &lda, &ipvt[0], &info);          \
    return size_type(abs(info));                                              \
  }

  getrf_interface(sgetrf_, BLAS_S)
  getrf_interface(dgetrf_, BLAS_D)
  getrf_interface(cgetrf_, BLAS_C)
  getrf_interface(zgetrf_, BLAS_Z)

  /* ********************************************************************* */
  /* LU solve.                                                             */
  /* ********************************************************************* */

# define getrs_interface(f_name, NorT, lapack_name, base_type)             \
  inline void                                                              \
  f_name(const dense_matrix<base_type> &A, const lapack_ipvt &ipvt,        \
         std::vector<base_type> &x, const std::vector<base_type> &b) {     \
    GMMLAPACK_TRACE("getrs_interface");                                    \
    const char t=NorT; const BLAS_INT n=BLAS_INT(mat_nrows(A)), nrhs(1);   \
    BLAS_INT info(0); gmm::copy(b, x);                                     \
    if (n) lapack_name(&t, &n, &nrhs, &A(0,0), &n, &ipvt[0], &x[0], &n,    \
                       &info);                                             \
  }

  getrs_interface(lu_solve, 'N', sgetrs_, BLAS_S)
  getrs_interface(lu_solve, 'N', dgetrs_, BLAS_D)
  getrs_interface(lu_solve, 'N', cgetrs_, BLAS_C)
  getrs_interface(lu_solve, 'N', zgetrs_, BLAS_Z)
  getrs_interface(lu_solve_transposed, 'T', sgetrs_, BLAS_S)
  getrs_interface(lu_solve_transposed, 'T', sgetrs_, BLAS_D)
  getrs_interface(lu_solve_transposed, 'T', sgetrs_, BLAS_C)
  getrs_interface(lu_solve_transposed, 'T', sgetrs_, BLAS_Z)

  /* ********************************************************************* */
  /* LU inverse.                                                           */
  /* ********************************************************************* */

# define getri_interface(lapack_name, base_type)                           \
  inline void lu_inverse(const dense_matrix<base_type> &LU,                \
                         const lapack_ipvt &ipvt,                          \
                         dense_matrix<base_type> &A) {                     \
    GMMLAPACK_TRACE("getri_interface");                                    \
    const BLAS_INT n=BLAS_INT(mat_nrows(A));                               \
    BLAS_INT info(0), lwork(-1); base_type work1;                          \
    if (n) {                                                               \
      gmm::copy(LU, A);                                                    \
      lapack_name(&n, &A(0,0), &n, &ipvt[0], &work1, &lwork, &info);       \
      const size_type worksize=size_type(gmm::real(work1));                \
      std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);     \
      lapack_name(&n, &A(0,0), &n, &ipvt[0], &work[0], &lwork, &info);     \
    }                                                                      \
  }

  getri_interface(sgetri_, BLAS_S)
  getri_interface(dgetri_, BLAS_D)
  getri_interface(cgetri_, BLAS_C)
  getri_interface(zgetri_, BLAS_Z)

  /* ********************************************************************** */
  /* QR factorization.                                                      */
  /* ********************************************************************** */

# define geqrf_interface(lapack_name, base_type)                           \
  inline void qr_factor(dense_matrix<base_type> &A) {                      \
    GMMLAPACK_TRACE("geqrf_interface");                                    \
    const size_type mm(mat_nrows(A)), nn(mat_ncols(A));                    \
    if (mm && nn) {                                                        \
      const BLAS_INT m=BLAS_INT(mm), n=BLAS_INT(nn);                       \
      BLAS_INT info(0), lwork(-1);                                         \
      std::vector<base_type> tau(nn); base_type work1;                     \
      lapack_name(&m, &n, &A(0,0), &m, &tau[0], &work1, &lwork, &info);    \
      const size_type worksize=size_type(gmm::real(work1));                \
      std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);     \
      lapack_name(&m, &n, &A(0,0), &m, &tau[0], &work[0], &lwork, &info);  \
      GMM_ASSERT1(!info, "QR factorization failed");                       \
    }                                                                      \
  }

  geqrf_interface(sgeqrf_, BLAS_S)
  geqrf_interface(dgeqrf_, BLAS_D)
    // For complex values, housholder vectors are not the same as in
    // gmm::lu_factor. Impossible to interface for the moment.
    //  geqrf_interface(cgeqrf_, BLAS_C)
    //  geqrf_interface(zgeqrf_, BLAS_Z)

# define geqrf_interface2(lapack_name1, lapack_name2, base_type) inline    \
  void qr_factor(const dense_matrix<base_type> &A,                         \
                 dense_matrix<base_type> &Q, dense_matrix<base_type> &R) { \
    GMMLAPACK_TRACE("geqrf_interface2");                                   \
    const size_type mm(mat_nrows(A)), nn(mat_ncols(A));                    \
    if (mm && nn) {                                                        \
      const BLAS_INT m=BLAS_INT(mm), n=BLAS_INT(nn);                       \
      BLAS_INT info(0), lwork(-1);                                         \
      std::copy(A.begin(), A.end(), Q.begin());                            \
      std::vector<base_type> tau(nn); base_type work1;                     \
      lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work1, &lwork, &info);   \
      const size_type worksize=size_type(gmm::real(work1));                \
      std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);     \
      lapack_name1(&m, &n, &Q(0,0), &m, &tau[0], &work[0], &lwork, &info); \
      GMM_ASSERT1(!info, "QR factorization failed");                       \
      base_type *p = &R(0,0), *q = &Q(0,0);                                \
      for (BLAS_INT j = 0; j < n; ++j, q += m-n)                           \
        for (BLAS_INT i = 0; i < n; ++i, ++p, ++q)                         \
          *p = (j < i) ? base_type(0) : *q;                                \
      lapack_name2(&m, &n, &n, &Q(0,0), &m,&tau[0],&work[0],&lwork,&info); \
    }                                                                      \
    else gmm::clear(Q);                                                    \
  }

  geqrf_interface2(sgeqrf_, sorgqr_, BLAS_S)
  geqrf_interface2(dgeqrf_, dorgqr_, BLAS_D)
  geqrf_interface2(cgeqrf_, cungqr_, BLAS_C)
  geqrf_interface2(zgeqrf_, zungqr_, BLAS_Z)

  /* ********************************************************************** */
  /* QR algorithm for eigenvalues search.                                   */
  /* ********************************************************************** */

# define gees_interface(lapack_name, base_type)                            \
  template <typename VECT> inline void implicit_qr_algorithm(              \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q,                                       \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    GMMLAPACK_TRACE("gees_interface");                                     \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    const size_type nn(mat_nrows(A));                                      \
    if (!nn) return;                                                       \
    dense_matrix<base_type> H(nn,nn); gmm::copy(A, H);                     \
    const char jobvs=(compvect ? 'V' : 'N'), sort='N';                     \
    const BLAS_INT n=BLAS_INT(nn);                                         \
    std::vector<double> rwork(nn), eigv1(nn), eigv2(nn);                   \
    BLAS_INT info(0), lwork(-1), sdim; base_type work1;                    \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
                &eigv2[0], &Q(0,0), &n, &work1, &lwork, &rwork[0], &info); \
    const size_type worksize=size_type(gmm::real(work1));                  \
    std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);       \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigv1[0],       \
                &eigv2[0], &Q(0,0), &n, &work[0], &lwork, &rwork[0],&info);\
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    extract_eig(H, eigval_, tol);                                          \
  }

# define gees_interface_cplx(lapack_name, base_type)                       \
  template <typename VECT> inline void implicit_qr_algorithm(              \
         const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q,                                       \
         double tol=gmm::default_tol(base_type()), bool compvect = true) { \
    GMMLAPACK_TRACE("gees_interface2");                                    \
    typedef bool (*L_fp)(...);  L_fp p = 0;                                \
    const size_type nn(mat_nrows(A));                                      \
    if (!nn) return;                                                       \
    dense_matrix<base_type> H(nn,nn); gmm::copy(A, H);                     \
    const char jobvs=(compvect ? 'V' : 'N'), sort='N';                     \
    const BLAS_INT n=BLAS_INT(nn);                                         \
    std::vector<double> rwork(nn), eigvv(nn*2);                            \
    BLAS_INT info(0), lwork(-1), sdim; base_type work1;                    \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work1, &lwork, &rwork[0], &rwork[0], &info); \
    const size_type worksize=size_type(gmm::real(work1));                  \
    std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);       \
    lapack_name(&jobvs, &sort, p, &n, &H(0,0), &n, &sdim, &eigvv[0],       \
                &Q(0,0), &n, &work[0], &lwork, &rwork[0], &rwork[0],&info);\
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    extract_eig(H, eigval_, tol);                                          \
  }

  gees_interface(sgees_, BLAS_S)
  gees_interface(dgees_, BLAS_D)
  gees_interface_cplx(cgees_, BLAS_C)
  gees_interface_cplx(zgees_, BLAS_Z)


# define geev_interface(f_name, lapack_name, lapack_name_cplx,             \
                        base_type, cplx_type, _jobvl, _jobvr)              \
  template <typename VECT> inline void                                     \
  f_name(const dense_matrix<base_type> &A, VECT &eigval_,                  \
         dense_matrix<base_type> &Q) {                                     \
    GMMLAPACK_TRACE("geev_interface");                                     \
    const size_type nn(mat_nrows(A));                                      \
    if (!nn) return;                                                       \
    dense_matrix<base_type> H(nn,nn); gmm::copy(A, H);                     \
    const char jobvl=_jobvl, jobvr=_jobvl;                                 \
    const BLAS_INT n=BLAS_INT(nn); BLAS_INT info(0), lwork(-1);            \
    std::vector<base_type> eigvr(nn), eigvi(nn);                           \
    base_type work1;                                                       \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigvr[0], &eigvi[0],     \
                &Q(0,0), &n, &Q(0,0), &n, &work1, &lwork, &info);          \
    const size_type worksize=size_type(gmm::real(work1));                  \
    std::vector<base_type> work(worksize); lwork=BLAS_INT(worksize);       \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigvr[0], &eigvi[0],     \
                &Q(0,0), &n, &Q(0,0), &n, &work[0], &lwork, &info);        \
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    gmm::copy(eigvr, gmm::real_part(eigval_));                             \
    gmm::copy(eigvi, gmm::imag_part(eigval_));                             \
  }                                                                        \
  template <typename VECT> inline void                                     \
  f_name(const dense_matrix<cplx_type> &A, VECT &eigval_,                  \
         dense_matrix<cplx_type> &Q) {                                     \
    GMMLAPACK_TRACE("geev_interface");                                     \
    const size_type nn(mat_nrows(A));                                      \
    if (!nn) return;                                                       \
    dense_matrix<cplx_type> H(nn,nn); gmm::copy(A, H);                     \
    const char jobvl=_jobvl, jobvr=_jobvl;                                 \
    const BLAS_INT n=BLAS_INT(nn); BLAS_INT info(0), lwork(-1);            \
    std::vector<cplx_type> eigv(nn); std::vector<base_type> rwork(2*nn);   \
    cplx_type work1;                                                       \
    lapack_name(&jobvl, &jobvr, &n, &H(0,0), &n, &eigv[0], &Q(0,0), &n,    \
                &Q(0,0), &n, &work1, &lwork, &rwork[0], &info);            \
    const size_type worksize=size_type(gmm::real(work1));                  \
    std::vector<cplx_type> work(worksize); lwork=BLAS_INT(worksize);       \
    lapack_name_cplx(&jobvl, &jobvr, &n, &H(0,0), &n, &eigv[0], &Q(0,0),   \
                     &n, &Q(0,0), &n, &work[0], &lwork,  &rwork[0], &info);\
    GMM_ASSERT1(!info, "QR algorithm failed");                             \
    gmm::copy(eigv, eigval_);                                              \
  }

  geev_interface(geev_interface_right, sgeev_, cgeev_, BLAS_S, BLAS_C, 'N', 'V')
  geev_interface(geev_interface_right, dgeev_, zgeev_, BLAS_D, BLAS_Z, 'N', 'V')
  geev_interface(geev_interface_left, sgeev_, cgeev_, BLAS_S, BLAS_C, 'V', 'N')
  geev_interface(geev_interface_left, dgeev_, zgeev_, BLAS_D, BLAS_Z, 'V', 'N')

  /* ********************************************************************** */
  /* SCHUR algorithm:                                                       */
  /*  A = Q*S*(Q^T), with Q orthogonal and S upper quasi-triangula          */
  /* ********************************************************************** */

# define geesx_interface(lapack_name, lapack_name_cplx,                 \
                         base_type, cplx_type)                          \
  inline void schur(const dense_matrix<base_type> &A,                   \
                    dense_matrix<base_type> &S,                         \
                    dense_matrix<base_type> &Q) {                       \
    GMMLAPACK_TRACE("geesx_interface");                                 \
    const size_type mm(mat_nrows(A)), nn(mat_ncols(A)), worksize(8*nn); \
    GMM_ASSERT1(mm==nn, "Schur decomposition requires square matrix");  \
    resize(S, nn, nn); copy(A, S); resize(Q, nn, nn);                   \
    const char jobvs='V', sort='N', sense='N'; const bool select=false; \
    const BLAS_INT n=BLAS_INT(nn), lwork=BLAS_INT(worksize), liwork(1); \
    BLAS_INT info(0), sdim(0); base_type rconde(0), rcondv(0);          \
    std::vector<base_type> work(worksize), wr(nn), wi(nn);              \
    std::vector<BLAS_INT> iwork(1), bwork(1);                           \
    lapack_name(&jobvs, &sort, &select, &sense, &n, &S(0,0), &n,        \
                &sdim, &wr[0], &wi[0], &Q(0,0), &n, &rconde, &rcondv,   \
                &work[0], &lwork, &iwork[0], &liwork, &bwork[0], &info);\
    GMM_ASSERT1(!info, "SCHUR algorithm failed");                       \
  }                                                                     \
  inline void schur(const dense_matrix<cplx_type> &A,                   \
                    dense_matrix<cplx_type> &S,                         \
                    dense_matrix<cplx_type> &Q) {                       \
    GMMLAPACK_TRACE("geesx_interface");                                 \
    const size_type mm(mat_nrows(A)), nn(mat_ncols(A)), worksize(8*nn); \
    GMM_ASSERT1(mm==nn, "Schur decomposition requires square matrix");  \
    resize(S, nn, nn); copy(A, S); resize(Q, nn, nn);                   \
    const char jobvs='V', sort='N', sense='N'; const bool select=false; \
    const BLAS_INT n=BLAS_INT(nn), lwork=BLAS_INT(worksize);            \
    BLAS_INT info(0), sdim(0); base_type rconde(0), rcondv(0);          \
    std::vector<cplx_type> work(worksize), w(nn);                       \
    std::vector<base_type> rwork(worksize);                             \
    std::vector<BLAS_INT> bwork(1);                                     \
    lapack_name_cplx(&jobvs, &sort, &select, &sense, &n, &S(0,0), &n,   \
                     &sdim, &w[0], &Q(0,0), &n, &rconde, &rcondv,       \
                     &work[0], &lwork, &rwork[0], &bwork[0], &info);    \
    GMM_ASSERT1(!info, "SCHUR algorithm failed");                       \
  }

  geesx_interface(sgeesx_, cgeesx_, BLAS_S, BLAS_C)
  geesx_interface(dgeesx_, zgeesx_, BLAS_D, BLAS_Z)


  /* ********************************************************************** */
  /* Interface to SVD. Does not correspond to a Gmm++ functionality.        */
  /* Author : Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>        */
  /* ********************************************************************** */

# define gesvd_interface(lapack_name, lapack_name_cplx,                   \
                         base_type, cplx_type)                            \
  inline void svd(const dense_matrix<base_type> &X,                       \
                  dense_matrix<base_type> &U,                             \
                  dense_matrix<base_type> &Vtransposed,                   \
                  std::vector<base_type> &sigma) {                        \
    GMMLAPACK_TRACE("gesvd_interface");                                   \
    const size_type mm(mat_nrows(X)), nn(mat_ncols(X)),                   \
                    mn_min(mm < nn ? mm : nn), worksize(15*mn_min);       \
    sigma.resize(mn_min); resize(U, mm, mm); resize(Vtransposed, nn, nn); \
    const char job='A'; const BLAS_INT m=BLAS_INT(mm), n=BLAS_INT(nn),    \
                                       lwork=BLAS_INT(worksize);          \
    std::vector<base_type> work(worksize); BLAS_INT info(0);              \
    lapack_name(&job, &job, &m, &n, &X(0,0), &m, &sigma[0], &U(0,0),      \
                &m, &Vtransposed(0,0), &n, &work[0], &lwork, &info);      \
  }                                                                       \
  inline void svd(const dense_matrix<cplx_type> &X,                       \
                  dense_matrix<cplx_type> &U,                             \
                  dense_matrix<cplx_type> &Vtransposed,                   \
                  std::vector<base_type> &sigma) {                        \
    GMMLAPACK_TRACE("gesvd_interface");                                   \
    const size_type mm(mat_nrows(X)), nn(mat_ncols(X)),                   \
                    mn_min(mm < nn ? mm : nn), worksize(15*mn_min);       \
    sigma.resize(mn_min); resize(U, mm, mm); resize(Vtransposed, nn, nn); \
    const char job='A'; const BLAS_INT m=BLAS_INT(mm), n=BLAS_INT(nn),    \
                                       lwork=BLAS_INT(worksize);          \
    std::vector<cplx_type> work(worksize);                                \
    std::vector<base_type> rwork(5*mn_min); BLAS_INT info(0);             \
    lapack_name_cplx(&job, &job, &m, &n, &X(0,0), &m, &sigma[0], &U(0,0), \
                     &m, &Vtransposed(0,0), &n, &work[0], &lwork,         \
                     &rwork[0], &info);                                   \
  }

  gesvd_interface(sgesvd_, cgesvd_, BLAS_S, BLAS_C)
  gesvd_interface(dgesvd_, zgesvd_, BLAS_D, BLAS_Z)

}

#else

namespace gmm
{
template <typename MAT>
void schur(const MAT &, MAT &, MAT &)
{
  GMM_ASSERT1(false, "Use of function schur(A,S,Q) requires GetFEM "
                     "to be built with Lapack");
}

template <typename BLAS_TYPE>
inline void svd(dense_matrix<BLAS_TYPE> &, dense_matrix<BLAS_TYPE> &,
         dense_matrix<BLAS_TYPE> &, std::vector<BLAS_TYPE> &)
{
  GMM_ASSERT1(false, "Use of function svd(X,U,Vtransposed,sigma) requires GetFEM "
                     "to be built with Lapack");
}

}// namespace gmm

#endif // GMM_USES_LAPACK

#endif // GMM_LAPACK_INTERFACE_H
