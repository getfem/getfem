/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_superlu_interface.h : interface with superlu,            */
/*            LU factorization and solve for sparse matrices.              */
/*     									   */
/* Date : October 17, 2003.                                                */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
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


#if defined(GMM_USES_SUPERLU)

#ifndef GMM_SUPERLU_INTERFACE_H
#define GMM_SUPERLU_INTERFACE_H

#include <gmm_kernel.h>

typedef int int_t;
#include "SRC/Cnames.h"
#include "SRC/supermatrix.h"
#include "SRC/util.h"

namespace SuperLU_S {
#include "SRC/ssp_defs.h"
}
namespace SuperLU_D {
#include "SRC/dsp_defs.h"
}
namespace SuperLU_C {
#include "SRC/csp_defs.h"
}
namespace SuperLU_Z {
#include "SRC/zsp_defs.h" 
}



namespace gmm {

  /*  interface for Create_CompCol_Matrix */

  void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     float *a, int *ir, int *jc) {
    SuperLU_S::sCreate_CompCol_Matrix(A, m, n, nnz, a, ir, jc,
				      SLU_NC, SLU_S, SLU_GE);
  }
  
  void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     double *a, int *ir, int *jc) {
    SuperLU_D::dCreate_CompCol_Matrix(A, m, n, nnz, a, ir, jc,
				      SLU_NC, SLU_D, SLU_GE);
  }
  
  void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     std::complex<float> *a, int *ir, int *jc) {
    SuperLU_C::cCreate_CompCol_Matrix(A, m, n, nnz, (SuperLU_C::complex *)(a),
				      ir, jc, SLU_NC, SLU_C, SLU_GE);
  }
  
  void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     std::complex<double> *a, int *ir, int *jc) {
    SuperLU_Z::zCreate_CompCol_Matrix(A, m, n, nnz,
				      (SuperLU_Z::doublecomplex *)(a), ir, jc,
				      SLU_NC, SLU_Z, SLU_GE);
  }

  /*  interface for Create_Dense_Matrix */

  void Create_Dense_Matrix(SuperMatrix *A, int m, int n, float *a, int k)
  { SuperLU_S::sCreate_Dense_Matrix(A, m, n, a, k, SLU_DN, SLU_S, SLU_GE); }
  void Create_Dense_Matrix(SuperMatrix *A, int m, int n, double *a, int k)
  { SuperLU_D::dCreate_Dense_Matrix(A, m, n, a, k, SLU_DN, SLU_D, SLU_GE); }
  void Create_Dense_Matrix(SuperMatrix *A, int m, int n,
			   std::complex<float> *a, int k) {
    SuperLU_C::cCreate_Dense_Matrix(A, m, n, (SuperLU_C::complex *)(a),
				    k, SLU_DN, SLU_C, SLU_GE);
  }
  void Create_Dense_Matrix(SuperMatrix *A, int m, int n, 
			   std::complex<double> *a, int k) {
    SuperLU_Z::zCreate_Dense_Matrix(A, m, n, (SuperLU_Z::doublecomplex *)(a),
				    k, SLU_DN, SLU_Z, SLU_GE);
  }

  /*  interface for gssv */

#define DECL_GSSV(NAMESPACE,FNAME,FLOATTYPE,KEYTYPE) \
  void SuperLU_gssv(superlu_options_t *options, SuperMatrix *A, int *p, \
  int *q, SuperMatrix *L, SuperMatrix *U, SuperMatrix *B,               \
  SuperLUStat_t *stats, int *info, KEYTYPE) {                           \
  NAMESPACE::FNAME(options, A, p, q, L, U, B, stats, info);             \
  }

  DECL_GSSV(SuperLU_S,sgssv,float,float)
  DECL_GSSV(SuperLU_C,cgssv,float,std::complex<float>)
  DECL_GSSV(SuperLU_D,dgssv,double,double)
  DECL_GSSV(SuperLU_Z,zgssv,double,std::complex<double>)


#define DECL_GSSVX(NAMESPACE,FNAME,FLOATTYPE,KEYTYPE) \
  void SuperLU_gssvx(superlu_options_t *options, SuperMatrix *A,         \
		     int *perm_c, int *perm_r, int *etree, char *equed,  \
		     FLOATTYPE *R, FLOATTYPE *C, SuperMatrix *L,         \
		     SuperMatrix *U, void *work, int lwork,              \
		     SuperMatrix *B, SuperMatrix *X,                     \
		     FLOATTYPE *recip_pivot_growth,                      \
		     FLOATTYPE *rcond, FLOATTYPE *ferr, FLOATTYPE *berr, \
		     SuperLUStat_t *stats, int *info, KEYTYPE) {         \
    NAMESPACE::mem_usage_t mem_usage;                                    \
    NAMESPACE::FNAME(options, A, perm_c, perm_r, etree, equed, R, C, L,  \
		     U, work, lwork, B, X, recip_pivot_growth, rcond,    \
		     ferr, berr, &mem_usage, stats, info);               \
  }

  DECL_GSSVX(SuperLU_S,sgssvx,float,float)
  DECL_GSSVX(SuperLU_C,cgssvx,float,std::complex<float>)
  DECL_GSSVX(SuperLU_D,dgssvx,double,double)
  DECL_GSSVX(SuperLU_Z,zgssvx,double,std::complex<double>)

  /* ********************************************************************* */
  /*   SuperLU solve                                                       */
  /* ********************************************************************* */

  template <typename MAT, typename VECTX, typename VECTB>
  void SuperLU_solve(const MAT &A, const VECTX &X_, const VECTB &B,
		     int permc_spec = 1) {
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     *   permc_spec = 3: use approximate minimum degree column ordering
     */
    VECTX &X = const_cast<VECTX &>(X_);
    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    int m = mat_nrows(A), n = mat_ncols(A), nrhs = 1, info;

    csc_matrix<T> csc_A(m, n); gmm::copy(A, csc_A);
    std::vector<T> rhs(m);
    gmm::copy(B, rhs);

    int nz = nnz(csc_A);
    if ((2 * nz / n) >= m)
      DAL_WARNING(2, "CAUTION : it seems that SuperLU has a problem"
		  " for nearly dense sparse matrices");

    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = NATURAL;
    switch (permc_spec) {
    case 1 : options.ColPerm = MMD_ATA; break;
    case 2 : options.ColPerm = MMD_AT_PLUS_A; break;
    case 3 : options.ColPerm = COLAMD; break;
    }
    SuperLUStat_t stat;
    StatInit(&stat);

    SuperMatrix SA, SL, SU, SB; // SuperLU format.
    Create_CompCol_Matrix(&SA, m, n, nz, csc_A.pr,
			  (int *)(csc_A.ir), (int *)(csc_A.jc));
    Create_Dense_Matrix(&SB, m, nrhs, &rhs[0], m);
    
    std::vector<int> perm_r(m), perm_c(3*n);

    SuperLU_gssv(&options, &SA, &perm_c[n], &perm_r[0], &SL, &SU, &SB,
		 &stat, &info, T());
    if (info != 0)
      DAL_THROW(failure_error, "SuperLU solve failed: info=" << info);
    gmm::copy(rhs, X);
    Destroy_SuperMatrix_Store(&SB);
    Destroy_SuperMatrix_Store(&SA);
    Destroy_SuperNode_Matrix(&SL);
    Destroy_CompCol_Matrix(&SU);
  }

  template <typename MAT, typename VECTX, typename VECTB>
  void SuperLU_solve(const MAT &A, const VECTX &X_, const VECTB &B,
		     double& rcond_, int permc_spec = 1) {
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     *   permc_spec = 3: use approximate minimum degree column ordering
     */
    VECTX &X = const_cast<VECTX &>(X_);
    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    int m = mat_nrows(A), n = mat_ncols(A), nrhs = 1, info=0;

    csc_matrix<T> csc_A(m, n); gmm::copy(A, csc_A);
    std::vector<T> rhs(m), sol(m);
    gmm::copy(B, rhs);

    int nz = nnz(csc_A);
    if ((2 * nz / n) >= m)
      DAL_WARNING(2, "CAUTION : it seems that SuperLU has a problem"
		  " for nearly dense sparse matrices");

    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = NATURAL;
    switch (permc_spec) {
    case 1 : options.ColPerm = MMD_ATA; break;
    case 2 : options.ColPerm = MMD_AT_PLUS_A; break;
    case 3 : options.ColPerm = COLAMD; break;
    }
    SuperLUStat_t stat;
    StatInit(&stat);


    SuperMatrix SA, SL, SU, SB, SX; // SuperLU format.
    Create_CompCol_Matrix(&SA, m, n, nz, csc_A.pr,
			  (int *)(csc_A.ir), (int *)(csc_A.jc));
    Create_Dense_Matrix(&SB, m, nrhs, &rhs[0], m);
    Create_Dense_Matrix(&SX, m, nrhs, &sol[0], m);

    std::vector<int> etree(n);
    char equed[] = "B";
    // fact[] = "E", refact[] = "N", istrans[] = "N";
    std::vector<R> Rscale(m),Cscale(n); // row scale factors
    std::vector<R> ferr(nrhs), berr(nrhs);
    R recip_pivot_gross, rcond;
    // perm_c oversized to turn around a bug (?) of superlu 
    // with almost full matrices.
    std::vector<int> perm_r(m), perm_c(3*n);
    
    SuperLU_gssvx(&options, &SA, &perm_c[n], &perm_r[0], 
		  &etree[0] /* output */, equed /* output         */, 
		  &Rscale[0] /* row scale factors (output)        */, 
		  &Cscale[0] /* col scale factors (output)        */,
		  &SL /* fact L (output)*/, &SU /* fact U (output)*/, 
		  NULL /* work                                    */, 
		  0 /* lwork: superlu auto allocates (input)      */, 
		  &SB /* rhs */, &SX /* solution                  */,
		  &recip_pivot_gross /* reciprocal pivot growth   */
		  /* factor max_j( norm(A_j)/norm(U_j) ).         */,  
		  &rcond /*estimate of the reciprocal condition   */
		  /* number of the matrix A after equilibration   */,
		  &ferr[0] /* estimated forward error             */,
		  &berr[0] /* relative backward error             */,
		  &stat, &info, T());
    if (info != 0)
      DAL_THROW(failure_error, "SuperLU solve failed: info=" << info);
    gmm::copy(sol, X);
    rcond_ = rcond;
    Destroy_SuperMatrix_Store(&SB);
    Destroy_SuperMatrix_Store(&SX);
    Destroy_SuperMatrix_Store(&SA);
    Destroy_SuperNode_Matrix(&SL);
    Destroy_CompCol_Matrix(&SU);
  }
}

  
#endif // GMM_SUPERLU_INTERFACE_H

#endif // GMM_USES_SUPERLU
