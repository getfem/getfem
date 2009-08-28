// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2008 Julien Pommier
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "getfem/getfem_superlu.h"

typedef int int_t;

/* because SRC/util.h defines TRUE and FALSE ... */
#ifdef TRUE
# undef TRUE
#endif
#ifdef FALSE
# undef FALSE
#endif

#include "superlu/slu_Cnames.h"
#include "superlu/supermatrix.h"
#include "superlu/slu_util.h"

namespace SuperLU_S {
#include "superlu/slu_sdefs.h"
}
namespace SuperLU_D {
#include "superlu/slu_ddefs.h"
}
namespace SuperLU_C {
#include "superlu/slu_cdefs.h"
}
namespace SuperLU_Z {
#include "superlu/slu_zdefs.h" 
}



namespace gmm {

  /*  interface for Create_CompCol_Matrix */

  inline void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     float *a, int *ir, int *jc) {
    SuperLU_S::sCreate_CompCol_Matrix(A, m, n, nnz, a, ir, jc,
				      SLU_NC, SLU_S, SLU_GE);
  }
  
  inline void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     double *a, int *ir, int *jc) {
    SuperLU_D::dCreate_CompCol_Matrix(A, m, n, nnz, a, ir, jc,
				      SLU_NC, SLU_D, SLU_GE);
  }
  
  inline void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     std::complex<float> *a, int *ir, int *jc) {
    SuperLU_C::cCreate_CompCol_Matrix(A, m, n, nnz, (SuperLU_C::complex *)(a),
				      ir, jc, SLU_NC, SLU_C, SLU_GE);
  }
  
  inline void Create_CompCol_Matrix(SuperMatrix *A, int m, int n, int nnz,
			     std::complex<double> *a, int *ir, int *jc) {
    SuperLU_Z::zCreate_CompCol_Matrix(A, m, n, nnz,
				      (SuperLU_Z::doublecomplex *)(a), ir, jc,
				      SLU_NC, SLU_Z, SLU_GE);
  }

  /*  interface for Create_Dense_Matrix */

  inline void Create_Dense_Matrix(SuperMatrix *A, int m, int n, float *a, int k)
  { SuperLU_S::sCreate_Dense_Matrix(A, m, n, a, k, SLU_DN, SLU_S, SLU_GE); }
  inline void Create_Dense_Matrix(SuperMatrix *A, int m, int n, double *a, int k)
  { SuperLU_D::dCreate_Dense_Matrix(A, m, n, a, k, SLU_DN, SLU_D, SLU_GE); }
  inline void Create_Dense_Matrix(SuperMatrix *A, int m, int n,
			   std::complex<float> *a, int k) {
    SuperLU_C::cCreate_Dense_Matrix(A, m, n, (SuperLU_C::complex *)(a),
				    k, SLU_DN, SLU_C, SLU_GE);
  }
  inline void Create_Dense_Matrix(SuperMatrix *A, int m, int n, 
			   std::complex<double> *a, int k) {
    SuperLU_Z::zCreate_Dense_Matrix(A, m, n, (SuperLU_Z::doublecomplex *)(a),
				    k, SLU_DN, SLU_Z, SLU_GE);
  }

  /*  interface for gssv */

#define DECL_GSSV(NAMESPACE,FNAME,FLOATTYPE,KEYTYPE) \
  inline void SuperLU_gssv(superlu_options_t *options, SuperMatrix *A, int *p, \
  int *q, SuperMatrix *L, SuperMatrix *U, SuperMatrix *B,               \
  SuperLUStat_t *stats, int *info, KEYTYPE) {                           \
  NAMESPACE::FNAME(options, A, p, q, L, U, B, stats, info);             \
  }

  DECL_GSSV(SuperLU_S,sgssv,float,float)
  DECL_GSSV(SuperLU_C,cgssv,float,std::complex<float>)
  DECL_GSSV(SuperLU_D,dgssv,double,double)
  DECL_GSSV(SuperLU_Z,zgssv,double,std::complex<double>)

  /*  interface for gssvx */

#define DECL_GSSVX(NAMESPACE,FNAME,FLOATTYPE,KEYTYPE) \
    inline float SuperLU_gssvx(superlu_options_t *options, SuperMatrix *A,	\
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
    return mem_usage.for_lu; /* bytes used by the factor storage */     \
  }

  DECL_GSSVX(SuperLU_S,sgssvx,float,float)
  DECL_GSSVX(SuperLU_C,cgssvx,float,std::complex<float>)
  DECL_GSSVX(SuperLU_D,dgssvx,double,double)
  DECL_GSSVX(SuperLU_Z,zgssvx,double,std::complex<double>)

  /* ********************************************************************* */
  /*   SuperLU solve interface                                             */
  /* ********************************************************************* */

  template<typename T>
  void SuperLU_solve(const gmm::csc_matrix<T> &csc_A, T *sol, T *rhs,
		     double& rcond_, int permc_spec) {
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     *   permc_spec = 3: use approximate minimum degree column ordering
     */
    typedef typename gmm::number_traits<T>::magnitude_type R;

    int m = int(mat_nrows(csc_A)), n = int(mat_ncols(csc_A));
    int nrhs = 1, info = 0, nz = int(nnz(csc_A));

    GMM_ASSERT1(nz != 0, "Cannot factor a matrix full of zeros!");
    GMM_ASSERT1(n == m, "Cannot factor a non-square matrix");

    if ((2 * nz / n) >= m)
      GMM_WARNING2("CAUTION : it seems that SuperLU has a problem"
		  " for nearly dense sparse matrices");

    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = NATURAL;
    options.PrintStat = NO;
    options.ConditionNumber = YES;
    switch (permc_spec) {
    case 1 : options.ColPerm = MMD_ATA; break;
    case 2 : options.ColPerm = MMD_AT_PLUS_A; break;
    case 3 : options.ColPerm = COLAMD; break;
    }
    SuperLUStat_t stat;
    StatInit(&stat);

    SuperMatrix SA, SL, SU, SB, SX; // SuperLU format.
    Create_CompCol_Matrix(&SA, m, n, nz, const_cast<T*>(&csc_A.pr[0]),
			  (int *)(&csc_A.ir[0]), (int *)(&csc_A.jc[0]));
    Create_Dense_Matrix(&SB, m, nrhs, &rhs[0], m);
    Create_Dense_Matrix(&SX, m, nrhs, &sol[0], m);

    memset(&SL,0,sizeof SL);
    memset(&SU,0,sizeof SU);

    std::vector<int> etree(n);
    char equed[] = "B";
    std::vector<R> Rscale(m),Cscale(n); // row scale factors
    std::vector<R> ferr(nrhs), berr(nrhs);
    R recip_pivot_gross, rcond;
    std::vector<int> perm_r(m), perm_c(n);

    SuperLU_gssvx(&options, &SA, &perm_c[0], &perm_r[0], 
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
    rcond_ = rcond;
    if (SB.Store) Destroy_SuperMatrix_Store(&SB);
    if (SX.Store) Destroy_SuperMatrix_Store(&SX);
    if (SA.Store) Destroy_SuperMatrix_Store(&SA);
    if (SL.Store) Destroy_SuperNode_Matrix(&SL);
    if (SU.Store) Destroy_CompCol_Matrix(&SU);
    StatFree(&stat);
    GMM_ASSERT1(info != -333333333, "SuperLU was cancelled."); // user interruption (for matlab interface)
    GMM_ASSERT1(info == 0, "SuperLU solve failed: info=" << info);
  }

  template void SuperLU_solve(const gmm::csc_matrix<float> &csc_A, float *sol, float *rhs, double& rcond_, int permc_spec);
  template void SuperLU_solve(const gmm::csc_matrix<double> &csc_A, double *sol, double *rhs, double& rcond_, int permc_spec);
  template void SuperLU_solve(const gmm::csc_matrix<std::complex<float> > &csc_A, std::complex<float> *sol, std::complex<float> *rhs, double& rcond_, int permc_spec);
  template void SuperLU_solve(const gmm::csc_matrix<std::complex<double> > &csc_A, std::complex<double> *sol, std::complex<double> *rhs, double& rcond_, int permc_spec);

  struct SuperLU_factor_impl_common {
    mutable SuperMatrix SA, SL, SB, SU, SX;
    mutable SuperLUStat_t stat;
    mutable superlu_options_t options;
    float memory_used;
    mutable bool is_init;
    mutable char equed;
    void free_supermatrix() {
      if (is_init) {
	if (SB.Store) Destroy_SuperMatrix_Store(&SB);
	if (SX.Store) Destroy_SuperMatrix_Store(&SX);
	if (SA.Store) Destroy_SuperMatrix_Store(&SA);
	if (SL.Store) Destroy_SuperNode_Matrix(&SL);
	if (SU.Store) Destroy_CompCol_Matrix(&SU);
      }
    }
    SuperLU_factor_impl_common() : is_init(false) {}
    virtual ~SuperLU_factor_impl_common() { free_supermatrix(); }
  };
  
  template <typename T> struct SuperLU_factor_impl : public SuperLU_factor_impl_common {
    typedef typename gmm::number_traits<T>::magnitude_type R;

    std::vector<int> etree, perm_r, perm_c;
    std::vector<R> Rscale, Cscale;
    std::vector<R> ferr, berr;
    std::vector<T> rhs;
    std::vector<T> sol;
    void build_with(const gmm::csc_matrix<T> &A, int permc_spec);
    void solve(int transp);
  };

  template <typename T>
  void SuperLU_factor_impl<T>::build_with(const gmm::csc_matrix<T> &A, int permc_spec) {
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     *   permc_spec = 3: use approximate minimum degree column ordering
     */
    free_supermatrix();
    int n = int(mat_nrows(A)), m = int(mat_ncols(A)), info = 0;

    rhs.resize(m); sol.resize(m);
    gmm::clear(rhs);
    int nz = int(nnz(A));

    GMM_ASSERT1(nz != 0, "Cannot factor a matrix full of zeros!");
    GMM_ASSERT1(n == m, "Cannot factor a non-square matrix");
    
    set_default_options(&options);
    options.ColPerm = NATURAL;
    options.PrintStat = NO;
    options.ConditionNumber = NO;
    switch (permc_spec) {
      case 1 : options.ColPerm = MMD_ATA; break;
      case 2 : options.ColPerm = MMD_AT_PLUS_A; break;
      case 3 : options.ColPerm = COLAMD; break;
    }
    StatInit(&stat);
    
    Create_CompCol_Matrix(&SA, m, n, nz, const_cast<T*>(&A.pr[0]),
                          (int *)(&A.ir[0]), (int *)(&A.jc[0]));
    Create_Dense_Matrix(&SB, m, 0, &rhs[0], m);
    Create_Dense_Matrix(&SX, m, 0, &sol[0], m);
    memset(&SL,0,sizeof SL);
    memset(&SU,0,sizeof SU);


    equed = 'B';
    Rscale.resize(m); Cscale.resize(n); etree.resize(n);
    ferr.resize(1); berr.resize(1);
    R recip_pivot_gross, rcond;
    perm_r.resize(m); perm_c.resize(n);
    memory_used = SuperLU_gssvx(&options, &SA, &perm_c[0], &perm_r[0], 
                                &etree[0] /* output */, &equed /* output        */, 
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
    
    Destroy_SuperMatrix_Store(&SB);
    Destroy_SuperMatrix_Store(&SX);
    Create_Dense_Matrix(&SB, m, 1, &rhs[0], m);
    Create_Dense_Matrix(&SX, m, 1, &sol[0], m);
    StatFree(&stat);
    
    GMM_ASSERT1(info != -333333333, "SuperLU was cancelled.");
    GMM_ASSERT1(info == 0, "SuperLU solve failed: info=" << info);
    is_init = true;
  }

  template <typename T> 
  void SuperLU_factor_impl<T>::solve(int transp) {
    options.Fact = FACTORED;
    options.IterRefine = NOREFINE;
    switch (transp) {
      case SuperLU_factor<T>::LU_NOTRANSP: options.Trans = NOTRANS; break;
      case SuperLU_factor<T>::LU_TRANSP: options.Trans = TRANS; break;
      case SuperLU_factor<T>::LU_CONJUGATED: options.Trans = CONJ; break;
      default: GMM_ASSERT1(false, "invalid value for transposition option");
    }
    StatInit(&stat);
    int info = 0;
    R recip_pivot_gross, rcond;
    SuperLU_gssvx(&options, &SA, &perm_c[0], &perm_r[0], 
                  &etree[0] /* output */, &equed /* output        */, 
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
    StatFree(&stat);
    GMM_ASSERT1(info == 0, "SuperLU solve failed: info=" << info);
  }
   
  template<typename T> void 
  SuperLU_factor<T>::build_with(const gmm::csc_matrix<T> &A, int permc_spec) {
    ((SuperLU_factor_impl<T>*)impl)->build_with(A,permc_spec);
  }

  template<typename T> void
  SuperLU_factor<T>::solve(int transp) const {
    ((SuperLU_factor_impl<T>*)impl)->solve(transp);
  }

  template<typename T> std::vector<T> &
  SuperLU_factor<T>::sol() const {
    return ((SuperLU_factor_impl<T>*)impl)->sol;
  }

  template<typename T> std::vector<T> &
  SuperLU_factor<T>::rhs() const {
    return ((SuperLU_factor_impl<T>*)impl)->rhs;
  }

  template<typename T> 
  SuperLU_factor<T>::SuperLU_factor() {
    impl = new SuperLU_factor_impl<T>();
  }

  template<typename T> 
  SuperLU_factor<T>::~SuperLU_factor() {
    delete impl;
  }

  template<typename T> 
  SuperLU_factor<T>::SuperLU_factor(const SuperLU_factor& other) {
    impl = new SuperLU_factor_impl<T>();
    GMM_ASSERT1(!(other.impl->is_init),
		"copy of initialized SuperLU_factor is forbidden");
    other.impl->is_init = false;
  }

  template<typename T> SuperLU_factor<T>&  
  SuperLU_factor<T>::operator=(const SuperLU_factor& other) {
    GMM_ASSERT1(!(other.impl->is_init) && !(impl->is_init),
		"assignment of initialized SuperLU_factor is forbidden");
    return *this;
  }

  template<typename T> float 
  SuperLU_factor<T>::memsize() const {
    return impl->memory_used;
  }

  /*  void force_instantiation() {
    SuperLU_factor<float> a;
    SuperLU_factor<double> b;
    SuperLU_factor<std::complex<float> > c;
    SuperLU_factor<std::complex<double> > d;
    //a = 0; b = 0; c = 0; d = 0;
  } 
  */ 
}

template class gmm::SuperLU_factor<float>;
template class gmm::SuperLU_factor<double>;
template class gmm::SuperLU_factor<std::complex<float> >;
template class gmm::SuperLU_factor<std::complex<double> >;

static int (*superlu_callback)();

/* this one is called from superlu (see dcolumn_bmod) */
extern "C" int handle_getfem_callback() {
  if (superlu_callback) return superlu_callback();
  else return 0;
}

extern "C" void set_superlu_callback(int (*cb)()) {
  superlu_callback = cb;
}
