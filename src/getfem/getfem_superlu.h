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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file getfem_superlu.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date August 2004.
   @brief SuperLU interface for getfem
   
   We do not use gmm_superlu_interface.h for a good reason. This file
   does not include any of the superlu headers, hence when getfem is
   installed, it does not need to install the superlu headers.
*/

#ifndef GETFEM_SUPERLU
#define GETFEM_SUPERLU
#ifndef GMM_USES_SUPERLU
#define GMM_USES_SUPERLU
#endif
#include "getfem_config.h"
#include "gmm/gmm_kernel.h"

namespace gmm {

  template<typename T>
  void SuperLU_solve(const gmm::csc_matrix<T> &A, T *X_, T *B, double& rcond_, int permc_spec = 3);
  /** solve a sparse linear system AX=B (float, double, complex<float>
      or complex<double>) via SuperLU.

      @param A the matrix (a copy is made if A is not a gmm::csc_matrix)
      @param X the solution.
      @param B the right hand side.
      @param rcond_ contains on output an estimate of the condition number of A.
      @param permc_spec specify the kind of renumbering than SuperLU should do.
  */
  template<typename MAT, typename V1, typename V2>
  void SuperLU_solve(const MAT &A, const V1& X, const V2& B, double& rcond_, int permc_spec = 3) {
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    
    int m = int(mat_nrows(A)), n = int(mat_ncols(A));
    gmm::csc_matrix<T> csc_A(m,n); 
    gmm::copy(A,csc_A);
    std::vector<T> rhs(m), sol(m);
    gmm::copy(B, rhs);
    
    SuperLU_solve(csc_A, &sol[0], &rhs[0], rcond_, permc_spec);
    gmm::copy(sol, const_cast<V1 &>(X));
  }
  
  struct SuperLU_factor_impl_common;

  /** Factorization of a sparse matrix with SuperLU.
      
      This class can be used as a preconditioner for gmm iterative solvers.
  */
  template <class T> class SuperLU_factor {
    SuperLU_factor_impl_common *impl;
  public :
    enum { LU_NOTRANSP, LU_TRANSP, LU_CONJUGATED };

    /** Do the factorization of the supplied sparse matrix. */
    template <class MAT> void build_with(const MAT &A,  int permc_spec = 3) {
      int m = int(mat_nrows(A)), n = int(mat_ncols(A));
      gmm::csc_matrix<T> csc_A(m,n); 
      gmm::copy(A,csc_A);
      build_with(csc_A, permc_spec);
    }
    void build_with(const gmm::csc_matrix<T> &A, int permc_spec = 3);
    template <typename VECTX, typename VECTB> 
    /** After factorization, do the triangular solves.
       transp = LU_NOTRANSP   -> solves Ax = B
       transp = LU_TRANSP     -> solves A'x = B
       transp = LU_CONJUGATED -> solves conj(A)X = B 
    */
    void solve(const VECTX &X, const VECTB &B, int transp=LU_NOTRANSP) const {
      gmm::copy(B, rhs());
      solve(transp);
      gmm::copy(sol(),const_cast<VECTX &>(X));
    }
    void solve(int transp=LU_NOTRANSP) const;
    std::vector<T> &sol() const;
    std::vector<T> &rhs() const;
    SuperLU_factor();
    ~SuperLU_factor();
    float memsize() const;
    SuperLU_factor(const SuperLU_factor& other);
    SuperLU_factor& operator=(const SuperLU_factor& other);
  };

  template <typename T, typename V1, typename V2> inline
  void mult(const SuperLU_factor<T>& P, const V1 &v1, const V2 &v2) {
    P.solve(v2,v1);
  }

  template <typename T, typename V1, typename V2> inline
  void transposed_mult(const SuperLU_factor<T>& P,const V1 &v1,const V2 &v2) {
    P.solve(v2, v1, SuperLU_factor<T>::LU_TRANSP);
  }
}

extern "C" void set_superlu_callback(int (*cb)());

#endif
