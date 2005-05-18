// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (Getfem)
// File    : getfem_superlu.h : .
//           
// Date    : August 2004.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Julien Pommier
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#ifndef GETFEM_SUPERLU
#define GETFEM_SUPERLU
#ifndef GMM_USES_SUPERLU
#define GMM_USES_SUPERLU
#endif
#include <gmm_kernel.h>
#include <getfem_config.h>

namespace gmm {

  template<typename T>
  void SuperLU_solve(const gmm::csc_matrix<T> &A, T *X_, T *B, double& rcond_, int permc_spec = 3);
  
  template<typename MAT, typename V1, typename V2>
  void SuperLU_solve(const MAT &A, const V1& X, const V2& B, double& rcond_, int permc_spec = 3) {
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    
    int m = mat_nrows(A), n = mat_ncols(A);
    gmm::csc_matrix<T> csc_A(m,n); 
    gmm::copy(A,csc_A);
    std::vector<T> rhs(m), sol(m);
    gmm::copy(B, rhs);
    
    SuperLU_solve(csc_A, &sol[0], &rhs[0], rcond_, permc_spec);
    gmm::copy(sol, const_cast<V1 &>(X));
  }
  
  class SuperLU_factor_impl_common;

  template <class T> class SuperLU_factor {
    SuperLU_factor_impl_common *impl;
  public :
    enum { LU_NOTRANSP, LU_TRANSP, LU_CONJUGATED };
    template <class MAT> void build_with(const MAT &A,  int permc_spec = 3) {
      int m = mat_nrows(A), n = mat_ncols(A);
      gmm::csc_matrix<T> csc_A(m,n); 
      gmm::copy(A,csc_A);
      build_with(csc_A, permc_spec);
    }
    void build_with(const gmm::csc_matrix<T> &A, int permc_spec = 3);
    template <typename VECTX, typename VECTB> 
    /* transp = LU_NOTRANSP   -> solves Ax = B
       transp = LU_TRANSP     -> solves A'x = B
       transp = LU_CONJUGATED -> solves conj(A)X = B */
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
  private:
    /*    SuperLU_factor(const SuperLU_factor&);
    SuperLU_factor& operator=(const SuperLU_factor& other);
    */
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

#endif
