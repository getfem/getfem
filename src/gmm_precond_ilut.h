/* -*- c++ -*- (enables emacs c++ mode)                                    */
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_ilut.h : modified version from I.T.L.            */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
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
#ifndef GMM_PRECOND_ILUT_H
#define GMM_PRECOND_ILUT_H

//: ILUT:  Incomplete LU with threshold and K fill-in Preconditioner.
//  The algorithm of ILUT(A, 0, 1.0e-6) is slower than ILU(A). If No fill-in 
//  is arrowed, you can use ILU instead of ILUT.
//
// Notes: The idea under a concrete Preconditioner such 
//        as ilut is to create a Preconditioner
//        object to use in iterative methods. 
//

/*
  Performane comparing for SSOR, ILU and ILUT based on sherman 5 matrix 
  in Harwell-Boeing collection on Sun Ultra 30 UPA/PCI (UltraSPARC-II 296MHz)
  Preconditioner & Factorization time  &  Number of Iteration \\ \hline
  SSOR        &   0.010577  & 41 \\
  ILU         &   0.019336  & 32 \\
  ILUT with 0 fill-in and threshold of 1.0e-6 & 0.343612 &  23 \\
  ILUT with 5 fill-in and threshold of 1.0e-6 & 0.343612 &  18 \\ \hline
*/

#include <gmm_solvers.h>
#include <gmm_tri_solve.h>
#include <gmm_interface.h>

namespace gmm {

  template <class Matrix>
  class ilut_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;
    typedef row_matrix<svector> LU_Matrix;

    bool invert;
    LU_Matrix L, U;

  protected:
    int K;
    double eps;    

    template<class M> void do_ilut(const M&, row_major, int = 0);
    void do_ilut(const Matrix&, col_major);

  public:
    ilut_precond(const Matrix& A, int k_, double eps_) 
      : invert(false), L(mat_nrows(A), mat_ncols(A)),
	U(mat_nrows(A), mat_ncols(A)), K(k_), eps(eps_) {
   
      if (!is_sparse(A))
	DAL_THROW(failure_error,
		  "Matrix should be sparse for incomplete ilu");
      do_ilut(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    
    ilut_precond(void) {}

  };

  template<class Matrix> template<class M> 
  void ilut_precond<Matrix>::do_ilut(const M& A, row_major, int _try) {
    std::vector<value_type> indiag(mat_nrows(A));
    svector w(mat_ncols(A));
    value_type tmp;
    gmm::clear(U); gmm::clear(L);

    for (size_type i = 0; i < mat_nrows(A); ++i) {
      gmm::copy(mat_const_row(A, i), w);
      double norm_row = gmm::vect_norm2(w);

      size_type nL = 0;
      typename linalg_traits<svector>::iterator it = vect_begin(w);
      for (; it != vect_end(w); ++it) if (i > it.index()) nL++;
      size_type nU = w.nb_stored() - nL - 1;

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	tmp = (wk->e) * indiag[k];
	if (dal::abs(tmp) < eps * norm_row) { w.sup(k); --krow; } 
	else { wk->e += tmp; gmm::add(scaled(mat_row(U, k), -tmp), w); }
      }
      
      if ((tmp = w[i]) == value_type(0)) {
	DAL_WARNING(2, "pivot " << i << " is zero");
	if (_try > 10)
	  tmp = 1.0;
	else {
	  ++K; eps /= value_type(2);
	  DAL_WARNING(2, "trying with " << K
		      << " additional elements and threshold " << eps);
	  do_ilut(A, row_major(), ++_try);
	  return;
	}
      }

      indiag[i] = 1.0 / tmp;
      U(i,i) = tmp; gmm::clean(w, eps * norm_row); w[i] = 0.0;
      std::sort(w.begin(), w.end(), _elt_rsvector_value_less<value_type>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      size_type nnl = 0, nnu = 0;
      for (; wit != wite; ++wit)
	if (wit->c < i) { if (nnl < nL+K) L(i, wit->c) = wit->e; ++nnl; }
	else            { if (nnu < nU+K) U(i, wit->c) = wit->e; ++nnu; }
    }
  }

  template<class Matrix> 
  void ilut_precond<Matrix>::do_ilut(const Matrix& A, col_major) {
    do_ilut(gmm::transposed(A), row_major());
    invert = true;
  }

  template <class Matrix, class V1, class V2> inline
  void mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    if (P.invert) {
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      gmm::lower_tri_solve(P.L, v2, true);
      gmm::upper_tri_solve(P.U, v2, false);
    }
  }

  template <class Matrix, class V1, class V2> inline
  void transposed_mult(const ilut_precond<Matrix>& P,const V1 &v1,V2 &v2) {
    gmm::copy(v1, v2);
    if (P.invert) {
      gmm::lower_tri_solve(P.L, v2, true);
      gmm::upper_tri_solve(P.U, v2, false);
    }
    else {
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }

  template <class Matrix, class V1, class V2> inline
  void left_mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    else gmm::lower_tri_solve(P.L, v2, true);
  }

  template <class Matrix, class V1, class V2> inline
  void right_mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    else gmm::upper_tri_solve(P.U, v2, false);
  }

  template <class Matrix, class V1, class V2> inline
  void transposed_left_mult(const ilut_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(P.U, v2, false);
    else gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
  }

  template <class Matrix, class V1, class V2> inline
  void transposed_right_mult(const ilut_precond<Matrix>& P, const V1 &v1,
			     V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(P.L, v2, true);
    else gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
  }

}

#endif 

