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
/* This file is a part of GETFEM++                                         */
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

  template<class T> struct entry1 {
    size_type  index;
    T value;
    entry1(void) {}
    entry1(size_type i) : index(i) {}
    bool operator <(const entry1& e) const { return index < e.index; }
    bool operator ==(const entry1& e) const { return index == e.index; }
    bool operator !=(const entry1& e) const { return index != e.index; }
  };

  template<class T> struct entry1_value_less {
    inline bool operator()(const entry1<T>& a, const entry1<T>& b) const
    { return (dal::abs(a.value) < dal::abs(b.value)); }
  };
  
  template <class Matrix>
  class ilut_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef row_matrix<rsvector<value_type> > LU_Matrix;
    typedef std::vector<entry1<value_type> > entry_vec;
    typedef std::vector<value_type> value_type_vec;

    bool invert;
    LU_Matrix L, U;

  protected:
    int K, dropped;
    double eps;    

    void do_ilut(const Matrix&, row_major);
    void do_ilut(const Matrix&, col_major);

  public:
    ilut_precond(const Matrix& A, int k_, double eps_) 
      : invert(false), L(mat_nrows(A), mat_ncols(A)),
	U(mat_nrows(A), mat_ncols(A)), K(k_), dropped(0), eps(eps_) {
   
      if (!is_sparse(A))
	DAL_THROW(failure_error,
		  "Matrix should be sparse for incomplete ilu");
      do_ilut(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
  };

  template<class Matrix> 
  void ilut_precond<Matrix>::do_ilut(const Matrix& A, row_major) {
    value_type_vec indiag(mat_nrows(A));
    
    entry_vec* w = new entry_vec;
    entry_vec* wswap = new entry_vec;

    for (size_type i = 0; i < mat_nrows(A); ++i) {

      typedef typename linalg_traits<Matrix>::const_sub_row_type row_type;
      row_type row = mat_const_row(A, i);
      typename linalg_traits<row_type>::const_iterator
	it = vect_const_begin(row), ite = vect_const_end(row);
      size_type ninrow = 0;
      for (; it != ite; ++it) ++ninrow;
      it = vect_const_begin(row);
      
      w->resize(ninrow);
      double norm_row = 0.;
      size_type nL = 0, nU = 0, inrow = 0;
      for (; it != ite; ++it) {
	(*w)[inrow].value = *it;
	(*w)[inrow].index = it.index();
	inrow++;
	norm_row += dal::sqr(dal::abs(*it));
	if (i > it.index()) nL++;
      }
      std::sort(w->begin(), w->end());
      
      norm_row = sqrt(norm_row);
      nU = ninrow - nL - 1;
      norm_row = 1./norm_row;
      
      typename entry_vec::iterator wk= w->begin(), wkend = w->end();
      entry_vec tmp_w;
      size_type krow = 0;
      while( wk != wkend ) {
	size_type k = (*wk).index;
	if ( k >= i ) break;
	value_type tmp = (*wk).value;
	tmp = tmp * indiag[k];  
	if ( dal::abs(tmp) * norm_row < eps ) {
	  w->erase(wk); 
	  dropped++;
	} else {
	  (*wk).value = tmp;

	  typedef typename linalg_traits<LU_Matrix>::sub_row_type LU_row_type;
	  LU_row_type U_k = mat_row(U, k);
	  typename linalg_traits<LU_row_type>::const_iterator
	    U_kj = vect_const_begin(U_k), U_kend = vect_const_end(U_k);
	  size_type U_knnz = 0;
	  for (; U_kj != U_kend; ++U_kj) ++U_knnz;
	  U_kj = vect_const_begin(U_k);
	  ++U_kj;
	  tmp_w.resize((U_knnz == 0) ? 0 : U_knnz - 1);
	  
	  size_type irow = 0;
	  for (; U_kj != U_kend; ++U_kj) {
	    tmp_w[irow].value = -tmp * (*U_kj);
	    tmp_w[irow].index = U_kj.index();
	    irow++;
	  }

	  wswap->resize(w->size() + tmp_w.size());
	  typename entry_vec::iterator wj = w->begin(), tmp_wj = tmp_w.begin();
	  size_type j = 0;
	  for(;;) {
	    if (wj == w->end()) {
	      for (; tmp_wj != tmp_w.end(); ++tmp_wj, ++j)
		(*wswap)[j] = *tmp_wj;
	      break;
	    }
	    
	    if (tmp_wj == tmp_w.end()) {
	      for (; wj != w->end(); ++wj, ++j) (*wswap)[j] = *wj;
	      break;
	    }
	    
	    if (*wj == *tmp_wj) {
	      (*wswap)[j] = *wj; 
	      (*wswap)[j].value += (*tmp_wj).value;
	      ++tmp_wj; ++wj; ++j;
	      continue;
	    }
	    
	    if ( *wj < *tmp_wj )
	      { (*wswap)[j] = *wj; ++wj; ++j; continue; }
	    
	    (*wswap)[j] = *tmp_wj; 
	    ++tmp_wj;  ++j;
	  }
	  
	  wswap->resize(j);
	  entry_vec* tswap = wswap;
	  wswap = w;
	  w = tswap;
	  
	  krow++; 
	}
	wk = w->begin()+krow;
	wkend = w->end();
      }
      
      if (i) {
	typename entry_vec::iterator wi = w->begin(), wend = w->end();
	for (; wi != wend; ) {
	  if ((*wi).index != i && dal::abs((*wi).value)*norm_row < eps)
	    {  w->erase(wi); dropped++; }
	  else ++wi;
	}
	
	typename entry_vec::iterator diag = 
	  std::find(w->begin(), w->end(), entry1<value_type>(i));
	size_type m = diag-w->begin();
	std::make_heap(w->begin(), diag, entry1_value_less<value_type>());
	size_type jmax = std::min(nL+K, m);
	dropped += m-jmax;
	for (size_type j=m; j>m-jmax; --j) {   
	  typename entry_vec::iterator first = w->begin();
	  L(i, (*first).index) = (*first).value;
	  std::pop_heap(first, first+j, entry1_value_less<value_type>());
	}
	
	U(i, i) = (*diag).value;
	
	m = w->end() - diag - 1;
	jmax = std::min(nU+K, m);
	if (m > 0) {
	  dropped += m-jmax;
	  std::make_heap(diag+1, w->end(), entry1_value_less<value_type>());
	  
	  for (size_type j=m; j>m-jmax; --j) {   
	    typename entry_vec::iterator first = diag+1;
	    U(i, (*first).index) = (*first).value;
	    std::pop_heap(first, first+j, entry1_value_less<value_type>());
	  }
	}
      }
      else { 
	for (typename entry_vec::iterator wj= w->begin(); wj != w->end(); ++wj)
	  U(i, (*wj).index) = (*wj).value;
      }
      indiag[i] = 1.0 / U(i,i);
    }

    delete w;
    delete wswap;
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

