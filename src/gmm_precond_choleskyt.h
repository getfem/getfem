/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_choleskyt.h : Incomplete LDLT factorisation      */
/*                                      with fill-in and threshold.        */
/*					                                   */
/*     									   */
/* Date : June 30, 2003.                                                   */
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
#ifndef GMM_PRECOND_CHOLESKYT_H
#define GMM_PRECOND_CHOLESKYT_H

#include <gmm_solvers.h>
#include <gmm_tri_solve.h>
#include <gmm_interface.h>

namespace gmm {

#define ildltt_precond choleskyt_precond

  template <typename Matrix>
  class choleskyt_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;

    row_matrix<svector> U;
    std::vector<value_type> indiag;

  protected:
    int K;
    double eps;    

    template<typename M> void do_choleskyt(const M&, row_major, int = 0);
    void do_choleskyt(const Matrix&, col_major);

  public:
    choleskyt_precond(const Matrix& A, int k_, double eps_) 
      : U(mat_nrows(A),mat_ncols(A)),
	indiag(std::min(mat_nrows(A), mat_ncols(A))),
	K(k_), eps(eps_)
    {
      if (!is_sparse(A))
	DAL_THROW(failure_error,
		  "Matrix should be sparse for incomplete cholesky");
      do_choleskyt(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    choleskyt_precond(void) {}
  };

  template<typename Matrix> template<typename M> 
  void choleskyt_precond<Matrix>::do_choleskyt(const M& A,row_major,int _try) {
    svector w(mat_ncols(A));
    value_type tmp;

    gmm::clear(U);
    for (size_type i = 0; i < mat_nrows(A); ++i) {
      gmm::copy(mat_const_row(A, i), w);
      double norm_row = gmm::vect_norm2(w);

      size_type nU = 0;
      typename linalg_traits<svector>::iterator it = vect_begin(w);
      for (; it != vect_end(w); ++it) if (i < it.index()) nU++;

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	tmp = (wk->e);
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
	  do_choleskyt(A, row_major(), ++_try);
	  return;
	}
      }

      indiag[i] = 1.0 / tmp;
      gmm::clean(w, eps * norm_row);
      gmm::scale(w, indiag[i]);
      std::sort(w.begin(), w.end(), _elt_rsvector_value_less<value_type>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      for (size_type nnu = 0; wit != wite; ++wit)
	if (wit->c > i) { if (nnu < nU+K) U(i, wit->c) = wit->e; ++nnu; }
    }
  }

  template<typename Matrix> 
  void choleskyt_precond<Matrix>::do_choleskyt(const Matrix& A, col_major)
  { do_choleskyt(gmm::transposed(A), row_major()); }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    gmm::lower_tri_solve(gmm::transposed(P.U), v2, true);
    for (size_type i = 0; i < P.indiag.size(); ++i) v2[i] *= P.indiag[i];
    gmm::upper_tri_solve(P.U, v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const choleskyt_precond<Matrix>& P,const V1 &v1, V2 &v2)
  { mult(P, v1, v2); }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    gmm::lower_tri_solve(gmm::transposed(P.U), v2, true);
    for (size_type i = 0; i < P.indiag.size(); ++i) v2[i] *= P.indiag[i];
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.U, v2, true); }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const choleskyt_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    gmm::upper_tri_solve(P.U, v2, true);
    for (size_type i = 0; i < P.indiag.size(); ++i) v2[i] *= P.indiag[i];
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const choleskyt_precond<Matrix>& P, const V1 &v1,
			     V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::transposed(P.U), v2, true); }

}

#endif 

