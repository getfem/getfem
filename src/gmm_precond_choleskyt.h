/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_choleskyt.h : Incomplete Cholesky with fill-in   */
/*                                      and threshlod.                     */
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
#ifndef GMM_PRECOND_CHOLESKYT_H
#define GMM_PRECOND_CHOLESKYT_H

#include <gmm_solvers.h>
#include <gmm_tri_solve.h>
#include <gmm_interface.h>

namespace gmm {

  template <class Matrix>
  class choleskyt_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;

    row_matrix<svector> U;

  protected:
    int K;
    double eps;    

    template<class M> void do_choleskyt(const M&, row_major);
    void do_choleskyt(const Matrix&, col_major);

  public:
    choleskyt_precond(const Matrix& A, int k_, double eps_) 
      : U(mat_nrows(A), mat_ncols(A)), K(k_), eps(eps_) {
      if (!is_sparse(A))
	DAL_THROW(failure_error,
		  "Matrix should be sparse for incomplete ilu");
      do_choleskyt(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
  };

  template<class Matrix> template<class M> 
  void choleskyt_precond<Matrix>::do_choleskyt(const M& A, row_major) {
    std::vector<value_type> indiag(mat_nrows(A));
    svector w(mat_ncols(A));
    value_type a;

    for (size_type i = 0; i < mat_nrows(A); ++i) {
      gmm::copy(mat_const_row(A, i), w);
      double norm_row = gmm::vect_norm2(w);

      size_type nU = 0;
      typename linalg_traits<svector>::iterator it = vect_begin(w);
      for (; it != vect_end(w); ++it) if (i < it.index()) nU++;

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	a = (wk->e) / indiag[k];
	if (dal::abs(a) < eps * norm_row) { w.sup(k); --krow; } 
	else { wk->e += a; gmm::add(scaled(mat_row(U, k), -a), w); }
      }

      if ((a = w[i]) <= 0) DAL_THROW(failure_error, "negative value found");
      a = sqrt(a);

      U(i,i) = a; gmm::clean(w, eps * norm_row); gmm::scale(w, 1.0 / a);
      std::sort(w.begin(), w.end(), _elt_rsvector_value_less<value_type>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      size_type nnu = 0;
      for (; wit != wite; ++wit)
	if (wit->c > i) { if (nnu < nU+K) U(i, wit->c) = wit->e; ++nnu; }
      
      indiag[i] = U(i,i);
    }
  }

  template<class Matrix> 
  void choleskyt_precond<Matrix>::do_choleskyt(const Matrix& A, col_major)
  { do_choleskyt(gmm::transposed(A), row_major()); }

  template <class Matrix, class V1, class V2> inline
  void mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    gmm::lower_tri_solve(gmm::transposed(P.U), v2);
    gmm::upper_tri_solve(P.U, v2);
  }

  template <class Matrix, class V1, class V2> inline
  void transposed_mult(const choleskyt_precond<Matrix>& P,const V1 &v1,V2 &v2)
  { mult(P, v1, v2); }

  template <class Matrix, class V1, class V2> inline
  void left_mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::transposed(P.U), v2); }

  template <class Matrix, class V1, class V2> inline
  void right_mult(const choleskyt_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.U, v2); }

  template <class Matrix, class V1, class V2> inline
  void transposed_left_mult(const choleskyt_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) 
  { copy(v1, v2); gmm::upper_tri_solve(P.U, v2); }

  template <class Matrix, class V1, class V2> inline
  void transposed_right_mult(const choleskyt_precond<Matrix>& P, const V1 &v1,
			     V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::transposed(P.U), v2); }

}

#endif 

