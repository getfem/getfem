/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_idgmres.h : implicitly deflated GMRES             */
/*     									   */
/* Date : October 6, 2003.                                                 */
/* Authors :  Caroline Lecalvez, Caroline.Lecalvez@gmm.insa-tlse.fr        */
/*            Yves Renard, Yves.Renard@gmm.insa-tlse.fr                    */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Caroline Lecalvez, Yves Renard.                     */
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

#ifndef GMM_IDGMRES_H
#define GMM_IDGMRES_H

namespace gmm {

  // Implicitly restarted and deflated Generalized Minimum Residual
  //
  //   See: C. Le Calvez, B. Molina, Implicitly restarted and deflated
  //        FOM and GMRES, numerical applied mathematics,
  //        (30) 2-3 (1999) pp191-212.
  //
  //
  // A : Real or complex unsymmetric matrix.
  // x : initial guess vector and final result.
  // b : right hand side
  // M : preconditionner
  // restart : size of the subspace between two restarts
  // p : number of converged ritz values seeked
  // k : size of the remaining Krylov subspace when the p ritz values
  //      have not yet converged 0 <= p <= k < restart.
  // tol_vp : tolerance on the ritz values.

  template < class Mat, class Vec, class VecB, class Precond, class Basis >
  void idgmres(const Mat &A, Vec &x, const VecB &b, const Precond &M,
	     int restart, int p, int k, double tol_vp,
	     iteration &outer, Basis& KS) {

    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    R a, beta;

    std::vector<T> w(vect_size(x)), r(vect_size(x)), u(vect_size(x));
    std::vector<T> c_rot(restart+1), s_rot(restart+1), s(restart+1);
    gmm::dense_matrix<T> H(restart+1, restart), Hess(restart+1, restart),
      Hess1(restart+1, restart);
    std::vector<T> em1(restart+1);

    gmm::clear(H);
    gmm::clear(em1);
    em1[restart+1] = T(1);

    outer.set_rhsnorm(gmm::vect_norm2(b));
    if (outer.get_rhsnorm() == 0.0) { clear(x); return; }
    
    mult(A, scaled(x, -1.0), b, w);
    mult(M, w, r);
    beta = gmm::vect_norm2(r);

    iteration inner = outer;
    inner.reduce_noisy();
    inner.set_maxiter(restart);
    inner.set_name("GMRes inner iter");
    
    while (! outer.finished(beta)) {
      
      gmm::copy(gmm::scaled(r, 1.0/beta), KS[0]);
      gmm::clear(s);
      s[0] = beta;
      
      size_type i = 0; inner.init();
      
      do {
	gmm::mult(A, KS[i], u);
	gmm::mult(M, u, KS[i+1]);
	orthogonalize_with_refinment(KS, mat_col(H, i), i);
	H(i+1, i) = a = gmm::vect_norm2(KS[i+1]);
	gmm::scale(KS[i+1], R(1) / a);

	gmm::copy(mat_col(H, i), mat_col(Hess, i));
	gmm::copy(mat_col(H, i), mat_col(Hess1, i));
	

	for (size_type k = 0; k < i; k++)
	  Apply_Givens_rotation_left(H(k,i), H(k+1,i), c_rot[k], s_rot[k]);
	
	Givens_rotation(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	Apply_Givens_rotation_left(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	H(i+1, i) = T(0); // utile ?
	Apply_Givens_rotation_left(s[i], s[i+1], c_rot[i], s_rot[i]);
	
	++inner, ++outer, ++i;
      } while (! inner.finished(dal::abs(s[i])));

      gmm::upper_tri_solve(H, s, i, false);
      gmm::combine(KS, s, x, i);
      gmm::mult(A, gmm::scaled(x, T(-1)), b, w);
      gmm::mult(M, w, r);
      beta = gmm::vect_norm2(r);
    }
  }


  template < class Mat, class Vec, class VecB, class Precond >
  void idgmres(const Mat &A, Vec &x, const VecB &b,
	     const Precond &M, int restart, iteration& outer) {
    typedef typename linalg_traits<Mat>::value_type T;
    modified_gram_schmidt<T> orth(restart, vect_size(x));
    gmres(A, x, b, M, restart, outer, orth); 
  }

}

#endif
