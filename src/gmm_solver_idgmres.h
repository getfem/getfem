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
  // m : size of the subspace between two restarts
  // p : number of converged ritz values seeked
  // k : size of the remaining Krylov subspace when the p ritz values
  //      have not yet converged 0 <= p <= k < m.
  // tol_vp : tolerance on the ritz values.


  template <class T> compare_vp {
    bool operator()(const std::pair<T, size_type> &a,
		    const std::pair<T, size_type> &b) const
    { return (dal::abs(a.first) > dal::abs(b.first)); }
  }


  template < class Mat, class Vec, class VecB, class Precond, class Basis >
  void idgmres(const Mat &A, Vec &x, const VecB &b, const Precond &M,
	     int m, int p, int k, double tol_vp,
	     iteration &outer, Basis& KS) {

    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    R a, beta;

    std::vector<T> w(vect_size(x)), r(vect_size(x)), u(vect_size(x));
    std::vector<T> c_rot(m+1), s_rot(m+1), s(m+1);
    std::vector<T> y(m+1), ztest(m+1), gam(m+1);
    std::vector<T> gamma(m+1);
    gmm::dense_matrix<T> H(m+1, m), Hess(m+1, m),
      Hess1(m+1, m), W(vect_size(x), m+1);

    gmm::clear(H);
   
    int tb_def = 0, tb_deb = 1;

    outer.set_rhsnorm(gmm::vect_norm2(b));
    if (outer.get_rhsnorm() == 0.0) { clear(x); return; }
    
    mult(A, scaled(x, -1.0), b, w);
    mult(M, w, r);
    beta = gmm::vect_norm2(r);

    iteration inner = outer;
    inner.reduce_noisy();
    inner.set_maxiter(m);
    inner.set_name("GMRes inner iter");
    
    while (! outer.finished(beta)) {
      
      gmm::copy(gmm::scaled(r, 1.0/beta), KS[0]);
      gmm::clear(s);
      s[0] = beta;
      gmm::copy(s, gamma);

      inner.set_maxiter(m - tb_deb + 1);
      size_type i = tb_deb - 1; inner.init();
      
      do {
	mult(A, KS[i], u);
	mult(M, u, KS[i+1]);
	orthogonalize_with_refinment(KS, mat_col(H, i), i);
	H(i+1, i) = a = gmm::vect_norm2(KS[i+1]);
	gmm::scale(KS[i+1], R(1) / a);

	gmm::copy(mat_col(H, i), mat_col(Hess, i));
	gmm::copy(mat_col(H, i), mat_col(Hess1, i));
	

	for (size_type l = 0; l < i; ++l)
	  Apply_Givens_rotation_left(H(l,i), H(l+1,i), c_rot[l], s_rot[l]);
	
	Givens_rotation(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	Apply_Givens_rotation_left(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	H(i+1, i) = T(0); 
	Apply_Givens_rotation_left(s[i], s[i+1], c_rot[i], s_rot[i]);
	
	++inner, ++outer, ++i;
      } while (! inner.finished(dal::abs(s[i])));

      gmm::copy(s, y);
      upper_tri_solve(H, y, i, false);
      combine(KS, y, x, i);
      mult(A, gmm::scaled(x, T(-1)), b, w);
      mult(M, w, r);
      beta = gmm::vect_norm2(r);

      gmm::clear(gam); gam[m] = s[i];
      for (size_type l = m; l > 0; --l)
	Apply_Givens_rotation_left(gam[l-1], gam[l], c_rot[l-1], s_rot[l-1]);

      mult(KS.mat(), gam, r);
      beta = gmm::vect_norm2(r);
      
      mult(Hess, scaled(y, T(-1)), gamma, ztest);
      // En fait, d'après Caroline qui s'y connait ztest et gam devrait
      // être confondus

      if (tb_def < p) {
	ns = dal::sgn(ztest[m]);
	gmm::copy(V, W); gmm::copy(scaled(r, ns / beta), mat_col(W, m));
	
	// Computation of the oblique matrix
	sub_interval SUBI(0, m);
	add(scaled(sub_vector(ztest, SUBI), -Hess1(m, m-1) / ztest[m]),
	    sub_vector(mat_col(Hess1, m-1), SUBI));
	Hess1(m, m-1) *= ns * beta / ztest[m]; 

	/********************************************************************/
	/*  Locking                                                         */
	/********************************************************************/

	// Computation of the Ritz eigenpairs.

	dense_matrix<T> evect(m-tb_def, m-tb_def);
	std::vector<std::complex<R> > eval(m);
	std::vector<R> ritznew(m, T(-1));
	
	// dense_matrix<T> evect_lock(tb_def, tb_def);
	

	sub_interval SUB1(tb_def, m-tb_def);
	implicit_qr_algorithm(sub_matrix(Hess1, SUB1, SUB1),
			      sub_vector(eval, SUB1), evect);
	sub_interval SUB2(0, tb_def);
	implicit_qr_algorithm(sub_matrix(Hess1, SUB2, SUB2),
			      sub_vector(eval, SUB2), /* evect_lock */);

	for (size_type l = tb_def; l < m; ++l)
	  ritznew[l] = dal::abs(evect(m-tb_def-1, l-tb_def) * Hess1(m, m-1));
	
	std::vector< std::pair<T, size_type> > eval_sort(m);
	for (size_type l = 0; l < m; ++l)
	  eval_sort[l] = std::pair<T, size_type>(eval[l], l);
	std::sort(eval_sort.begin(), eval_sort.end(), compare_vp());

	std::vector<bool> kept(m, false);
	std::fill(kept.begin(), kept.begin()+tb_def, true);
	
	
	//	Which are the eigenvalues that converged ?
	//
	//	nb_want is the number of eigenvalues of 
	//	Hess(tb_def+1:n,tb_def+1:n) that converged and are WANTED
	//
	//	nb_unwant is the number of eigenvalues of 
	//	Hess(tb_def+1:n,tb_def+1:n) that converged and are UNWANTED
	//
	//	nb_nolong is the number of eigenvalues of 
	//	Hess(1:tb_def,1:tb_def) that are NO LONGER WANTED. 
	//
	//	tb_deftot is the number of the deflated eigenvalues
	//	that is tb_def + nb_want + nb_unwant
	//
	//	tb_defwant is the number of the wanted deflated eigenvalues
	//	that is tb_def + nb_want - nb_nolong

	dense_matrix<T> YB(m-tb_def, m-tb_def);
	std::vector<T> pure(m-tb_def, T(0));
	
	size_type nb_want = 0, nb_unwant = 0, nb_nolong = 0, j;

	for (j = 0, ind = 0; j < m-p; ++j) {
	  if (ritznew[eval_sort[j].second] == R(-1)) {
	    if (std::imag(eval_sort[j].first) != R(0)) {
	      nb_nolong += 2; ++j;
	    } 
	    else nb_nolong++;
	  }
	  else {
	    if (ritznew[eval_sort[j].second]
		< tol_vp * dal::abs(eval_sort[j].first)) {
	      
	      for (size_type l = 0, l < m-tb_def; ++l)
		YB(l, ind) = std::real(evect(l, eval_sort[j].second));
	      kept[eval_sort[j].second] = true;
	      ++j; ++nb_unwant; ind++;
	      
	      if (std::imag(eval_sort[j].first) != R(0)) {
		for (size_type l = 0, l < m-tb_def; ++l)
		  YB(l, ind) = std::imag(evect(l, eval_sort[j].second));
		pure[ind-1] = T(1);
		pure[ind] = T(2);
		
		kept[eval_sort[j].second] = true;
		
		nb_unwant++;
		++ind;
	      }
	    }
	  }
	}


	for (; j < m; ++j) {
	  if (ritznew[eval_sort[j].second] != R(-1)) {

	      for (size_type l = 0, l < m-tb_def; ++l)
		YB(l, ind) = std::real(evect(l, eval_sort[j].second));
	      pure[ind] = T(1);
	      ++ind;
	      kept[eval_sort[j].second] = true;
	      ++nb_want;

	    if (ritznew[eval_sort[j].second]
		< tol_vp * dal::abs(eval_sort[j].first)) {
		for (size_type l = 0, l < m-tb_def; ++l)
		  YB(l, ind) = std::imag(evect(l, eval_sort[j].second));
		pure[ind] = T(2);
		
		j++;
		kept[eval_sort[j].second] = true;
		
		nb_want++;
		++ind;	      
	    }
	  }
	}
      

      }
      
    }
  }


  template < class Mat, class Vec, class VecB, class Precond >
  void idgmres(const Mat &A, Vec &x, const VecB &b,
	     const Precond &M, int m, iteration& outer) {
    typedef typename linalg_traits<Mat>::value_type T;
    modified_gram_schmidt<T> orth(m, vect_size(x));
    gmres(A, x, b, M, m, outer, orth); 
  }

}

#endif
