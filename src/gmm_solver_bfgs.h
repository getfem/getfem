/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_bfgs.h : a bfgs algorithm                         */
/*     									   */
/* Date : October 14 2004.                                                 */
/* Author : Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
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
#ifndef GMM_BFGS_H
#define GMM_BFGS_H

#include <gmm_kernel.h>
#include <gmm_iter.h>

namespace gmm {

  // BFGS algorithm (Broyden, Fletcher, Goldfarb, Shanno)
  // Quasi Newton method for optimization problems.


  // delta[k] = x[k+1] - x[k]
  // gamma[k] = grad f(x[k+1]) - grad f(x[k])
  // H[0] = I
  // Hgamma[k] = H[k] gamma[k]
  // alpha[k] = (1 + (gamma[k]^T H[k] gamma[k]) / delta[k]^T gamma[k])
  //           / delta[k]^T gamma[k]
  // beta[k] = gamma[k]^T gamma[k]
  // H[k+1] = H[k] + alpha[k] delta[k] delta[k]^T
  //           - (gamma[k] gamma[k]^T H[k] + Hgamma[k] delta[k]^T) / beta[k]

  template <typename VECTOR> struct bfgs_invhessian {
    
    typedef typename linalg_traits<VECTOR>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    std::vector<VECTOR> delta, gamma, Hgamma;
    std::vector<T> alpha, beta;

    void hmult(const VECTOR &X, VECTOR &Y) {
      copy(X, Y);
      for (k = 0 ; k < delta.size(); ++k) {
	add(scaled(gamma[k], -gmm::vect_sp(Y, gamma[k])/beta[k]), Y);
	double xgamma = vect_sp(X, delta[k]);
	add(scaled(delta[k], xgamma*alpha[k]), Y);
	add(scaled(Hgamma[k], -xgamma / beta[k]), Y);
      }
    }
    
    void restart(void) {
      delta = std::vector<VECTOR>();
      gamma = std::vector<VECTOR>();
      Hgamma = std::vector<VECTOR>();
      alpha.resize(0); beta.resize(0);
    }
    
    template<typename VECT1, typename VECT2>
    void add(const VECT1 &deltak, const VECT2 &gammak) {
      size_type N = vect_size(deltak), k = delta.size();
      VECTOR Y(N);
      hmult(gammak, Y);
      delta.resize(k+1); gamma.resize(k+1); Hgamma.resize(k+1);
      alpha.resize(k+1); beta.resize(k+1);
      resize(delta[k], N); resize(gamma[k], N); resize(Hgamma[k], N); 
      gmm::copy(deltak, delta[k]);
      gmm::copy(gammak, gamma[k]);
      gmm::copy(Y, Hgamma[k]);
      T a = vect_sp(deltak, gammak);
      T alphak = (T(1) + vect_sp(gammak, Y) / a) / a;
      T betak = vect_sp(gammak, gammak);
      alpha[k] = alphak;
      beta[k] = betak;
    }

    
  }


  template <typename FUNCTION, typename DERIVATIVE, typename VECTOR> 
  void bfgs(const FUNCTION &f, const DERIVATIVE &d, VECTOR &x,
	    int restart, iteration& iter) {

    typedef typename linalg_traits<VECTOR>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    
    bfgs_invhessian<VECTOR> invhessian;
    VECTOR r(vect_size(x)), d(vect_size(x)), y(vect_size(x)), r2(vect_size(x));
    d(x, r);
    T lambda;
    R m1 = 0.3, m2 = 0.6;
    
    while (! iter.finished_vect(r)) {

      invhessian.mult(r, d); gmm::scale(d, T(-1));
      

      // Wolfe Line search
      T lambda_min = 0, lambda_max, val = f(x);
      T derivative = gmm::vect_sp(r, d);
      bool unbounded = true;
      lambda = 1;
      
      for(;;) {

	add(x, scaled(lambda, d), y);
	
	if (f(y) <= val + m1 * lambda * derivative) {
	  d(y, r2);
	  T derivative2 = gmm::vect_sp(r2, d);
	  if (derivative2 >= m2 derivative) break;
	  lambda_min = lambda;
	}
	else {
	  lambda_max = lambda;
	  unbounded = false;
	}
	if (unbounded) lambda *= T(2);
	else  lambda = (lambda_max + lambda_min) / T(2);
      }

      // Rank two update
      gmm::add(scaled(r2, T(-1)), r);
      if (iter.get_iteration() % restart == 0)
	invhessian.restart();
      else
	invhessian.add(gmm::scaled(d, lambda), gmm::scaled(r, T(-1)));
      copy(r2, r); copy(y, x);
    }

  }


}

#endif 

