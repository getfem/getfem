// -*- c++ -*-
//
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
//                                                                           
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_qmr.h : from I.T.L.                               */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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
#ifndef GMM_QMR_H
#define GMM_QMR_H

#include <gmm_solvers.h>

namespace gmm {

  //: Quasi-Minimal Residual
  //
  //  This routine solves the unsymmetric linear system Ax = b using the
  //  Quasi-Minimal Residual method.
  // 
  //See: R. W. Freund and N. M. Nachtigal, A quasi-minimal residual method for 
  //non-Hermitian linear systems, Numerical Math., 60(1991), pp. 315-339
  //
  // Preconditioner -  Incomplete LU, Incomplete LU with threshold,
  //                   SSOR or identity_preconditioner.

  template <class Matrix, class Vector, class VectorB, class Precond1>
  void qmr(const Matrix &A, Vector &x, const VectorB &b, const Precond1 &M1,
	   iteration& iter)
  {
    typedef typename linalg_traits<Vector>::value_type value_type;
    value_type delta(0), ep(0), beta(0), rho_1(0), gamma_1(0), theta_1(0);

    typedef typename temporary_vector<Vector>::vector_type TmpVec;
    size_t nn = vect_size(x);
    TmpVec r(nn), v_tld(nn), y(nn), w_tld(nn), z(nn), v(nn), w(nn);
    TmpVec y_tld(nn), z_tld(nn), p(nn), q(nn), p_tld(nn), d(nn), s(nn);

    iter.set_rhsnorm(gmm::vect_norm2(b));
    if (iter.get_rhsnorm() == 0.0) { clear(x); return; }

    gmm::mult(A, gmm::scaled(x, -1.0), b, r);
    gmm::copy(r, v_tld);

    gmm::mult_left(M1, v_tld, y);
    value_type rho = gmm::vect_norm2(y);

    gmm::copy(r, w_tld);
    gmm::transposed_mult_right(M1, w_tld, z);
    value_type xi = gmm::vect_norm2(z);
  
    value_type gamma = 1.0, eta = -1.0, theta = 0.0;
  
    while (! iter.finished_vect(r)) {
    
      if (rho == 0.0 || xi == 0.0)
	DAL_THROW(failure_error, "QMR failed to converge");

      gmm::copy(gmm::scaled(v_tld, 1./rho), v);
      gmm::scale(y, 1./rho);

      gmm::copy(gmm::scaled(w_tld, 1./xi), w);
      gmm::scale(z, 1./xi);

      delta = gmm::vect_sp(z, y);
      if (delta == 0.0) DAL_THROW(failure_error, "QMR failed to converge");

      gmm::mult_right(M1, y, y_tld);		
      gmm::transposed_mult_left(M1, z, z_tld);

      if (iter.first()) {
	gmm::copy(y_tld, p);
	gmm::copy(z_tld, q);
      } else {
	gmm::add(y_tld, gmm::scaled(p, -(xi  * delta / ep)), p);
	gmm::add(z_tld, gmm::scaled(q, -(rho * delta / ep)), q);
      }
    
      gmm::mult(A, p, p_tld);

      ep = gmm::vect_sp(q, p_tld);
      if (ep == 0.0) DAL_THROW(failure_error, "QMR failed to converge");

      beta = ep / delta;
      if (beta == 0.0) DAL_THROW(failure_error, "QMR failed to converge");

      gmm::add(p_tld, gmm::scaled(v, -beta), v_tld);
      gmm::mult_left(M1, v_tld, y);

      rho_1 = rho;
      rho = gmm::vect_norm2(y);

      gmm::mult(gmm::transposed(A), q, w_tld);
      gmm::add(w_tld, gmm::scaled(w, -beta), w_tld);
      gmm::transposed_mult_right(M1, w_tld, z);

      xi = gmm::vect_norm2(z);

      gamma_1 = gamma;
      theta_1 = theta;

      theta = rho / (gamma_1 * beta);
      gamma = 1.0 / sqrt(1.0 + theta * theta);

      if (gamma == 0.0) DAL_THROW(failure_error, "QMR failed to converge");

      eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);

      if (iter.first()) {
	gmm::copy(gmm::scaled(p, eta), d);
	gmm::copy(gmm::scaled(p_tld, eta), s);
      } else {
	value_type tmp = (theta_1 * theta_1 * gamma * gamma);
	gmm::add(gmm::scaled(p, eta), gmm::scaled(d, tmp), d);
	gmm::add(gmm::scaled(p_tld, eta), gmm::scaled(s, tmp), s);
      }
      gmm::add(d, x);
      gmm::add(gmm::scaled(s, -1.), r);

      ++iter;
    }
  }


}

#endif 

