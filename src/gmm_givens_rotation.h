// -*- c++ -*-
//
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
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
/* File    :  gmm_given_rotation.h : from I.T.L.                           */
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

#ifndef GMM_GIVENS_ROTATION_H
#define GMM_GIVENS_ROTATION_H

#include <gmm_solvers.h>


namespace gmm {

template <class T>
class givens_rotation {
public:

  //: Default constructor
  inline givens_rotation() : a_(0), b_(0), c_(0), s_(0) { }

  //: Givens Plane Rotation Constructor
  inline givens_rotation(const T& a_in, const T& b_in) {
    T roe;
    if (dal::abs(a_in) > dal::abs(b_in))
      roe = a_in;
    else
      roe = b_in;
    
    T scal = dal::abs(a_in) + dal::abs(b_in);
    T r, z;
    if (scal != T(0)) {
      T a_scl = a_in / scal;
      T b_scl = b_in / scal;
      r = scal * ::sqrt(a_scl * a_scl + b_scl * b_scl);
      if (roe < T(0)) r *= -1;
      c_ = a_in / r;
      s_ = b_in / r;
      z = 1;
      if (dal::abs(a_in) > dal::abs(b_in))
        z = s_;
      else if (dal::abs(b_in) >= dal::abs(a_in) && c_ != T(0))
        z = T(1) / c_;
    } else {
      c_ = 1; s_ = 0; r = 0; z = 0;      
    }
    a_ = r;
    b_ = z;
  }

  inline void set_cs(const T& cin, const T& sin) { c_ = cin; s_ = sin; }

  //: Apply plane rotation to two real scalars. (name change a VC++ workaround)
  inline void scalar_apply(T& x, T& y) {
    T tmp = c_ * x + s_ * y;
    y = c_ * y - s_ * x;
    x = tmp;
  }

  //: Apply plane rotation to two vectors.
  template <class VecX, class VecY>
  inline void apply(const VecX& x_, const VecY& y_) {
    VecX& x = const_cast<VecX&>(x_);
    VecY& y = const_cast<VecY&>(y_);

    typename linalg_traits<VecX>::iterator xi = x.begin(), xend = x.end();
    typename linalg_traits<VecY>::iterator yi = y.begin();

    while ( xi != xend ) { scalar_apply(*xi, *yi); ++xi; ++yi; }
  }

  inline T a() { return a_; }
  inline T b() { return b_; }
  inline T c() { return c_; }
  inline T s() { return s_; }
protected:
  T sign(const T& t) { T ret = 1; if ( t < T(0) ) ret = -1; return ret; }
  T a_, b_;
  T c_, s_;
};



//:  The specialization for complex numbers.

  /*
  Using a and b to represent elements of an input complex vector, the CROTG
  and ZROTG functions calculate the elements real c and complex s of an
  orthogonal matrix such that:

              c*a + s*b = r
  -conjugate(s)*a + c*b = 0
  */


//!category: functors
//!component: type
template <class T>
class givens_rotation < std::complex<T> > {
  typedef std::complex<T> C;
public:
  //:
  inline givens_rotation() : cs(0), sn(0) { }
  
  inline T abs_sq(C t) 
  { return std::real(t) * std::real(t) + std::imag(t) * std::imag(t); }
    
  //:
  inline givens_rotation(const C& a_in, const C& b_in) {

    T a = dal::abs(a_in), b = dal::abs(b_in);
    if ( a == T(0) ) {
      cs = T(0);
      sn = C(1.);
    } else {
      T scale = a + b;
      T norm = ::sqrt(abs_sq(a_in/scale)+abs_sq(b_in/scale)) * scale;
    
      cs = a / norm;
      sn = a_in/a * std::conj(b_in)/norm;
      
      //in zrotg there is an assignment for ca, what is that for? 
    }
  }
  //:  Apply plane rotation to two vectors.
  template <class VecX, class VecY>
  inline void apply(const VecX& x_, const VecY& y_) {
    VecX& x = const_cast<VecX&>(x_);
    VecY& y = const_cast<VecY&>(y_);
    
    typename linalg_traits<VecX>::iterator xi = x.begin(), xend = x.end();
    typename linalg_traits<VecY>::iterator yi = y.begin();
    
    while (xi != xend ) { scalar_apply(*xi, *yi); ++xi; ++yi; }
  }
  //: Apply plane rotation to two complex scalars.
  inline void scalar_apply(C& x, C& y) {
    //complex<T> temp  =  std::conj(cs) * x + std::conj(sn) * y;
    //y = cs * y - sn * x; 
    std::complex<T> temp  =  cs * x + sn * y;
    y = cs * y - std::conj(sn) * x;     
    x = temp;
  }
  inline void set_cs(const T& cs_, const C& sn_) {
    cs = cs_; sn = sn_;
  }

  inline T c() { return cs; }
  inline C s() { return sn; }

protected:
  T cs;
  C sn;
};
}

#endif
