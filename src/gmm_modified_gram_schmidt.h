// -*- c++ -*-
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
/* File    :  gmm_modified_gram_schmidt.h : from I.T.L.                    */
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

#ifndef GMM_MODIFIED_GRAM_SCHMIDT_ORTHOGONALIZATION_H
#define GMM_MODIFIED_GRAM_SCHMIDT_ORTHOGONALIZATION_H

#include <gmm_solvers.h>

namespace gmm {

  template <class Vec>
  class modified_gram_schmidt {
  public:
    typedef typename linalg_traits<Vec>::size_type size_type; 

    template <class Size>
    modified_gram_schmidt(int restart, const Size& s) 
      : V(restart+1) { do_initialize(s); }
    
    const Vec& operator[](size_type i) const { return V[i]; }

    Vec& operator[](size_type i) { return V[i]; }
    
  protected:
  
    template <class Size>
    void do_initialize(const Size& s)
    { std::fill(V.begin(), V.end(), Vec(s)); }
    
    std::vector<Vec> V;
  };

  template <class Vec, class VecHi, class Size>
  void orthogonalize(modified_gram_schmidt<Vec>& V, const VecHi& _Hi, Size i) {
    VecHi& Hi = const_cast<VecHi&>(_Hi);
    
    for (Size k = 0; k <= i; k++) {
      Hi[k] = gmm::vect_sp(V[i+1], V[k]);
      gmm::add(gmm::scaled(V[k], -Hi[k]), V[i+1]);
    }
  }
  
  template <class Vec, class VecS, class VecX, class Size>
  void combine(modified_gram_schmidt<Vec>& V, const VecS& s, VecX& x, Size i) {
    for (Size j = 0; j < i; ++j)
      gmm::add(gmm::scaled(V[j], s[j]), x);
  }
}

#endif
