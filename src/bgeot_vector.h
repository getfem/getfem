// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_vector.h : plain vectors.
//           
// Date    : June 01, 1995.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1995-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file bgeot_vector.h
   @brief some (old) vector class, and some typedefs.

   vsvector should be replaced by gmm calls.
*/

#ifndef BGEOT_VECTOR_H__
#define BGEOT_VECTOR_H__

#include <bgeot_config.h>
#include <gmm_kernel.h>
#include <bgeot_small_vector.h>

namespace bgeot {

  typedef std::vector<scalar_type> base_vector;
  typedef small_vector<scalar_type> base_small_vector;
  typedef base_small_vector base_node;
  typedef gmm::dense_matrix<scalar_type> base_matrix;

  template <class VEC_CONT>
  void vectors_to_base_matrix(base_matrix &G, const VEC_CONT &a) {
    size_type P = (*(a.begin())).size(), NP = a.end() - a.begin();
    G.resize(P, NP);
    typename VEC_CONT::const_iterator it = a.begin(), ite = a.end();
    base_matrix::iterator itm = G.begin();
    for (; it != ite; ++it, itm += P)
      std::copy((*it).begin(), (*it).end(), itm);
  }

}  /* end of namespace bgeot.                                           */

namespace std {
  inline void swap(bgeot::base_node& a, bgeot::base_node& b) { a.swap(b); }
}

#include <gmm_interface_bgeot.h>

#endif  /* BGEOT_VECTOR_H__ */
