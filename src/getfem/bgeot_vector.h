/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 1995-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file bgeot_vector.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 01, 1995.
   @brief some (old) vector class, and some typedefs.

   vsvector should be replaced by gmm calls.
*/

#ifndef BGEOT_VECTOR_H__
#define BGEOT_VECTOR_H__

#include "bgeot_config.h"
#include "gmm/gmm_kernel.h"
#include "bgeot_small_vector.h"

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

#include "gmm/gmm_interface_bgeot.h"

#endif  /* BGEOT_VECTOR_H__ */
