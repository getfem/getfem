/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2015 Julien Pommier
 
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
#include <vector>
#include "bgeot_config.h"

namespace bgeot {
  /**
     generation of permutations, and ranking/unranking of these.

     based on algorithms detailed in "Ranking and Unranking Permutations in linear time", W. Myrvold, F. Ruskey 
     ( http://www.csr.uvic.ca/~fruskey/Publications/RankPerm.html )
     note that this is not lexigraphical order, and to_rank(0) != {0,1,2,3,...} (it is {1,2,3,...,n,0})

     however, the reset(), finished(), and ++ operator are based on the lexicagraphical ordering
   */
  class permutation : public std::vector<dim_type> {
    size_type remaining;
  public:
    permutation(size_type n) : std::vector<dim_type>(n) { reset(); }
    size_type nb_permutations() { return permutation::nb_permutations(size()); }
    static size_type nb_permutations(size_type n)
    { size_type f=1; for (; n>1; --n) f *= n; return f; }
    void reset() {
      remaining = 1;
      for (size_type i=0; i < size(); ++i)
	{ (*this)[i] = dim_type(i); remaining *= (i+1); }
    }
    permutation& to_rank(size_type r);
    permutation inversed() const {
      permutation pinv(*this);
      for (size_type i=0; i < size(); ++i) pinv[(*this)[i]] = dim_type(i);
      return pinv;
    }
    size_type rank() const;
    bool finished() const { return remaining == 0; }
    /* increment in lexicographical order (not the best, but it is simple) */
    const permutation &operator ++();
    template < typename CONT1, typename CONT2 > void apply_to(const CONT1& src, CONT2& dest) 
    { for (size_type i=0; i < size(); ++i) dest[i] = src[(*this)[i]]; }
  };
  inline permutation& permutation::to_rank(size_type r) {
    reset();
    for (size_type n = size(); n; --n) {
      std::swap((*this)[n-1], (*this)[r % n]);
      r /= n;
    }
    return (*this);
  }
  inline size_type permutation::rank() const { 
    permutation p(*this);
    permutation pinv(p.inversed());
    size_type mul=1, r=0;
    for (size_type n=size(); n>1; --n) {
      dim_type s = p[n-1];
      std::swap(p[n-1], p[pinv[n-1]]);
      std::swap(pinv[s],pinv[n-1]);
      r += s*mul; mul*=n;
    }
    return r;
  }  
  inline const permutation &permutation::operator ++() {
    if (--remaining == 0) return (*this);
    size_type i = size()-2, j=size()-1;
    while ((*this)[i] > (*this)[i+1]) i--;      
    while ((*this)[i] > (*this)[j]) j--;
    std::swap((*this)[i], (*this)[j]);
    for (size_type r = size()-1, s=i+1; r>s; --r, ++s) std::swap((*this)[r],(*this)[s]);
    return (*this);
  }
}
