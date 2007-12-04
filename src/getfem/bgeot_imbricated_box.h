// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

/**@file bgeot_imbricated_box.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date January 26, 1999.
   @brief A comparison function for bgeot::base_node
*/
#ifndef BGEOT_IMBRICATED_BOX
#define BGEOT_IMBRICATED_BOX

#include "bgeot_vector.h"

namespace bgeot {
  inline scalar_type sfloor(scalar_type x)
  { return (x >= 0) ? floor(x) : -floor(-x); }

  /// A comparison function for bgeot::base_node
  struct imbricated_box_less
    : public std::binary_function<base_node, base_node, int>
  { 
    mutable int exp_max, exp_min;
    mutable scalar_type c_max;
    unsigned base;

    /// comparaison function
    int operator()(const base_node &x, const base_node &y) const;
    
    imbricated_box_less(unsigned ba = 10, int emi = -15, int ema = -2) {
      base = ba; exp_max = ema; exp_min = emi;
      c_max = pow(double(base), double(-exp_max));
    }
  };
}
#endif
