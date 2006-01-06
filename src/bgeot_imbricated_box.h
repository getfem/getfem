// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_imbricated_box.h : particular point sort.
//           
// Date    : January 26, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2006 Yves Renard
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

/**@file bgeot_imbricated_box.h
   @brief A comparison function for bgeot::base_node
*/
#ifndef BGEOT_IMBRICATED_BOX
#define BGEOT_IMBRICATED_BOX

#include <bgeot_vector.h>

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
