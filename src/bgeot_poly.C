/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_poly.C : plain polynomials with several variables.     */
/*     									   */
/*                                                                         */
/* Date : December 01, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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



#include <bgeot_poly.h>
#include <bgeot_matrix.h>

namespace bgeot
{

  size_type alpha(short_type n, short_type d)
  {
    #define STORED 150
    static bgeot::fsmatrix<size_type, STORED> M;
    static bool init = false;
    if (!init)
    {
      for (short_type i = 0; i < STORED; ++i)
      { 
	M(i, 0) = M(0, i) = 1;
	for (short_type j = 1; j <= i; ++j)
	  M(i,j) = M(j,i) = (M(i, j-1) * (i+j)) / j;
      }
      init = true;
    }
    if (n >= STORED || d >= STORED)
      DAL_THROW(internal_error,
		"alpha called with n = " << n << " and d = " << d);

    return M(d, n);
  }

  const power_index &power_index::operator ++()
  { 
    short_type n = size(), l;
    if (n > 0) {
      iterator it = begin() + (n-2);
      for (l = n-2; l != short_type(-1); --l, --it)
	if (*it != 0) break;
      short_type a = (*this)[n-1]; (*this)[n-1] = 0;
      (*this)[short_type(l+1)] = a + 1;
      if (l != short_type(-1)) ((*this)[l])--;
    }
    return *this;
  }
  
  const power_index &power_index::operator --()
  {
    short_type n = size(), l;
    if (n > 0) {
      iterator it = begin() + (n-1);
      for (l = n-1; l != short_type(-1); --l, --it)
	if (*it != 0) break;
      if (l != short_type(-1)) {
	short_type a = (*this)[l]; (*this)[l] = 0; (*this)[n-1] = a - 1;
	if (l > 0) ((*this)[l-1])++;
      }
    }
    return *this;
  }
  
  short_type power_index::degree(void) const
  {
    short_type d = 0; const_iterator it = begin(), ite = end();
    for ( ; it != ite; ++it) d += *it;
    return d;
  }

  size_type power_index::global_index(void) const
  {
    short_type d = degree(), n = size();
    size_type l = 0;
    const_iterator it = begin(), ite = end();
    for ( ; it != ite && d > 0; ++it)
    { l += alpha(n, d-1); d -= *it; --n; }
    return l;
  }

  power_index::power_index(short_type nn) : std::vector<short_type>(nn)
  { std::fill(begin(), end(), short_type(0)); }

}  /* end of namespace bgeot.                                             */
