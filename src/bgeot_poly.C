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
#include <bgeot_vector.h>
#include <numeric>
namespace bgeot
{
    #define STORED 150
  static gmm::dense_matrix<size_type> _alpha_M(STORED, STORED);
  static void _alpha_init() {
    static bool init = false;
    if (!init)
    {
      for (short_type i = 0; i < STORED; ++i)
      { 
	_alpha_M(i, 0) = _alpha_M(0, i) = 1;
	for (short_type j = 1; j <= i; ++j)
	  _alpha_M(i,j) = _alpha_M(j,i) = (_alpha_M(i, j-1) * (i+j)) / j;
      }
      init = true;
    }
  }
  static inline size_type _alpha(short_type n, short_type d) { return _alpha_M(d,n); }

  size_type alpha(short_type n, short_type d)
  {
    _alpha_init();
    if (n >= STORED || d >= STORED)
      DAL_THROW(internal_error,
		"alpha called with n = " << n << " and d = " << d);
    return _alpha(n,d);
  }

  const power_index &power_index::operator ++()
  { 
    short_type n = size(), l;
    if (n > 0) {
      size_type g_idx = _global_index; short_type deg = _degree;
      iterator it = begin() + (n-2);
      for (l = n-2; l != short_type(-1); --l, --it)
	if (*it != 0) break;
      short_type a = (*this)[n-1]; (*this)[n-1] = 0;
      (*this)[short_type(l+1)] = a + 1;
      if (l != short_type(-1)) ((*this)[l])--;
      else if (deg+1) _degree = deg+1;
      if (g_idx+1) _global_index = g_idx+1;
      //_degree = short_type(-1);
    }
    return *this;
  }
  
  const power_index &power_index::operator --()
  {
    short_type n = size(), l;
    if (n > 0) {
      size_type g_idx = _global_index; short_type deg = _degree;
      iterator it = begin() + (n-1);
      for (l = n-1; l != short_type(-1); --l, --it)
	if (*it != 0) break;
      if (l != short_type(-1)) {
	short_type a = (*this)[l]; (*this)[l] = 0; (*this)[n-1] = a - 1;
	if (l > 0) ((*this)[l-1])++; 
        else if (deg+1) _degree = deg-1;
      }
      if (g_idx+1) _global_index = g_idx-1;
    }
    return *this;
  }
  
  short_type power_index::degree() const
  { 
    if (_degree != short_type(-1)) return _degree;
    _degree = std::accumulate(begin(), end(), 0); 
    return _degree;
  }

  size_type power_index::global_index(void) const
  {
    if (_global_index != size_type(-1)) return _global_index;
    short_type d = degree(), n = size();
    _global_index = 0;
    const_iterator it = begin(), ite = end();
    for ( ; it != ite && d > 0; ++it)
    { _global_index += _alpha(n, d-1); d -= *it; --n; }
    return _global_index;
  }
  /*
  size_type power_index::global_index_shift(size_type idx, short_type k) {
    _degree = degree() + k - v[idx];
    
  }

  size_type power_index::set(size_type idx, short_type v) {
    short_type deg = degree(); size_type pi = power_index();
    
  }
  */
  power_index::power_index(short_type nn) : v(nn), _degree(0), _global_index(0)
  { std::fill(begin(), end(), short_type(0)); _alpha_init(); }

}  /* end of namespace bgeot.                                             */
