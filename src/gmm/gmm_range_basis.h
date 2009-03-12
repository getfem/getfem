// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard
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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file gmm_range_basis.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date March 10, 2009.
   @brief Extract a basis of the range of a (sparse) matrix from the columns of
          this matrix.
*/
#ifndef GMM_RANGE_BASIS_H
#define GMM_RANGE_BASIS_H

#include "gmm_kernel.h"
#include "gmm_iter.h"

namespace gmm {

  /** Range Basis :     
    Extract a basis of the range of a (sparse) matrix from the columns of
    this matrix.
  */
  template <typename Mat>
  void range_basis(const Mat &B, std::vector<bool> &columns,
		   double EPS, col_major) {
   
    size_type nc = mat_ncols(B), nc_r = nc;
    columns.resize(nc); gmm::fill(columns, true);

    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    R norm_max = R(0);
    for (size_type i = 0; i < nc; ++i)
      norm_max=std::max(norm_max, vect_norm2(mat_col(B, i)));

    // Is there a better strategy ?
    for (size_type i = 0; i < nc; ++i) {
      // cout << "i : " << vect_norm2(mat_col(B, i)); 
      if (vect_norm2(mat_col(B, i)) < norm_max*EPS)
	{ columns[i] = false; --nc_r; /* cout << " elminated" << endl; */ }
      // else  cout << " kept" << endl;
    }

    std::vector<T> w(mat_nrows(B));
    
    while (nc_r) {

      // cout << "nc_r = " << nc_r << endl;
      std::vector<T> v(nc_r), v0(nc_r);

      // Spectral radius of B^* B
      
      R rho = R(0), rho2 = R(0);

      gmm::fill_random(v);
      for (size_type i = 0; i < 1000000; ++i) {
	R rho_old = rho;
	gmm::clear(w);
	gmm::copy(v, v0);
	size_type k = 0;
	for (size_type j = 0; j < nc; ++j)
	  if (columns[j]) { add(scaled(mat_col(B, j), v[k]), w); ++k; }
	
	k = 0;
	for (size_type j = 0; j < nc; ++j)
	  if (columns[j]) { v[k] = vect_hp(w, mat_col(B, j)); ++k; }

	rho = vect_hp(v, v0) / vect_hp(v0, v0); // Rayleigh quotient

	if (gmm::abs(rho_old-rho) <= rho*1E-3) break;
	  
	gmm::scale(v, T(1)/vect_norm2(v));
      }

      rho *= R(4)/R(7);

      // Computing an element of the null space of de B^* B

      gmm::fill_random(v);
      for (size_type i = 0; i < 1000000; ++i) {
	R rho_old = rho2;
	gmm::clear(w);
	gmm::copy(v, v0);
	size_type k = 0;
	for (size_type j = 0; j < nc; ++j)
	  if (columns[j]) { add(scaled(mat_col(B, j), v[k]), w); ++k; }
	
	k = 0;
	for (size_type j = 0; j < nc; ++j)
	  if (columns[j]) {
	    v[k] = vect_hp(w, mat_col(B, j)) - v[k]*rho;
	    ++k;
	  }
	GMM_ASSERT1(k == nc_r, "Internal error");

	rho2 = vect_hp(v, v0) / vect_hp(v0, v0); // Rayleigh quotient
	if (gmm::abs(rho_old-rho2) <= rho*EPS/R(1000)) break;

	gmm::scale(v, T(1)/vect_norm2(v));
      }
      
      // cout << "rho = " << rho << " rho2 = " << rho2 << " rho+rho2 " << (rho+rho2)/rho << endl;
    
      if (gmm::abs(rho+rho2) < EPS*rho) {
	size_type k = 0, j_max = size_type(-1);
	R val_max = R(0);
	for (size_type j = 1; j < nc; ++j)
	  if (columns[j]) {
	    if (gmm::abs(v[k]) > val_max) {
	      val_max = gmm::abs(v[k]);
	      j_max = j;
	    }
	    ++k;
	  }
	GMM_ASSERT1(j_max != size_type(-1), "Internal error");
	cout << "Eliminating " << j_max << endl;
	cout << "v = " << v << endl;
	columns[j_max] = false; --nc_r;
      }
      else break;
    }

  }

  template <typename Mat>
  void range_basis(const Mat &B, std::vector<bool> &columns,
		   double EPS, row_major) {
    typedef typename  linalg_traits<Mat>::value_type T;
    gmm::col_matrix< rsvector<T> > BB(mat_nrows(B), mat_ncols(B));
    GMM_WARNING3("A copy of a row matrix is done into a column matrix "
		 "for range basis algorithm.");
    gmm::copy(B, BB);
    range_basis(BB, columns, EPS);
  }

  template <typename Mat>
  void range_basis(const Mat &B, std::vector<bool> &columns,
		   double EPS=1E-12) {
    range_basis(B, columns, EPS,
		typename principal_orientation_type
		<typename linalg_traits<Mat>::sub_orientation>::potype());
  }

}

#endif
