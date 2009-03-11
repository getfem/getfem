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
  template <typename Mat> range_basis(const Mat &B, double EPS,
				      col_major) {
    size_type nc = mat_ncols(B), nc_r = nc;
    std::vector<bool> b(nc, true);

    typedef typename linalg_traits<V>::value_type T
    typedef typename number_traits<T>::magnitude_type R;

    R norm_max = R(0);
    for (size_type i = 0; i < nc; ++i)
      norm_max=std::max(norm_max, vect_norm2(mat_col(B, i)));

    // Is there a better strategy ?
    for (size_type i = 0; i < nc; ++i)
      if (vect_norm2(mat_col(B, i)) < norm_max*EPS)
	{ b[i] = false; nc--; }

    std::vector<T> w(mat_nrows(B))
    
    while (nc_r) {

      std::vector<T> v(nc_r);

      // recherche du rayon spectral de B^* B
      
      R rho = 0;

      gmm::fill_random(v);
      for (size_type i = 0; i < 1000 /* ...*/ ; ++i) {
	gmm::clear(w);
	size_type k = 0;
	for (size_type j = 0; j < nc)
	  if (b[j]) { gmm::add(gmm::scaled(mat_col(B, j), v[k]), w); k++; }
	
	k = 0;
	for (size_type j = 0; j < nc)
	  if (b[j]) { v[k] = gmm::vect_hp(w, mat_col(B, j)); k++; }
	  
	gmm::scale(v, T(1)/vect_norm2(v));

	// + calcul coefficient de rayleight

      }



      // recherche d'un élément du noyau de B^* B

      gmm::fill_random(v);
      for (size_type i = 0; i < 1000 /* ...*/ ; ++i) {
	gmm::clear(w);
	size_type k = 0;
	for (size_type j = 0; j < nc)
	  if (b[j]) { gmm::add(gmm::scaled(mat_col(B, j), v[k]), w); k++; }
	
	k = 0;
	for (size_type j = 0; j < nc)
	  if (b[j]) {
	    v[k] = gmm::vect_hp(w, mat_col(B, j)) - v[k]*rho*R(4)/R(7);
	    k++;
	  }
	  
	gmm::scale(v, T(1)/vect_norm2(v));

	// + calcul coefficient de rayleight

      }
      

      // + selection du terme maximal dans v si la valeur propre est proche
      // de rho*R(4)/R(7) sinon break




    }


  }

  template <typename Mat> range_basis(const Mat &B, double EPS,
                                      row_major) {
    col_matrix< rs_vector<typename linalg_traits<Mat>::value_type >
      BB(mat_nrows(B), mat_ncols(B));
    GMM_WARNING3("A copy of a row matrix is done into a column matrix "
		 "for range basis algorithm.");
    gmm::copy(B, BB);
    range_basis(BB, EPS);
  }

  template <typename Mat> range_basis(const Mat &B, double EPS=1E-14) {
    range_basis(B, EPS,
		typename principal_orientation_type
		<typename linalg_traits<Mat>::sub_orientation>::potype());
  }
  


}

#endif
