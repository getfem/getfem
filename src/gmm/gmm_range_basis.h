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
   @brief Extract a basis of the range of a (large sparse) matrix from the 
          columns of this matrix.
*/
#ifndef GMM_RANGE_BASIS_H
#define GMM_RANGE_BASIS_H

#include "gmm_kernel.h"
#include "gmm_iter.h"
#include <set>


namespace gmm {

  template <typename Mat>
  void range_basis_eff(const Mat &B, std::set<size_type> &columns,
		       double EPS) {
   
    typedef std::set<size_type> TAB;
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc_r = columns.size(), k;
    double SQRTEPS = gmm::sqrt(EPS);

    std::vector<T> w(mat_nrows(B));
    
    while (nc_r) {

      std::vector<T> v(nc_r), v0(nc_r);

      // Spectral radius of B^* B
      
      R rho = R(0), rho2 = R(0);

      gmm::fill_random(v);
      for (size_type i = 0; i < 1000000; ++i) {
	R rho_old = rho;
	gmm::clear(w);
	gmm::copy(v, v0);
	k = 0;
	for (TAB::iterator it = columns.begin(); it!=columns.end(); ++it, ++k)
	  add(scaled(mat_col(B, *it), v[k]), w);
	
	k = 0;
	for (TAB::iterator it = columns.begin(); it!=columns.end(); ++it, ++k)
	  v[k] = vect_hp(w, mat_col(B, *it));

	rho = gmm::abs(vect_hp(v, v0) / vect_hp(v0, v0)); // Rayleigh quotient

	if (gmm::abs(rho_old-rho) <= rho*1E-5) break;
	  
	gmm::scale(v, T(1)/vect_norm2(v));
      }

      rho *= R(8)/R(15);

      // Computing an element of the null space of de B^* B
      gmm::fill_random(v);
      for (size_type i = 0; i < 1000000; ++i) {
	R rho_old = rho2;
	gmm::clear(w);
	gmm::copy(v, v0);
	k = 0;
	for (TAB::iterator it = columns.begin(); it!=columns.end(); ++it, ++k)
	  add(scaled(mat_col(B, *it), v[k]), w);
	
	k = 0;
	for (TAB::iterator it = columns.begin(); it!=columns.end(); ++it, ++k)
	  v[k] = vect_hp(w, mat_col(B, *it)) - v[k]*rho;

	rho2 = gmm::abs(vect_hp(v, v0) / vect_hp(v0, v0)); // Rayleigh quotient
	if (gmm::abs(rho_old-rho2) <= rho*EPS/R(50)) break;
	if (gmm::abs(rho_old-rho2) <= rho*SQRTEPS
	    && gmm::abs(rho2 - rho)*SQRTEPS > gmm::abs(rho_old-rho2))
	  break;

	gmm::scale(v, T(1)/vect_norm2(v));
      }
      
      if (gmm::abs(rho-rho2) < EPS*rho) {
	size_type j_max = size_type(-1);
	R val_max = R(0);

	k = 0;
	for (TAB::iterator it=columns.begin(); it!=columns.end(); ++it) {
	  if (gmm::abs(v[k]) > val_max) {
	    val_max = gmm::abs(v[k]);
	    j_max = *it;
	  }
	  ++k;
	}
	GMM_ASSERT1(j_max != size_type(-1), "Internal error");
	columns.erase(j_max); nc_r = columns.size();
      }
      else break;
    }
  }


  template <typename Mat>
  void range_basis_eff_dense(const Mat &B, std::set<size_type> &columns,
			     double EPS) {
   
    typedef std::set<size_type> TAB;
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc_r = columns.size(), nc = mat_ncols(B), nr = mat_nrows(B);


    gmm::row_matrix< gmm::wsvector<T> > H(nc, nc_r), BB(nr, nc_r);
    
    size_type i = 0;
    for (TAB::iterator it=columns.begin(); it!=columns.end(); ++it, ++i)
      H(*it, i) = T(1);

    gmm::mult(B, H, BB);
    gmm::dense_matrix<T> M(nc_r, nc_r);
    gmm::mult(gmm::transposed(BB), BB, M);

    std::vector<size_type> ipvt(nc_r);
    gmm::lu_factor(M, ipvt);

    R elt_max = R(0);
    for (i = 0; i < nc_r; ++i) {
      elt_max = std::max(elt_max, gmm::abs(M(i,i)));
      // cout << gmm::abs(M(i,i)) << "  ";
    }
    // cout << endl << endl;

    i = 0;
    std::set<size_type> c = columns;
    for (TAB::iterator it = c.begin(); it != c.end(); ++it, ++i)
      if (gmm::abs(M(i,i)) <= EPS*elt_max) columns.erase(*it);   
  }





//   template <typename Mat>
//   void range_basis_rec(const Mat &B, std::set<size_type> &columns,
// 		   double EPS) {
   
//     typedef typename linalg_traits<Mat>::value_type T;
//     typedef typename number_traits<T>::magnitude_type R;

//     size_type nc_r = columns.size();
//     if (nc_r > 100) {
//       std::set<size_type> c1, c2, c3, c4, c5, c6, c7, c8;
//       size_type k = 0;
//       for (std::set<size_type>::iterator it = columns.begin();
// 	   it != columns.end(); ++it, ++k) 
// 	if (k < nc_r/8) c1.insert(*it);
// 	else if (k < 2*nc_r/8) c2.insert(*it);
// 	else if (k < 3*nc_r/8) c3.insert(*it);
// 	else if (k < 4*nc_r/8) c4.insert(*it);
// 	else if (k < 5*nc_r/8) c5.insert(*it);
// 	else if (k < 6*nc_r/8) c6.insert(*it);
// 	else if (k < 7*nc_r/8) c7.insert(*it);
// 	else c8.insert(*it);
//       size_type c1size = c1.size();
//       range_basis_rec(B, c1, EPS);
//       if (c1.size() < c1size) {
// 	range_basis_rec(B, c2, EPS);
// 	range_basis_rec(B, c3, EPS);
// 	range_basis_rec(B, c4, EPS);
// 	range_basis_rec(B, c5, EPS);
// 	range_basis_rec(B, c6, EPS);
// 	range_basis_rec(B, c7, EPS);
// 	range_basis_rec(B, c8, EPS);
// 	columns = c1;
// 	for (std::set<size_type>::iterator it = c2.begin(); it != c2.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c3.begin(); it != c3.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c4.begin(); it != c4.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c5.begin(); it != c5.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c6.begin(); it != c6.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c7.begin(); it != c7.end(); ++it)
// 	  columns.insert(*it);
// 	for (std::set<size_type>::iterator it = c8.begin(); it != c8.end(); ++it)
// 	  columns.insert(*it);
//       }
//     }
//     range_basis_eff(B, columns, EPS);
//   }


  template <typename Mat>
  void range_basis_rec(const Mat &B, std::set<size_type> &columns,
		   double EPS) {
   
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc_r = columns.size();
    if (nc_r > 50) {
      cout << "begin small range basis" << endl;
      std::set<size_type> c1, cres;
      for (std::set<size_type>::iterator it = columns.begin();
	   it != columns.end(); ++it) {
	c1.insert(*it);
	if (c1.size() >= 250) {
	  size_type c1size = c1.size();
	  range_basis_eff_dense(B, c1, EPS);
	  for (std::set<size_type>::iterator it2=c1.begin(); it2 != c1.end();
	       ++it2) cres.insert(*it2);

	  if (c1.size() == c1size) {
	    for (size_type i = 0; it != columns.end() && i < 1000; ++it, ++i)
	      cres.insert(*it);
	  }
	  c1.clear(); 
	}
      }
      
      if (c1.size() > 10) {
	range_basis_eff_dense(B, c1, EPS);
	for (std::set<size_type>::iterator it = c1.begin(); it != c1.end(); ++it)
	  cres.insert(*it); 
      }
      columns = cres;
    }
    cout << "begin global range basis for " << columns.size() << " columns " << endl;
    range_basis_eff(B, columns, EPS);
    cout << "begin global dense range basis for " << columns.size() << " columns " << endl;
    range_basis_eff_dense(B, columns, EPS);
  }



  template <typename Mat>
  void range_basis(const Mat &B, std::set<size_type> &columns,
		   double EPS, col_major) {
   
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;
  
    R norm_max = R(0);
    for (size_type i = 0; i < mat_ncols(B); ++i)
      norm_max=std::max(norm_max, vect_norm2(mat_col(B, i)));
    
    columns.clear();
    for (size_type i = 0; i < mat_ncols(B); ++i)
      if (vect_norm2(mat_col(B, i)) >= norm_max*EPS)
	columns.insert(i);

    range_basis_rec(B, columns, EPS);
  }


  template <typename Mat>
  void range_basis(const Mat &B, std::set<size_type> &columns,
		   double EPS, row_major) {
    typedef typename  linalg_traits<Mat>::value_type T;
    gmm::col_matrix< rsvector<T> > BB(mat_nrows(B), mat_ncols(B));
    GMM_WARNING3("A copy of a row matrix is done into a column matrix "
		 "for range basis algorithm.");
    gmm::copy(B, BB);
    range_basis(BB, columns, EPS);
  }

  /** Range Basis :     
    Extract a basis of the range of a (large sparse) matrix from the columns of
    this matrix.
  */
  template <typename Mat>
  void range_basis(const Mat &B, std::set<size_type> &columns,
		   double EPS=1E-13) {
    range_basis(B, columns, EPS,
		typename principal_orientation_type
		<typename linalg_traits<Mat>::sub_orientation>::potype());
  }

}

#endif
