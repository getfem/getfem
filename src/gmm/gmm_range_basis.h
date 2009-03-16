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

  // Range basis with the power method
  // Complex version not verified
  template <typename Mat>
  void range_basis_eff_power(const Mat &B, std::set<size_type> &columns,
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
	if (gmm::abs(rho_old-rho2) <= rho*EPS) break;
	if (gmm::abs(rho_old-rho2) <= rho*SQRTEPS
	    && gmm::abs(rho2 - rho)*SQRTEPS > gmm::abs(rho_old-rho2))
	  break;

	gmm::scale(v, T(1)/vect_norm2(v));
      }
      
      if (gmm::abs(rho-rho2) < EPS*rho*1000) {
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
	cout << "column " << j_max << " supressed, remaining "
	     << columns.size() << endl;
      }
      else break;
    }
  }

  // Range basis with LU decomposition. Not stable from a numerical viewpoint.
  // Complex version not verified
  template <typename Mat>
  void range_basis_eff_lu(const Mat &B, std::set<size_type> &columns,
			     std::vector<bool> &c_ortho, double EPS) {
   
    typedef std::set<size_type> TAB;
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc_r = 0, nc_o = 0, nc = mat_ncols(B), nr = mat_nrows(B), i, j;

    for (TAB::iterator it=columns.begin(); it!=columns.end(); ++it)
      if (!(c_ortho[*it])) ++nc_r; else nc_o++;

    if (nc_r > 0) {
 
      gmm::row_matrix< gmm::rsvector<T> > Hr(nc, nc_r), Ho(nc, nc_o);
      gmm::row_matrix< gmm::rsvector<T> > BBr(nr, nc_r), BBo(nr, nc_o);
      
      i = j = 0;
      for (TAB::iterator it=columns.begin(); it!=columns.end(); ++it)
	if (!(c_ortho[*it]))
	  { Hr(*it, i) = T(1)/ vect_norminf(mat_col(B, *it)); ++i; }
	else
	  { Ho(*it, j) = T(1)/ vect_norm2(mat_col(B, *it)); ++j; }
      
      gmm::mult(B, Hr, BBr);
      gmm::mult(B, Ho, BBo);
      gmm::dense_matrix<T> M(nc_r, nc_r), BBB(nc_r, nc_o), MM(nc_r, nc_r);
      gmm::mult(gmm::conjugated(BBr), BBr, M);
      gmm::mult(gmm::conjugated(BBr), BBo, BBB);
      gmm::mult(BBB, gmm::conjugated(BBB), MM);
      gmm::add(gmm::scaled(MM, T(-1)), M);
      
      std::vector<int> ipvt(nc_r);
      gmm::lu_factor(M, ipvt);
            
      R emax = R(0);
      for (i = 0; i < nc_r; ++i) emax = std::max(emax, gmm::abs(M(i,i)));

      i = 0;
      std::set<size_type> c = columns;
      for (TAB::iterator it = c.begin(); it != c.end(); ++it)
	if (!(c_ortho[*it])) {
	  if (gmm::abs(M(i,i)) <= EPS*emax) columns.erase(*it);
	  ++i;
	}
    }
  }


  // Range basis with Gram-Schmidt orthogonalization (sparse version)
  // Complex version not verified
  template <typename Mat>
  void range_basis_eff_gram_schmidt_sparse(const Mat &BB,
					   std::set<size_type> &columns,
					   std::vector<bool> &c_ortho,
					   double EPS) {
   
    typedef std::set<size_type> TAB;
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc = mat_ncols(BB), nr = mat_nrows(BB);

    gmm::col_matrix< rsvector<T> > B(nr, nc);
    for (std::set<size_type>::iterator it = columns.begin();
	 it != columns.end(); ++it) {
      gmm::copy(mat_col(BB, *it), mat_col(B, *it));
      gmm::scale(mat_col(B, *it), T(1)/vect_norm2(mat_col(B, *it)));
    }

    std::set<size_type> c = columns, rc = columns;
    // cout << "debut ortho groupee" << endl;
    
    for (std::set<size_type>::iterator it = c.begin(); it != c.end(); ++it)
      if (c_ortho[*it]) {
	for (std::set<size_type>::iterator it2 = rc.begin();
	     it2 != rc.end(); ++it2)
	  if (!(c_ortho[*it2])) {
	    T r = -vect_hp(mat_col(B, *it2), mat_col(B, *it));
	    add(scaled(mat_col(B, *it), r), mat_col(B, *it2));
	  }
	rc.erase(*it);
      }
    
    // cout << "debut ortho un par un" << endl;
    
    while (rc.size()) {
      R nmax = R(0); size_type cmax = size_type(-1);
      for (std::set<size_type>::iterator it=rc.begin(); it!=rc.end(); ++it) {
	R n = vect_norm2(mat_col(B, *it));
	if (nmax < n) { nmax = n; cmax = *it; }
	if (n < EPS) { columns.erase(*it); rc.erase(*it); }
      }
      
      if (nmax < EPS) break;
      
//       cout << "selecting " << cmax << " nmax = " << nmax
// 	   << " nmin = " << nmin << endl;
      
      gmm::scale(mat_col(B, cmax), T(1)/vect_norm2(mat_col(B, cmax)));
      rc.erase(cmax);
      for (std::set<size_type>::iterator it=rc.begin(); it!=rc.end(); ++it) {
	T r = -vect_hp(mat_col(B, *it), mat_col(B, cmax));
	add(scaled(mat_col(B, cmax), r), mat_col(B, *it));
      }
    }
    
    for (std::set<size_type>::iterator it=rc.begin(); it!=rc.end(); ++it)
      columns.erase(*it);

  }


  // Range basis with Gram-Schmidt orthogonalization (dense version)
  // Complex version not verified
  template <typename Mat>
  void range_basis_eff_gram_schmidt_dense(const Mat &B,
					  std::set<size_type> &columns,
					  std::vector<bool> &c_ortho,
					  double EPS) {
    
    typedef std::set<size_type> TAB;
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc_r = columns.size(), nc = mat_ncols(B), nr = mat_nrows(B), i;
    std::set<size_type> rc;
 
    row_matrix< gmm::rsvector<T> > H(nc, nc_r), BB(nr, nc_r);
    std::vector<T> v(nc_r);
    std::vector<size_type> ind(nc_r);
      
    i = 0;
    for (TAB::iterator it = columns.begin(); it != columns.end(); ++it, ++i)
      H(*it, i) = T(1) / vect_norm2(mat_col(B, *it));
     
    mult(B, H, BB);
    dense_matrix<T> M(nc_r, nc_r);
    mult(gmm::conjugated(BB), BB, M);
    
    i = 0;
    for (TAB::iterator it = columns.begin(); it != columns.end(); ++it, ++i)
      if (c_ortho[*it]) {
	gmm::copy(mat_row(M, i), v);
	rank_one_update(M, scaled(v, T(-1)), v);
	M(i, i) = T(1);
      }
      else { rc.insert(i); ind[i] = *it; }

    while (rc.size() > 0) {

      // Next pivot
      R nmax = R(0); size_type imax = size_type(-1);
      for (TAB::iterator it = rc.begin(); it != rc.end(); ++it) {
	R a = M(*it, *it);
	if (a > nmax) { nmax = a; imax = *it; }
	if (a < EPS) { rc.erase(*it); columns.erase(ind[*it]); }
      }

      // cout << "nmax = " << nmax << endl;
      if (nmax < EPS) break;

      // Normalization
      gmm::scale(mat_row(M, imax), T(1) / sqrt(nmax));
      gmm::scale(mat_col(M, imax), T(1) / sqrt(nmax));
      
      // orthogonalization
      copy(mat_row(M, imax), v);
      rank_one_update(M, scaled(v, T(-1)), v);
      M(imax, imax) = T(1);

      rc.erase(imax);
      
    }
    
    // cout << "remaining = " << rc.size() << endl; getchar();
    for (std::set<size_type>::iterator it=rc.begin(); it!=rc.end(); ++it)
      columns.erase(ind[*it]);
  }



  template <typename L> size_type nnz_eps(const L& l, double eps) { 
    typename linalg_traits<L>::const_iterator it = vect_const_begin(l),
      ite = vect_const_end(l);
    size_type res(0);
    for (; it != ite; ++it) if (gmm::abs(*it) >= eps) ++res;
    return res;
  }

  template <typename L>
  bool reserve__rb(const L& l, std::vector<bool> &b, double eps) {
    typename linalg_traits<L>::const_iterator it = vect_const_begin(l),
      ite = vect_const_end(l);
    bool ok = true;
    for (; it != ite; ++it)
      if (gmm::abs(*it) >= eps && b[it.index()]) ok = false;
    if (ok) {
      for (it = vect_const_begin(l); it != ite; ++it)
	if (gmm::abs(*it) >= eps) b[it.index()] = true;
    }
    return ok;
  }

  template <typename Mat>
  void range_basis(const Mat &B, std::set<size_type> &columns,
		   double EPS, col_major) {
   
    typedef typename linalg_traits<Mat>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type nc = mat_ncols(B), nr = mat_nrows(B);

    std::vector<R> norms(nc);
    std::vector<bool> c_ortho(nc), booked(nr);
    std::vector< std::set<size_type> > nnzs(mat_nrows(B));

    R norm_max = R(0);
    for (size_type i = 0; i < nc; ++i) {
      norms[i] = vect_norminf(mat_col(B, i));
      norm_max = std::max(norm_max, norms[i]);
    }

    columns.clear();
    for (size_type i = 0; i < nc; ++i)
      if (norms[i] >= norm_max*EPS) { 
	columns.insert(i);
	nnzs[nnz_eps(mat_col(B, i), EPS * norms[i])].insert(i);
      }
    
    for (size_type i = 1; i < nr; ++i)
      for (std::set<size_type>::iterator it = nnzs[i].begin();
	   it != nnzs[i].end(); ++it)
	if (reserve__rb(mat_col(B, *it), booked, EPS * norms[*it]))
	  c_ortho[*it] = true;

    size_type sizesm[5] = {50, 125, 200, 350, 450};
    for (int k = 0; columns.size() > sizesm[k] && k < 4; ++k) {
      size_type nc_r = columns.size();
      cout << "begin small range basis with " << columns.size()
 	   << " columns, sizesm =  " << sizesm[k] <<  endl;
      std::set<size_type> c1, cres;
      for (std::set<size_type>::iterator it = columns.begin();
	   it != columns.end(); ++it) {
	c1.insert(*it);
	if (c1.size() >= sizesm[k]) {
	  // sizesm = 100 + size_type(gmm::random() * 100);
	  size_type c1size = c1.size();
	  // range_basis_eff_lu(B, c1, c_ortho, EPS);
	  range_basis_eff_gram_schmidt_dense(B, c1, c_ortho, EPS);
	  for (std::set<size_type>::iterator it2=c1.begin(); it2 != c1.end();
	       ++it2) cres.insert(*it2);

	  if (c1.size() == c1size && false) { // a supprimer ?
	    for (size_type i = 0; it != columns.end() && i < 1000; ++it, ++i)
	      cres.insert(*it);
	    if (it != columns.end()) cres.insert(*it);
	  }
	  
	  c1.clear(); 
	}
      }
      // if (c1.size() > 10) range_basis_eff_lu(B, c1, c_ortho, EPS);
      if (c1.size() > 10)
	range_basis_eff_gram_schmidt_dense(B, c1, c_ortho, EPS);
      for (std::set<size_type>::iterator it = c1.begin(); it != c1.end(); ++it)
	cres.insert(*it);
      columns = cres;
      if (columns.size() == nc_r) break;

    }
    cout << "begin global dense range basis for " << columns.size()
 	 << " columns " << endl;

    if (columns.size() > 500)
      range_basis_eff_power(B, columns, EPS);
    else
      range_basis_eff_gram_schmidt_dense(B, columns, c_ortho, EPS);
    // range_basis_eff_lu(B, columns, c_ortho, EPS);

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
    Extract a basis of the range of a (large sparse) matrix selecting some
    column vectors of this matrix. This is in particular usefull to select
    an independent set of linear constraints.

    The algorithm is optimized for two cases :
       - when the (not trivial) kernel is small. An iterativ algorithm
         based on Lanczos method is applied
       - when the (not trivial) kernel is large and most of the dependencies
         can be detected locally. An block Gram-Schmidt is applied first then
         the Lanczos method when the remaining kernel is greatly smaller.
    The LU decomposition has been tested fr local elimination but gives bad
    results : the algorithm is unstable and do not permit to give the right
    number of vector at the end of the process. Moreover, the number of final
    vector depend greatly on the number of vectors in a block of the local
    analysis.
  */
  template <typename Mat>
  void range_basis(const Mat &B, std::set<size_type> &columns,
		   double EPS=1E-12) {
    range_basis(B, columns, EPS,
		typename principal_orientation_type
		<typename linalg_traits<Mat>::sub_orientation>::potype());
  }

}

#endif
