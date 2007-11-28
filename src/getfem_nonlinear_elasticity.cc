// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
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
//===========================================================================


#include "getfem/getfem_nonlinear_elasticity.h"

namespace getfem {


  int check_symmetry(const base_tensor &t) {
    int flags = 7; size_type N = 3;
    for (size_type n = 0; n < N; ++n)
      for (size_type m = 0; m < N; ++m)
	for (size_type l = 0; l < N; ++l)
	  for (size_type k = 0; k < N; ++k) {
	    if (gmm::abs(t(n,m,l,k) - t(l,k,n,m))>1e-10) flags &= (~1); 
	    if (gmm::abs(t(n,m,l,k) - t(m,n,l,k))>1e-10) flags &= (~2); 
	    if (gmm::abs(t(n,m,l,k) - t(n,m,k,l))>1e-10) flags &= (~4);
	  }
    return flags;
  }

  void abstract_hyperelastic_law::random_E(base_matrix &E) {
    size_type N = gmm::mat_nrows(E);
    base_matrix Phi(N,N); gmm::fill_random(Phi);
    gmm::mult(gmm::transposed(Phi),Phi,E);
    gmm::scale(E,-1.); gmm::add(gmm::identity_matrix(),E); 
    gmm::scale(E,-0.5);
  }

  void abstract_hyperelastic_law::test_derivatives
  (size_type N, scalar_type h, const base_vector& param) const {
    base_matrix E(N,N), E2(N,N), DE(N,N); 
    random_E(E); random_E(DE);
    gmm::scale(DE,h);
    gmm::add(E,DE,E2);
    
    base_matrix sigma1(N,N), sigma2(N,N);
    getfem::base_tensor tdsigma(N,N,N,N);
    base_matrix dsigma(N,N);
    gmm::copy(E,E2); gmm::add(DE,E2);
    sigma(E, sigma1, param); sigma(E2, sigma2, param);
    
    scalar_type d = strain_energy(E2, param) - strain_energy(E, param);
    scalar_type d2 = 0;
    for (size_type i=0; i < N; ++i) 
      for (size_type j=0; j < N; ++j) d2 += sigma1(i,j)*DE(i,j);
    if (gmm::abs(d-d2) > h*1e-5) 
      cout << "wrong derivative of strain_energy, d=" << d
	   << ", d2=" << d2 << "\n";
    
    grad_sigma(E,tdsigma,param);
    for (size_type i=0; i < N; ++i) {
      for (size_type j=0; j < N; ++j) {
	dsigma(i,j) = 0;
	for (size_type k=0; k < N; ++k) {
	  for (size_type m=0; m < N; ++m) {
	    dsigma(i,j) += tdsigma(i,j,k,m)*DE(k,m);
	  }
	}
	sigma2(i,j) -= sigma1(i,j);
	if (gmm::abs(dsigma(i,j) - sigma2(i,j)) > h*1e-5) {
	  cout << "wrong derivative of sigma, i=" << i << ", j=" 
	       << j << ", dsigma=" << dsigma(i,j) << ", var sigma = " 
	       << sigma2(i,j) << "\n";
	}
      }
    }
  }
    
  scalar_type SaintVenant_Kirchhoff_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    return gmm::sqr(gmm::mat_trace(E)) * params[0] / scalar_type(2)
      + gmm::mat_euclidean_norm_sqr(E) * params[1];
  }
  
  void SaintVenant_Kirchhoff_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, params[0] * gmm::mat_trace(E));
    gmm::add(gmm::scaled(E, 2 * params[1]), result);
  }
  void SaintVenant_Kirchhoff_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    std::fill(result.begin(), result.end(), scalar_type(0));
    size_type N = gmm::mat_nrows(E);
    for (size_type i = 0; i < N; ++i)
      for (size_type l = 0; l < N; ++l) {
	result(i, i, l, l) = params[0];
	result(i, l, i, l) += params[1];
	  result(i, l, l, i) += params[1];
      }
  }


  scalar_type Mooney_Rivlin_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    scalar_type C1 = params[0], C2 = params[1];
    return scalar_type(2) *
      (gmm::mat_trace(E) * (C1 + scalar_type(2)*C2)
       + C2*(gmm::sqr(gmm::mat_trace(E)) - gmm::mat_euclidean_norm_sqr(E)));
  }

  void Mooney_Rivlin_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    scalar_type C12 = scalar_type(2) * params[0];
    scalar_type C24 = scalar_type(4) * params[1];
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, C24*(gmm::mat_trace(E)+scalar_type(1)) + C12);
    gmm::add(gmm::scaled(E, -C24), result);
  }
  void Mooney_Rivlin_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    scalar_type C22 = scalar_type(2) * params[1];
    std::fill(result.begin(), result.end(), scalar_type(0));
    size_type N = gmm::mat_nrows(E);
    for (size_type i = 0; i < N; ++i)
      for (size_type l = 0; l < N; ++l) {
	result(i, i, l, l) = scalar_type(2) * C22;
	result(i, l, i, l) -= C22;
	result(i, l, l, i) -= C22;
      }
  }
 
  scalar_type Ciarlet_Geymonat_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type a = params[1] + params[2] / scalar_type(2);
    scalar_type b = -(params[1] + params[2]) / scalar_type(2);
    scalar_type c = params[0]/scalar_type(4)  - b;
    scalar_type d = params[0]/scalar_type(2) + params[1];
    //scalar_type d = params[0] - scalar_type(2)*params[2] - scalar_type(4)*b;
    scalar_type e = -(scalar_type(3)*(a+b) + c);
    base_matrix C(N, N);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    scalar_type det = gmm::lu_det(C);
    return a * gmm::mat_trace(C)
      + b * (gmm::sqr(gmm::mat_trace(C)) - 
	     gmm::mat_euclidean_norm_sqr(C))/scalar_type(2)
      + c * det - d * log(det) / scalar_type(2) + e;
  }
  void Ciarlet_Geymonat_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type a = params[1] + params[2] / scalar_type(2);
    scalar_type b = -(params[1] + params[2]) / scalar_type(2);
    scalar_type c = params[0]/scalar_type(4)  - b;
    scalar_type d = params[0]/scalar_type(2) + params[1]; 
    //d=params[0] - scalar_type(2)*params[2] - scalar_type(4)*b;
    base_matrix C(N, N);
    assert(gmm::abs(2*a+4*b+2*c-d)<1e-5);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, scalar_type(2) * (a + b * gmm::mat_trace(C)));
    gmm::add(gmm::scaled(C, -scalar_type(2) * b), result);
    scalar_type det = gmm::lu_inverse(C);
    gmm::add(gmm::scaled(C, scalar_type(2) * c * det - d), result);
  }
  void Ciarlet_Geymonat_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type b2 = -(params[1] + params[2]); // b * 2
    scalar_type c = (params[0]  - 2*b2) / scalar_type(4);
    //scalar_type d = params[0] - scalar_type(2)*params[2] - 2*b2;
    scalar_type d = params[0]/scalar_type(2) + params[1]; 
    base_matrix C(N, N);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    scalar_type det = gmm::lu_inverse(C);
    std::fill(result.begin(), result.end(), scalar_type(0));
    for (size_type i = 0; i < N; ++i)
      for (size_type j = 0; j < N; ++j) {
	result(i, i, j, j) += 2*b2;
	result(i, j, i, j) -= b2;
	result(i, j, j, i) -= b2;
	for (size_type  k = 0; k < N; ++k)
	  for (size_type  l = 0; l < N; ++l)
	    result(i, j, k, l) += 
	      (C(i, k)*C(l, j) + C(i, l)*C(k, j)) * (d-scalar_type(2)*det*c)
	      + (C(i, j) * C(k, l)) * det*c*scalar_type(4);
      }
  }

}  /* end of namespace getfem.                                             */

