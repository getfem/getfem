/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/
#include <gmm.h>
#include <gmm.h>

// scalar product working also for matrices (to be done in GMM++ ...
template<class VAR> 
typename gmm::linalg_traits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y)
{ return local_sp(X, Y, typename gmm::linalg_traits<VAR>::linalg_type()); }

template<class VAR> 
typename gmm::linalg_traits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y, gmm::abstract_vector)
{ return gmm::vect_sp(X, Y); }

template<class VAR> 
typename gmm::linalg_traits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y, gmm::abstract_matrix) {
  typename gmm::linalg_traits<VAR>::value_type res(0);
  for (gmm::size_type i = 0; i < gmm::mat_nrows(X); ++i) 
    for (gmm::size_type j = 0; j < gmm::mat_ncols(X); ++j)
      res += X(i, j) * Y(i, j);
  return res;
}


// Make a test of the gradient around X.
template <class FUNC, class GRAD, class VAR> 
void test_grad_at(FUNC f, GRAD grad, const VAR &X) {
  
  typedef typename gmm::linalg_traits<VAR>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  VAR Y(X), Z(X), G(X);
  
  grad(X, G);
  T valx = f(X);

  R eps(1), max_ratio(1), ecart, ecart_old, min_ecart(1);
  gmm::fill_random(Z);
  T derdir = local_sp(G, Z), estimate_derdir;
  for (int i = 0; i < 10; ++i, eps /= R(10)) {
    gmm::add(gmm::scaled(Z, eps), X, Y);
    estimate_derdir = (f(Y) - valx) / eps;
    ecart = gmm::abs(derdir - estimate_derdir);
    min_ecart = std::min(ecart, min_ecart);
    // The goal is of course to obtain a clear decreasing sequence
    cout << " " << ecart;
    if (i >= 1)
      if (ecart != T(0)) max_ratio = std::max(max_ratio, ecart_old / ecart);
      else max_ratio = R(10);
    ecart_old = ecart;
  }
  cout << endl;
  if (max_ratio < R(9) && min_ecart > 1E-9) {
    cout << "ERROR, The gradient does not seem to be ok !! max_ratio = "
	 << max_ratio << "\n";
    exit(1);
  }
}

template <class FUNC, class GRAD, class VAR> 
void test_grad(FUNC f, GRAD grad, const VAR &X) {
  VAR Y(X);
  for (long i = 0; i < 10000; ++i) {
    gmm::fill_random(Y);
    // gmm::scale(Y, rand() / 1000 + 1);
    cout << "Expe " << i+1 << " X = " << Y;
    test_grad_at(f, grad, Y);
    cout << endl;
  }
  cout << "The gradient seems to be ok !!\n";
}

//
// Gradient of the Frobenius condition number
//

template <typename MAT, typename MAT2> void
squared_Frobenius_condition_number_gradient(const MAT& M, MAT2& G) { 
  typedef typename gmm::linalg_traits<MAT>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  
  gmm::size_type n = gmm::mat_ncols(M);
  gmm::dense_matrix<T> B(n,n), C(n,n);
  gmm::mult(gmm::transposed(M), M, B);
  R trB = gmm::mat_trace(B);
  gmm::lu_inverse(B);
  R trBinv = gmm::mat_trace(B);
  gmm::mult(B,B,C);
  gmm::mult(gmm::scaled(M, T(-2)*trB), C, G);
  gmm::add(gmm::scaled(M, T(2)*trBinv), G);
}


typedef gmm::dense_matrix<double> DM;

struct func {
  double operator()(const DM &M) { return Frobenius_condition_number_sqr(M); } 
};

struct grad {
  void operator()(const DM &M, DM &G)
    { squared_Frobenius_condition_number_gradient(M, G); }
};

//
// Signed distance for the torus
//

typedef std::vector<double> base_node;
typedef double scalar_type;

struct func2 {
  scalar_type operator()(const base_node &P) const {
    scalar_type R = 2.0, r = 0.5;

    scalar_type x = P[0], y = P[1], z = P[2];
    scalar_type c = sqrt(x*x + y*y);
    if (c == 0.) return R - r;
    return sqrt(gmm::sqr(c-R) + z*z) - r;
  }
};


struct grad2 {
  void operator()(const base_node &P, base_node &G) const {
    gmm::clear(G); 
    scalar_type R = 2.0, r = 0.5;
    scalar_type x = P[0], y = P[1], z = P[2];
    scalar_type c = sqrt(x*x + y*y);
    if (c == 0.) return;
    scalar_type w = 1. - R / c;
    scalar_type e = sqrt(gmm::sqr(c-R) + z*z);
    if (e == 0) return;
    G[0] = x * w / e;
    G[1] = y * w / e;
    G[2] = z / e;
  }
};

int main(void) {

  test_grad(func(), grad(), DM(5, 5));
  test_grad(func2(), grad2(), base_node(3));

  return 0;
}
