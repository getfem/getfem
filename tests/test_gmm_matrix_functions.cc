/*===========================================================================
 
 Copyright (C) 2014-2015 Konstantinos Poulios.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/


#include "gmm/gmm.h"
#include "gmm/gmm_dense_matrix_functions.h"
#include "gmm/gmm_inoutput.h"


using gmm::size_type;


template <class T> void test_sqrtm(T) {

  size_type n = 4;

  gmm::dense_matrix<T> X(n, n), SQRTMX(n, n), Z(n, n);
  gmm::copy(gmm::identity_matrix(), X);
  gmm::copy(gmm::identity_matrix(), Z);
  X(1,1) = X(2,2) = pow(T(2),-26);
  X(0,3) = T(1);
  Z(1,1) = Z(2,2) = pow(T(2),-13);
  Z(0,3) = T(0.5);

  gmm::sqrtm(X, SQRTMX);

  std::cout << "X = " << X << std::endl;
  std::cout << "sqrtm(X) = " << SQRTMX << std::endl;
  gmm::add(gmm::scaled(Z, -T(1)), SQRTMX);
  std::cout << "err(sqrtm(X)) = " << SQRTMX << std::endl;
  assert(gmm::mat_norm1(SQRTMX) < gmm::default_tol(T()));
}


template <class T> void test_logm(T) {
  typedef typename gmm::number_traits<T>::magnitude_type R;

  size_type n = 2;
  gmm::dense_matrix<T> X(n, n), LOGMX(n, n), Z(n, n);
  gmm::copy(gmm::identity_matrix(), X);
  X(0,1) = -T(1);
  Z(0,1) = -T(1);
  gmm::logm(X, LOGMX);
  std::cout << "X = " << X << std::endl;
  std::cout << "logm(X) = " << LOGMX << std::endl;
  gmm::add(gmm::scaled(Z, -T(1)), LOGMX);
  std::cout << "err(logm(X)) = " << LOGMX << std::endl;
  assert(gmm::mat_norm1(LOGMX) < R(3*n)*gmm::default_tol(T()));
  //assert (norm (logm ([1 -1;0 1]) - [0 -1; 0 0]) < 1e-5)

  //assert (norm (expm (logm ([-1 2 ; 4 -1])) - [-1 2 ; 4 -1]) < 1e-5)

  n = 3;
  gmm::resize(X, n, n); gmm::resize(LOGMX, n, n); gmm::resize(Z, n, n);
  gmm::copy(gmm::identity_matrix(), X);
  gmm::clear(LOGMX); gmm::clear(Z);
  X(0,1) = X(0,2) = X(1,2) = -T(1);
  Z(0,1) = Z(1,2) = -T(1);
  Z(0,2) = -T(1.5);
  gmm::logm(X, LOGMX);
  std::cout << "X = " << X << std::endl;
  std::cout << "logm(X) = " << LOGMX << std::endl;
  gmm::add(gmm::scaled(Z, -T(1)), LOGMX);
  std::cout << "err(logm(X)) = " << LOGMX << std::endl;
  assert(gmm::mat_norm1(LOGMX) < R(3*n)*gmm::default_tol(T()));
  //assert (logm ([1 -1 -1;0 1 -1; 0 0 1]), [0 -1 -1.5; 0 0 -1; 0 0 0], 1e-5)

  //assert (logm (10), log (10))
  //assert (full (logm (eye (3))), logm (full (eye (3))))
  //assert (full (logm (10*eye (3))), logm (full (10*eye (3))), 8*eps)
  //assert (logm (expm ([0 1i; -1i 0])), [0 1i; -1i 0], 10 * eps)
}


int main(void)
{
  srand(1459);

  try {

    test_sqrtm(float());
    test_sqrtm(double());
    test_sqrtm(std::complex<float>());
    test_sqrtm(std::complex<double>());
    
    test_logm(float());
    test_logm(double());
    test_logm(std::complex<float>());
    test_logm(std::complex<double>());
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
