/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2003  Yves Renard.                                   */
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

// à compiler avec la ligne de commande pour lapack/blas
// g++ -I ../../src -O3 ../../tests/test_gmm_lapack.C -o test_gmm_lapack -llapack -lblas -lg2c

// à compiler avec la ligne de commande pour atlas
// g++ -I ../../src -O3 ../../tests/test_gmm_lapack.C -o test_gmm_lapack /usr/lib/atlas/liblapack.a /usr/lib/atlas/libblas.a -latlas  -lg2c

// options d'optimisations avec g++ :
//  -funroll-all-loops -ffast-math -fstrict-aliasing -fomit-frame-pointer

// #define GMM_USES_LAPACK
#include <ftool.h>
#include <gmm.h>
#include <gmm_inoutput.h>


using gmm::size_type;

template<class MAT> void my_mult(const MAT &A, const MAT &B, MAT &C) {
  gmm::mult(gmm::conjugated(A), gmm::conjugated(B), C);
}

template <class T> void test_with(T) {
  size_type n = 5;

  gmm::dense_matrix<T> A(n, n), B(n, n), C(n, n);
  std::vector<T> x(n), y(n), z(n);
  
  gmm::fill_random(A);
  gmm::fill_random(B);
  gmm::fill_random(x);
  gmm::fill_random(y);
  
  
  double exectime = ftool::uclock_sec();
  implicit_qr_algorithm(A, x, C);
  

  cout << "A = " << A << endl;
  cout << "x = " << x << endl;
  cout << "C = " << C << endl;
  // my_mult(A, B, C);
  cout << "cpu time = " << ftool::uclock_sec() - exectime << endl;
  // cout << "col(B,2) = " << gmm::mat_const_col(B,2) << endl;
  // cout << "col(C,2) = " << gmm::mat_const_col(C,2) << endl;
  
}

int main(void)
{
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb);

  srand(1459);

# if defined(GMM_USES_LAPACK) || defined(GMM_USES_ATLAS)
  cout << "Trying using Lapack\n";
# else
  cout << "Not using Lapack\n";
# endif

  try {
    
    //    test_with(float());
    //    test_with(double());
    //    test_with(std::complex<float>());
    test_with(std::complex<float>());
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
