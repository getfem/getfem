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
// g++ -I ../../src -O3 ../../tests/test_superlu.C -o test_superlu superlu.a -I ~/source++/ -DGMM_USES_SUPERLU


// options d'optimisations avec g++ :
//  -funroll-all-loops -ffast-math -fstrict-aliasing -fomit-frame-pointer

// #define GMM_USES_LAPACK
#include <ftool.h>
#include <gmm.h>

using gmm::size_type;

template <class T> void test_with(T) {
  size_type n = 50;

  gmm::row_matrix<gmm::wsvector<T> > A(n, n), B(n, n), C(n, n);
  std::vector<T> x(n), y(n), z(n);
  
  gmm::copy(gmm::identity_matrix(), A);
  gmm::fill_random(A, 1.0);
  gmm::fill_random(B);
  gmm::fill_random(x);
  gmm::fill_random(y);

  A(0,1) = 0;
  A(1,2) = 0;
  A(2,4) = 0;
  A(3,0) = 0;
  A(4,1) = 0;

  gmm::SuperLU_solve(A, x, y);

  // gmm::lu_solve(A, z, y);

  cout << "y = " << y << endl;
  cout << "x = " << x << endl;
  // cout << "z = " << z << endl;
  gmm::mult(A, x, y);
  cout << "Ax = " << y << endl;
  // gmm::mult(A, z, y);
  // cout << "Az = " << y << endl;
}

int main(void)
{
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb);

  srand(1459);

# if defined(GMM_USES_SUPERLU)
  cout << "Trying using SuperLU\n";
# else
  cout << "Not using SuperLU\n";
# endif

  try {

    cout << "sizeof(int) = " << sizeof(int)
	 << " sizeof(long) = " << sizeof(long) << endl;
    
    // test_with(float());
    test_with(double());
    // test_with(std::complex<float>());
    // test_with(std::complex<double>());
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
