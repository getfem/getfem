/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard.
 
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

// à compiler avec la ligne de commande pour lapack/blas
// g++ -I ../../src -O3 ../../tests/test_gmm_lapack.C -o test_gmm_lapack -llapack -lblas -lg2c

// à compiler avec la ligne de commande pour atlas
// g++ -I ../../src -O3 ../../tests/test_gmm_lapack.C -o test_gmm_lapack /usr/lib/atlas/liblapack.a /usr/lib/atlas/libblas.a -latlas  -lg2c

// options d'optimisations avec g++ :
//  -funroll-all-loops -ffast-math -fstrict-aliasing -fomit-frame-pointer

// pour qd ou dd :
// /home/gmmpc15/renard/usr/pc_g++/lib/libqd.a 
// #define NO_INLINE
// #include <qd.h> 
// #include <dd.h>
// #include <x86.h>

#define NO_INLINE
#include <dd.h>
#include <qd.h>
#include <x86.h>

// #define GMM_USES_LAPACK
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"


using gmm::size_type;

template<class MAT> void my_mult(const MAT &A, const MAT &B, MAT &C) {
  gmm::mult(gmm::conjugated(A), gmm::conjugated(B), C);
  
}

template <class T> void test_with(T) {
  size_type n = 7;

  gmm::dense_matrix<T> A(n, n), B(n, n), C(n, n);
  std::vector<T> x(n), y(n), z(n);
  
  gmm::fill_random(A);
  gmm::fill_random(B);
  gmm::fill_random(x);
  gmm::fill_random(y);

  gmm::lu_solve(A, x, y);
  gmm::mult(A, x, gmm::scaled(y, T(-1)), z);
  cout << "z = " << z << endl;
  
  
  double exectime = dal::uclock_sec();
  implicit_qr_algorithm(A, x, C, 1E-10);

  // gmm::mult(A, x, gmm::scaled(y, T(-1)), z);
  // cout << "z = " << z << endl;

  cout << "A = " << A << endl;
  cout << "x = " << x << endl;
  cout << "C = " << C << endl;
  // my_mult(A, B, C);

  gmm::mult(C, conjugated(C), B);
  cout << "B = " << B << endl;
  cout << "cpu time = " << dal::uclock_sec() - exectime << endl;
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

    unsigned short old_cw;
    x86_fix_start(&old_cw);
    
    // test_with(float());
    // test_with(double());
    // test_with(std::complex<float>());
    // test_with(std::complex<double>());

    dd_real a = "1.23456789012345678901234567890123456789";

    cout << "a dd-real : " << a << endl;

    test_with(qd_real());
    test_with(std::complex<double>());
    test_with(std::complex<qd_real>());

    x86_fix_end(&old_cw);
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
