/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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
#include <bgeot_matrix.h>
#include <gmm.h>

int main(void)
{
  try {

    cout.precision(16);

    cout << "/***********************************************************/\n";
    cout << "/*                   Test of fsmatrix                      */\n";
    cout << "/***********************************************************/\n";

    bgeot::fsmatrix<double, 10> m;
    bgeot::fsvector<double, 10> y, x, b;
    gmm::copy(gmm::identity_matrix(), m);
    int j = 5, k = 12;
    for (int i = 0; i < 10; ++i) { 
      m(i, j) = double(k);
      j = (j + 6) % 10; k = (k + 15) % 31; 
      m(i, j) = double(k);
      j = (j + 6) % 10; k = (k + 15) % 31;
      x[i] = 10.0 * double(i - 5 + ((i >= 5) ? 1 : 0));
    }
    cout << "m = " << m << endl;
    cout << "x = " << x << endl;

    gmm::mult(m, x, b);
    cout << "b = " << b << endl;
    gmm::clear(y);
    gmm::bicgstab(m, y, b, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y = " << y << endl;
    gmm::add(gmm::scaled(x, double(-1.0)), y);
    cout << "y = " << y << endl;
    double error = gmm::vect_norm2(y);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);
    
    cout << "/***********************************************************/\n";
    cout << "/*                   Test of vsmatrix                      */\n";
    cout << "/***********************************************************/\n";

    bgeot::vsmatrix<double> m2(10,10);
    gmm::copy(m, m2);
    cout << "m2 = " << m2 << endl;
    std::vector<double> y2(10), b2(10);
    gmm::copy(b, b2);
    gmm::clear(y2);
    gmm::bicgstab(m2, y2, b2, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);
    
    cout << "transposed(m2) = " << gmm::transposed(m2) << endl;
    cout << "transposed(transposed(m2)) = "
	 << gmm::transposed(gmm::transposed(m2)) << endl;
    
    gmm::transposed(gmm::transposed(m2))(3, 1) = 5.0;
    if (m2(3, 1) != 5.0) 
      DAL_THROW(dal::failure_error, "write error on matrix m2.");
    gmm::transposed(gmm::transposed(m2))(3, 1) = 0.0;

    cout << "/***********************************************************/\n";
    cout << "/*             Test of row_matrix<wsvector>                */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<gmm::wsvector<double> > m3(10, 10);
    gmm::copy(m, m3);
    cout << "m3 = " << m3 << endl;
    cout << "transposed(m3) = " << gmm::transposed(m3) << endl;
    gmm::clear(y2);
    gmm::bicgstab(m3, y2, b, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);
    
    cout << "/***********************************************************/\n";
    cout << "/*         Test of row_matrix<std::vector>                 */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<std::vector<double> > m4(10, 10);
    gmm::copy(m, m4);
    cout << "m4 = " << m4 << endl;
    cout << "transposed(m4) = " << gmm::transposed(m4) << endl;
    gmm::clear(y2);
    gmm::bicgstab(m4, y2, b, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);
    
    cout << "/***********************************************************/\n";
    cout << "/*             Test of col_matrix<wsvector>                */\n";
    cout << "/***********************************************************/\n";

    gmm::col_matrix<gmm::wsvector<double> > m5(10, 10);
    gmm::copy(m3, m5);
    cout << "m5 = " << m5 << endl;
    cout << "transposed(m5) = " << gmm::transposed(m5) << endl;
    gmm::clear(y2);
    gmm::bicgstab(m5, y2, b, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    cout << "/***********************************************************/\n";
    cout << "/*             Test of row_matrix<slvector>                */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<gmm::slvector<double> > m5bis(10, 10);
    gmm::copy(m3, m5bis);
    cout << "m5bis = " << m5bis << endl;
    cout << "transposed(m5bis) = " << gmm::transposed(m5bis) << endl;
    gmm::clear(y2);
    gmm::bicgstab(m5bis, y2, b, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    cout << "/***********************************************************/\n";
    cout << "/*         Test of sub_matrices of fsmatrix                */\n";
    cout << "/***********************************************************/\n";

    gmm::sub_interval sint1(10, 10), sint2(15, 10);
    gmm::sub_slice ssli1(10, 10, 3), ssli2(15, 10, 2);
    std::vector<gmm::size_type> ind1(10), ind2(10);
    ind1[0] = 1; ind1[1] = 5; ind1[2] = 3; ind1[3] = 25; ind1[4] = 15; 
    ind1[5] = 2; ind1[6] = 8; ind1[7] = 6; ind1[8] = 12; ind1[9] = 16; 
    ind2[0] = 0; ind2[1] = 2; ind2[2] = 5; ind2[3] = 17; ind2[4] = 24; 
    ind2[5] = 1; ind2[6] = 3; ind2[7] = 6; ind2[8] = 10; ind2[9] = 12; 
    gmm::sub_index sind1(ind1), sind2(ind2);

    bgeot::fsmatrix<double, 38> m6;
    gmm::clear(m6);
    gmm::copy(m, gmm::sub_matrix(m6, sint1, sint2));
    cout << "m6 = " << m6 << endl;
    cout << "gmm::sub_matrix(m6, sint1, sint2)  = "
//           << gmm::transposed(gmm::transposed
// 			     (gmm::sub_matrix(m6,sint2, sint1))) << endl;
	 << gmm::sub_matrix(m6,sint1, sint2) << endl;
    gmm::clear(y2);
    gmm::bicgstab(gmm::sub_matrix(m6, sint1, sint2), y2, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    gmm::copy(m3, gmm::sub_matrix(m6, ssli1, ssli2));
    cout << "m6 = " << m6 << endl;
    cout << "gmm::sub_matrix(m6, ssli1, ssli2)  = "
	 << gmm::sub_matrix(m6, ssli1, ssli2) << endl;
    gmm::clear(y2);
    gmm::bicgstab(gmm::sub_matrix(m6, ssli1, ssli2), y2, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    gmm::copy(m3, gmm::sub_matrix(m6, sint1, ssli2));
    cout << "m6 = " << m6 << endl;
    cout << "gmm::sub_matrix(m6, sint1, ssli2)  = "
	 << gmm::sub_matrix(m6, sint1, ssli2) << endl;
    gmm::clear(y2);
    gmm::bicgstab(gmm::sub_matrix(m6, sint1, ssli2), y2, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y2 = " << y2 << endl;
    gmm::add(gmm::scaled(x, -1.0), y2);
    error = gmm::vect_norm2(y2);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    gmm::copy(m4, gmm::sub_matrix(m6, sind1, sind2));
    gmm::wsvector<double> y3(10);
    cout << "m6 = " << m6 << endl;
    cout << "gmm::sub_matrix(m6, sind1, sind2)  = "
	 << gmm::sub_matrix(m6, sind1, sind2) << endl;
    gmm::clear(y3);
    gmm::bicgstab(gmm::sub_matrix(m6, sind1, sind2), y3, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y3 = " << y3 << endl;
    gmm::add(gmm::scaled(x, -1.0), y3);
    error = gmm::vect_norm2(y3);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    cout << "/***********************************************************/\n";
    cout << "/*      Test of sub_matrices of row_matrix<wsvector>       */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<gmm::wsvector<double> > m7(38, 40);
    gmm::clear(m7);
    gmm::copy(m3, gmm::sub_matrix(m7, sind1, sind2));
    cout << "m7 = " << m7 << endl;
    cout << "gmm::transposed(gmm::sub_matrix(m7, sind1, sind2))  = "
	 << gmm::transposed(gmm::sub_matrix(m7, sind1, sind2)) << endl;
    gmm::clear(y);
    gmm::bicgstab(gmm::sub_matrix(m7, sind1, sind2), y, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y = " << y << endl;
    gmm::add(gmm::scaled(x, -1.0), y);
    error = gmm::vect_norm2(y);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    gmm::copy(m, gmm::sub_matrix(m7, sint1, ssli2));
    cout << "m7 = " << m7 << endl;
    cout <<
      "gmm::transposed(gmm::transposed(gmm::sub_matrix(m7, sint1, ssli2))  =\n"
	 << gmm::transposed(gmm::transposed(gmm::sub_matrix(m7, sint1, ssli2)))
	 << endl;
    gmm::clear(y3);
    gmm::bicgstab(gmm::sub_matrix(m7, sint1, ssli2), y3, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y3 = " << y3 << endl;
    gmm::add(gmm::scaled(x, -1.0), y3);
    error = gmm::vect_norm2(y3);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    cout << "/***********************************************************/\n";
    cout << "/*      Test of sub_matrices of row_matrix<rsvector>       */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<gmm::rsvector<double> > m8(40, 40);
    gmm::clear(m8);
    gmm::copy(m3, gmm::sub_matrix(m8, sind1, sind2)); //essayer m3 par la suite
    cout << "m8 = " << m8 << endl;
    cout << "gmm::transposed(gmm::sub_matrix(m8, sind1, sind2))  = "
	 << gmm::transposed(gmm::sub_matrix(m8, sind1, sind2)) << endl;
    gmm::clear(y);
    gmm::bicgstab(gmm::sub_matrix(m8, sind1, sind2), y, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y = " << y << endl;
    gmm::add(gmm::scaled(x, -1.0), y);
    error = gmm::vect_norm2(y);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    gmm::copy(m, gmm::sub_matrix(m8, sint1, ssli2));
    cout << "m8 = " << m8 << endl;
    cout <<
      "gmm::transposed(gmm::transposed(gmm::sub_matrix(m8, sint1, ssli2))  =\n"
	 << gmm::transposed(gmm::transposed(gmm::sub_matrix(m8, sint1, ssli2)))
	 << endl;
    gmm::clear(y3);
    gmm::bicgstab(gmm::sub_matrix(m8, sint1, ssli2), y3, b,
		  gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y3 = " << y3 << endl;
    gmm::add(gmm::scaled(x, -1.0), y3);
    error = gmm::vect_norm2(y3);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);

    cout << "/***********************************************************/\n";
    cout << "/*      matrix-matrix multiplication                       */\n";
    cout << "/***********************************************************/\n";

    gmm::row_matrix<gmm::rsvector<double> > m9(10, 10);
    gmm::mult(gmm::transposed(gmm::sub_matrix(m8, sint1, ssli2)), 
	      gmm::sub_matrix(m8, sint1, ssli2), m9);
    cout << "m9 = " << m9 << endl;
    gmm::mult(gmm::transposed(gmm::sub_matrix(m8, sint1, ssli2)), b, y);
    gmm::cg(m9, y3, y, gmm::identity_matrix(), 1000, 1E-16, 0);
    cout << "y3 = " << y3 << endl;
    gmm::add(gmm::scaled(x, -1.0), y3);
    error = gmm::vect_norm2(y3);
    cout << "Error : " << error << endl;
    if (error > 1.0E-10)
      DAL_THROW(dal::failure_error, "computation error too large : " << error);
    

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
