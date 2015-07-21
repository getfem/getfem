/*===========================================================================
 
 Copyright (C) 2006-2015 Yves Renard, Julien Pommier.
 
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
// Little program which computes the base functions of the Argyris element
// on the reference element.

#include <iostream>
#include <gmm.h>
#include <getfem_config.h>

/// Print P to the output stream o. for instance cout << P;
template<typename T> void poly_cpp_display(std::ostream &o,
					   const bgeot::polynomial<T>& P) { 
  bool first = true; unsigned n = 0;
  typename bgeot::polynomial<T>::const_iterator it = P.begin(), ite = P.end();
  bgeot::power_index mi(P.dim());
  if (it != ite && *it != T(0))
    { o << *it; first = false; ++it; ++n; ++mi; }
  for ( ; it != ite ; ++it, ++mi ) {
    if (*it != T(0)) {
      if (!first) { if (*it < T(0)) o << " - "; else o << " + "; }
      else if (*it < T(0)) o << "-";
      if (gmm::abs(*it)!=T(1)) o << gmm::abs(*it);
      for (unsigned j = 0; j < P.dim(); ++j)
	if (mi[j] != 0) {
	  if (j != 0 || gmm::abs(*it) != T(1)) o << "*";
	  for (unsigned k=0; k < mi[j]; ++k) {
	    if (k) o << "*"; o << "xyz"[j];
	  }
	}
      first = false; ++n;
    }
  }
  if (n == 0) o << "0";
}


int main(void) {

  try {
    bgeot::base_poly one(3, 0), x(3, 1, 0), y(3, 1, 1), z(3, 1, 2); one.one();
    bgeot::base_poly base[20];
    // base for P5
    base[ 0] = one;
    base[ 1] = x;
    base[ 2] = y;
    base[ 3] = z;
    base[ 4] = x*x;
    base[ 5] = x*y;
    base[ 6] = x*z;
    base[ 7] = y*y;
    base[ 8] = y*z;
    base[ 9] = z*z;
    base[10] = x*x*x;
    base[11] = x*x*y;
    base[12] = x*x*z;
    base[13] = x*y*y;
    base[14] = x*y*z;
    base[15] = x*z*z;
    base[16] = y*y*y;
    base[17] = y*y*z;
    base[18] = y*z*z;
    base[19] = z*z*z;
    

    bgeot::base_matrix M(20, 20);
    
    for (int i = 0; i < 20; ++i) {
      bgeot::base_poly p = base[i], q;
      M( 0, i) = p.eval(bgeot::base_node(0.0, 0.0, 0.0).begin());
      M( 1, i) = p.eval(bgeot::base_node(1.0, 0.0, 0.0).begin());
      M( 2, i) = p.eval(bgeot::base_node(0.0, 1.0, 0.0).begin());
      M( 3, i) = p.eval(bgeot::base_node(0.0, 0.0, 1.0).begin());
      
      double u_3 = 1.0 / 3.0;
      
      M( 4, i) = p.eval(bgeot::base_node(u_3, u_3, u_3).begin());
      M( 5, i) = p.eval(bgeot::base_node(0.0, u_3, u_3).begin());
      M( 6, i) = p.eval(bgeot::base_node(u_3, 0.0, u_3).begin());
      M( 7, i) = p.eval(bgeot::base_node(u_3, u_3, 0.0).begin());
      
      q = p; q.derivative(0);
      M( 8, i) = q.eval(bgeot::base_node(0.0, 0.0, 0.0).begin());
      M( 9, i) = q.eval(bgeot::base_node(1.0, 0.0, 0.0).begin());
      M(10, i) = q.eval(bgeot::base_node(0.0, 1.0, 0.0).begin());
      M(11, i) = q.eval(bgeot::base_node(0.0, 0.0, 1.0).begin());
      
      q = p; q.derivative(1);
      M(12, i) = q.eval(bgeot::base_node(0.0, 0.0, 0.0).begin());
      M(13, i) = q.eval(bgeot::base_node(1.0, 0.0, 0.0).begin());
      M(14, i) = q.eval(bgeot::base_node(0.0, 1.0, 0.0).begin());
      M(15, i) = q.eval(bgeot::base_node(0.0, 0.0, 1.0).begin());
      
      q = p; q.derivative(2);
      M(16, i) = q.eval(bgeot::base_node(0.0, 0.0, 0.0).begin());
      M(17, i) = q.eval(bgeot::base_node(1.0, 0.0, 0.0).begin());
      M(18, i) = q.eval(bgeot::base_node(0.0, 1.0, 0.0).begin());
      M(19, i) = q.eval(bgeot::base_node(0.0, 0.0, 1.0).begin());
    }
    
    gmm::clean(M, 1E-10);
    cout << "M = " << M << endl;
    
    gmm::lu_inverse(M);
    
    gmm::clean(M, 1E-10);
    cout << "inv M = " << M << endl;
    
    cout.precision(13);
    
    for (int i = 0; i < 20; ++i) {
      bgeot::base_poly p(3,3);
      for (int j = 0; j < 20; ++j)
	if (gmm::abs(M(j, i)) > 1E-8) p += base[j]*M(j, i);
      
      cout << "\\hat{\\varphi}_{" << i << "}(x,y) = " << p << ",\\\\" << endl;
    }
    
    for (int i = 0; i < 20; ++i) {
      bgeot::base_poly p(3,3);
      for (int j = 0; j < 20; ++j)
	if (gmm::abs(M(j, i)) > 1E-8) p += base[j]*M(j, i);
      
      cout << "base_[" << i << "]=";
      poly_cpp_display(cout, p);
      cout << ";\n";
    }

  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
