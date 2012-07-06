/*===========================================================================
 
 Copyright (C) 2006-2012 Yves Renard, Julien Pommier.
 
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

using bgeot::size_type;

template<typename T> bool recognize_frac(T a, int &i, int &j) 
 {
  for (i=1; i < 100; ++i) {
    for (j=1; j < 100; ++j) {
      if (gmm::abs(a - double(i)/j)<1e-12) {
	return true;
      }
    }
  }
  return false;
}

template<typename T> void print_const(std::ostream &o, T a) {
  if (a < 0) o << "-";
  a = gmm::abs(a);

  if (gmm::abs(a - int(a)) < 1e-12) { o << a; return; }

  int ii, jj;
  for (unsigned k=1; k < 100; ++k) {
    if (recognize_frac(a/sqrt(k), ii, jj)) {
      bool m=false;
      if (k != 1) { 
	o << "sqrt(" << k << ")"; 
	if (ii != 1 || jj != 1) o << "*"; else return;
      }
      o << ii;
      if (jj != 1) o << "/" << jj;
      return;
    }
  }
  o << a;
}

template<typename T> void spec_print(std::ostream &o,
				     const bgeot::polynomial<T>& P) { 
  bool first = true; size_type n = 0;
  typename bgeot::polynomial<T>::const_iterator it = P.begin(), ite = P.end();
  bgeot::power_index mi(P.dim());
  if (it != ite && *it != T(0))
    {  print_const(o, *it); first = false; ++it; ++n; ++mi; }
  for ( ; it != ite ; ++it, ++mi ) {
    if (*it != T(0)) {
      bool first_var = true;
      if (!first) { if (*it < T(0)) o << " - "; else o << " + "; }
      else if (*it < T(0)) o << "-";
      if (gmm::abs(gmm::abs(*it) - 1) > 1E-14) {
	print_const(o, gmm::abs(*it));
	first_var = false;
      }
      for (size_type j = 0; j < P.dim(); ++j)
	if (mi[j] != 0) {
	  if (!first_var) o << "*"; first_var = false;
	  if (P.dim() <= 7) o << "xyzwvut"[j];
	  else o << "x_" << j; 
	  if (mi[j] > 1) o << "^" << mi[j];
	}
      first = false; ++n;
    }
  }
  if (n == 0) o << "0";
}



int main(void) {

  try {
    bgeot::base_poly one(2, 0), x(2, 1, 0), y(2, 1, 1); one.one();
    bgeot::base_poly base[20];
    // base for P5
    base[0] = one;
    base[1] = x;
    base[2] = y;
    base[3] = x*x;
    base[4] = x*y;
    base[5] = y*y;
    base[6] = x*x*x;
    base[7] = x*x*y;
    base[8] = x*y*y;
    base[9] = y*y*y;
    

    bgeot::base_matrix M(30, 30);
    
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 10; ++i) {
	bgeot::base_poly p = base[i], q, q2;
	if (j == 1)
	  M( 0, i+10*j) = p.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M( 1, i+10*j) = p.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 1)
	  M( 2, i+10*j) = p.eval(bgeot::base_node(0.0, 1.0).begin());

	q = p; q.derivative(0);
	if (j == 1)
	  M( 3, i+10*j) = q.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M( 4, i+10*j) = q.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 1)
	  M( 5, i+10*j) = q.eval(bgeot::base_node(0.0, 1.0).begin());

	q = p; q.derivative(1);
	if (j == 1)
	  M( 6, i+10*j) = q.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M( 7, i+10*j) = q.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 1)
	  M( 8, i+10*j) = q.eval(bgeot::base_node(0.0, 1.0).begin());

	//
	// relations internes
	//

	// derivées normales en P1
	q = p; q.derivative(1); q2 = p; q2.derivative(0); q += q2;
	q /= sqrt(2);
	if (j == 0)
	  M( 9, i+10*j) = q.eval(bgeot::base_node(0.5, 0.5).begin())
	    - (q.eval(bgeot::base_node(0.0, 1.0).begin())
	       + q.eval(bgeot::base_node(1.0, 0.0).begin())) / 2.0;
	
	q = p; q.derivative(0);
	if (j == 1)
	  M(10, i+10*j) = q.eval(bgeot::base_node(0.0, 0.5).begin())
	    - (q.eval(bgeot::base_node(0.0, 0.0).begin())
	       + q.eval(bgeot::base_node(0.0, 1.0).begin())) / 2.0;
	
	q = p; q.derivative(1);
	if (j == 2)
	  M(11, i+10*j) = q.eval(bgeot::base_node(0.5, 0.0).begin())
	    - (q.eval(bgeot::base_node(0.0, 0.0).begin())
	       + q.eval(bgeot::base_node(1.0, 0.0).begin())) / 2.0;

	// raccord en (0.0, 0.0)
	if (j == 1)
	  M(12, i+10*j) = p.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M(12, i+10*j) = -p.eval(bgeot::base_node(0.0, 0.0).begin());
	q = p; q.derivative(0);
	if (j == 1)
	  M(13, i+10*j) = q.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M(13, i+10*j) = -q.eval(bgeot::base_node(0.0, 0.0).begin());
	q = p; q.derivative(1);
	if (j == 1)
	  M(14, i+10*j) = q.eval(bgeot::base_node(0.0, 0.0).begin());
	if (j == 2)
	  M(14, i+10*j) = -q.eval(bgeot::base_node(0.0, 0.0).begin());

	// raccord en (1.0, 0.0)
	if (j == 0)
	  M(15, i+10*j) = p.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 2)
	  M(15, i+10*j) = -p.eval(bgeot::base_node(1.0, 0.0).begin());
	q = p; q.derivative(0);
	if (j == 0)
	  M(16, i+10*j) = q.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 2)
	  M(16, i+10*j) = -q.eval(bgeot::base_node(1.0, 0.0).begin());
	q = p; q.derivative(1);
	if (j == 0)
	  M(17, i+10*j) = q.eval(bgeot::base_node(1.0, 0.0).begin());
	if (j == 2)
	  M(17, i+10*j) = -q.eval(bgeot::base_node(1.0, 0.0).begin());

	// raccord en (0.0, 1.0)
	if (j == 0)
	  M(18, i+10*j) = p.eval(bgeot::base_node(0.0, 1.0).begin());
	if (j == 1)
	  M(18, i+10*j) = -p.eval(bgeot::base_node(0.0, 1.0).begin());
	q = p; q.derivative(0);
	if (j == 0)
	  M(19, i+10*j) = q.eval(bgeot::base_node(0.0, 1.0).begin());
	if (j == 1)
	  M(19, i+10*j) = -q.eval(bgeot::base_node(0.0, 1.0).begin());
	q = p; q.derivative(1);
	if (j == 0)
	  M(20, i+10*j) = q.eval(bgeot::base_node(0.0, 1.0).begin());
	if (j == 1)
	  M(20, i+10*j) = -q.eval(bgeot::base_node(0.0, 1.0).begin());

	// raccord en (1/3, 1/3)
	double u_3 = 1.0 / 3.0;
	if (j == 0)
	  M(21, i+10*j) =  p.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 1)
	  M(21, i+10*j) = -p.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 0)
	  M(22, i+10*j) =  p.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 2)
	  M(22, i+10*j) = -p.eval(bgeot::base_node(u_3, u_3).begin());
	q = p; q.derivative(0);
	if (j == 0)
	  M(23, i+10*j) =  q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 1)
	  M(23, i+10*j) = -q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 0)
	  M(24, i+10*j) =  q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 2)
	  M(24, i+10*j) = -q.eval(bgeot::base_node(u_3, u_3).begin());
	q = p; q.derivative(1);
	if (j == 0)
	  M(25, i+10*j) =  q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 1)
	  M(25, i+10*j) = -q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 0)
	  M(26, i+10*j) =  q.eval(bgeot::base_node(u_3, u_3).begin());
	if (j == 2)
	  M(26, i+10*j) = -q.eval(bgeot::base_node(u_3, u_3).begin());

	// raccord en (1/6, 1/6)
	double u_6 = 1.0 / 6.0;
	q = p; q.derivative(0); q2 = p; q2.derivative(1); q -= q2;
	if (j == 1)
	  M(27, i+10*j) =  q.eval(bgeot::base_node(u_6, u_6).begin());
	if (j == 2)
	  M(27, i+10*j) = -q.eval(bgeot::base_node(u_6, u_6).begin());
	
	// raccord en (1/6, 2/3)
	double u_2 = 2.0 / 3.0;
	q = p; q.derivative(0); q2 = p;
	q2.derivative(1); q2 *= 2.0; q += q2;
	if (j == 1)
	  M(28, i+10*j) =  q.eval(bgeot::base_node(u_6, u_2).begin());
	if (j == 0)
	  M(28, i+10*j) = -q.eval(bgeot::base_node(u_6, u_2).begin());
	
	// raccord en (2/3, 1/6)
	q = p; q.derivative(0); q *= 2.0; q2 = p;
	q2.derivative(0); q += q2;
	if (j == 2)
	  M(29, i+10*j) =  q.eval(bgeot::base_node(u_2, u_6).begin());
	if (j == 0)
	  M(29, i+10*j) = -q.eval(bgeot::base_node(u_2, u_6).begin());
    }
    
    gmm::clean(M, 1E-10);
    cout << "M = " << M << endl;
    
    gmm::lu_inverse(M);
    
    gmm::clean(M, 1E-10);
    cout << "inv M = " << M << endl;
    
    cout.precision(11);
    
    bool latex = false;
    
    for (int i = 0; i < 9; ++i)
      for (int j = 0; j < 3; ++j) {
	bgeot::base_poly p(2,3);
	for (int k = 0; k < 10; ++k)
	  if (gmm::abs(M(k+10*j, i)) > 1E-8) p += base[k]*M(k+10*j, i);
      
	if (latex)
	  cout << "\\hat{\\varphi}_{" << i << "}^{" << j << "}(x,y) = ";
	else 
	  cout << "    \"";
	spec_print(cout, p);
	if (latex)
	  cout << ",\\\\" << endl;
	else
	  cout << ";\"\n";
    }

  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
