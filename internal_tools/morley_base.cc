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
      if (gmm::abs(a - double(i)/j)<1e-10) {
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

void morley2(void) {
  bgeot::base_poly one(2, 0), x(2, 1, 0), y(2, 1, 1); one.one();
  bgeot::base_poly base[6];
  // base for P5
  base[ 0] = one;
  base[ 1] = x;
  base[ 2] = y;
  base[ 3] = x*x;
  base[ 4] = x*y;
  base[ 5] = y*y;

  bgeot::base_matrix M(6, 6);

  for (int i = 0; i < 6; ++i) {
    bgeot::base_poly p = base[i], q;
    M(0, i) = p.eval(bgeot::base_node(0.0, 0.0).begin());
    M(1, i) = p.eval(bgeot::base_node(1.0, 0.0).begin());
    M(2, i) = p.eval(bgeot::base_node(0.0, 1.0).begin());

    q = p; q.derivative(0);  bgeot::base_poly r = p; r.derivative(1);
    M(3, i) = (q.eval(bgeot::base_node(0.5, 0.5).begin())
		+ r.eval(bgeot::base_node(0.5, 0.5).begin())) / ::sqrt(2.0);

    q = p; q.derivative(0);
    M(4, i) = -q.eval(bgeot::base_node(0.0, 0.5).begin());

    q = p; q.derivative(1);
    M(5, i) = -q.eval(bgeot::base_node(0.5, 0.0).begin());

  }

  gmm::clean(M, 1E-10);
  cout << "M = " << M << endl;

  gmm::lu_inverse(M);

  gmm::clean(M, 1E-10);
  cout << "inv M = " << M << endl;

  cout.precision(13);
 
  cout << "Morley in dimension 2: \n";
 
  for (int i = 0; i < 6; ++i) {
    bgeot::base_poly p(2,5);
    for (int j = 0; j < 6; ++j)
      if (gmm::abs(M(j, i)) > 1E-8) p += base[j]*M(j, i);

    cout << "base_[" << i << "]="; spec_print(cout, p);
    cout << ";\n";
  }

 
}








int main(void) {
  morley2();
  return 0;
}
