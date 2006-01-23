// -*- c++ -*- (enables emacs c++ mode)
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

	q = p; q.derivative(1); q2 = p; q2.derivative(0); q += q2;
	if (j == 0)
	  M( 9, i+10*j) = q.eval(bgeot::base_node(0.5, 0.5).begin());
	
	q = p; q.derivative(0);
	if (j == 1)
	  M(10, i+10*j) = -q.eval(bgeot::base_node(0.0, 0.5).begin());
	
	q = p; q.derivative(1);
	if (j == 2)
	  M(11, i+10*j) = -q.eval(bgeot::base_node(0.5, 0.0).begin());
	
	//
	// raccord internes
	//

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
	q = p; q.derivative(1); q2 = p; q2.derivative(0); q -= q2;
	if (j == 1)
	  M(27, i+10*j) =  q.eval(bgeot::base_node(u_6, u_6).begin());
	if (j == 2)
	  M(27, i+10*j) = -q.eval(bgeot::base_node(u_6, u_6).begin());
	
	// raccord en (1/6, 2/3)
	double u_2 = 2.0 / 3.0;
	q = p; q.derivative(1); q *= u_3; q2 = p;
	q2.derivative(0); q2 *= u_2; q += q2;
	if (j == 1)
	  M(28, i+10*j) =  q.eval(bgeot::base_node(u_6, u_2).begin());
	if (j == 0)
	  M(28, i+10*j) = -q.eval(bgeot::base_node(u_6, u_2).begin());
	
	// raccord en (2/3, 1/6)
	q = p; q.derivative(1); q *= u_2; q2 = p;
	q2.derivative(0); q2 *= u_3; q += q2;
	if (j == 2)
	  M(28, i+10*j) =  q.eval(bgeot::base_node(u_2, u_6).begin());
	if (j == 0)
	  M(28, i+10*j) = -q.eval(bgeot::base_node(u_2, u_6).begin());
    }
    
    gmm::clean(M, 1E-10);
    cout << "M = " << M << endl;
    
    gmm::lu_inverse(M);
    
    gmm::clean(M, 1E-10);
    cout << "inv M = " << M << endl;
    
    cout.precision(11);
    
    bool latex = true;
    
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 3; ++j) {
	bgeot::base_poly p(2,3);
	for (int k = 0; k < 10; ++k)
	  if (gmm::abs(M(k+10*j, i)) > 1E-8) p += base[k]*M(k+10*j, i);
      
	if (latex)
	  cout << "\\hat{\\varphi}_{" << i << "}^{" << j << "}(x,y) = ";
	else 
	  cout << "    \"";
	cout << p;
	if (latex)
	  cout << ",\\\\" << endl;
	else
	  cout << ";\"\n";
    }

  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
