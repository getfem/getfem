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
  bgeot::base_poly one(2, 0), x(2, 1, 0), y(2, 1, 1); one.one();
  bgeot::base_poly base[21];
  // base for P5
  base[ 0] = one;
  base[ 1] = x;
  base[ 2] = y;
  base[ 3] = x*x;
  base[ 4] = x*y;
  base[ 5] = y*y;
  base[ 6] = x*x*x;
  base[ 7] = x*x*y;
  base[ 8] = x*y*y;
  base[ 9] = y*y*y;
  base[10] = x*x*x*x;
  base[11] = x*x*x*y;
  base[12] = x*x*y*y;
  base[13] = x*y*y*y;
  base[14] = y*y*y*y;
  base[15] = x*x*x*x*x;
  base[16] = x*x*x*x*y;
  base[17] = x*x*x*y*y;
  base[18] = x*x*y*y*y;
  base[19] = x*y*y*y*y;
  base[20] = y*y*y*y*y;


  bgeot::base_matrix M(21, 21);

  for (int i = 0; i < 21; ++i) {
    bgeot::base_poly p = base[i];
    M(0, i) = p.eval(bgeot::base_node(0.0, 0.0).begin());
    M(1, i) = p.eval(bgeot::base_node(1.0, 0.0).begin());
    M(2, i) = p.eval(bgeot::base_node(0.0, 1.0).begin());
  
    bgeot::base_poly q = p; q.derivative(0);
    M(3, i) = q.eval(bgeot::base_node(0.0, 0.0).begin());
    M(4, i) = q.eval(bgeot::base_node(1.0, 0.0).begin());
    M(5, i) = q.eval(bgeot::base_node(0.0, 1.0).begin());
  
    q = p; q.derivative(1);
    M(6, i) = q.eval(bgeot::base_node(0.0, 0.0).begin());
    M(7, i) = q.eval(bgeot::base_node(1.0, 0.0).begin());
    M(8, i) = q.eval(bgeot::base_node(0.0, 1.0).begin());

    q = p; q.derivative(0); q.derivative(0);
    M(9, i) = q.eval(bgeot::base_node(0.0, 0.0).begin());
    M(10, i) = q.eval(bgeot::base_node(1.0, 0.0).begin());
    M(11, i) = q.eval(bgeot::base_node(0.0, 1.0).begin());
    
    q = p; q.derivative(1); q.derivative(0);
    M(12, i) = q.eval(bgeot::base_node(0.0, 0.0).begin());
    M(13, i) = q.eval(bgeot::base_node(1.0, 0.0).begin());
    M(14, i) = q.eval(bgeot::base_node(0.0, 1.0).begin());
    
    q = p; q.derivative(1); q.derivative(1);
    M(15, i) = q.eval(bgeot::base_node(0.0, 0.0).begin());
    M(16, i) = q.eval(bgeot::base_node(1.0, 0.0).begin());
    M(17, i) = q.eval(bgeot::base_node(0.0, 1.0).begin());
    

    q = p; q.derivative(0);  bgeot::base_poly r = p; r.derivative(1);
    M(18, i) = (q.eval(bgeot::base_node(0.5, 0.5).begin())
		+ r.eval(bgeot::base_node(0.5, 0.5).begin())) / ::sqrt(2.0);

    q = p; q.derivative(0);
    M(19, i) = -q.eval(bgeot::base_node(0.0, 0.5).begin());

    q = p; q.derivative(1);
    M(20, i) = -q.eval(bgeot::base_node(0.5, 0.0).begin());

  }

  gmm::clean(M, 1E-10);
  cout << "M = " << M << endl;

  gmm::lu_inverse(M);

  gmm::clean(M, 1E-10);
  cout << "inv M = " << M << endl;

  cout.precision(13);
  
  for (int i = 0; i < 21; ++i) {
    bgeot::base_poly p(2,5);
    for (int j = 0; j < 21; ++j)
      if (gmm::abs(M(j, i)) > 1E-8) p += base[j]*M(j, i);

    cout << "\\bar{\\varphi}_{" << i << "}(x,y) = " << p << ",\\\\" << endl;
  }

  for (int i = 0; i < 21; ++i) {
    bgeot::base_poly p(2,5);
    for (int j = 0; j < 21; ++j)
      if (gmm::abs(M(j, i)) > 1E-8) p += base[j]*M(j, i);

    cout << "base_[" << i << "]=";
    poly_cpp_display(cout, p);
    cout << ";\n";
  }

  return 0;
}
