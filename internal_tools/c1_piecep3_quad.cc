// -*- c++ -*- (enables emacs c++ mode)
// Little program which computes the base functions of the Argyris element
// on the reference element.

#include <iostream>
#include <gmm.h>
#include <getfem_config.h>

#define REDUCED 1

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
    bgeot::base_poly base[10];
    bgeot::base_matrix M(40, 40);
    int nbbase = 16;
    
    // base for P3
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
    
    for (int i = 0; i < 10; ++i) {
      bgeot::base_poly p = base[i], px = p, py = p;
      px.derivative(0); py.derivative(1);

      M( 0, i+10*1) =  p.eval(bgeot::base_node(0.0, 0.0).begin());
      M( 1, i+10*1) = px.eval(bgeot::base_node(0.0, 0.0).begin());
      M( 2, i+10*1) = py.eval(bgeot::base_node(0.0, 0.0).begin());
      M( 3, i+10*3) =  p.eval(bgeot::base_node(1.0, 0.0).begin());
      M( 4, i+10*3) = px.eval(bgeot::base_node(1.0, 0.0).begin());
      M( 5, i+10*3) = py.eval(bgeot::base_node(1.0, 0.0).begin());
      M( 6, i+10*1) =  p.eval(bgeot::base_node(0.0, 1.0).begin());
      M( 7, i+10*1) = px.eval(bgeot::base_node(0.0, 1.0).begin());
      M( 8, i+10*1) = py.eval(bgeot::base_node(0.0, 1.0).begin());
      M( 9, i+10*2) =  p.eval(bgeot::base_node(1.0, 1.0).begin());
      M(10, i+10*2) = px.eval(bgeot::base_node(1.0, 1.0).begin());
      M(11, i+10*2) = py.eval(bgeot::base_node(1.0, 1.0).begin());
	
#if !REDUCED
      
      nbbase = 16;
      M(12, i+10*0) = px.eval(bgeot::base_node(1.0, 0.5).begin());
      M(13, i+10*1) = -px.eval(bgeot::base_node(0.0, 0.5).begin());
      M(14, i+10*2) = py.eval(bgeot::base_node(0.5, 1.0).begin());
      M(15, i+10*3) = -py.eval(bgeot::base_node(0.5, 0.0).begin());

#else

      nbbase = 12;
      M(12, i+10*0) = px.eval(bgeot::base_node(1.0, 0.5).begin())
	- (px.eval(bgeot::base_node(1.0, 0.0).begin())
	   + px.eval(bgeot::base_node(1.0, 1.0).begin())) / 2.0;
      M(13, i+10*1) = px.eval(bgeot::base_node(0.0, 0.5).begin())
	- (px.eval(bgeot::base_node(0.0, 0.0).begin())
	   + px.eval(bgeot::base_node(0.0, 1.0).begin())) / 2.0;
      M(14, i+10*2) = py.eval(bgeot::base_node(0.5, 1.0).begin())
	- (py.eval(bgeot::base_node(0.0, 1.0).begin())
	   + py.eval(bgeot::base_node(1.0, 1.0).begin())) / 2.0;
      M(15, i+10*3) = py.eval(bgeot::base_node(0.5, 0.0).begin())
	- (py.eval(bgeot::base_node(0.0, 0.0).begin())
	   + py.eval(bgeot::base_node(1.0, 0.0).begin())) / 2.0;
#endif	

      //
      // raccord internes
      //
      
      // raccord en (0.0, 0.0)
      M(16, i+10*1) = p.eval(bgeot::base_node(0.0, 0.0).begin());
      M(16, i+10*3) = -p.eval(bgeot::base_node(0.0, 0.0).begin());
      M(17, i+10*1) = px.eval(bgeot::base_node(0.0, 0.0).begin());
      M(17, i+10*3) = -px.eval(bgeot::base_node(0.0, 0.0).begin());
      M(18, i+10*1) = py.eval(bgeot::base_node(0.0, 0.0).begin());
      M(18, i+10*3) = -py.eval(bgeot::base_node(0.0, 0.0).begin());

      // raccord en (1.0, 0.0)
      M(19, i+10*0) = p.eval(bgeot::base_node(1.0, 0.0).begin());
      M(19, i+10*3) = -p.eval(bgeot::base_node(1.0, 0.0).begin());
      M(20, i+10*0) = px.eval(bgeot::base_node(1.0, 0.0).begin());
      M(20, i+10*3) = -px.eval(bgeot::base_node(1.0, 0.0).begin());
      M(21, i+10*0) = py.eval(bgeot::base_node(1.0, 0.0).begin());
      M(21, i+10*3) = -py.eval(bgeot::base_node(1.0, 0.0).begin());
      
      // raccord en (0.0, 1.0)
      M(22, i+10*2) = p.eval(bgeot::base_node(0.0, 1.0).begin());
      M(22, i+10*1) = -p.eval(bgeot::base_node(0.0, 1.0).begin());
      M(23, i+10*2) = px.eval(bgeot::base_node(0.0, 1.0).begin());
      M(23, i+10*1) = -px.eval(bgeot::base_node(0.0, 1.0).begin());
      M(24, i+10*2) = py.eval(bgeot::base_node(0.0, 1.0).begin());
      M(24, i+10*1) = -py.eval(bgeot::base_node(0.0, 1.0).begin());
      
      // raccord en (1.0, 1.0)
      M(25, i+10*0) = p.eval(bgeot::base_node(1.0, 1.0).begin());
      M(25, i+10*2) = -p.eval(bgeot::base_node(1.0, 1.0).begin());
      M(26, i+10*0) = px.eval(bgeot::base_node(1.0, 1.0).begin());
      M(26, i+10*2) = -px.eval(bgeot::base_node(1.0, 1.0).begin());
      M(27, i+10*0) = py.eval(bgeot::base_node(1.0, 1.0).begin());
      M(27, i+10*2) = -py.eval(bgeot::base_node(1.0, 1.0).begin());
      
      
      // raccord en (0.5, 0.5)
      M(28, i+10*0) = p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(28, i+10*1) = -p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(29, i+10*0) = px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(29, i+10*1) = -px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(30, i+10*0) = py.eval(bgeot::base_node(0.5, 0.5).begin());
      M(30, i+10*1) = -py.eval(bgeot::base_node(0.5, 0.5).begin());
      M(31, i+10*0) = p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(31, i+10*2) = -p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(32, i+10*0) = px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(32, i+10*2) = -px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(33, i+10*0) = py.eval(bgeot::base_node(0.5, 0.5).begin());
      M(33, i+10*2) = -py.eval(bgeot::base_node(0.5, 0.5).begin());
      M(34, i+10*0) = p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(34, i+10*3) = -p.eval(bgeot::base_node(0.5, 0.5).begin());
      M(35, i+10*0) = px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(35, i+10*3) = -px.eval(bgeot::base_node(0.5, 0.5).begin());
      M(36, i+10*0) = py.eval(bgeot::base_node(0.5, 0.5).begin());
      M(36, i+10*3) = -py.eval(bgeot::base_node(0.5, 0.5).begin());
      
      // raccord en (0.25, 0.25)
      M(37, i+10*1) = (px-py).eval(bgeot::base_node(0.25, 0.25).begin());
      M(37, i+10*3) = -(px-py).eval(bgeot::base_node(0.25, 0.25).begin());
      
      // raccord en (0.75, 0.75)
      M(38, i+10*0) = (px-py).eval(bgeot::base_node(0.75, 0.75).begin());
      M(38, i+10*2) = -(px-py).eval(bgeot::base_node(0.75, 0.75).begin());
      
      // raccord en (0.25, 0.75)
      M(39, i+10*1) = (px+py).eval(bgeot::base_node(0.25, 0.75).begin());
      M(39, i+10*2) = -(px+py).eval(bgeot::base_node(0.25, 0.75).begin());
      
//    // raccord en (0.75, 0.25) non nécessaire
//    M(40, i+10*0) = (px+py).eval(bgeot::base_node(0.75, 0.25).begin());
//    M(40, i+10*3) = -(px+py).eval(bgeot::base_node(0.75, 0.25).begin());
    }

    gmm::clean(M, 1E-13);
    cout << "M = " << M << endl;
    
    double det = gmm::lu_det(M);
    cout << "det = " << det << endl;
    
    if (gmm::abs(det) < 1e-15) {
      cout << "Non invertible matrix, non-unisolvant finite element" << endl;
      bgeot::base_matrix MM = M, Q = M;
      std::vector<std::complex<double> > eigval(gmm::mat_ncols(M));
      gmm::implicit_qr_algorithm(MM, eigval, Q);
      cout << "eigval : " << eigval << endl;
      exit(1);
    }

    gmm::lu_inverse(M);
    gmm::clean(M, 1E-10);
    cout << "inv M = " << M << endl;
    
    cout.precision(13);
    
    bool latex = false;
    
    for (int i = 0; i < nbbase; ++i) {
      for (int j = 0; j < 4; ++j) {
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
      cout << endl;
    }

  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
