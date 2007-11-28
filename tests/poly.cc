// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================
#include "getfem/bgeot_poly.h"

std::string horner_print(bgeot::short_type degree, bgeot::power_index &mi,
			 bgeot::short_type k, bgeot::short_type de) {
  char s[1024];
  const char *xyz = "xyzabcdefghijklmnop";
  if (k == 0) {
    sprintf(s, "P[%d]", int(mi.global_index()));
    return s;
  } else {
    std::string str;
    //T v = (*(it+k-1)), res = T(0);
    for (mi[k-1] = degree-de; mi[k-1] != bgeot::short_type(-1); (mi[k-1])--) {
      //res = horner(mi, k-1, de + mi[k-1], it) + v * res;
      if (str.size())
	sprintf(s, "%s + %c*(%s)",
		horner_print(degree, mi,k-1,de+mi[k-1]).c_str(), xyz[k-1],
		str.c_str());
      else 
	sprintf(s, "%s", horner_print(degree, mi,k-1,de+mi[k-1]).c_str());
      str = s;
    }
    mi[k-1] = 0;
    return str;
  }
}


void dump_poly_eval() {
  for (unsigned dim = 1; dim <= 3; ++dim) {
    cout << "case " << dim << ": {\n";
    for (unsigned k=0; k < dim; ++k) {
      cout << "T " << "xyzZ"[k] << " = it[" << k << "];\n";
    }
    for (unsigned dg=2; dg <= 6; ++dg) {
      cout << "  if (deg == " << dg << ") ";
      bgeot::power_index mi(dim);
      cout << "    return " << horner_print(dg, mi, dim, 0) << ";\n";
    }
    cout << "} break;\n";
  }
}

int main(void)
{
  try {

    bgeot::base_poly W, Z; W[0] = 1.0;
    Z[0] = 2.0;
    cout << "rd = " << W.real_degree() << endl;
    cout << "W = " << W << endl;
    W.direct_product(Z);
    cout << "W = " << W << endl;

    bgeot::base_poly P(2,2), Q(2,2);
    P[2] = 1.0;
    Q[3] = 2.0;
    cout << "Le nombre de monomes de P est " << P.size() << endl;
    cout << "P = " << P << endl;
    cout << "Q = " << Q << endl;
    cout << "P + Q = " << P + Q << endl;
    cout << "P * 2.0 * Q = " << P * 2.0 * Q << endl;
    bgeot::base_poly R = P * Q;
    cout << "Le nombre de monomes de R est " << R.size() << endl;
    cout << "Le degre de R est " << R.degree() << endl;
    P.direct_product(Q);
    cout << "Produit direct de P et Q : " << P << endl;
    
    Z.direct_product(P);
    cout << "Produit direct de P et Z : " << Z << endl;    

    P = bgeot::base_poly(3,1,1);
    P *= 3.0; P *= bgeot::base_poly(3,1,0);
    P += bgeot::base_poly(3,1,1);
    P *= bgeot::base_poly(3,1,2);
    P += bgeot::base_poly(3,1,0);
    cout << "P = " << P << " : degree=" << P.degree() << endl;
    
    
    bgeot::opt_long_scalar_type tab[3];
    tab[0] = 1.0; tab[1] = 2.0; tab[2] = -1.0;
    
    cout << "P(1.0, 2.0) = " << P.eval(&(tab[0])) << endl;

    for (unsigned dg=0; dg <= 6; ++dg) {
      for (unsigned dim=0; dim <= 3; ++dim) {
	bgeot::base_poly PP(dim, dg);
	for (unsigned i=0; i < PP.size(); ++i) 
	  PP[i] = bgeot::opt_long_scalar_type(rand())
	    / bgeot::opt_long_scalar_type(RAND_MAX);
	std::vector<bgeot::opt_long_scalar_type> X(dim); 
	for (unsigned i=0; i < dim; ++i) X[i] = 
	  bgeot::opt_long_scalar_type(rand())
	  / bgeot::opt_long_scalar_type(RAND_MAX);
	bgeot::opt_long_scalar_type a = PP.eval(X.begin());
	bgeot::power_index mi(dim);
	bgeot::opt_long_scalar_type b = PP.horner(mi,dim,0,X.begin());

	cout << "[d=" << dim << ", dg=" << PP.degree() 
	  //<< ", P=" << PP << " -> " 
	     << a << " == " << b 
	     << "?\n";
	assert(gmm::abs(a-b) < 1e-14);
	
	//cout << "Horner: " << PP.horner_print(mi,dim,0) << "\n";
      }
    }
    cout << "\n--------------------------------------------------------\n";
    dump_poly_eval();
    cout << "\n--------------------------------------------------------\n";

    
    Q = P;
    P *= Q;
    cout << "PP = " << P << "\n";
    P.derivative(0);
    cout << "PP.derivative(0)=" << P << "\n";

    bgeot::power_index p(3);
    for (int i=0; i < 20; ++i, ++p) {
      cout << "i=" << i << ", p=";
      for (unsigned k=0; k < p.size(); ++k) cout << p[k] << " ";
      cout << "degree=" << p.degree() << ", global_index(p)="
	   << p.global_index() << "\n";      
    }

    bgeot::base_poly S(1,2); S[0] = -2; S[1] = 3; S[2] = 1;
    cout << "P=" << P << ", S=" << S << " \n";
    cout << "P(S,x)=" << bgeot::poly_substitute_var(P,S,0) << "\n";
    bgeot::opt_long_scalar_type t0 = gmm::uclock_sec();
    std::vector<bgeot::opt_long_scalar_type> v(3);
    for (unsigned i=0; i < 100000; ++i) {
      for (unsigned k=0; k < v.size(); ++k)
	v[k] =  bgeot::opt_long_scalar_type(rand())
	  / bgeot::opt_long_scalar_type(RAND_MAX);
      P.eval(v.begin());
    }
    cout << "poly eval : " << gmm::uclock_sec() - t0 << "sec \n";
    bgeot::base_poly QQ(P); QQ.derivative(1); QQ.derivative(2);
    cout << "QQ=" << QQ << "\n";
    for (unsigned i=0; i < 100000; ++i) {
      QQ.eval(v.begin());
    }
    cout << "poly eval : " << gmm::uclock_sec() - t0 << "sec \n";

    t0 = gmm::uclock_sec();
    bgeot::opt_long_scalar_type z=0;
    for (unsigned i=0; i < 100000; ++i) {
      bgeot::base_poly P2(P);
      for (unsigned k=0; k < P.dim(); ++k) { 
        P2.derivative(k); z += P2[0];
      }
    }
    cout << "poly derivative : " << gmm::uclock_sec() - t0 << "sec\n";

    bgeot::base_poly P2;
    std::stringstream ss; ss << P;
    P2 = bgeot::read_base_poly(P.dim(), ss);
    cout << "P=" << P << "\nread_base_poly=" << P2 << "\n";
    assert(P == P2);
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;
}
