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
#include <bgeot_poly.h>

typedef bgeot::polynomial<double> base_poly;

int main(void)
{
  try {

    bgeot::polynomial<double> W, Z; W[0] = 1.0;
    Z[0] = 2.0;
    cout << "rd = " << W.real_degree() << endl;
    cout << "W = " << W << endl;
    W.direct_product(Z);
    cout << "W = " << W << endl;

    bgeot::polynomial<double> P(2,2), Q(2,2);
    P[2] = 1.0;
    Q[3] = 2.0;
    cout << "Le nombre de monomes de P est " << P.size() << endl;
    cout << "P = " << P << endl;
    cout << "Q = " << Q << endl;
    cout << "P + Q = " << P + Q << endl;
    cout << "P * 2.0 * Q = " << P * 2.0 * Q << endl;
    bgeot::polynomial<double> R = P * Q;
    cout << "Le nombre de monomes de R est " << R.size() << endl;
    cout << "Le degre de R est " << R.degree() << endl;
    P.direct_product(Q);
    cout << "Produit direct de P et Q : " << P << endl;
    
    Z.direct_product(P);
    cout << "Produit direct de P et Z : " << Z << endl;    

    P = bgeot::polynomial<double>(3,1,1);
    P *= 3.0; P *= bgeot::polynomial<double>(3,1,0);
    P += bgeot::polynomial<double>(3,1,1);
    P *= bgeot::polynomial<double>(3,1,2);
    P += bgeot::polynomial<double>(3,1,0);
    cout << "P = " << P << endl;
    
    
    double tab[3]; tab[0] = 1.0; tab[1] = 2.0; tab[2] = -1.0;
    
    cout << "P(1.0, 2.0) = " << P.eval(&(tab[0])) << endl;
    
    Q = P;
    P *= Q;
    cout << "PP = " << P << "\n";
    P.derivative(0);
    cout << "PP.derivative(0)=" << P << "\n";

    bgeot::power_index p(3);
    for (int i=0; i < 20; ++i, ++p) {
      cout << "i=" << i << ", p=";
      for (unsigned k=0; k < p.size(); ++k) cout << p[k] << " ";
      cout << "degree=" << p.degree() << ", global_index(p)=" << p.global_index() << "\n";      
    }

    base_poly S(1,2); S[0] = -2; S[1] = 3; S[2] = 1;
    cout << "P=" << P << ", S=" << S << " \n";
    cout << "P(S,x)=" << bgeot::poly_substitute_var(P,S,0) << "\n";
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
